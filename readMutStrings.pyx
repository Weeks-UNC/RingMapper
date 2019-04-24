#cython: boundscheck=False, wraparound=False


########################################################
#   This code  contains optimized cython functions for reading
#   classified mutations string from new and old shapemapper
#
#   Part of the ringmapper pipeline
#
#   This file is licensed under the terms of the MIT license
#
#   Lead developer: Anthony Mustoe
#   Contributors: Nicole Lama, Steve Busan
#   
#   Version 4.3
#   February 2019
#
########################################################


from libc.stdio cimport FILE, fopen, fclose, getline
from libc.stdlib cimport atoi, free
from libc.string cimport strsep, strtok, strcmp, strdup

import numpy as np
cimport numpy as np



###########################################################################
# Define READ object

cdef struct READ:
    int start
    int stop
    char* read
    char* muts
    char* qual


###########################################################################
#
# Fill mutation matrices
#
###########################################################################


def fillMatrices(str inputFile, int[:,::1] read_arr, int[:, ::1] comut_arr, int[:, ::1] inotj_arr,
                 int window, int mincoverage, int fileformat, int undersample):
    """ inputFile = path to the classified mutations files
    read_arr 
    comut_arr 
    inotj_arr    =  NxN numpy matrices to fill
    window       =  window >=1 over which to compute correlations
    mincoverage  =  exclude reads with coverage less than this 
                    (i.e. must have at least this many positions as 1 in read string)
    fileformat   =  Shapemapper mut string format (1,2,3)
    undersample  =  undersample this number of reads (without replacement) from read file
                    Default is -1, which reads all lines

    returns (totreads, readreads)
    """ 
    
    cdef int maxindex = read_arr.shape[0]

    # working arrays storing read events
    # first position contains element counter
    cdef np.int32_t[:] readnts = np.zeros(maxindex+1, dtype=np.int32)
    cdef np.int32_t[:] mutnts  = np.zeros(maxindex+1, dtype=np.int32)
    

    # mem alloc for reading the file using c
    cdef FILE* cfile
    cdef size_t lsize = 0
    cdef ssize_t endfile
    cdef char* line = NULL
    
    # mem alloc for read parsing and for loops   
    cdef READ r
    cdef int i, j, i_index
    
    # mem alloc for undersampling index array 
    cdef np.int32_t[:] readidxArray = np.zeros(max(0,undersample), dtype=np.int32)
    cdef int ridx = 0
    cdef int sufficientreads = 0


    # initialize undersample array if undersampling specified
    if undersample > 0:
        
        sufficientreads = undersampleIndices(inputFile, fileformat, undersample, readidxArray)
        
        if sufficientreads < 1:
            print('WARNING: Cannot undersample file {0}, which has only has {1} non-empty reads'.format(inputFile, -1*sufficientreads))
            undersample = -1 # turn off undersampling



    # open the file
    cfile = fopen(inputFile, "r")
    if cfile == NULL:
        raise IOError(2, "No such file or directory: '{0}'".format(inputFile))


    cdef int linenum = -1
    cdef int readreads = 0
    cdef int skipped_reads = 0

    # iterate through lines
    while True:
        
        linenum += 1 
        
        # get the line of text using c-function
        endfile = getline(&line, &lsize, cfile)
        if endfile == -1:
            break
        

        # logic to randomly undersample reads (if specified)
        if undersample > 0:
            if ridx >= undersample:
                break
            elif readidxArray[ridx] == linenum:
                ridx += 1
            else:
                continue
            


        # parse line into individual values, reset and fill readnts/mutnts
        try:
            r = parseLine(line, fileformat)
            if r.read == NULL:
                raise IndexError()
            elif r.stop >= maxindex:
                print "Skipping line {0} with out-of-array-bounds = ({1}, {2})".format(linenum, r.start, r.stop)
                continue

            fillReadMut(readnts, mutnts, r, window, mincoverage)
 
        except:
            skipped_reads += 1
            #print "Skipping incorrectly formatted line {0}".format(linenum)
            continue
        
        # check if read was read, and if so, increment counter
        if readnts[0]>1:
            readreads += 1
        

        # fill the read count matrix; first element of readnts is +1 of last index
        for i in xrange(1, readnts[0]):
            i_index = readnts[i]
            for j in xrange(i, readnts[0]):
                read_arr[ i_index, readnts[j] ] += 1
        
        
        # fill in the mut matrices; first element of mutnts is +1 of last index
        for i in xrange(1, mutnts[0]):
        
            i_index = mutnts[i]
            
            for j in xrange(i, mutnts[0]):
                comut_arr[ i_index, mutnts[j] ] += 1


            # fill in inotj
            # Note this loop overcounts for j=mutated
            # Diagonal is not used, so don't worry about i=j case
            for j in xrange(1, readnts[0]): 
                inotj_arr[ i_index, readnts[j] ] += 1

            # correct for the over addition in inotj in above loop
            for j in xrange(1, mutnts[0]):
                inotj_arr[ i_index, mutnts[j] ] -= 1
        

    fclose(cfile)
    
    if skipped_reads > 0:
        print("skipped {} lines".format(skipped_reads))


    return linenum, readreads
 


       
##################################################################################

def fillIndependentProbArrays(str inputFile, int[:,::1] countArray,
                              int window, int mincoverage, int fileformat, int undersample):

    """inputFile = path to the classified mutations files
    totalCount   =  Nx8 array containing modification countss
    window       =  window >=1 over which to compute correlations
    mincoverage  =  exclude reads with coverage less than this 
                    (i.e. must have at least this many positions as 1 in read string)
    
    fileformat   =  Shapemapper mut string format (1,2,3)
    undersample  =  undersample this number of reads (without replacement) from read file
                    Default is -1, which reads all lines
    
    returns (totreads, readreads)
    """
    
    # establish counting dimensions and initialize internal countArray    
    cdef int maxindex = countArray.shape[0]
    cdef int maxcount = countArray.shape[1]/2-1
    cdef int[:,:] workingCountArray = np.zeros((maxindex, 2*(maxcount+2)), dtype=np.int32)


    # working arrays storing read events
    # first position contains element counter
    cdef np.int32_t[:] readnts = np.zeros(maxindex+1, dtype=np.int32)
    cdef np.int32_t[:] mutnts  = np.zeros(maxindex+1, dtype=np.int32)
    

    # mem alloc for reading the file using c
    cdef FILE* cfile
    cdef size_t lsize = 0
    cdef ssize_t endfile
    cdef char* line = NULL
    
    # mem alloc for read parsing and for loops   
    cdef READ r
    cdef int i, j, i_index, totalmuts
    

    # mem alloc for undersampling index array 
    cdef np.int32_t[:] readidxArray = np.zeros(max(0,undersample), dtype=np.int32)
    cdef int ridx = 0
    cdef int sufficientreads = 0


    # initialize undersample array if undersampling specified
    if undersample > 0:
        
        sufficientreads = undersampleIndices(inputFile, fileformat, undersample, readidxArray)
        
        if sufficientreads < 1:
            print('WARNING: Cannot undersample file {0}, which has only has {1} non-empty reads'.format(inputFile, -1*sufficientreads))
            undersample = -1 # turn off undersampling



    # open the file
    cfile = fopen(inputFile, "r")
    if cfile == NULL:
        raise IOError(2, "No such file or directory: '{0}'".format(inputFile))


    cdef int linenum = -1
    cdef int readreads = 0
    cdef int skipped_reads = 0

    # iterate through lines
    while True:
        
        linenum += 1 
        
        # get the line of text using c-function
        endfile = getline(&line, &lsize, cfile)
        if endfile == -1:
            break
        

        # logic to randomly undersample reads (if specified)
        if undersample > 0:
            if ridx >= undersample:
                break
            elif readidxArray[ridx] == linenum:
                ridx += 1
            else:
                continue
            


        # parse line into individual values, reset and fill readnts/mutnts
        try:
            r = parseLine(line, fileformat)
            if r.read == NULL:
                raise IndexError()
            elif r.stop >= maxindex:
                print "Skipping line {0} with out-of-array-bounds = ({1}, {2})".format(linenum, r.start, r.stop)
                continue
            
            fillReadMut(readnts, mutnts, r, window, mincoverage)
 
        except:
            skipped_reads += 1
            #print "Skipping incorrectly formatted line {0}".format(linenum)
            continue
        

        # check if read was read, and if so, increment counter
        if readnts[0]>1:
            readreads += 1
        
        
        # compute totalmuts
        totalmuts = 0
        for i in xrange( r.stop - r.start + 1):
            totalmuts += r.muts[i]-48 # will be 0 if no mut, 1 if mut
        

        if totalmuts > maxcount+1:
            totalmuts = maxcount+1
        

        for i in xrange(1, readnts[0]):
            workingCountArray[readnts[i], totalmuts] +=1

        for i in xrange(1, mutnts[0]):
        
            # need to subtract since mutated
            workingCountArray[mutnts[i], totalmuts] -= 1
            
            # now need to add to mutated
            workingCountArray[mutnts[i], maxcount+1+totalmuts] += 1
            

    fclose(cfile)


    if skipped_reads > 0:
        print("skipped {} lines".format(skipped_reads))

    
    # transfer workingCountArray to countArray
    countArray[:,:maxcount+1] = workingCountArray[:,:maxcount+1]
    for i in xrange(maxindex+1):
        countArray[i,maxcount] += workingCountArray[i,maxcount+1]
        
    countArray[:,maxcount+1:] = workingCountArray[:,maxcount+2:2*maxcount+3]
    # note don't need to transfer the last column, since will always be 0


    return linenum, readreads
 





##################################################################################

cdef READ parseLine(char* line, int fileformat):
    # accessory function for fillMatrices
    # uses c functions to read through line; avoids having to use python split
    # which outputs python array and slows things down

    cdef READ r
    cdef int i

    cdef int leadersize
    if fileformat==3:
        leadersize=2
    else:
        leadersize=fileformat
    

    # init token values
    cdef char* running = line
    cdef char* token = strsep(&running, " \t")
    cdef int tokindex = 1
    cdef int valueset = 0
    

    # pop off leader information so that token = r.start
    for i in xrange(leadersize):
        token = strsep(&running, " \t")
    
   
    while token:
        
        if tokindex == 1:
            r.start = atoi(token)

        elif tokindex == 2:
            r.stop = atoi(token)
        
        elif fileformat==3:
            if tokindex==3 and strcmp(token, "INCLUDED") != 0:
                break
            elif tokindex == 6:
                r.read = token
            elif tokindex == 7:
                r.muts = token
                valueset = 1
                break
        else:
            if tokindex == 3:
                r.read = token
            elif tokindex == 4:
                r.muts = token
                valueset = 1
                break
        
        token = strsep(&running, " \t")
        tokindex+=1
    

    if valueset != 1 or r.start >= r.stop:
        r.start = 0
        r.stop = 0
        r.read = NULL
    

    return r



##################################################################################


cdef int undersampleIndices(str inputFile, int fileformat, int undersample, np.int32_t[:] readidxArray):
    """This function will fill readidxArray with valid reads from inputFile
    
    inputFile = mutstring file
    fileformat   =  Shapemapper mut string format (1,2,3)
    undersample  =  undersample this number of reads (without replacement) from read file
                    Default is -1, which reads all lines
    readidxArray = preallocated array with size = undersample

    inputFile is checked to make sure there are sufficient valid reads to sample    
    Return 1 if sufficient reads (and undersampling performed)
    Return <0 if not sufficient reads (actual value is -1 * validReads)
    """


    # mem alloc for reading the file using c
    cdef FILE* cfile
    cdef size_t lsize = 0
    cdef ssize_t endfile
    cdef char* line = NULL
    cdef READ r
    cdef int i    

    cdef int linenum = -1    
    cdef int validlines_idx = 0 
 
    # init the numpy array to hold valid line idxes    
    cdef int filelength = countReadLines(inputFile)

    validlines = np.zeros(filelength, dtype=np.int32)


    cfile = fopen(inputFile, "r")
    if cfile == NULL:
        raise IOError(2, "No such file or directory: '{0}'".format(inputFile))
    

    # iterate through lines
    while True:
        
        linenum += 1 
        
        # get the line of text using c-function
        endfile = getline(&line, &lsize, cfile)
        if endfile == -1:
            break
        
        # parse line and check whether it is valid       
        r = parseLine(line, fileformat)
        if r.read != NULL:
            validlines[validlines_idx] = linenum
            validlines_idx += 1                      
    
    # check to make sure there are sufficient valid reads to undersample
    if validlines_idx < undersample:
        return -validlines_idx
    
    # resize the validlines array to exact size
    validlines = validlines[:validlines_idx]   
    
    # randomly undersample without replacement and sort
    idxArray = np.random.choice(validlines, undersample, replace=False)
    idxArray.sort()

    # assign readidxArray    
    for i in range(undersample):    
        readidxArray[i] = idxArray[i]
    
    return 1    




##################################################################################

cdef int countReadLines(str inputFile):
    # Return the number of lines in inputFile
    
    cdef str line
    cdef int nlines = 0
    
    with open(inputFile) as inp:
        for line in inp:
            nlines += 1
    
    return nlines


##################################################################################


cdef int fillReadMut(np.int32_t[:] readnts, np.int32_t[:] mutnts, READ r, int window, int mincoverage) except -1:
    # accessory function for fillMatrices
    # parse the read for mut events
    # readnts and mutnts contain seq indices of valid reads and muts, respectively
    # Element 0 of readnts and mutnts is used as a counter...
    # Reads with coverage less than mincoverage are ignored
    #
    #
    # Note that there are a few ways we could treat the data for windowed correlations
    # 1) Require full coverage across the window
    #    -this is unrealistic for modifications, since we frequently expect modifications
    #     to be a del or multimismathch that results in 'noread' events 5' of the modification
    #
    # 2) Require partial coverage across the window
    #    -Problem of what counts as partial? 
    #    -Again have issue of muts likely introducing 'noread' events
    #
    # 3) Require only single position to be covered in window
    #    -fairest for modifications
    #    -fairest for regions  5' of highly modified nts that have a lot of 'noread' events
    #    -in general, windows should have low mutation rates, so assuming 1 'read' event
    #     in absence of other data indicates 'no mutation' is most likely
    #    -Note that the best way to treat this would be imputing based on mutation rates
    #
    # --> This code uses option 3

    
    cdef int step = window-1
    
    cdef int endseq = len(readnts)-1-window

    cdef int i,j

    cdef int readindex = 1 
    cdef int mutindex = 1 
    
    cdef int rcount = 0
    cdef int mcount = 0
    cdef int validpos = 0

    # reset index of arrays
    mutnts[0] = 0
    readnts[0] = 0
    
    # if read is too short, abort   
    if r.stop-r.start+1 < mincoverage:
        return 1
    

    # Handle the beginning of the read
    # For window>1, need to consider that mutation/read events at the beginning
    # of the read contribute to mut/reads at upstream positions
    # In other words, if read starts at nt 20 and window=3, then we need to increment
    # nt 18 and 19, which partially overlap 20
    # For window=1 this loop is not executed ( range(0,0,-1) = [] )
    for i in xrange(step,0,-1):
        
        rcount = 0
        mcount = 0

        if r.start-i<0:
            continue
        
        # tabulate reads/muts for the window extending into the read
        for j in xrange(window-i):
            rcount += r.read[j]-48 # '0' = unicode 48; '1' = 49
            mcount += r.muts[j]-48
        
        # add positions
        if rcount and mcount:
            mutnts[mutindex] = r.start-i
            mutindex+=1
            readnts[readindex] = r.start-i
            readindex+=1
        elif rcount:
            readnts[readindex] = r.start-i
            readindex+=1


    # Start main loop
    # This is executed for all window sizes
    # Reset rcount and mcount and start accumulating the sum for the first position
    rcount = 0
    mcount = 0
    for i in xrange(step):
        rcount += r.read[i]-48  # '0' = unicode 48; '1' = 49
        mcount += r.muts[i]-48
        validpos += r.read[i]-48
    

    for i in xrange( r.stop - r.start + 1 - step):
        
        if i>0:
            rcount -= r.read[i-1]-48
            mcount -= r.muts[i-1]-48

        rcount += r.read[i+step]-48
        mcount += r.muts[i+step]-48
        validpos += r.read[i+step]-48
        
        # this will count any window with at least one mutation
        if rcount and mcount:
            mutnts[mutindex] = r.start+i
            readnts[readindex] = r.start+i
            mutindex += 1
            readindex += 1
        
        # original: count only positions with window number of reads
        # CHANGE 5/2/18: count windows with at least one 'read' event
        elif rcount:
            readnts[readindex] = r.start+i
            readindex += 1



    # handle the end of the read: muts/reads contribute the downstream window
    # this loop is not executed for window = 1 ( range(1,1) = [] )
    # i is carried over from end of above for loop
    for j in xrange(1, window):
        
        if r.start+i+j > endseq:
            break

        rcount -= r.read[i+j-1]-48
        mcount -= r.muts[i+j-1]-48

        if rcount and mcount:
            mutnts[mutindex] = r.start+i+j
            readnts[readindex] = r.start+i+j
            mutindex += 1
            readindex += 1
        elif rcount:
            readnts[readindex] = r.start+i+j
            readindex += 1



    # if too few read positions, ignore this read by setting index to 0
    if validpos < mincoverage:
        mutindex=0
        readindex=0


    mutnts[0] = mutindex  # mutindex is +1 of last filled array position
    readnts[0] = readindex
    
    return 1


    

###########################################################################
###########################################################################
#
# Functions for reading old ShapeMapper mutation strings
#
###########################################################################
###########################################################################

def fillMatrices_Old(str inputFile, int[:,::1] read_arr, int[:, ::1] comut_arr, int[:, ::1] inotj_arr,
                     int window, int phred_cutoff, char* accepted_events, int mutseparation, int maxdel):
    
    """ inputFile = path to oldsytle mutation string file
    read_arr 
    comut_arr 
    inotj_arr       =  NxN numpy matrices to fill
    window          =  window >=1 over which to compute correlations
    phred_cutoff    =  phred value required to count as valid mutation
    accepted_events =  Mutation characters to count as valid mutations (e.g. 'AGTC-')
    mutseparation   =  Minimum seperation required between mut events
    maxdel          =  Maximum deletion size allowed
    """ 
    


    cdef int maxindex = read_arr.shape[0]

    # working arrays storing read events
    # first position contains element counter
    cdef int[:] readnts = np.zeros(maxindex+1, dtype=np.int32)
    cdef int[:] mutnts  = np.zeros(maxindex+1, dtype=np.int32)
    
    # mem alloc for reading the file using c
    cdef FILE* cfile
    cdef size_t lsize = 0
    cdef ssize_t endfile
    cdef char* line = NULL
    
    
    cdef READ r
    cdef int i, j, i_index

    
    # open the file
    cfile = fopen(inputFile, "r")
    if cfile == NULL:
        raise IOError(2, "No such file or directory: '{0}'".format(inputFile))

    
    cdef int linenum = -1
    cdef int readreads = 0

    # iterate through lines
    while True:
        
        linenum += 1 
        # get the the line of text using c-function
        endfile = getline(&line, &lsize, cfile)
        if endfile == -1:
            break
        

        # parse line into individual values, reset and fill readnts/mutnts
        try:
            r = parseLine_Old(line)
            if r.read == NULL:
                raise IndexError()
            elif r.stop > maxindex:
                print "Line {0} outside array bounds :: read bounds = ({1}, {2})".format(linenum, r.start, r.stop)
                continue
            
            fillReadMut_Old(readnts, mutnts, r, window, phred_cutoff, accepted_events, mutseparation, maxdel)
 
        except:
            print "Skipping incorrectly formatted line {0}".format(linenum)
            continue
        
        
        # check if read was read, and if so, increment counter
        if readnts[0]>0:
            readreads += 1


        # fill the read count matrix. Want to fill upper-right of matrix,
        # so traverse through arrays in reverse
        for i in xrange(readnts[0], 0, -1):
        
            i_index = readnts[i]
            read_arr[ i_index, i_index ] += 1
            
            for j in xrange(i-1, 0, -1):
                read_arr[ i_index, readnts[j] ] += 1
               
        
        
        # fill in the mut matrices
        for i in xrange(mutnts[0], 0, -1):
        
            i_index = mutnts[i]
            comut_arr[ i_index, i_index ] += 1
            
            for j in xrange(i-1, 0, -1):
                comut_arr[i_index, mutnts[j] ] += 1


            # fill in inotj
            # Note this loop overcounts for j=mutated
            # Diagnol is not used, so don't worry about i=j case
            for j in xrange(readnts[0], 0, -1): 
                inotj_arr[ i_index, readnts[j] ] += 1

            # correct for the over addition in inotj in above loop
            for j in xrange(mutnts[0], 0, -1):
                inotj_arr[ i_index, mutnts[j] ] -= 1
        


    fclose(cfile)

    return linenum, readreads



###########################################################################



cdef int fillReadMut_Old(int[:] readnts, int[:] mutnts, READ r, int window, int phred_cut, char* accepted_events, int mutdist_cut, int maxdelsize) except -1:
    # Accessory function for fillMatrices_Old
    # parse the old mutation string read for mut events
    # readnts and mutnts contain seq indices of valid reads and muts, respectively
    # Element 0 of readnts and mutnts is used as a counter...
    # Mutations within mutdist_cut of eachother are ignored
    # Reads with deletions larger than maxdelsize are ignored

    cdef char mutcode
    cdef int i,j, seqpos, qscore, counter, end

    cdef int readindex = 1 
    cdef int mutindex = 1 
    cdef int mutdist = mutdist_cut
    
    cdef int readlen = len(r.read) - 1
    cdef int lastmut = readlen+1

    cdef int delsize = 0
    cdef int boolval = 0

    cdef int[:] rmut = np.zeros(readlen+1, dtype=np.int32)
    cdef int[:] rread = np.zeros(readlen+1, dtype=np.int32)
    

    # reset readnts/mutnts arrays
    readnts[0]=0
    mutnts[0]=0
    
    # first quality filter read 3' -> 5'
    for i, mutcode in enumerate(r.read[::-1]):
        
        seqpos = readlen-i
        qscore = r.qual[seqpos] - 33 
        
        # Prob of match actually being mismatch is very low at any qscore,
        # thus use permissive cutoff
        if mutcode == '|' and qscore >= 10:
            mutdist += 1
            delsize = 0
            rread[seqpos] = 1
        
        elif mutcode == '~':
            mutdist = 0
            delsize += 1
            lastmut = seqpos

            # ignore reads that have too long mutations
            if delsize > maxdelsize:
                boolval = 1
                break

        elif mutcode in accepted_events:
                        
            if mutdist >= mutdist_cut and qscore >= phred_cut:
                rmut[seqpos] = 1
                rread[seqpos] = 1
            elif mutdist < mutdist_cut:
                rread[seqpos:lastmut] = 0
            
            lastmut = seqpos
            mutdist = 0 
    
    # exited loop above early -- read rejected, so return (readnts/mutnts were reset at top)
    if boolval:
        return 1    
    
    end = -1
    if r.start < 1:
        end = -r.start -1 # -1 here corrects for 1-indexing correction in parseLine
    
    for i in xrange(readlen-window+1, end, -1):
        # if there is mutation in the window, assume its good
        # regardless if some nucs are missing data...
        
        boolval = 0
        counter = 0
        for j in xrange(i, i+window):
            if rmut[j]:
                boolval = 1
                break
            if rread[j]:
                counter+=1

        if boolval:
            mutnts[mutindex] = r.start+i
            readnts[readindex] = r.start+i
            mutindex += 1
            readindex += 1

        elif counter == window:
            readnts[readindex] = r.start+i
            readindex += 1
            

    # adjust index so that it points to last element, rather than forward
    # note this for looping in read_data, since we are looping in reverse
    mutnts[0] = mutindex-1
    readnts[0] = readindex-1

    return 1


###########################################################################


cdef READ parseLine_Old(char* line):
    # accessory function for fillMatrices_Old
    # uses c functions to read through line; avoids having to use python split
    # which outputs python array and slows things down

    cdef READ r

    cdef char* token = strtok(line, " \t");
    token = strtok(NULL, " \t") #pop off the read name
    
    cdef int tokindex = 1

    while token:
        
        if tokindex == 1:
            r.start = atoi(token) - 1 # correct for 1-indexing 

        elif tokindex == 2:
            r.stop = atoi(token) - 1

        elif tokindex == 3:
            r.read = token

        elif tokindex == 5:
            r.qual = token
            break
        
        token = strtok(NULL, " \t")
        tokindex+=1
    
    # perform quaility checks...
    if tokindex != 5 or r.start > r.stop:  
        r.read = NULL

    return r











# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
