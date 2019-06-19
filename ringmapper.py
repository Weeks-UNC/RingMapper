#!/usr/bin/env python

###############################################################################
#
#  Primary code to perform RING-MaP analysis on ShapeMapper output files
#   
#  See README for further details
#  Run with -h flag for arguments
#
#  Lead developer:  Anthony Mustoe
#  Contributors: Nicole Lama, Steve Busan
#
#  This file is licensed under the terms of the MIT license  
#
#  Version 4.2
#  December 2018
#
###############################################################################



import sys, argparse, itertools, math
import numpy as np

import readMutStrings  # cython code containing I/O funcs



class RINGexperiment(object):
    """RINGexperiment objects contain matrices and methods for computing correlated mutations
    from RINGexperiments parsed by ShapeMapper
    """

    def __init__(self, fasta = None, exfile=None, bgfile=None, arraysize=1000, 
                 corrtype = 'g', verbal=False, **kwargs):
        """
        fasta = fasta file of the sequence being analyzed
        exfile = datafile containing experiment data
        bgfile = datafile containing bg data
        arraysize = optional size of the arrays, if fasta not provided
        corrtype = type of correlation 
        kwargs are passed to initDataMatrices
        """

        if fasta is not None:
            self.sequence = self.readFasta(fasta, verbal=verbal)
            # set arraysize to sequence length, plus a little padding to guard against indexErrors
            self.arraysize = len(self.sequence)
        else:
            self.sequence = None
            self.arraysize = arraysize
        

        # initialize matrix values
        # ex matrices hold raw experiment data
        self.ex_readarr = None
        self.ex_comutarr = None
        self.ex_inotjarr = None
        
        # bg matrices hold raw bg data
        self.bg_readarr = None
        self.bg_comutarr = None
        self.bg_inotjarr = None
        
        # correlation matrices hold correlations
        self.ex_correlations = None
        self.bg_correlations = None 
        self.ex_zscores = None

        self.window = None
        
        self.setCorrType(corrtype, verbal=verbal)

        if exfile:
            self.initDataMatrices('ex', exfile, verbal=verbal, **kwargs)
        if bgfile:
            self.initDataMatrices('bg', bgfile, verbal=verbal, **kwargs)
        


    def setCorrType(self, corrtype, verbal=False):
        """Set correlation type
        Valid options are 'phi' and 'mi'
        """
        
        # convert to lower case
        corrtype = corrtype.lower()

        if corrtype == 'chi':
            self.correlationfunc = self._phiyates
            #self.significancefunc = self._phiyates
            if verbal: print("Using Yates-corrected Chi2 based correlation metric")

        elif corrtype == 'g':
            self.correlationfunc = self._mistatistic
            #self.significancefunc = self._mistatistic
            if verbal: print("Using G-test correlation metric")
        
        elif corrtype == 'apc':
            self.correlationfunc = self._mistatistic
            #self.significancefunc = self._mistatistic
            if verbal: print("Using APC corrected G-test correlation metric")
        
        elif corrtype == 'mi':
            self.correlationfunc = self._mutualinformation
            if verbal: print("Using normalized MI correlation metric")


        elif corrtype == 'nmi':
            self.correlationfunc = self._norm_mutualinformation
            if verbal: print("Using normalized MI correlation metric")

        else:
            raise ValueError("Unrecognized correlation metric : {0}".format(corrtype))

        self.corrtype = corrtype


    def readFasta(self, fasta, verbal=False):
        """Read the sequence in the provided sequence file"""
        
        with open(fasta) as inp:
            inp.readline()
            
            seq = ''
            for line in inp:
                if line[0] == '>':
                    break
                seq += line.strip()

        
        if verbal:
            print("Sequence length={0} read from {1}".format(len(seq), fasta))

        return seq    



    def initDataMatrices(self, prefix, datafile, window=1, verbal=False, **kwargs):
        """initialize and fill the read, comut, inotj matrices
        prefix    = ex/bg, indicating the type of the data
        datafile  = new/old mutation string file
        **kwargs are passed onto appropriate fillMatrices function
        """
        
        if self.window is None:
            self.window = window
            if verbal: print("Computing correlations using window={0}".format(self.window))

        elif window != self.window:
            raise ValueError("Window already set to {0}; passed window={1}".format(self.window, window))

    
        # initialize the matrices
        read = np.zeros( (self.arraysize, self.arraysize), dtype=np.int32)
        comut = np.zeros( (self.arraysize, self.arraysize), dtype=np.int32)
        inotj = np.zeros( (self.arraysize, self.arraysize), dtype=np.int32)
        
        # determine whether new or old mutstring format
        filetype = self._filetype(datafile)
        
        if filetype > 0:
            if verbal: print("Filling {0} arrays from {1}".format(prefix, datafile))
            self._fillMatrices(datafile, read, comut, inotj, filetype, verbal=verbal, **kwargs)
        else:
            if verbal: print("Filling {0} arrays from OLD format file {1}".format(prefix, datafile))
            self._fillMatrices_Old(datafile, read, comut, inotj, verbal=verbal, **kwargs)


        # assign the matrices
        setattr(self, prefix+'_readarr', read)
        setattr(self, prefix+'_comutarr', comut)
        setattr(self, prefix+'_inotjarr', inotj)
 



    def _filetype(self, datafile):  
        """Determine format of datafile:
             return 0 if old format
             return 1 if ShapeMapper 2.1.1 format
             return 2 if ShapeMapper 2.1.2 format
             return 3 if ShapeMapper 2.1.4-rc or higher format
        """
        
        try:
            fileformat = 999
            with open(datafile) as inp:
                
                line = inp.readline()
                spl = line.split()

                if '|' in spl[3]:
                    fileformat = 0
                elif spl[0][:4] in ('MERG', 'PAIR','UNPA', 'UNSP'):
                    if not spl[4].isdigit():
                        fileformat = 3
                    else:
                        fileformat = 2
                else:
                    fileformat = 1
                
            return fileformat


        except:
            raise IOError('{0} has unrecognized format'.format(datafile))



    def _fillMatrices(self, datafile, read, comut, inotj, fileformat, mincoverage=0, undersample=-1, verbal=False, **kwargs):
        """Call the cython fillMatrices function for new classified mutation file format
        datafile    = New classified mutations file to read
        read
        comut
        inotj       = NxN matrices to fill
        fileformat  = parsed mutation file code from _filetype
        mincoverage = Minimum number of valid 'read' positions required per read
        undersample = Maximum number of reads to read; default of -1 means read all reads
        """
        
        if 0<mincoverage<1:
            validpos = sum([x.isupper() for x in self.sequence])
            mincoverage *= validpos
        
        if verbal and mincoverage>0:
            print("Read length filtering ON\n\tMatch threshold = {0}".format(mincoverage))

        fillstats = readMutStrings.fillMatrices(datafile, read, comut, inotj, self.window, mincoverage, fileformat, undersample)
    
        if verbal:
            print("Input summary:")
            print("\t{0} reads in {1}".format(fillstats[0], datafile))
            print("\t{0} reads passed filtering".format(fillstats[1]))
            


    def _fillMatrices_Old(self, datafile, read, comut, inotj, phred_cut=30,
                          accepted_events = 'AGCT-', mutseparation=5, maxdel=1000, verbal=False, **kwargs):
        """Call the cython fillMatrices_Old function
        datafile        = Old mutation string file to read
        read
        comut
        inotj           = NxN matrices to fill
        phred_cut       = Minimum phred value required for valid mutations
        accepted_events = Accepted valid mutations events
        mutseparation   = Separation distance required between valid mutations
        maxdel          = maximum deletion/no-data region allowed for valid reads
        """
        
        if verbal:
            print("Post-processing old ShapeMapper called mutations:")
            print("\tPhred cutoff = {0}".format(phred_cut))
            print("\tMut. event separation = {0}".format(mutseparation))
            print("\tMaximum deletion cutoff = {0}".format(maxdel))
            print("\tAccepted mut. events = {0}".format(accepted_events))


        fillstats = readMutStrings.fillMatrices_Old(datafile, read, comut, inotj, self.window,
                                                    phred_cut, accepted_events, mutseparation, maxdel)


        if verbal:
            print("Input summary:")
            print("\t{0} reads in {1}".format(fillstats[0], datafile))
            print("\t{0} reads passed filtering".format(fillstats[1]))
           
    


    def getMaxArrayIndex(self, prefix='ex'):
        """Return index of the last non-zero diagonal element. Equal to sequence length if set.
        Otherweise, determine length of molecule after read matrices are filled"""
        
        try:
            return self.maxarrayindex
        
        except AttributeError:
            
            if self.sequence is not None:
                self.maxarrayindex = len(self.sequence)

            else:
                arr = getattr(self, prefix+'_readarr')
            
                last = 0
                for i in xrange(arr.shape[0]):
                    if arr[i,i] != 0:
                        last = i
            
                self.maxarrayindex = last+1
            
            return self.maxarrayindex



    def getReactiveNts(self, ratecut, prefix='ex'):
        """Return indices of nts with mutation rates above ratecut"""
        
        readarr = getattr(self, prefix+'_readarr')
        comutarr = getattr(self, prefix+'_comutarr')

        ntlist = []
        for i in xrange( self.getMaxArrayIndex() ):
        
            if readarr[i,i] == 0:
                continue

            mutrate = float(comutarr[i,i])/readarr[i,i]

            if mutrate > ratecut:
                ntlist.append(i)

        return ntlist

    

    def getUnreactiveNts(self, ratecut, prefix='ex'):
        """Return indices of nts with mutation rates below ratecut"""

        readarr = getattr(self, prefix+'_readarr')
        comutarr = getattr(self, prefix+'_comutarr')

        ntlist = []
        for i in xrange( self.getMaxArrayIndex() ):
        
            if readarr[i,i] == 0:
                continue

            mutrate = float(comutarr[i,i])/readarr[i,i]

            if mutrate < ratecut:
                ntlist.append(i)

        return ntlist



    def _phistatistic(self, phi, n):
        """convert phi coefficient to chi2 statistic
        phi = phi coeff
        n   = total number observations
        """
        if np.isnan(phi):
            return 0

        return n*phi**2

    
    def _phiyates(self, n,b,c,d):
        """Compute yates chi2 from contigency table values"""

        af = float(n-b-c-d)
        bf = float(b)

        # multiply floats, which should avoid overflow errs
        bot = (af+bf)*(c+d)*(af+c)*(bf+d)
        if bot < 1:
            return 0
        
        # multiply floats which should avoid overflow errs
        top = n*(abs(af*d - bf*c) - 0.5*n)**2
        return top/bot 


    def _phi(self, n,b,c,d):
        """ Return Phi
        n = a+b+c+d
        - a, b, c, d correspond to values in the 2 x 2
        contingency table tested for nucs i and j:
         
               i
            0      1
          -----------
         0|  a     b
          |
         1|  c     d
        """
        
        # convert to float for multiplication to avoid overflow
        af = float(n-b-c-d)
        bf = float(b)
    
        #return (af*d - bf*c) / min( (af+bf)*(bf+d), (af+c)*(c+d) )
        bot = (af+bf)*(c+d)*(af+c)*(bf+d)
        if bot < 1:
            return 0

        return (af*d - bf*c)/np.sqrt(bot)
        

    def _mutualinformation(self, n, b,c,d):
        """Compute Mutual Information for a given nt pair
        n = a+b+c+d
        - a, b, c, d correspond to values in the 2 x 2
        contingency table tested for nucs i and j:
         
               i
            0      1
          -----------
         0|  a     b
          |
         1|  c     d
        """
        
        bf = float(b)
        df = float(d)
        a = n-bf-c-df

        if min(a,b,c,d) <  1:
            return 0 

        mi = a*np.log(a) + bf*np.log(bf) + c*np.log(c) + df*np.log(df)
        mi += n*np.log(n)
        mi -= (a+c)*np.log(a+c) + (a+bf)*np.log(a+bf) + (bf+df)*np.log(bf+df) + (c+df)*np.log(c+df)
        mi /= n
        
        return mi
    
 

    def _mistatistic(self, n, b,c,d):
        """convert mutual information value to g statistic
        n  = total number observations
        """
        return 2*n*self._mutualinformation(n,b,c,d)

    
    def _norm_mutualinformation(self, n, b, c, d):

        mi = self._mutualinformation(n,b,c,d)
    
        bf = float(b)
        df = float(d)
        cf = float(c)
        af = n-bf-c-df

        hx = -1*( (af+bf)*np.log(af+bf) + (cf+df)*np.log(cf+df) - n*np.log(n) ) / n 
        hy = -1*( (af+cf)*np.log(af+cf) + (bf+df)*np.log(bf+df) - n*np.log(n) ) / n

        return mi / np.sqrt(hx*hy)



    def correlationsign(self,i,j, prefix='ex'):
        
        arr = getattr(self, prefix+'_readarr')
        n = float(arr[i,j])
        arr = getattr(self, prefix+'_inotjarr')
        b = float(arr[i,j])
        c = float(arr[j,i])
        arr = getattr(self, prefix+'_comutarr')
        d = float(arr[i,j])
        
        if n==0 or b+d==0 or c+d==0:
            return 0
        if (n*d)/((b+d)*(c+d))<1:
            return -1
        else:
            return 1
        

    def apcCorrection(self, prefix, mindefined=10):
        """ Perform APC correction to correlation matrix
        prefix = ex/bg
        """
        
        mimatrix = getattr(self, prefix+'_correlations')
        
        # Note that mimatrix is symmetric masked array.
        # mean function computes mean only over valid entries
        # Diagonal + buffer entries are masked out (invalid), as well as other 
        # entries masked because of low counts, etc.
        

        # determine nts that have low valid counts and mask out
        counts = mimatrix.count(axis=0)
        for i,c in enumerate(counts):
            if 0<c<mindefined:  # if c==0 everything is already masked, so need to redo
                mimatrix[i,:] = np.ma.masked
                mimatrix[:,i] = np.ma.masked
                

        # compute mean over all MI entries
        xyBar = mimatrix.mean() 
        
        if xyBar == 0:
            return 0  

        #compute MI(i,Xbar), mean along each column
        xBar = mimatrix.mean(axis=0) #(nx1) array
        
        yBar = xBar.reshape((-1, 1)) #(1xn) array
        
        apc = xBar*yBar/xyBar # (n x n) array

        miP = mimatrix - apc
        
        setattr(self, prefix+'_correlations', miP)
    


    def _correlationMatrix(self, prefix, corrbuffer, mindepth, mincount):
        """ Calculate the correlation for every position
        prefix        = data matrix to use (ex, bg)
        corrbuffer    = distance between correlations
        mindepth      = minimum required read depth
        mincount     = minimum required count 

        returns masked array
        """
               
        read = getattr(self, prefix+'_readarr')
        inotj = getattr(self, prefix+'_inotjarr')
        comut = getattr(self, prefix+'_comutarr')
        
        seqlen = self.getMaxArrayIndex()
        
        # initialize the matrix
        cmat = np.empty((seqlen, seqlen), dtype=np.float32)
        cmat[:] = np.nan
        
        for i in xrange(seqlen):
            for j in xrange(i+corrbuffer, seqlen):               
                if read[i,j]>=mindepth and min(inotj[i,j], inotj[j,i], comut[i,j])>=mincount:
                    cmat[i,j] = self.correlationfunc(read[i,j], inotj[i,j], inotj[j,i], comut[i,j])
                    cmat[j,i] = cmat[i,j] 
        

        # set correlations matrix to masked array with non-filled (nan) values masked out
        setattr(self, prefix+'_correlations', np.ma.masked_invalid(cmat))

    


    def computeCorrelationMatrix(self, corrbuffer=5, mindepth=10000, 
                                 mincount=50, ignorents = [], highbgrate=0.02, 
                                 highbgcorr=10.83, verbal=True):
        """Compute the correlation matrices and mask invalid entries
        corrbuffer     = buffer to keep between correlations (i.e. minimum correlation distance)
        mindepth       = minimum pairwise read depth
        mincount       = minimum number of counts allowed in contigency table entry
        ignorents      = array of nts to ignore in correlation calcs
        highbgrate     = maximum bg rate (will be scaled by window size)
        highbgcorr     = Chi2 cutoff indicating significant bg correlation
        """
        
                
        # adjust buffer to account for window size.
        self.corrbuffer = corrbuffer+self.window
        
        # scale bgrate by window
        self.highbgrate = highbgrate*self.window

        # compute correlation matrix
        self._correlationMatrix('ex', self.corrbuffer, mindepth, mincount)
        
            
        # mask out user specified values
        for i in ignorents:
            self.ex_correlations[i,:] = np.ma.masked
            self.ex_correlations[:,i] = np.ma.masked
            if verbal: print ("Nt {0} ignored: specified by user".format(i+1))
        
        # set containing masked out nts
        allinvalid = set(ignorents)
        
        
        # Mask out nt columns that have high bg mutation rates
        if self.bg_readarr is not None:

            highbgnts = self.getReactiveNts(self.highbgrate, prefix='bg')

            for i in highbgnts:
                self.ex_correlations[i,:] = np.ma.masked
                self.ex_correlations[:,i] = np.ma.masked
                
                if verbal and i not in allinvalid:
                    e = float(self.bg_comutarr[i,i])/self.bg_readarr[i,i]
                    print( "Nt {0} ignored: bg_rate={1:.3f}".format(i+1,e) )
            
            # add highbgnts to invalid nts
            allinvalid.update(highbgnts)
       


        # perform apc correction
        if self.corrtype == 'apc':
            self.apcCorrection('ex')
        

        # compute z-scores of ex matrix
        self.computeZscores()


        # now cross-reference and remove bg-correlated pairs
        self.maskBGcorrelated(highbgcorr=highbgcorr, invalid=allinvalid, verbal=verbal)





    def maskBGcorrelated(self, highbgcorr=10.83, invalid=[], verbal=False):
        """mask positions in ex_correlations and ex_zscores that are correlated
        in the bg sample and which have higher mi"""


        if self.bg_readarr is None:
            return

        # compute bg correlations and get pairs that are significantly correlated
        self._correlationMatrix('bg', self.corrbuffer, 10000, 20)
        bgcorrs = self.significantCorrelations('bg', highbgcorr)
         

        # search through significant correlations
        for i,j in bgcorrs:

            # compute mi of ex and bg samples
            exmi = self._mutualinformation(self.ex_readarr[i,j], self.ex_inotjarr[i,j],
                                           self.ex_inotjarr[j,i], self.ex_comutarr[i,j])
            bgmi = self._mutualinformation(self.bg_readarr[i,j], self.bg_inotjarr[i,j],
                                           self.bg_inotjarr[j,i], self.bg_comutarr[i,j])
            
            if 5*bgmi > exmi:
                
                # look to see if i,j was significant and if so print
                excorr = self.ex_correlations[i,j]
                if verbal and i not in invalid and j not in invalid and excorr>=20:
                    outstr = 'Correlated pair ({0},{1}) w/ chi2={2:.1f} ignored'.format(i+1, j+1, excorr)
                    outstr += ': correlated in bg w/ chi2={0:.1f}'.format(self.bg_correlations[i,j])
                    print(outstr)        
                
                # mask out values
                self.ex_correlations[i,j] = np.ma.masked
                self.ex_correlations[j,i] = np.ma.masked
                self.ex_zscores[i,j] = np.ma.masked
                self.ex_zscores[j,i] = np.ma.masked
                



    def significantCorrelations(self, prefix, chi2cut, sign=-1):
        """ Calculate the correlation for every position
        prefix    = data matrix to use (ex/bg)
        chi2cut   = significance cutoff
        sign      = option to filter correlation by requiring them 
                    to be positively correlated (set to 0/1 to enable)

        returns list of (i,j) pairs
        """

        seqlen = self.getMaxArrayIndex()
       
        corrmat = getattr(self, prefix+'_correlations')
        
        corrs = []
        for i in xrange(seqlen):
            for j in xrange(i, seqlen):
                # this will automatically skip masked out values
                if corrmat[i,j] >= chi2cut and self.correlationsign(i,j,prefix) >= sign:
                    corrs.append((i,j))

        
        return corrs



    def computeZscores(self):
        
        if self.ex_correlations is None:
            raise AttributeError("ex_correlations has not been initialized!")
        

        # initialize the matrix
        zscores = np.empty( self.ex_correlations.shape )
        zscores[:] = np.nan

        corrmat = self.ex_correlations

        # compute means and std for z-score calculation
        corrmean = corrmat.mean(axis=0)
        corrstd = corrmat.std(axis=0)
        
        seqlen = self.getMaxArrayIndex()
        
        for i in xrange(seqlen):
            for j in xrange(i+1, seqlen):
                if not corrmat.mask[i,j]:
                    zscores[i,j] = (corrmat[i,j]-corrmean[i])/corrstd[i]
                    zscores[j,i] = (corrmat[i,j]-corrmean[j])/corrstd[j]


        self.ex_zscores = np.ma.masked_invalid( zscores )

    

    def getMeanZ(self, i, j):
        """Return the mean zscore at i,j"""
        return (self.ex_zscores[i,j]+self.ex_zscores[j,i])/2



    def writeCorrelations(self, outfile, chi2cut=20, sign=-1):
        """Write out correlations in a human readable format that also conforms
        to pairing probability (dotplot) file format for subsequent plotting.
    
        -Outfile is the output file path
        -chi2cut is the significance cutoff
        -if sign 0/1 than only write out positive correlations
        """
 
        if self.ex_zscores is None:
            self.computeZscores()

        corrs = self.significantCorrelations('ex', chi2cut, sign=sign)
        

        with open(outfile,'w') as OUT:
        
            OUT.write("{0}\tWindow={1}\tMetric={2}\n".format(self.getMaxArrayIndex(), self.window, self.corrtype.upper()))
            
            OUT.write("i\tj\tStatistic\t+/-\tZij\tZi\tZj\tMod_Depth\tMod_Comuts\tUnt_Depth\tUnt_Comuts\n")

            for i,j in corrs:
                
                OUT.write("{0}\t{1}\t".format(i+1, j+1))
                
                OUT.write("{0:.2f}\t{1}\t".format(self.ex_correlations[i,j], self.correlationsign(i,j, 'ex')))
                
                OUT.write("{0:.2f}\t".format(self.getMeanZ(i,j)))

                OUT.write("{0:.2f}\t{1:.2f}\t".format(self.ex_zscores[i,j], self.ex_zscores[j,i]))
                OUT.write("{0}\t{1}\t".format(self.ex_readarr[i,j],self.ex_comutarr[i,j]))
                
                if self.bg_readarr is not None:
                    OUT.write("{0}\t{1}".format(self.bg_readarr[i,j],self.bg_comutarr[i,j]))
                
                OUT.write("\n")

     

    def writeDataMatrices(self, prefix, outprefix):
        """Save the data matrices to file
        prefix    = ex/bg matrix prefix
        outprefix = prefix of output file name
        """

        for a in ('readarr', 'comutarr', 'inotjarr'):
            name = '{0}_{1}'.format(prefix, a)
            matrix= getattr(self, name)
            np.savetxt('{0}_{1}.mat'.format(outprefix, name), matrix, delimiter=" ", fmt="%.6e")
 
        
    def readDataMatrices(self, prefix):
        
        for a in ('readarr', 'comutarr', 'inotjarr'):
            setattr(self, 'ex_'+a, np.loadtxt('{0}_ex_{1}.mat'.format(prefix, a)))

        try:
            for a in ('readarr', 'comutarr', 'inotjarr'):
                setattr(self, 'bg_'+a, np.loadtxt('{0}_bg_{1}.mat'.format(prefix, a)))
        except IOError:
            print('WARNING: no bg matrices found for {0}'.format(prefix))


        


    ###############################################################################



def parseArguments():

    parser = argparse.ArgumentParser(description = "Compute correlations from parsed mutations")
    parser.add_argument('inputFile', help='Path to mutation string file (can be new or old ShapeMapper format)')
    parser.add_argument('outputFile', help='Path for correlation output file')
    
    parser.add_argument('--fasta', help='Path to fasta sequence file')
    
    parser.add_argument('--untreated', help='Path to untreated (bg) mutation file. Used to remove high bg positions and bg correlations')
    
    parser.add_argument('--window', type=int, default=1, help="Nt window over which to compute correlations (default = 1)")
    parser.add_argument('--chisq_cut', type=float, default=20.0, help="Set chisq cutoff (default = 20)")
    parser.add_argument('--mindepth', type=int, default=10000, help='Minimum pairwise read depth allowed for calculating correlations (default = 10000)')
    parser.add_argument('--mincount', type=int, default=50, help="""Minimum required count in contigency table 
                        (default = 50). Nt pairs with fewer than this number of comutations are ignored""")
    
    parser.add_argument('--metric', type=str, default='apc', help="""Metric to use for computing correlations. 
                        options are chi/g/apc (Chi, G-test, or APC corrected G-test). (default = apc)""")
    
    parser.add_argument('--mincorrdistance', type=int, default=5, help="""Minimum distance allowed between correlations (default=5)""")

    parser.add_argument('--mincoverage', type=float, default=0, help="""Quality filter reads by requiring a minimum 
                        number of positional matches to reference sequence.
                        This parameter can be set as a fraction of the total molecule length 
                        (i.e. if 0.8 is passed, reads will be required to have at least 0.8*len(fasta) valid matches)/
                        Alternatively, an integer value can be passed 
                        (i.e. if 150 is passed, reads will be required to have at least 150 valid matches).
                        By default, this filter is disabled.""")

    parser.add_argument('--highbg_rate', type=float, default = 0.02, help="""Ignore nts with bg reactivity above
                        this value (default = 0.02). Value is multipled by window (so default=0.06 for window=3)""")
    
    # 15.13 or 10.83
    parser.add_argument('--highbg_corr', type=float, default = 10.83, help="""Ignore nt pairs correlated in the bg 
                        sample, with correlation determined via this significance value (default=10.83 --> P=1e-3)""")

    parser.add_argument('--ignorents', help="""A list of comma-separated (1-indexed) nts to ignore (e.g. 10,12,35)""")

    parser.add_argument('--molsize', type=int, default=1000, help="""Size of molecule arrays (default = 1000). 
                        Value must be larger than max read index. Only used if fasta file not provided.""") 
    
    parser.add_argument('--undersample', type=int, default=-1, help="""Randomly undersample specified number of reads 
                         from inputFile (default=-1 [disabled]).""") 

    parser.add_argument('--writematrixfile', help="Write mutation matrices to file (provide prefix)")
    
    

    args = parser.parse_args()   
        
    # check to make sure mincoverage is correct
    if 0<args.mincoverage<1:
        assert args.fasta is not None, "fasta file must be provided when using fraction mincoverage filter"


    # parse ignorents argument
    if args.ignorents:
        spl = args.ignorents.split(',')
        ig = []
        try:
            for x in spl:
                if ':' in x:
                    xspl = x.split(':')
                    ig.extend(range(int(xspl[0])-1, int(xspl[1])))
                else:
                    ig.append(int(x)-1)  # correct for 1-indexing of input 
        except ValueError:
            raise ValueError('Invalid ignorents option: {0}'.format(args.ignorents))

        args.ignorents = ig
    
    else:
        args.ignorents = []


    return args




###############################################################################



if __name__ == '__main__':

    args = parseArguments()

    verbal = True
    #verbal = False

    # initialize the object and read in matrices
    ringexp = RINGexperiment(fasta=args.fasta, 
                             exfile = args.inputFile,
                             bgfile = args.untreated,
                             arraysize = args.molsize,
                             window = args.window,
                             corrtype = args.metric, 
                             mincoverage = args.mincoverage,
                             undersample = args.undersample,   
                             verbal = verbal)
    

    # compute the correlation matrices
    ringexp.computeCorrelationMatrix(corrbuffer = args.mincorrdistance, 
                                     mindepth = args.mindepth,
                                     mincount = args.mincount,
                                     ignorents = args.ignorents,
                                     highbgrate = args.highbg_rate,
                                     highbgcorr = args.highbg_corr, 
                                     verbal = verbal)
    
    ringexp.writeCorrelations(args.outputFile, chi2cut=args.chisq_cut)
 
    
    
    if args.writematrixfile:
        ringexp.writeDataMatrices('ex', args.writematrixfile)
        
        if ringexp.bg_readarr is not None:
            ringexp.writeDataMatrices('bg', args.writematrixfile)
            






# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
