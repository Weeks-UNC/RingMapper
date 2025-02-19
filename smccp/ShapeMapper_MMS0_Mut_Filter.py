######################################################################
#  ShapeMapper MMS 0 Parsed Mutation File Filter
#  Version 1.0
#  (c) 2024 David Mitchell III
######################################################################

import itertools
import os
import argparse

######################################################################

## ADVANCED FILTER 10 ##

def scanLines(inl7, inl8):
    line_len = len(inl8)
    winsize = 10
    step = 1

    ## Special case for first 50 nt in a read ##
    testlist7 = list(map(int, inl7))            # Read string
    testlist8 = list(map(int, inl8))            # Mutation string
    testin8 = testlist8[0:50]
    testn = testin8.count(1)
    if testn / winsize >= 0.2:
        testlist7[0:50] = [0]*50
        testlist8[0:50] = [0]*50


    ## Case for the remainder of the read ##
    for w in range(50, line_len, step):
        Bstate = 0
        testin7 = testlist7[w:w+winsize]             # Read string
        testin8 = testlist8[w:w+winsize]             # Mutation string

        # Determine density of mutations in window #
        testn = testin8.count(1)
        testdens = testn / winsize
        if testdens >= 0.3 and testdens < 0.4:
            Bstate = 1
        elif testdens >= 0.4:
            Bstate = 2

        # Determine number of consecutive mutations in window #
        testnum = 0
        if Bstate >= 1:
            testnum = max(len(list(v)) for g,v in itertools.groupby(testin8, lambda x: x == 1) if g)

        # Identifiy windows with excessive mutations #
        #   Excessive mutations defined here as >= 3 consecutive mutations (a '1' in the window vs '0' for unmodified)
        #   Identifiy position of right-most mutation within the consecutive string and keep as '1'
        #   Change all mutations to the left of this mutation to '0' up to the end of the consecutive string
        if Bstate == 1 and testnum >= 3:
            a = 0
            b = a + testnum
            while a <= winsize-testnum:
                if 0 not in testin8[a:b]:
                    switchlen = testnum-1
                    # Change "1" to "0" in mut string (testlist8)
                    switchin8 = [0]*switchlen
                    testin8[a:b-1] = switchin8
                    # Change "1" to "0" in read string (testlist7)
                    switchin7 = [0]*switchlen
                    testin7[a:b-1] = switchin7
                    mvr = testnum
                else:
                    mvr = 1
                a+=mvr
                b+=mvr
        elif Bstate == 2:
            testin7 = [0]*winsize
            testin8 = [0]*winsize
        
        testlist7[w:winsize+w] = testin7
        testlist8[w:winsize+w] = testin8

    # Create output #
    outl7 = ''.join(map(str, testlist7))
    outl8 = ''.join(map(str, testlist8))

    return outl7, outl8

######################################################################

def filterFiles(datafile, inname, outDir):

    try:
        os.mkdir(outDir)
        print('Creating directory for filtered .mut files at ' + outDir + '.')
    except:
        print('Directory ' + outDir + ' already exists. Adding filtered .mut files here.')

    outFileName = outDir + '/MMS0Filter_' + inname
    outfile = open(outFileName, 'w')

    print('Input file name:', inname, '\nOutput file name:', outFileName)

    with open(datafile) as inp:

        line = inp.readlines()
        numlines = len(line)

        for b in range(numlines):
            spl = line[b].split()
            num = len(spl)

            if num <= 9 and spl[4] != 'OFF_TARGET':
                outfile.write(line[b])

            if num > 9 and spl[4] != 'OFF_TARGET':
                in_line7 = spl[7] 
                in_line8 = spl[8]
                out_line7, out_line8 = scanLines(in_line7, in_line8)
                for a in range(7):
                    outfile.write(spl[a] + '\t')
                outfile.write(out_line7 + '\t' + out_line8 + '\t')
                for c in range(9,num,1):
                    outfile.write(spl[c] + ' ')
                outfile.write('\n')

            if b % 50000 == 0:
                print('Working on line...', b+1)

    outfile.close()

    print('\n/// Completed MMS 0 Parsed Mutation File Filtering ///')


######################################################################


if __name__=='__main__':
    
    parser=argparse.ArgumentParser(description = "Program to filter ShapeMapper MMS 0 parsed mutation files prior to running PairMapper")
    parser.add_argument('datafile', help="Parsed mutation file (*.mut) generated from ShapeMapper with data generated using minimum mutation separation 0 (--min-mutation-separation 0)")
    parser.add_argument('--output_name', type=str, help="Output name for new filtered parsed mutation files created by this program, prepended by 'MMS0Filter_'. If unspecified, it uses the original input.mut name")
    parser.add_argument('--output_dir', type=str, default='MMS0_Filter_Output', help="Output directory for new filtered parsed mutation files. If unspecified, the default directory is 'MMS0_Filter_Output'")
    
    args=parser.parse_args()

    n = args.datafile
    i = n.split('/')
    inm = i[-1]

    # Check if the output name is specified #
    if args.output_name:
        inm = args.output_name
        
        # Check if the specified output name has a .mut extension #
        tmp = inm.split('.')
        ext = tmp[-1]
        if ext != 'mut':
            inm = inm + '.mut'
            

    ## Perform the Filtration ##
    filterFiles(args.datafile, inm, args.output_dir)
