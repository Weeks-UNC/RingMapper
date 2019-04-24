#!/usr/bin/env python

#########################################################################
#
#  Code to detect base pairs from correlated chemical probing data
#   
#  See README for further details
#  Run with -h flag for arguments
#
#  Lead developer:  Anthony Mustoe
#  Contributors: Nicole Lama, Steve Busan
#
#  This file is licensed under the terms of the MIT license  
#
#  Version 1
#  December 2018
#
#########################################################################



import argparse, sys, itertools
import numpy as np


from ringmapper import RINGexperiment
from ReactivityProfile import ReactivityProfile
#from arcPlot import ArcPlot



class PairMapper(object):
    """ class to hold matrix of pairs """


    def __init__(self, ringexp, profile, 
                       chi2cut=20, 
                       primary_reactivity=0.2, primary_zscore=2.0,
                       secondary_reactivity=0.5, secondary_zscore=2.0,
                       maxGU=None, maxNC=0):
        """Initialize PairMapper from ringexp object"""
        
        self.parent = ringexp
        self.profile = profile

        # first check whether filters are passed
        self.passfilter = self.checkMutationRates(primary_reactivity = primary_reactivity)
        
        self.complementarycorrs = None
        self.primary = None
        self.secondary = None

        # compute list of complementary correlated pairs 
        self.complementarycorrs = self.computeComplementaryCorrs(chi2cut, maxGU=maxGU, maxNC=maxNC)
       
        self.computePrimaryCorrs(primary_reactivity, primary_zscore)
        self.computeSecondaryCorrs(secondary_reactivity, secondary_zscore)



    def computePrimaryCorrs(self, primary_reactivity, primary_zscore):
        """Compute primary pair signals"""

        corrs = self.getStrongest()

        # now need to filter by reactivity 
        corrs = self.filterReactive(corrs, primary_reactivity)

        # filter by zave
        corrs = self.filterZscore(corrs, primary_zscore)
    
        self.primary = corrs

    
    def computeSecondaryCorrs(self, secondary_reactivity, secondary_zscore):
        """Compute secondary pair signals"""
            
        # start from all complementary corrs
        corrs = self.complementarycorrs
        
        # filter by reactivity 
        corrs = self.filterReactive(corrs, secondary_reactivity)

        # filter by zave
        corrs = self.filterZscore(corrs, secondary_zscore)
    
        # eliminate corrs that are primary
        corrs = list( set(corrs)-set(self.primary) )
        corrs.sort()
        
        self.secondary = corrs



    def filterZscore(self, corrlist, cutoff):
        """Filter the set of correlations in corrlist by average zscore"""

        filtered = []
        for i,j in corrlist:
            if self.parent.getMeanZ(i,j) >= cutoff:
                filtered.append((i,j))

        return filtered



    def filterReactive(self, corrlist, r_cutoff):
        """Filter the set of correlations in corrlist by reactivity """

        filtered = []
        for c in corrlist:
            r0 = safeReactivityMean(self.profile.normprofile[c[0]:c[0]+self.parent.window])
            r1 = safeReactivityMean(self.profile.normprofile[c[1]:c[1]+self.parent.window])
            
            if r0<=r_cutoff and r1<=r_cutoff:
                filtered.append(c)
        
        return filtered



    def getStrongest(self):
        """Return non-conflicting strongest correlations"""

        # get list of strongest correlation at each nt
        corrdict = {}
        for i,j in self.complementarycorrs:
            if i not in corrdict: 
                corrdict[i] = (i,j, self.parent.ex_correlations[i,j])
            elif self.parent.ex_correlations[i,j] > corrdict[i][2]:
                corrdict[i] = (i,j, self.parent.ex_correlations[i,j])
            
            if j not in corrdict: 
                corrdict[j] = (i,j, self.parent.ex_correlations[i,j])
            elif self.parent.ex_correlations[i,j] > corrdict[j][2]:
                corrdict[j] = (i,j, self.parent.ex_correlations[i,j])
       

        # eliminate correlations that are not mutually strongest
        ntkeys = corrdict.keys()
        for k in ntkeys:
            # need to check since we are deleting things
            if k not in corrdict:
                continue

            i,j,c = corrdict[k]
            if k == i and (j not in corrdict or corrdict[j] != corrdict[k]):
                del corrdict[k]
            elif k == j and (i not in corrdict or corrdict[i] != corrdict[k]):
                del corrdict[k]


        # function to test if two correlations are parallel
        def parallel(c1, c2):
            if c1[0]-c2[0] == -1*(c1[1]-c2[1]):
                return True
            return False


        # build final list while filtering out conflicting windows
        output = []
        for k in corrdict:
            
            # only consider/add each correlation once
            if k != corrdict[k][0]:
                continue
            
            # scan up and downstream of each correlation leg to see if there is 
            # a stronger conflicting correlation
            keep = True
            for i in range(1-self.parent.window,0)+range(1,self.parent.window):
                
                if k+i in corrdict and not parallel(corrdict[k], corrdict[k+i]) \
                                   and corrdict[k+i][2] > corrdict[k][2]:
                    keep = False
                    break
                
                # repeat for other side of the correlation
                j = corrdict[k][1]
                if j+i in corrdict and not parallel(corrdict[k], corrdict[j+i]) \
                                   and corrdict[j+i][2] > corrdict[k][2]:
                    keep = False
                    break

            if keep:
                output.append(corrdict[k][:2])
        
        output.sort()

        return output





    def computeComplementaryCorrs(self, chi2cut, maxGU=None, maxNC=0):
        """Return list of complementary correlations"""

        def isComplement(seq1, seq2, maxGU, maxNC):
            counts = [0,0,0]

            # Not considered complementary if different lengths
            if len(seq1) != len(seq2):
                return False
            
            for i in range(len(seq1)):
                pair = seq1[i] + seq2[-(i+1)]

                if pair in ('AU', 'UA', 'TA','AT', 'GC', 'CG'):        
                    counts[0]+=1
                elif pair in ('GU', 'UG', 'GT', 'TG'):
                    counts[1]+=1
                else:
                    counts[2]+=1
            
            if counts[2] > maxNC or counts[1] > maxGU:
                return False
            
            return True
        
        
        # default is all GU pairs are valid
        if maxGU is None:
            maxGU = self.parent.window
    

        # get all positive correlations 
        corrs = self.parent.significantCorrelations('ex', chi2cut, sign=1)
        
        compcorrs = []
        
        for i,j in corrs:
            
            if isComplement(self.parent.sequence[i:i+self.parent.window], 
                    self.parent.sequence[j:j+self.parent.window], maxGU=maxGU, maxNC=maxNC):
                compcorrs.append((i,j))

        return compcorrs



    def writePairs(self, outname):
        

        def writeline(OUT, i, j, code):
            #function to write a single line of output
            OUT.write("{0}\t{1}\t".format(i+1, j+1))
            OUT.write("{0:.2f}\t".format(self.parent.ex_correlations[i,j]))
            OUT.write("{0}\t".format(code))
            OUT.write("{0:.2f}\t".format(self.parent.getMeanZ(i,j)))
            OUT.write("{0:.2f}\t{1:.2f}\t".format(self.parent.ex_zscores[i,j], 
                                                  self.parent.ex_zscores[j,i]))
            OUT.write("{0}\t{1}\t".format(self.parent.ex_readarr[i,j],
                                          self.parent.ex_comutarr[i,j]))
            OUT.write("{0}\t{1}\n".format(self.parent.bg_readarr[i,j],
                                        self.parent.bg_comutarr[i,j]))


        # determine complementary correlations that aren't primary or secondary
        remainder = set(self.complementarycorrs)-set(self.primary)-set(self.secondary)
        remainder = list(remainder)
        remainder.sort()


        with open(outname,'w') as outf:
            
            outf.write("{0}\tWindow={1}\tMetric=PAIRMAPPER\t".format(len(self.parent.sequence), self.parent.window)) 
            outf.write("\n".format())

            outf.write("i\tj\tSig.\tClass\tZij\tZi\tZj\tMod_Depth\tMod_Comuts\tUnt_Depth\tUnt_Comuts\n")

            for i,j in self.primary:
                writeline(outf, i, j, 1)

            for i,j in self.secondary:
                writeline(outf, i, j, 2)

            for i,j in remainder:
                writeline(outf, i, j, 0)




    def writePairBonusFile(self, filepath, chi2cut, maxNC=1, scale=0.5, intercept=0):
        """Write matrix of pairing bonuses for use in RNAstructure folding"""

        seqlen = len(self.parent.sequence)
        pairmat = np.zeros((seqlen, seqlen))
 
        corrs = self.computeComplementaryCorrs(chi2cut, maxNC=maxNC)
        
        for i,j in corrs:
            
            meanz = self.parent.getMeanZ(i,j)  
            if meanz > 1:
                energy = -0.616*scale*np.log(meanz) - intercept
                
                # add to each base pair
                for w in range(self.parent.window):
                    pairmat[i+w, j+self.parent.window-1-w] += energy


        # make symmetric
        pairmat += pairmat.transpose()

        np.savetxt(filepath, pairmat, fmt='%.4f')



    def checkMutationRates(self, depthcut=10000, primary_reactivity = 0.2):
        """Check comutation rates to make sure they are high enough to accurately measure pairs"""
        
        seqlen = self.parent.ex_readarr.shape[0]
        
        # initiate the masked array
        comuts = np.ma.masked_array(self.parent.ex_comutarr, dtype=np.float64)
        
        # mask out the lower diagonal (which is void of info)
        comuts[np.tril_indices(seqlen)] = np.ma.masked
        

        # mask out nts with no mutation signal in the modifed
        nts = self.parent.getUnreactiveNts(5e-4, prefix='ex')
        for i in nts:
            comuts[i,:] = np.ma.masked
            comuts[:,i] = np.ma.masked
        
        # mask out high background positions
        nts = self.parent.getReactiveNts(self.parent.highbgrate, prefix='bg')
        for i in nts:
            comuts[i,:] = np.ma.masked
            comuts[:,i] = np.ma.masked
        

        # compute the comut rates 
        for i in range(seqlen):
            
            # mask out near diagonals
            for j in range(i+1, min(i+self.parent.corrbuffer, seqlen)):
                comuts[i,j] = np.ma.masked
            
            # now go through remaining (valid) nts
            for j in range(i+self.parent.corrbuffer, seqlen):

                # convert to rate
                if not comuts.mask[i,j] and self.parent.ex_readarr[i,j] > depthcut:
                    comuts[i,j] /= self.parent.ex_readarr[i,j] 
                    
                    # do background subtraction
                    if self.parent.bg_readarr[i,j] > depthcut:
                        comuts[i,j] -= float(self.parent.bg_comutarr[i,j])/self.parent.bg_readarr[i,j]
                    else:
                        comuts[i,j] = np.ma.masked

                else:
                    comuts[i,j] = np.ma.masked
 
        # now compute the median over valid (non-masked) values
        median_comut = np.ma.median(comuts)
        median_count = comuts.count()
        
        # now compute median comut rate of unreactive positions by masking out reactive nts
        for i in range(comuts.shape[0]):
            
            react = safeReactivityMean( self.profile.normprofile[i:i+self.parent.window] )

            if react is np.nan or react > primary_reactivity:
                comuts[i,:] = np.ma.masked
                comuts[:,i] = np.ma.masked


        median_unreact_comut = np.ma.median(comuts)
        median_unreact_count = comuts.count()
       

        if median_comut < 5e-4 or median_unreact_comut < 1e-4:
            sys.stdout.write('******************************\n')
            sys.stdout.write('WARNING: Low comutation rates!\n') 
            sys.stdout.write('\tMedian for all nts = {0:.2e}  ({1} pairs)\n'.format(median_comut, median_count))
            sys.stdout.write('\tMedian for unreactive nts = {0:.2e}  ({1} pairs)\n'.format(median_unreact_comut, median_unreact_count))
            sys.stdout.write('PAIR-MaP data likely untrustworthy\n') 

            return False

        else:
            sys.stdout.write("Reactivity rate quality checks passed\n")
            sys.stdout.write('\tMedian for all nts = {0:.2e}  ({1} pairs)\n'.format(median_comut, median_count))
            sys.stdout.write('\tMedian for unreactive nts = {0:.2e}  ({1} pairs)\n'.format(median_unreact_comut, median_unreact_count))
            return True



    
    def plot(self, output, minz=2, maxz=6):
        """Plot reactivity data and primary+secondary pairmap correlations"""
        
        plot = ArcPlot()

        plot.reactprofile = self.profile.normprofile
        plot.reactprofileType = 'DMS'
        plot.seq = self.profile.sequence

    
        corrs = [(i+1,j+1, self.parent.getMeanZ(i,j)) for i,j in self.secondary]
        
        plot.plotAlphaGradient(corrs, (30,194,255), (0.1,0.5), minz, maxz, 
                               window=self.parent.window, panel=-1)

        corrs = [(i+1,j+1, self.parent.getMeanZ(i,j)) for i,j in self.primary]
        
        plot.plotAlphaGradient(corrs, (0, 0, 243), (0.5,0.9), minz, maxz, 
                               window=self.parent.window, panel=-1)
        
        
        # check if we need to throw error message
        if self.comut_rates[0] < 5e-4 or self.comut_rate[1] < 1e-4:
            msg = 'WARNING! Low comutation rates!\n'
            msg += '\tMedian for all nts = {0:.2e}\n'.format(self.comut_rates[0])
            msg += '\tMedian for unreactive nts = {0:.2e}\n'.format(self.comut_rates[1])
            msg += 'PAIR-MaP data likely untrustworthy\n'

            msg_pos = (0.5, 0.7)
            fontsize = len(self.profile.normprofile)/10.
            msg_kwargs = {'fontsize':fontsize, 'bbox':dict(facecolor='red', alpha=0.3)}
        else:
            msg = 'Median for all nts = {0:.2e}\n'.format(self.comut_rates[0])
            msg += 'Median for unreactive nts = {0:.2e}\n'.format(self.comut_rates[1])
            msg_pos = (0,1)
            msg_kwargs = {'fontsize':12}
        

        plot.writePlot(output, msg=msg, msg_pos=msg_pos, msg_kwargs=msg_kwargs)






#########################################################################################

def safeReactivityMean(rarray):
    """Compute the mean reactivity of the provided array, setting negative values to 0"""
    
    c = 0
    n = 0

    for r in rarray:
        if r == r and r>-10: # make sure defined
            n += 1
            c += max(r, 0) # don't let negatives skew
    
    if n==0:
        return np.nan
    else:
        return c/n

   




def parseArguments():

    parser = argparse.ArgumentParser(description = "Detect base pairs from correlated DMS data")
    
    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required arguments')

    required.add_argument('--modified_parsed', help='Path to Modified parsed.mut')
    required.add_argument('--untreated_parsed', help='Path to Untreated parsed.mut')
    required.add_argument('--profile', help='Path to profile.txt')
    required.add_argument('--out', help='output prefix')
    
    ######################################################
    # optional arguments

    optional.add_argument('--undersample', type=int, default=-1, help="""Randomly undersample specified number of reads 
                                           from modified file (default=-1 [disabled]).""")

    optional.add_argument('--ignorents', help="""A list of comma-separated (1-indexed) nts to ignore (e.g. 10,12,35)""")
 

    optional.add_argument('--mincoverage', type=float, default=0, help="""Quality filter reads by requiring a minimum 
                                           number of positional matches to reference sequence.
                                           This parameter can be set as a fraction of the total molecule length 
                                           (i.e. if 0.8 is passed, reads will be required to have at least 0.8*len(fasta) 
                                           valid matches) Alternatively, an integer value can be passed 
                                           (i.e. if 150 is passed, reads will be required to have at least 150 valid matches).
                                           By default, this filter is disabled.""")

    optional.add_argument('--chisq_cut', type=float, default=20.0, help="Set chisq cutoff (default = 20)")


    optional.add_argument('--highbg_rate', type=float, default = 0.02, help="""Ignore nt position with bg reactivity above
                                           this value. For PAIR-MaP analysis, windows with bg 3X this rate will be ignored""")

    optional.add_argument('--highbg_corr', type=float, default = 10.83, help="""Ignore nt pairs correlated in the bg 
                                           sample, with correlation determined via this significance value 
                                           (default=10.83 --> P=1e-3)""")

    optional.add_argument('--mincorrdistance', type=int, default=5, help="""Minimum distance allowed between correlations 
                                               (default = 5)""")
    
    optional.add_argument('--mindepth', type=int, default=10000, help="""Minimum pairwise read depth allowed for calculating 
                                        correlations (default = 10000)""")

    optional.add_argument('--mincount', type=int, default=50, help="""Minimum required count in contigency table 
                                      (default = 50). Nt pairs with fewer than this number of comutations are ignored""")
    
    optional.add_argument('--primary_reactivity', type=float, default=0.2, help="""Reactivity cutoff for primary correlations
                                      (default = 0.2)""")
    optional.add_argument('--primary_zscore', type=float, default=2.0, help="""Zscore cutoff for primary correlations
                                      (default = 2.0)""")
    optional.add_argument('--secondary_reactivity', type=float, default=0.5, help="""Reactivity cutoff for secondary correlations
                                      (default = 0.5)""")
    optional.add_argument('--secondary_zscore', type=float, default=2.0, help="""Zscore cutoff for secondary correlations
                                      (default = 2.0)""")

    parser._action_groups.append(optional)


    args = parser.parse_args()
    

    # check that required args are defined
    if not args.modified_parsed or not args.untreated_parsed \
            or not args.profile or not args.out:
        
        sys.stderr.write("! Required arguments not provided !\n\n")
        parser.print_help(sys.stderr)
        sys.exit(1)

    
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







if __name__ == '__main__':
    
    verbal = True

    args = parseArguments()
    
    # read in the reactivity profile information
    profile = ReactivityProfile(args.profile, bg=args.highbg_rate, ignorents=args.ignorents)
    profile.normalize(DMS=True)
    # write normed reactivities to file
    profile.writeReactivity('{}.dms'.format(args.out)) 
    

    
    # initialize ringmapper and data structures
    ringexp = RINGexperiment(corrtype = 'apc', arraysize=len(profile.sequence))
    ringexp.sequence = profile.sequence
        
    ringexp.initDataMatrices('ex', args.modified_parsed, window=3, 
                             mincoverage=args.mincoverage, undersample = args.undersample,
                             verbal = verbal)

    ringexp.initDataMatrices('bg', args.untreated_parsed, window=3, 
                             mincoverage=args.mincoverage, undersample = args.undersample,
                             verbal = verbal)

    
    


    # compute the correlation matrices
    ringexp.computeCorrelationMatrix(corrbuffer = args.mincorrdistance,
                                     mindepth = args.mindepth,
                                     mincount = args.mincount,
                                     ignorents = args.ignorents,
                                     highbgrate = args.highbg_rate,
                                     highbgcorr = args.highbg_corr,
                                     verbal = verbal)
    
    # compute zscores
    ringexp.computeZscores()

    # at this point, all the major calculations have been done
    
    # write all correlations to file
    ringexp.writeCorrelations('{0}-allcorrs.txt'.format(args.out), chi2cut=args.chisq_cut)
    
    
    # filter out pairs
    pairs = PairMapper(ringexp, profile, chi2cut=args.chisq_cut,
                       primary_reactivity=args.primary_reactivity, 
                       primary_zscore=args.primary_zscore,
                       secondary_reactivity=args.secondary_reactivity, 
                       secondary_zscore=args.secondary_zscore)
    

    if pairs.passfilter:

        # write out pairmap data
        pairs.writePairs('{}-pairmap.txt'.format(args.out))
    
        # write out folding restraint matrix
        pairs.writePairBonusFile('{}.bp'.format(args.out), args.chisq_cut, maxNC=1)


    # make plot
    #try:
    #    pairs.plot('{}-pairmapper.pdf'.format(args.out))
    #except:
    #    print("Was unable to generate plot")


       


# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
