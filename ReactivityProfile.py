
###############################################################################
#
# Basic info for reading in profile files from shapemapper and 
#   performing various operations on them
#
# Anthony Mustoe
# Copyright 2018
# 
# This file is licensed under the terms of the MIT license 
#
# Change-log
#
#
###############################################################################



# import other packages
import sys
import numpy as np

ntorder = ('A','C','G','U')


class ReactivityProfile(object):
    """ Object containing reactivity profile data
    Can contain:
        seqeunce    = nt sequence
        nts         = nt numbering
        rawprofile  = rxn rate 
        rawerror    = estimated error of rxn rate
        subprofile  = background subtracted profile 
        suberror    = background subtracted errors
        normprofile = normalized profile
        normerror   = normed errors
    """

    def __init__(self, plusfile = None, minusfile = None, **kwargs):

        self.sequence = None
        self.nts = None
        self.rawprofile = None
        self.rawerror = None
        self.backprofile = None
        self.backerror = None
        self.normprofile = None
        self.normerror = None
        self.subprofile = None
        self.suberror = None
        

        # set the default normMethod
        self.normMethod = self.normBox

        # set the default profile
        self.defaultProfile = 'norm'

        
        if plusfile:
            self.readProfile(plusfile, **kwargs)
           
        if minusfile:
            self.backgroundSubtract(minusfile, **kwargs)



    def initArray(self, name, mollen=None):

        if mollen is None:
            mollen = len(self.profile())
        
        prof = np.empty(mollen)
        prof[:] = np.nan
        setattr(self, name+'profile', prof)

        err = np.empty(mollen)
        err[:] = np.nan
        setattr(self, name+'error', err)


    def copy(self):

        out = ReactivityProfile()

        out.sequence = np.copy(self.sequence)
        out.nts = np.copy(self.nts)
        
        out.rawprofile = np.copy(self.rawprofile)
        out.backprofile = np.copy(self.backprofile)
        out.subprofile = np.copy(self.subprofile)
        out.normprofile = np.copy(self.normprofile)
        
        out.rawerror = np.copy(self.rawerror)
        out.backerror = np.copy(self.backerror)
        out.suberror = np.copy(self.suberror)
        out.normerror = np.copy(self.normerror)

        out.normMethod = self.normMethod

        # set the default profile
        out.defaultProfile = self.defaultProfile

        return out

    
    def cutProfile(self, start=None, stop=None):
        """return new ReactivityProfile Object
        Start is the nt number (typically 1-based)
        Stop is the nt number (inclusive, typically 1-based)
        """

        out = self.copy()
        
        try:
            start = start-out.nts[0]
        except TypeError:
            assert start is None

        try:
            stop = stop+1-out.nts[0]
        except TypeErorr:
            assert stop is None
        
        sel = slice(start, stop)

        out.sequence = out.sequence[sel]
        out.nts = np.arange(1, len(out.sequence)+1)
        
        out.rawprofile = out.rawprofile[sel]
        out.backprofile = out.backprofile[sel]
        out.subprofile = out.subprofile[sel]
        out.normprofile = out.normprofile[sel]
        
        out.rawerror = out.rawerror[sel]
        out.backerror = out.backerror[sel]
        out.suberror = out.suberror[sel]
        out.normerror = out.normerror[sel]

        return out




    def profile(self, name = None, err=False):
        """return either the default or specified profile"""

        if name is None:
            name = self.defaultProfile

        if not err:
            return getattr(self, name+'profile')
        else:
            return getattr(self, name+'profile'), getattr(self, name+'error')

        

        
    def readProfile(self, filename, **kwargs):
        """determine from the file extension what type of data it has, and then read"""

        ext = filename.split('.')[-1].lower()
        
        if ext == 'txt':
            # assume it has profile format
            prof = self.readProfileFile(filename, **kwargs) #

        elif ext == 'csv':
            # test whether it is a normal or pivoted file
            with open(filename,'rU') as f:
                spl = f.readline().split(',')[2]
            
            if spl=='number':
                self.readMutationCSV(filename, **kwargs) #
            else:
                self.readMutationCSVpivot(filename, **kwargs) #

        elif ext == 'map':
            self.readMapFile(filename, **kwargs) #
        
        elif ext == 'tab':
            prof = self.readTabFile(filename, **kwargs)

        else:
            raise IOError('unrecognized profile file extension :: %s' % ext)

        
        self.convertSequence()
    
    def convertSequence(self):
        """Make sure the sequence contains 'Us' vs. T'"""
        
        mask = self.sequence == 'T'
        self.sequence[mask] = 'U'
        

    def readProfileFile(self, filepath, bg=0.02, depthcut=100, ignorents =[], **kwargs):
        """read in Profile file output by new shapemapper"""
        
        seq = []
        num = []
        plus = []
        minus = []
        pdepth = []
        mdepth = []
        shape = []
        shapeerr = []
        
        with open(filepath, 'rU') as f:
            
            header = f.readline().split()
            header = [x.lower() for x in header]
    
            nt_idx = header.index('nucleotide')
            seq_idx = header.index('sequence')
            pd_idx = header.index('modified_effective_depth')
            pr_idx = header.index('modified_rate')
            md_idx = header.index('untreated_effective_depth')
            mr_idx = header.index('untreated_rate')
            try:
                s_idx = header.index('norm_profile')
                se_idx = header.index('norm_stderr')
            except ValueError:
                s_idx = None


            try:
                for line in f:
                    spl=line.split()
                    num.append( spl[nt_idx] )
                    seq.append(spl[seq_idx].upper())
                    pdepth.append( spl[pd_idx] )
                    plus.append( spl[pr_idx] )
                    mdepth.append( spl[md_idx] )
                    minus.append( spl[mr_idx] )
                    if s_idx is not None:
                        shape.append( spl[s_idx] ) 
                        shapeerr.append( spl[se_idx] )
            except:
                raise IOError("Unrecognized profile file format")

        
        self.sequence = np.array(seq)
        self.nts = np.array(num, dtype=int)
        self.backprofile = np.array(minus, dtype=float)
        self.rawprofile = np.array(plus, dtype=float)
        
        if len(shape) > 0:
            self.normprofile = np.array(shape, dtype=float)
            self.normerror = np.array(shapeerr, dtype=float)
        else:
            self.normprofile = np.zeros(self.rawprofile.shape)
            self.normerror = np.array(self.rawprofile.shape)

        # compute rawerror
        arr = np.array(pdepth, dtype=float)
        self.rawprofile[arr<depthcut] = np.nan
        self.rawerror = np.sqrt(self.rawprofile/arr)

        arr = np.array(mdepth, dtype=float)
        self.backprofile[arr<depthcut] = np.nan
        self.backerror = np.sqrt(self.backprofile/arr)
        
        self.backgroundSubtract(normalize=False)

        with np.errstate(invalid='ignore'):
            mask = (self.backprofile>bg ) # | np.isnan(self.normprofile) | (self.normprofile<-10)
            self.subprofile[mask] = np.nan
            self.normprofile[mask] = np.nan
            
        self.maskignore(ignorents)

        return None
    

    def maskignore(self, ignorents):
        
        mask = np.zeros(self.nts.size, dtype=bool)
        
        for i in ignorents:
            mask = mask | (self.nts == i)

        self.subprofile[mask] = np.nan
        self.normprofile[mask] = np.nan
        self.backprofile[mask] = np.nan
        self.rawprofile[mask] = np.nan
            


    def readMutationCSV(self, filepath, exclude = [], **kwargs):
        """read in a mutation count"""
        
        with open(filepath,'rU') as f:
            
            for line in f:
                line = line.strip(', \n')
        
                spl = line.split(',')
                
                if spl[2] == 'number':
                    nucs = np.array(spl[3:], dtype=int)
                    mutations = np.zeros(len(nucs))
                    excludeCount = np.zeros(len(nucs))

                elif spl[2] == 'sequence':
                    seq = np.array(spl[3:])

                elif ('del' in spl[2] or '->' in spl[2]):
                    if not (spl[2][0:3] in exclude or spl[2] in exclude):
                        mutations+=np.array(spl[3:], dtype=float)
                    else:
                        # tally excluded events so that we can subtract them from depth
                        excludeCount += np.array(spl[3:], dtype=float)

                elif spl[2] == 'depth':
                    depth = np.array(spl[3:], dtype = float)
                    break

        depth -= excludeCount

        # find locations where depth is zero and set to nan
        depth[depth==0] = np.nan
        
        # find locations where mutation rate is zero and set to nan
        mutations[mutations==0] = np.nan
        
        mutations /= depth
        stderr = np.sqrt(mutations)/np.sqrt(depth)
        

        self.sequence = seq
        self.nts = nucs
        self.rawprofile = mutations
        self.rawerror = stderr
        
 

    def readMutationCSVpivot(self, filepath, exclude=[], **kwargs):
        """read in a pivoted mutation count file"""
        
        f = open(filepath,'rU')
        # pop off the first two lines
        for i in range(2):
            f.readline()
            
        dkey = f.readline().strip(', \n').split(',')
        data = [[] for i in dkey]
        
        for line in f:
            line = line.strip(', \n')
 
            if len(line)==0:
                break

            spl = line.split(',')
            for i,v in enumerate(spl):
                data[i].append(v)
        

        for i,k in enumerate(dkey):
            if k == 'number':
                nucs = np.array(data[i], dtype=int)
                excludeCount = np.zeros(len(nucs))
                mutations = np.zeros(len(nucs))

            elif k== 'sequence':
                seq = np.array(data[i])

        
            elif ('del' in k or '->' in k):
                if not (k[0:3] in exclude or k in exclude):
                    mutations += np.array(data[i], dtype=float) 
                else:
                    excludeCount += np.array(data[i], dtype=float)

            elif k == 'depth':
                depth = np.array(data[i], dtype=float)
                break

        depth -= excludeCount

        # find locations where depth is zero and set to nan
        depth[depth==0] = np.nan
        
        # find locations where mutation rate is zero and set to nan
        mutations[mutations==0] = np.nan
        
        mutations /= depth
        stderr = np.sqrt(mutations)/np.sqrt(depth)
  
        self.sequence = seq
        self.nts = nucs
        self.rawprofile = mutations
        self.rawerror = stderr
        

    def readMapFile(self, filepath):
        """convert the .map file to a list of tuples"""
        
        seq, nucs, prof, err = '', [], [], []
        
        with open(filepath) as f:
            for line in f:
                spl = line.split()

                nucs.append(int(spl[0]))
                prof.append(float(spl[1]))
                err.append(float(spl[2]))
                seq+=spl[3]

        seq = np.array(list(seq))
        nucs = np.array(nucs)
        prof = np.array(prof)

        # convert -999 to nan
        prof[prof <= -10] = np.nan

        err = np.array(err)
        # convert errors to nan as well
        err[np.isnan(prof)] = np.nan
        
        self.sequence = seq
        self.nts = nucs
        self.normprofile = prof
        self.normerror = err



    def readTabFile(self, filepath, bg=0.02, **kwargs):
        """read in tab file"""

        with open(filepath, 'rU') as f:
    
            header = f.readline().split('\t')
            header = [x.strip() for x in header] 

            data = [ [] for i in header ]
    
            for line in f:
                spl = line.split('\t')
                for i,v in enumerate(spl):
                    data[i].append(v)

        
        self.sequence = np.array(data[1])
        self.nts = np.array(data[0], dtype=int)
        
        idx = header.index('rx rate')
        arr = np.array(data[idx], dtype=float)
        arr[arr<-10] = np.nan
        self.rawprofile = arr
        
        idx = header.index('bg rate')
        arr = np.array(data[idx], dtype=float)
        arr[arr<-10] = np.nan
        self.backprofile = arr
        
        idx = header.index('Normalized Reactivity')
        arr = np.array(data[idx], dtype=float)
        arr[arr<-10] = np.nan
        self.normprofile = arr
        
        idx = header.index('rx depth')
        arr = np.array(data[idx], dtype=float)
        arr[arr<1000] = 0
        self.rawprofile[arr<1000] = np.nan
        self.rawerror = np.sqrt(self.rawprofile/arr)
        
        idx = header.index('bg depth')
        arr = np.array(data[idx], dtype=float)
        arr[arr<1000] = 0
        self.backprofile[arr<1000] = np.nan
        self.backerror = np.sqrt(self.backprofile/arr)
        
        self.backgroundSubtract(normalize=False)

        with np.errstate(invalid='ignore'):
            mask = (self.backprofile>bg) | np.isnan(self.normprofile) | (self.normprofile<-10)
            self.subprofile[mask] = np.nan
            self.normprofile[mask] = np.nan
                              


    def normalize(self, DMS=False, byNT=False, name=None, normfactors = None, errfactors = None, **kwargs):
        """normalize the profile; overwrites values in normprofile
        By default, normalization is done in a sequence agnostic way (SHAPE default)
        If byNT, nts are normalized independently
        If DMS, nts are normalized w/ A+C and U+G on the same scale
        If normfactors is passed, use these precomputed normalization factors
           (dict of w/ A,G,U,C as keys and norm factors as values)
        If normfactors is passed, use errfactors (optional)
        """



        if name is None and self.subprofile is not None:
            name = 'sub'
        elif name is None:
            name = 'raw'
        

        prof,err = self.profile(name, True)
        
        # initialize the profile and error array
        normprof = np.array(prof)
        
        # initialize the error arrays
        with np.errstate(invalid='ignore'):
            if err is not None:
                normerr = np.zeros(err.shape)
                mask = prof > 0
                normerr[mask] = (err[mask]/prof[mask])**2
            else:
                normerr = None
        
        
        if normfactors is None:
            normfactors = {}
            errfactors = {}

            if byNT:
                for i in ntorder:
                    nfac, nerr = self.normMethod(self.reactivityByNt(nts=i, name=name))
                    normfactors[i] = nfac
                    errfactors[i] = nerr


            elif DMS:
                mask = (self.sequence == 'A') | (self.sequence=='C')
                nfac, nerr = self.norm90( normprof[mask] )
                for i in ('A','C'):
                    normfactors[i] = nfac
                    errfactors[i] = nerr
                
                mask = (self.sequence == 'G') | (self.sequence=='U')
                nfac, nerr = self.norm90( normprof[mask] )
                for i in ('G','U'):
                    normfactors[i] = nfac
                    errfactors[i] = nerr
            
            else:
                nfac, nerr = self.normMethod(normprof)
                for i in ntorder:
                    normfactors[i] = nfac
                    errfactors[i] = nerr
 
 

        # normalize the data
        for i in ntorder:
            mask = (self.sequence == i)
            normprof[mask] /= normfactors[i]

        self.normprofile = normprof

        if DMS:
            print("Renormalized data using DMS mode")
        elif byNT:
            print("Renormalized data using byNT mode")
        else:
            print("Renormalized data using standard mode")


        if normerr is not None and errfactors is not None:
            for i in ntorder:
                mask = (self.sequence == i)
                normerr[mask] += (errfactors[i]/normfactors[i])**2
                normerr[mask] = np.abs( normprof[mask] ) * np.sqrt( normerr[mask] )
            
            self.normerror = normerr

 
        return normfactors



    def normalize_external(self, profilefiles = [], profileobjs = [], **kwargs):
        """normalize reactivities using a set of other data to compute normfactors
        profilefiles = list of profile.txt files to normalize against
        profileobjs = list of profile objects to normalize against
        """

        combined = ReactivityProfile()

        # need to manually set dtype so appending below works properly
        combined.subprofile = np.array([], dtype=float)
        combined.sequence = np.array([])

        for f in profilefiles:
            prof = ReactivityProfile(f)
            combined.subprofile = np.append(combined.subprofile, prof.subprofile)
            combined.sequence = np.append(combined.sequence, prof.sequence)


        for prof in profileobjs:

            combined.subprofile = np.append(combined.subprofile, prof.subprofile)
            combined.sequence = np.append(combined.sequence, prof.sequence)

        if len(combined.subprofile)==0:
            raise AttributeError('No external reactivity data provided')

        nfacs = combined.normalize(**kwargs)

        self.normalize(normfactors=nfacs)

        return nfacs




    def backgroundSubtract(self, normalize=True, filepath=None, **kwargs):
        """Set subprofile
        By default, will subtract 'back' from 'raw' profile
        Alternatively, if 'back' is not set than can be provided via filepath"""

        
        if self.rawprofile is None:
            raise ValueError('rawprofile is None')
        
        if filepath is None and self.backprofile is None:
            raise ValueError('backprofile is None and no alternative filepath is provided')
        
        
        if filepath is not None:
            
            if self.backprofile is not None:
                print('Overiding backprofile with profile from:{0}'.format(filepath))


            prof = ReactivityProfile(filepath, **kwargs)

            if ( ( self.sequence is not None and not np.array_equal(prof.sequence, self.sequence)) 
                or ( self.nts is not None and not np.array_equal(prof.nts, self.nts) ) ):
                sys.stderr.write('Background profile does not have the same sequence or nt numbering\n')
        
            elif self.sequence is None:
                self.sequence = prof.sequence
                self.nts = prof.nts
        
            self.backprofile, self.backerror = prof.profile('raw', True)


        self.subprofile, self.suberror = self.computeProfileDiff(self, compname='back', myname='raw')
        
        if normalize:
            self.normalize(**kwargs)



    def computeProfileDiff(self, profile, compname=None, myname=None):
        """compute the difference between two profiles
        profile is another ReactivityProfile object
        """
        
        p,e = self.profile(myname, True)
        p2, e2 = profile.profile(compname, True)
        
        diff = p-p2
        
        #for i in range(len(p)):
        #    print "{0} {1:.3f} {2:.3f} {3:.3f}".format(i, p[i], p2[i], diff[i])


        if e is not None and e2 is not None:
            error = np.sqrt(e**2 + e2**2)
        else:
            error = None

        return diff, error

    
    def assignComparisonProfile(self, profile, compname=None, myname=None, zfact = True, **kwargs):
        """Assign the comparison profile arrays
        If zfact is true, will filter by whether differeces are significant, assigning
        diffs to zero otherwise
        """
        
        diff,err = self.computeProfileDiff(profile, compname, myname)

        if zfact:
            zfact = self.computeZfactor(profile, compname=compname, myname=myname, **kwargs)

            # filter out negative zvalues by setting to zero
            diff[zfact < 0] = 0
            err[zfact < 0] = 0
        
        self.compareprofile = diff
        self.compareerror = err
        

    def computeZfactor(self, profile, compname=None, myname=None, siglevel = 1.96, **kwargs):
        
        p1,e1 = self.profile(myname,True)
        p2,e2 = profile.profile(compname, True)

        zfact = []
        for i,v1 in enumerate(p1):
            v2 = p2[i]

            if np.isnan(v1) or np.isnan(v2):
                zfact.append(np.nan)
            else:
                top = siglevel * (e1[i] + e2[i])
                bot = abs(v1-v2)
                zfact.append(1-top/bot)
                
        return np.array(zfact)

    
    def normBox(self, data):
        """ NOTE: This method varies slightly from normalization method used 
        in the SHAPEMapper pipeline. Shapemapper sets undefined values to 0, and 
        then uses these values when computing iqr and 90th percentile. Including these
        values can skew these result. This method excludes such nan values. 
        Other elements are the same"""

        data2 = data[np.isfinite(data)]
        iq = np.percentile(data2, [25., 75.])
        iqr = iq[1]-iq[0]
    
        # filter by iqr
        data3 = data2[ data2 < (1.5*iqr+iq[1]) ]
        
        d2len = float(len(data2))
        d3len = float(len(data3))

        # see if too much data is classified as outliers
        if d2len<100 and d3len/d2len < 0.95:
            data3 = data2[data2<=np.percentile(data2, 95.)]
        elif d2len>=100 and d3len/d2len < 0.9:
            data3 = data2[data2<np.percentile(data2, 90.)]
       
        p90 = np.percentile(data3, 90.)
        data4 = data3[data3>p90]

        return np.mean(data4), np.std(data4)/np.sqrt(len(data4))

    def zeroNeg(self, name=None):
        """ zero out negative values in the profile name"""

        prof = self.profile(name)
        tmp = np.nan_to_num(prof)
        prof[tmp<0] = 1e-5

    
    def normWinsor(self, data):

        finitedata = data[ np.isfinite(data) ]
        fac = np.percentile(finitedata, 95.)

        return fac, -1 
        



    def norm90(self, data):
        
        finitedata = data[ np.isfinite(data) ]
        
        bnds = np.percentile(finitedata, [90., 99.])
        mask = (finitedata>= bnds[0]) & (finitedata<=bnds[1])

        normset = finitedata[mask]
               
        ave = np.mean( normset )
        std = np.std( normset )

        return ave, std/np.sqrt(len(normset))
            

    def writeReactivity(self, writename, name=None):
        """ Write out reactivity profile to the file"""

        p = self.profile(name)

        with open(writename, 'w') as out:
            for i, v in enumerate(p):
                
                if np.isnan(v):
                    v = -999

                out.write("{0} {1:.4f}\n".format(self.nts[i], v))
    
    def writeRNAstructureSeq(self, writename, header=None):
        
        if header is None:
            header = writename

        with open(writename,'w') as out:
            out.write(';\n{0}\n'.format(header))
            out.write('{0}1'.format(''.join(self.sequence)))

    

    def writeRxnColors(self, writename, unnorm=False, name = None,
                       colors = [(0.3, 0.7,'255,164,26'), (0.7,100, '255,0,0')]):
        
        """Write out colors explicity based on norm reaction rate
        Will also write out T and G nts as gray if they have a high dRxn if 
        """

        prof = self.profile(name)

        if unnorm:
            p2 = self.profile('sub')
            if p2 is None:
                p2 = self.profile('raw')


        out = open(writename, 'w')

        for i,v in enumerate(prof):
        
            col = 'white'
            
            if np.isnan(v) or v<-10:
                col = '180,180,180'
        
            for c in colors:
                if c[0] < v < c[1]:
                    col = c[2]
            
            if unnorm:  
                if self.sequence[i] in ('T', 'G') and p2[i] > 0.005:
                    col = '180,180,180'

            out.write("%d %s %s\n" % (self.nts[i], col, col))
        
        out.close()
        


    def plotProfileBar(self, writename, name=None, ylabel='Reactivity', title=None,
                       colors = [(-1e10,-10,'0.7'),(-10, 0.4, 'black'), 
                                 (0.4, 0.85,(1,.643,.102)), (0.85,100, (1,0,0))]):

        """ This is a general method for creating bar graphs in matplotlib
            An optional input is colors, a list of ybounds by which to color bars"""



        pro, err = self.profile(name, True)
    
        if err is None:
            err = np.zeros(len(pro))
        
        dseries = [ [[],[],[]] for i in range(len(colors)) ]
        nodata = []

        # assign each data point to its correct color
        for i,d in enumerate(pro):
            
            if np.isnan(d): 
                nodata.append(self.nts[i])
                continue

            for j, col in enumerate(colors):
            
                if col[0] < d <= col[1]:
                    dseries[j][0].append(self.nts[i])
                    dseries[j][1].append(d)
                    dseries[j][2].append(err[i])
                    break

        # figure out the axes
        xmin,xmax = min(self.nts)-1, max(self.nts)+1

        pfil = pro[np.isfinite(pro)]         
        ymin, ymax = min(pfil)*1.1, max(pfil)*1.1

        # now create the plot
        absoluteBarWidth = 0.05
        fig = plot.figure(figsize=(absoluteBarWidth*(xmax-xmin), 5))
  
        for i,dset in enumerate(dseries):      
            # data includes error bars
            plot.bar(dset[0], dset[1], align="center", width = 1.05, color=colors[i][2],
                     edgecolor = colors[i][2], linewidth = 0.0,
                     yerr = dset[2], ecolor = "0.7", capsize = 1)
 
        # now plot no data
        plot.bar(nodata, [ymax-ymin for i in nodata], bottom=[ymin for i in nodata], 
                color='0.7', width=1.05, align='center', linewidth=0.0)

        plot.xlabel('Nucleotides')
        plot.ylabel(ylabel)
    
        if title is not None:
            plot.title(title)
    

        plot.xlim(xmin, xmax)
        plot.xticks(np.arange(10,xmax,10), rotation=45)

        plot.ylim(ymin, ymax)
        plot.yticks(np.arange(np.ceil(ymin), np.ceil(ymax), 1))

        plot.grid(which='major',axis='y', color='0.5', linewidth=1.0)
        

        plot.savefig(writename)


    def reactivityByNt(self, resnums = None, nts=None, name = None):
        """ return a list of reactivities for a given set of nts (array), or nt type"""

        pro = self.profile(name)

        mask = np.isfinite(pro)
        #with np.errstate(invalid='ignore'):
        #    mask = mask & (pro > -0.3) & (pro < 4.0)

        try:
            ntit = iter(nts)
            ntmask = (self.sequence == next(ntit))
            for n in ntit:
                ntmask = ntmask | (self.sequence == n)
            
            mask = mask & ntmask
        except TypeError:
            pass

        try:
            resnums = set(resnums)
            mask = mask & np.array([i in resnums for i in self.nts])
        except TypeError:
            pass

        return pro[mask]


    def mutHistogram(self, name = None, nts = None, resnums = None,
                     bins=20,axes = None, writename='mutHist.pdf', **kwargs):

        write = False
        if axes is None:
            fig, axes = plot.subplots()
            write = True

        rxn = self.reactivityByNt(name=name, nts=nts, resnums = resnums)
        
        hist = axes.hist(rxn, histtype = 'step', bins = bins, **kwargs)
        
        if write:
            plot.savefig(writename)


        return hist


    def ntMutHistogram(self, writename='mutHist.pdf', name=None):

        fig, axes = plot.subplots(nrows = 2, ncols =2)

        for i,n in enumerate(ntorder):
            k,j = divmod(i,2)
            self.mutHistogram(nts=n, name=name, axes=axes[k,j])
    
            axes[k,j].set_title(n)

            for tick in axes[k,j].get_xticklabels():
                tick.set_rotation(-30)

        plot.tight_layout()
        plot.savefig(writename)
        plot.close()


    def plotMutBox(self, name=None, writename='boxplot.pdf', axes = None, **kwargs):

        data = []
        for n in ntorder:
            data.append(self.reactivityByNt(nts = n, name=name))
        
        write = False
        if axes is None:
            fig,axes = plot.subplots()
            write = True

        axes.boxplot(data, **kwargs) 

        if write:
            plot.xticks(np.arange(1,len(ntorder)+1), ntorder)
            plot.savefig(writename)


    
    def plotMutationRates(self, writename, title=None):
        """ plot the mutation rate profile """

        
        # figure out the axes
        xmin = min(self.nts)-1
        xmax = max(self.nts)+1

        # now create the plot
        fig = plot.figure()
  
        absoluteBarWidth = 0.05
        fig = plot.figure(figsize=(absoluteBarWidth*(xmax-xmin), 5))

        # data includes error bars
        plot.step(self.nts, self.rawprofile, where="mid")            

        if self.backprofile is not None:
            plot.step(self.nts, self.backprofile, where="mid")            


        plot.xlabel('Nucleotides')
        plot.ylabel('Mutation Rate')
    
        if title is not None:
            plot.title(title)


        #aplot.bar(dset[0], dset[1], align="center", width = 1.05, color=colors[i][2],
        #        edgecolor = colors[i][2], linewidth = 0.0,
        #        yerr = dset[2], ecolor = "0.7", capsize = 1)
       

        plot.xlim(xmin, xmax)
        plot.xticks(np.arange(10,xmax,10), rotation=45)
        
        plot.grid(which='major',axis='y', color='0.5', linewidth=1.0)
        
        plot.savefig(writename)










    



# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4
