# RingMapper v1.3

Code for performing RING-MaP and PAIR-MaP analysis
(RingMapper & PairMapper)

*Copywrite 2019,2020 Anthony Mustoe*. This project is licensed under the terms of the MIT license

Contact: anthony.mustoe@bcm.edu

-------------------------------------------------------------


General Description
-------------------
RingMapper & PairMapper automate the detection of correlated modifcation 
events in mutational profiling (MaP) experiments.

RingMapper is a flexible code for computing correlations (RINGs)
between all nucleotides in an RNA. Users can choose to compute
correlations between single nucleotides (window=1; default) or between larger
nucleotide windows (e.g. window=3, which tests whether 3-nt windows are
correlated with each other). Run ringmapper.py with the --help flag for options.

PairMapper is a wrapper for RingMapper that performs PAIR-MaP analysis,
i.e. it searches for correlations indicative of secondary structure base pairs. 
Correlation calculations are done using a 3-nt window and then filtered by 
nt complementarity and reactivity. Run pairmapper.py with the --help flag for options.

Both RingMapper and PairMapper require read alignment and preprocessing 
by ShapeMapper2 using the --output-parsed flag. 
See the ShapeMapper2 documentation for further information regarding
alignment and processing options.

Following the release of ShapeMapper2.3 RingMapper may now be used to compute
N1/3-N7 and N7-N7 RINGs in addition to traditional N1/3-N1/3 RINGs. See the
section titled "RingMapper Usage for N7" below. Note - PairMapper not currently
tested for compatibility with N7 functionality.

RING and PAIR-MaP data can be visualized using the arcPlot plotting tool,
available at https://github.com/Weeks-UNC/arcPlot


Installation
------------
- Compile cython routines by running:
```
	python setup.py build_ext --inplace
```

- Run test to make sure everything is working appropriately
```    
    ./test.sh
```    
    Note that the test will take ~20 seconds and requires ~200 MB of working disk space


Dependencies
------------
- python 2.7, numpy
- cython (developed using v0.23)
- ShapeMapper version 2.2 for pre-processing data


Citation
--------
Please cite:

A.M. Mustoe, N.N. Lama, P.S. Irving, S.W. Olson, K.M. Weeks. RNA base-pairing complexity in living cells visualized by correlated chemical probing, PNAS (2019).




--------------------------------------------------------------------------

RingMapper Usage
-----------------
```
ringmapper.py <optional args> inputfile outputfile | --help
```

### Required inputs
```
inputfile       Parsed mutation file of modified sample
outputpath      Path where output Ringfile will be written
```

### Optional arguments
```
--fasta         Fasta used for shapemapper alignment. Used to compute molsize.

--molsize       Size of molecule (if fasta argument not used). Default = 1000.
                Note that molsize must be greater than or equal to the size of the RNA

--untreated     Parsed mutation file for untreated sample. Nts with high mutation rates
                or correlations in the untreated sample are ignored.

--window        Window size used for calculation correlations. Default = 1.

--metric        Metric to use for computing correlations. Default = APC G-test (apc).
                Other options include Yates Chi2 (chi) and uncorrected G (g).

A full list of optional arguments/parameters can be accessed by running with --help flag
```

### Outputs
```
Ringfile        List of correlations. See docs/ringfile-format.txt
```

--------------------------------------------------------------------------

RingMapper Usage for N7
-----------------
```
ringmapper.py <optional args> --parsedMutga parsed.mutga --fasta fasta.fa inputfile outputfile | --help
```

### Required inputs
```
inputfile       Parsed mutation file of modified sample
--parsedMutga   Parsed mutation file of modified sample containing N7 modifications
--fasta         Path to fasta sequence file
outputpath      Path where output Ringfiles will be written
```

### Optional arguments
```
--untreatedMutga   Path to mutation file for untreated sample containing N7 modifications.
                   Nts with high mutation rates or correlations in the untreated sample
                   are ignored (Required --untreated)

--keepconcat       Keep the concatenated mutation string files for use in other analyses

A full list of optional arguments/parameters can be accessed by running with --help flag
```

### Outputs
```
Ringfile_concatrings   List of N1/3-N1/3, N1/3-N7, and N7-N7 correlations. See docs/ringfile-format.txt
Ringfile_N1            List of N1/3-N1/3 correlations. See docs/ringfile-format.txt
Ringfile_N1N7          List of N1/3-N7 correlations. See docs/ringfile-format.txt
Ringfile_N7            List of N7-N7 correlations. See docs/ringfile-format.txt
```

--------------------------------------------------------------------------

PairMapper Usage
-----------------
```
pairmapper.py <optional args> --modified_parsed <modfile> --untreated_parsed <untfile> --profile <profile> --out <prefix> | --help
```

### Required inputs
```
--modified_parsed <modfile>     Parsed mutation file of modified sample 

--profile <profile>             Profile file output by shapemapper

--out <prefix>                  Prefix used for output files
```

### Optional arguments
```
--untreated_parsed <untfile>    Parsed mutation file of untreated sample

--override_qualcheck            Pairmapper will refrain from giving results if the RNA does not pass modification rate 
                                thresholds. This flag will override this internal quality check and produce results 
                                regardless of underlying data quality. Use caution when intrepretting such results.

--renormalize                   Renormalize DMS data using internal normalization as described in Mustoe et al, 2019. 
                                If using ShapeMapper2.2 --dms preprocessing, DO NOT use this option.

Additionally, most internal filtering parameters can be modified if desired. 
A full list of modifiable parameters can be accessed by running with --help flag.
Note that only default parameters have been benchmarked

```

### Outputs
```
<prefix>-allcorrs.txt     Ringfile output containing all correlations computed 
                          using a 3-nt window. 
                          See docs/ringfile-format.txt

<prefix>-pairmap.txt      Pairmapfile output containing PairMapper filtered correlations.
                          See docs/pairmapfile-format.txt

<prefix>.bp	              Base pair bonus file for PAIR-MaP structure modeling in RNAstructure.
                          See docs/pairingbonus-format.txt

<prefix>.dms              Normalized DMS reactivities. [Only output if using --renormalize option] 
                          See docs/dmsreactivity-format.txt


```


--------------------------------------------------------------------------

Note on Shapemapper minimum mutation separation 0 parsed mutation files
-----------------------------------------------------------------------
Parsed mutation files (*.mut) generated by Shapemapper when using minimum mutation
separation 0 require additional filtering to achieve optimal base pair
detection by PairMapper. A separate filter program, ShapeMapper_MMS0_Mut_Filter.py,
is installed with this distribution. Please be aware that additional disk
space equivalent to the size of the source *.mut files is required for
this filtering to occur. 

Both the modified and untreated parsed mutation files require filtering.

### Usage
```
python ShapeMapper_MMS0_Mut_Filter.py <input>.mut <optional args>
```

### Required input
```
<input>.mut               Parsed mutation file (*.mut) generated from 
                          ShapeMapper with data generated using minimum 
                          mutation separation 0 (--min-mutation-separation 0)
```

### Optional arguments
```
--output_name             Output name for new filtered parsed mutation 
                          files created by this program, prepended by 
                          'MMS0Filter_'. If unspecified, it uses the 
                          original <input>.mut name

--output_dir              Output directory for new filtered parsed 
                          mutation files. If unspecified, the default 
                          directory is 'MMS0_Filter_Output'
```


--------------------------------------------------------------------------

Note on reproducing legacy (v1.1 or below) results
--------------------------------------------------
The default filtering parameters have been updated in v1.2. 
To run ringmapper using original parameters, use the following flag:
--mincount 10

To run pairmapper using original parameters, use the following flag:
--renormalize --secondary_reactivity 0.5 --mincount 50


----------------------------------------------------------------

Timing/Performance
----------------
A rule of thumb is to expect between 1000-10000 reads/sec. Large RNAs will be 
on the slow side (1000-5000 reads/sec) whereas small RNAs will be on
the fast side. For example, for a 1,000 nt RNA with 10,000,000 aligned reads, 
RingMapper/PairMapper processing should take roughly 1 hour.


----------------------------------------------------------------

Note on experimental modification rates
---------------------------------------

PAIR-MaP analysis can only be performed on MaP datasets containing high levels of chemical modification. Datasets are automatically checked to make sure they exceed minimal modification thresholds, corresponding to median comodification rates >0.0001. PAIR-MaP output files will not be written if this quality check fails.

Standard RING-MaP analysis can be performed on any MaP dataset. However, sensitivity again depends strongly on modification rates.

As a general rule of thumb, sucessful RING-MaP and PAIR-MaP analysis requires experimental conditions be optimized to achieve 95th percentile modification rates = ~10%. 


----------------------------------------------------------------


Example workflow
----------------
(1) First process data using Shapemapper
	
    shapemapper --name example --target rna.fa --output-parsed \
      --modified --folder Modified_example --untreated --folder Untreated_example

    This will generate a number of files. You need the following:
  example_Untreated_rna_parsed.mut
	example_Modified_rna_parsed.mut
	example_rna_profile.txt

(1a) If using minimum mutation separation 0 in Shapemapper, filter the .mut files

    python ShapeMapper_MMS0_Mut_Filter.py example_Untreated_rna_parsed.mut
    python ShapeMapper_MMS0_Mut_Filter.py example_Modified_rna_parsed.mut

    This generates the following output .mut files:
  MMS0Filter_example_Untreated_rna_parsed.mut
	MMS0Filter_example_Modified_rna_parsed.mut



(2) Measure correlations using ringmapper (using default window-size = 1)
	
    ringmapper.py --fasta rna.fa --untreated example_Untreated_rna_parsed.mut \
      example_Modified_rna_parsed.mut example-corrs.txt

    Correlations will be written into example-corrs.txt 


(3) Measure PAIR-MaP signals
	
    pairmapper.py --profile example_rna_profile.txt \
      --untreated example_Untreated_rna_parsed.mut \
      --modified example_Modified_rna_parsed.mut \
      --out example

    This will generate the following files:
	example-allcorrs.txt     (All correlations for 3-nt window)
	example-pairmap.txt      (PAIR-MaP filtered data)
	example.dms              (Normalized DMS reactivities)
	example.bp	             (Base pair bonuses for RNAstructure modeling)
	

(4) Plot correlation data according to Z-score using arcPlot.py
 
    arcPlot.py --fasta rna.fa --ringz example-corrs.txt \
      --dmsprofile example.dms example-corrs.pdf


(5) Plot PAIR-MaP data using arcPlot.py
 
    arcPlot.py --fasta rna.fa --pairmap example-pairmap.txt \
      --dmsprofile example.dms example-pairmap.pdf



