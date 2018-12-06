# RingMapper
Code for performing RING-MaP and PAIR-MaP analysis

(RingMapper & PairMapper)

-------------------------------------------------------------
Copywrite 2018 Anthony Mustoe
This project is licensed under the terms of the MIT license

Contact: amustoe@unc.edu
-------------------------------------------------------------


Description
-----------
RingMapper & Pairmapper automate the detection of correlated modifcation 
events in mutational profiling (MaP) experiments.

RingMapper is a flexible code for computing correlations (RINGs)
between all nucleotides in an RNA. Users can choose to compute
correlations between single nucleotides (window=1; default) or between larger
nucleotide windows (e.g. window=3, which tests whether 3-nt windows are
correlated with each other). Run ringmapper.py --help for options.

PairMapper is a wrapper for RingMapper that performs PAIR-MaP analysis --
i.e. searches for correlations indicative of secondary structure base pairs. 
Correlation calculations are done using a 3-nt window and then filtered by 
nt complementarity and reactivity. Run pairmapper.py --help for options.

Both RingMapper and PairMapper require read alignment and preprocessing 
by ShapeMapper2 using the --output-parsed flag. 
See the ShapeMapper2 documentation for further information regarding
alignment and processing options.

Also included in this distribution is arcPlot.py, which can be used to plot
RINGs and/or PAIR-MaP files. Run arcPlot.py --help for options.



Installation
------------
- Compile cython routines by running:
	python setup.py build_ext --inplace



Dependencies
------------
- python 2.7, numpy
- cython (developed using v0.23)
- matplotlib (for arcPlot)
- RNAtools2



Example workflow
----------------
(1) First process data using Shapemapper
	
    shapemapper --name example --target rna.fa --output-parsed \
      --modified --folder Modified_example --untreated --folder Untreated_example

    This will generate a number of files. You need the following:
	example_Untreated_rna_parsed.mut
	example_Modified_rna_parsed.mut
	example_rna_profile.txt


(2) Measure correlations using ringmapper
	
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
	example.dms              (normalized DMS reactivities)
	example.bp	         (PAIR-MAP Energy restraints for structure
                                  modeling with RNAstructure)
	example-pairmapper.pdf	 (Figure with DMS reactivities and PAIR-MaP
                                  correlations. Note this figure will not be 
                                  generated if your system does not have
                                  matplotlib installed or sometimes if the 
                                  X server is not setup correctly)


(4) Plot correlation data according to Z-score using arcPlot.py
 
    arcPlot.py --fasta rna.fa --ringz example-corrs.txt \
      --dmsprofile example.dms example-corrs.pdf


(5) Plot PAIR-MaP data using arcPlot.py
 
    arcPlot.py --fasta rna.fa --pairmap example-pairmap.txt \
      --dmsprofile example.dms example-pairmap.pdf



Timing/Performance
----------------
A rule of thumb is to expect between 1000-10000 reads/sec. Large RNAs will be 
on the slow side (1000-5000 reads/sec) whereas small RNAs will be on
the fast side. For example, for a 1,000 nt RNA with 10,000,000 aligned reads, 
RingMapper/ShapeMapper processing should take roughly 1 hour.







