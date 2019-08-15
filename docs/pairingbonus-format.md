Base Pair Bonus Format
======================

Base pair bonuses may be written in one of two formats, both which 
are recognized by RNAstructure (v6.2 and higher)


Matrix Format
-------------
A square, symmetric (seqlength, seqlenth) matrix of energy bonuses

0.000 entry at i,j indicates no bonus 
-X.XXX entry at i,j indicates favorable bonus between nts i,j

Bonuses are reported in kcal/mol computed at 37° C 



List Format
-----------
A list of i,j positions with pairing bonuses

IMPORTANT: First row must be header beginning with ';' for proper recognition by RNAstructure

Following rows contain non-zero i,j, bonus entries

Bonuses are reported in kcal/mol computed at 37° C 




