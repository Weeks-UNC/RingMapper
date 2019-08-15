RING-MaP File Format
====================

The first row contains header metadeta:
<seqlength>  Window=<window>    Metric=<metric>

    <seqlength>     length of the RNA sequence inferred by ringmapper

    <window>        window sized used for computing correlations

    <metric>        metric used for computing correlations


The second row list labels for each data column, briefly described below:

i,j         Indexes of correlated positions, denoting the 5' start of the window.
            Indexes are 1-based and inclusive. For example, i=3 denotes the 3-nt
            window containing (3,4,5)

Statistic   The correlation statistic (Chi2/G/APC, depending on selected metric)

+/-         Sign of the correlation. (+) indicates that the joint probability (pij)
            of comutation is greater than the product of the marginal probabilities:
            pij > pi*pj. (-) indicates the converse.


Zij         Average of Zi and Zj

Zi, Zj      Z-score of the statistic computed relative to all other correlations 
            originating at position i and j, respectively

Mod_Depth   Number of times windows i,j were read together in the Modified sample

Mod_Comuts  Number of times windows i,j were mutated together in the Modified sample

Unt_Depth   Number of times windows i,j were read together in the Untreated sample

Unt_Comuts  Number of times windows i,j were mutated together in the Untreated sample



