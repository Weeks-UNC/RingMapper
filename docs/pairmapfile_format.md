PAIR-MaP File Format
====================

The first row contains header metadeta indicating the sequence length 
    and how correlations were computed, which enables recognition by
    other software tools


The second row provides labels for each column, briefly described below:

i,j         Indexes of correlated positions, denoting the 5' start of the window.
            Indexes are 1-based and inclusive. For example, i=3 denotes the 3-nt
            window containing (3,4,5)

Sig.        The statistical significance of the correlation (APC G-statistic; G_APC)

Class       1 indicates the correlation is "principal/primary"
            2 indicates the correlation is "minor/secondary"
            0 denotes complementary positions that are significantly 
              correlated (G_APC > 20) but which do not pass PAIR-MaP filters

Zij         Average of Zi and Zj

Zi, Zj      Z-score of the G_APC computed relative to all other correlations 
            originating at position i and j, respectively

Mod_Depth   Number of times windows i,j were read together in the Modified sample

Mod_Comuts  Number of times windows i,j were mutated together in the Modified sample

Unt_Depth   Number of times windows i,j were read together in the Untreated sample

Unt_Comuts  Number of times windows i,j were mutated together in the Untreated sample



