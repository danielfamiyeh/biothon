# Biothon
Biothon provides a library of functions and structures that facilitate the analysis of biological data, namely DNA, RNA and protein  sequences. Numerous functions for global and local pairwise sequence alignments allow for regions of similarity and conserved domains to be identified between pairs of  biosequences. Multiple sequence alignment is  integrated into the library, using a distance-based, neighbour-joining implementation of the ClustalW algorithm. This same idea of hierarchical cluster analysis also facilitates the construction of phylogenetic dendrograms from lists of biosequences allowing for potential homologies between closely and distantly-related taxa to be shown visually. The library is standalone, requiring only the Python Standard Library to use.
## Features
- [x] Basic Biosequence Manipulation (Transciption with slicing, Translation, G-C content and More)
- [x] Pairwise Sequence Alignment (Needleman-Wunsch, Smith-Waterman, Hirschberg's, FASTA and BLAST)
- [x] Multiple Sequence ALignment (ClustalW)
- [x] Clustering (Phylogenetic Tree Construction, k-Medoid)
