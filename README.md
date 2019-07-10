# GQ-corr

This project aims to determine any dependencies of PQS (putative quadruplex sequences) locations on non-transcriptional factors. Calculations are made for human genome (hg38), by chromosome. For the first version, PQS were only found in coding strand.

Python script PQS_type.py enumerates G-quadruplexes divided into categories by loop length. File PQS_g2corr.py is used to calculate pair correlation function (PCF) of the PQS locations. 

Folder expr contains gene expression levels of several chromosomes from GTEx RNA-seq project, downloaded from UCSC genome browser. Folder hist contains PQS population histograms by chromosome (x10^4 scale to the X axis). Folder g2 contains correlation graphs for a few chromosomes.
