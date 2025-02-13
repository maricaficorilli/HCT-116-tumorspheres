# Tailoring a novel colorectal cancer stem cell-targeted therapy by inhibiting the methyltransferase SMYD3

This script takes samples of HCT-116 tumorspheres as input to perform differential expression analysis (DEA) using the edgeR Bioconductor/R package (version 3.42.4). Following DEA, an enrichment analysis is conducted on the differentially expressed (DE) genes using the clusterProfiler Bioconductor/R package (version 4.8.3).

In DEA-Results folder, there is file on the DEA results of the comparison. In this case we used a log fold change (logFC) >= 1 or <= -1 and False Discovery Rate (FDR) < 0.05.

