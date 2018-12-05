# AID-ISA

The increasing availability of gene expression data has encouraged the development of tailored intelligent data analysis techniques. Grouping genes characterized by similar expression patterns is a widespread accepted (and often mandatory) analysis step. Despite the fact that a number of biclustering methods have been developed to discover clusters of genes exhibiting a similar expression profile under a subgroup of experimental conditions, approaches based on similarity measures computed only on expression profiles may lead to biologically meaningless groups. The integration of additional information, such as functional annotations, into biclustering algorithms can instead provide an effective support for identifying meaningful gene associations.

To solve this problem, we developed a new biclustering approach called Additional Information Driven Iterative Signature Algorithm, AID-ISA. It supports the extraction of biologically relevant biclusters by leveraging additional knowledge.
AID-ISA is based on a refinement process, called Additional Information-Driven (AID), embedded into the well known Iterative Signature Algorithm (ISA). AID-ISA takes both expression profiles and additional sources of information into account to discover biologically meaningful biclusters.

## Citation

Please cite AID-ISA as:

> Visconti A., Cordero F., and Pensa R.G., *Leveraging additional knowledge to support coherent bicluster discovery in gene expression data*, Intelligent Data Analysis (2014), [DOI: 10.3233/IDA-140671](https://content.iospress.com/articles/intelligent-data-analysis/ida00671) 

[Read the PDF here](https://pdfs.semanticscholar.org/a326/786cd7c243c20a64cdf02208597fd5789a6f.pdf)


## License

AID-ISA is licensed under GNU GPL v3.


