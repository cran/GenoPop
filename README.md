# GenoPop

GenoPop is a R package designed to assist with population genomic analyses of data sets from non-model organisms or with low sequencing quality. It's created with the intention to simplify and streamline the analysis of large genomic data sets in VCF (Variant Call Format) files in a efficient manner, while handling problems of missing data.

The GenoPop package can be divided into parts. One part it the genotype imputation method GenoPop-impute, which is described in this [preprint](https://doi.org/10.22541/au.172515591.10119928/v1). The second part contains several function to calculate commonly used population genomics metrics, like Fst, and Dxy.

This document will give an overview about GenoPops functions and usability. Starting with an overview of its methods, a guide how to install and get started with the package.

## GenoPop-Impute Overview

GenoPop-Impute is a genotype imputation algorithm specifically designed for whole genome data sets. It does not require a SNP reference panel and, therefore, can be used for non-model organisms. A key aspect of GenoPop-Impute is its approach to handling large genomic datasets. Recognizing that SNPs within a linkage block share the same evolutionary history, GenoPop-Impute employs the assumption that these SNPs exhibit more comparable patterns than those from different linkage blocks. This assumption justifies segmenting the dataset into smaller blocks for parallel processing. Essentially, GenoPop-Impute performs batch-based imputation, where each batch contains SNPs likely to be correlated due to their close proximity in the genome and linkage disequilibrium. Using the missForest algorithm (Stekhoven & Bühlmanm, 2012) for the imputation of each batch, this approach enhances the efficiency and accuracy of the imputation process. A guide on how to execute GenoPop-Impute can be found using the R help options.

## Population Genomics Metrics Overview

### General Functionality

Each function in this part of GenoPop is designed to calculate specific population genomics metrics directly from VCF (Variant Call Format) files. These functions are designed for efficiency and handle large genomic datasets by processing data in parallel. For this, there are two modi available: processing in batches of equal numbers of SNPs and processing in windows of a specific genomic size in base pairs. In batch mode, the entire VCF file is processed at once to provide a general overview. In window mode, the file is processed in genomic sections to identify specific regions of interest. These functions typically return single metrics for batch mode or data frames detailing metrics per window.

### Metrics Overview

- **FixedSites**: Counts the number of sites fixed for the alternative allele. It helps identify regions with a complete fixation of an allele, potentially indicating selective sweeps or other evolutionary pressures.

- **SegregatingSites**: Counts the number of polymorphic sites, which are not fixed for the alternative allele. It's a measure of genetic variability within the population.

- **SingletonSites**: Counts the number of singleton sites, where a minor allele occurs only once in the sample. It can be an indicator of recent mutations.

- **PrivateAlleles**: Calculates the number of private alleles in two populations. Private alleles are present in one population but absent in another, providing insight into population differentiation.

- **Heterozygosity Rate**: Calculates the observed heterozygosity for each variant. It's a measure of genetic diversity within a population.

- **NucleotideDiversity (Pi)**: Measures the average number of nucleotide differences per site between two sequences. It's a key indicator of genetic diversity within a population.

- **Tajima's D**: A neutrality test comparing the number of segregating sites to the average number of nucleotide differences. It can suggest population expansion, selection, or bottlenecks.

- **Watterson's Theta**: A measure of genetic diversity within a population, based on the number of segregating sites.

- **Average Nucleotide Differences (Dxy)**: Measures the average number of nucleotide differences per site between two populations. It's a measure of genetic differentiation.

- **Fst**: The fixation index, measuring genetic differentiation between populations. It ranges from 0 (no differentiation) to 1 (complete differentiation).

- **OneDimSFS**: Calculates a one-dimensional site frequency spectrum, either folded or unfolded. It provides insights into allele frequency distributions within a population.

- **TwoDimSFS**: Calculates a two-dimensional site frequency spectrum for two populations. It's used to infer demographic history and population relationships.

Please note that this summary provides an overview of the functions and  their purposes. For complete understanding and appropriate usage, refer  to the detailed documentation of each function.
