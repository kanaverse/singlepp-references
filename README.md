# Reference datasets for `singlepp`

## Overview

This provides a number of prebuilt references for easy consumption by the [**singlepp**](https://github.com/LTLA/singlepp) library.
It is primarily based on the datasets available in the [**celldex**](https://bioconductor.org/packages/celldex) package,
themselves derived from the original [**SingleR**](https://bioconductor.org/packages/SingleR) publication.
Each reference dataset is represented by four files that are described below.

## File specification

### `matrix.csv`

This is a CSV file containing a dense matrix of ranks for all samples.
Each row corresponds to a sample where ranks are computed within each sample.
We store the ranks to avoid the need to store and transmit floating-point values.
This file does not contain any row or column names.

### `genes.csv`

This is a CSV file containing information about each gene in `matrix.csv`.
The first column contains the Ensembl identifiers, the second column contains the gene symbols, and the third column contains Entrez IDs.
Missing identifiers in any field are represented as empty strings.
No row or column names should be listed here.

### `labels_*.csv`

This is a CSV file containing listing the label for each sample in `matrix.csv`.
Each row corresponds to a column in `matrix.csv` and contains the label for the corresponding sample.
No row or column names should be listed here.
Each `matrix.csv` may be associated with multiple label sets (e.g., broad or fine labels).

### `markers_*.gmt`

This is a tab-delimited GMT file containing the ranked markers for each pairwise comparison between labels.
On each line, the first entry is label `X` and the second entry is label `Y`.
The remaining entries represent the row indices of markers defining `X` relative to `Y`.
Row indices are 0-based and refer to the rows of `matrix.csv`.
Note that indices are ranked, i.e., the strongest markers are listed first.
