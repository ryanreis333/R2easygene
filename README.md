# R2easygene is a package intended to provide an easy and quick way to update gene names of a seurat object.

### It provides a fast and simple way to update gene names to the latest approved symbols, ensuring consistency across samples, especially in projects that compile data from multiple experiments. It addresses discrepancies in gene naming conventions that may cause batch effects when raw data is not available for realignment.

### The motivation for this package came from a finding in my own work in which I found that gene aliases were introducing significant batch effect, even after integration methods were implemented. I discovered that these alternate and outdated, as well as sample-specific gene names, were a small but significant source of noise in the data. I hope that this provides an easy and relatively quick way of removing this source of variation.

![Logo](images/HGNC_image.png)

### The package relies on the HGNChelper package from the HUGO Gene Nomenclature Committee, above.

### Please reach out or submit requests/pulls with any questions, comments, or if you can help. Thanks!
