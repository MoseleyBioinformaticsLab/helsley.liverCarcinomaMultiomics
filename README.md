## Code and Data for Multi-omic Analyses Reveal Novel Gene-Metabolite Relationships in Human Steatohepatitic-driven HCC

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18233569.svg)](https://doi.org/10.5281/zenodo.18233569)

This repository provides the data and code for the analyses reported in "Multi-omic Analyses Reveal Novel Gene-Metabolite Relationships in Human Steatohepatitic-driven HCC".

The permanent record is archived on [Zenodo](https://doi.org/10.5281/zenodo.18227400).

## Use

If the GitHub repo is still available, you can clone that repository directly using `git`:

```
git clone https://github.com/MoseleyBioinformaticsLab/helsley.liverCarcinomaMultiomics.git
```

If not, then you should download the repo from Zenodo:

```
wget 'https://zenodo.org/records/18233569/files/MoseleyBioinformaticsLab/helsley.liverCarcinomaMultiomics-v_1.1.zip?download=1' --output-document=helsley.liverCarcinomaMultiomics.zip
unzip helsley.liverCarcinomaMultiomics.zip
```

Then you should `cd` into the directory, and download the Zenodo repository with the `data` and `_targets` directories.

```
cd helsley.liverCarcinomaMultiomics
wget 'https://zenodo.org/records/18271120/files/helsley_hcc_targets_data.tgz?download=1' --output-document=hcc_targets_data.tgz
tar -xzvf hcc_targets_data.tgz
```

Start an R session within that directory, and then you should be able to install `renv` and the associated packages.

```r
install.packages("renv")
library(renv)
renv::restore()
```