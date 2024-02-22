# TopoWindows

TopoWindows offers R functions to calculate phylogenetic trees in genomic intervals (sliding windows or user-specified regions) based on vcf files. For now only Neighbor-Joining trees can be computed, but there is an option to write sequences in fasta format for external maximum-likelihood or bayesian analyses.

# Installation

To use TopoWindows in an R environment, simply run `source("/path/to/Topo_windows_v02.R")`.

Two dependencies must be installed beforehand:
- `ape` (https://cran.r-project.org/web/packages/ape/index.html)
- `vcfR` (https://cran.r-project.org/web/packages/vcfR/index.html)

To run TopoWindows from a unix command line environment (see below), `optparse` (https://cran.r-project.org/web/packages/optparse/index.html) must also be installed.

# Usage

There are two ways to run TopoWindows: either directly within an R environment, or from a unix command line environment (which is useful to run on a High Performance Cluster).

## Running in an R environment

Once TopoWindows is loaded into the environment, three main functions are available to calculate trees in sliding windows or in user-specified regions.

### Trees in sliding windows

Use either `topo.windows.sites(vcf, size, incr, phased, prefix, write.seq, nj, dna.dist)` or `topo.windows.coord(vcf, size, incr, phased, prefix, write.seq, nj, dna.dist)`. `topo.windows.sites` defines sliding windows based on a fixed number of SNPs, while `topo.windows.coord` defines them based on fixed coordinates. In both cases, the first window starts with the first position in the vcf file, and the last one ends with the end of the vcf.

*Arguments:*

- `vcf`: path to the vcf file to use. Can be gzipped. *WARNING: as of now, the script is not able to handle vcfs that contain several chromosomes/contigs, so these must be split beforehand and the script run on each independently. I am currently working on an update to account for that.*
- `size`: size of the windows (either SNPs or bp).
- `incr`: overlap between windows. incr=0 means no overlap.
- `phased`: boolean defining whether the vcf is phased. If `T`, then two sequences are used for each (diploid) individuals, corresponding to the two haplotypes. If `F`, a consensus sequence is called, using IUPAC code for heterozyguous sites.
- `prefix`: a prefix for the output files.
- `write.seq`: a boolean defining whether to write a sequence in fasta format for each window.
- `nj`: a boolean defining whether to calculate neighbor-joining trees. if `T`, the `njs` function of `ape` is used to accomodate for missing data in the distance matrix (see `ape` manual for more details). 
- `dna.dist`: the model used to calculate the distance matrices for the nj trees (e.g., "raw", "JC69"). See manual for function `dist.dna` in the `ape` package for more details. Note that a pairwise deletion approach is used. Some models may return missing values if two sequences are too different, which often happens when there are lots of missing data. In such cases it may be worth using `dna.dist="raw"`. Distance calculation will also fail if one or more sequence(s) contains only missing data. Such sequences are thus removed before calculating the tree (see log file below).

*Output:*

Three output files are written:
- `prefix_NJ_trees.trees`: the NJ trees, one per line in newick format (only if `nj=T`). If tree inference failed, NA is written instead.
- `prefix_stats.tsv`: a tab-separated table giving information about the windows analyzed. It contains the following columns: CHR=chromosome, CHR.START=starting coordinate, CHR.END=ending coordinate, WIN.SIZE=size of the window in bp, NSITES=number of SNPs in the window, PROP.PIS=proportion of parsimony informative sites in the window, PROP.MISS=proportion of missing genotypes in the window, TREE=whether a tree could be inferred (YES) or not (NA).
- `prefix.log`: a log file indicating windows for which some sequences were removed (if they contained only missing data), and windows for which tree inference failed.

In addition, if `write.seq=T`, a folder `prefix_sequences` will be created, containing a fasta sequence for each window.

*Exemples:*

`topo.windows.sites(vcf="/path/to/file.vcf", size=500, incr=0, phased=T, prefix="500_SNPs_no_overlap", write.seq=F, nj=T, dna.dist="JC69")` will calculate trees in 500 SNPs non-overlapping sliding windows based on a phased vcf. NJ trees will be calculated using Jukes-Cantor 69 distances, and fasta sequences will not be written.

`topo.windows.coord(vcf="/path/to/file.vcf.gz", size=50000, incr=10000, phased=F, prefix="50000kb_10kb_overlap", write.seq=T, nj=T, dna.dist="raw")` will calculate trees in 50kb sliding windows with a 10kb overlap, based on an unphased vcf. NJ trees will be calculated using raw distances, and fasta sequences will be written.

### Trees in user-specified regions

Use `tree.region(vcf, regions, phased, write.seq, nj, prefix)` to calculate trees in one or more pre-defined regions.

*Arguments:*

Same as before, except for:
- `regions`: path to a bed-like table with three columns defining the windows: chromosome, start and end. The table shouldn't have a header.

*Output:*

For each region, a neighbor-joining tree (if `nj=T`) and/or a fasta sequence (if `write.seq=T`).

*Exemple:*

`tree.region(vcf="/path/to/file.vcf", regions="/path/to/regions.txt", phased=T, write.seq=T, nj=T, prefix="regions")`

## Running in an unix command line environment

It may be useful to execute TopoWindows directly from a unix CL environment, especially when working on a HPC. To help with that I provide a wrapper script which can be executed as follow: `Rscript Topo_windows_v02_cl_wrapper.R --vcf --prefix  --phased --nj --ali --dist [--regions | --type --size --incr]`. The combination of arguments used will determine which of the three functions is run. Before running the script, one must edit the first line to provide the full path to `Topo_windows_v02.R`.

*Arguments:*

- `--regions`: path to a regions table for `tree.region`. If this argument is specified, `tree.region` is run, and `--type`, `--size` and `--incr` are not needed. (default: NULL)
- `--type`: the type of sliding windows to use (either `c` for coordinates or `s` for sites). If `--type c` is specified `topo.windows.coord` is run, if `--type s` is specified `topo.windows.sites` is run.

