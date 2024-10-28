source("/path/to/Topo_windows_v02.R")

library(optparse)

option_list <- list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="prefix for output file", metavar="character"),
  make_option(c("-v", "--vcf"), type="character", default=NULL, 
              help="path of input vcf file", metavar="character"),
  make_option(c("-w", "--window-type"), type="character", default=NULL, 
              help="type of windows to use (coordinates or sites)", metavar="character"),
  make_option(c("-s", "--size"), type="numeric", default=NULL, 
              help="size of windows", metavar="numeric"),
  make_option(c("-d", "--phased"), type="logical", default=NULL, 
              help="whether the vcf file is phased", metavar="logical"),
  make_option(c("-t", "--tree"), type="character", default=NULL, 
              help="whether to calculate neighbor-joining (NJ), Maximum-Likelihood (ML) or no (N) trees", metavar="character"),
  make_option(c("-a", "--ali"), type="logical", default=NULL, 
              help="whether to output alignments in fasta format for each window", metavar="logical"),
  make_option(c("-i", "--incr"), type="numeric", default=0, 
              help="overlap between windows", metavar="numeric"),
  make_option(c("-r", "--regions"), type="character", default=NULL, 
              help="A bed file specifying regions in which to calculate trees", metavar="character"),
  make_option(c("-m", "--dna-model"), type="character", default=NULL, 
              help="Substitution model used to calculate the phylogenies", metavar="character"),
  make_option(c("-f", "--force"), type="logical", default=NULL, 
              help="Whether to over-write pre-existing output files", metavar="logical"),
  make_option(c("-n", "--missingness"), type="numeric", default=NULL, 
              help="Maximum proportion of missing data allowed to keep a sequence in a given window", metavar="numeric")
)

opt_parser <- OptionParser(option_list = option_list)

opt <- parse_args(opt_parser)

PREF <- opt$prefix
VCF <- opt$vcf
TYPE <- opt$window-type
SIZE <- opt$size
PHASED <- opt$phased
TREE <- opt$tree
ALN <- opt$ali
INCR <- opt$incr
REG <- opt$regions
MODEL <- opt$dna-model
FORCE <- opt$force
MISS <- opt$missingness

if(!is.null(REG)){
  tree.region(vcf=VCF, regions=REG, phased=PHASED, write.seq=ALN, tree=TREE, dna.model=MODEL, missing.thresh=MISS, force=FORCE)
} else if( TYPE == "s"){
  topo.windows.sites(vcf=VCF, size=SIZE, incr=INCR, prefix=PREF, phased=PHASED, write.seq=ALN, tree=TREE, dna.model=MODEL, missing.thresh=MISS, force=FORCE)
} else if( TYPE == "c" ){
  topo.windows.coord(vcf=VCF, size=SIZE, incr=INCR, prefix=PREF, phased=PHASED, write.seq=ALN, tree=TREE, dna.model=MODEL, missing.thresh=MISS, force=FORCE)
}
