library(optparse)

option_list <- list(
  make_option(c("-p", "--prefix"), type="character", default=NULL, 
              help="prefix for output file", metavar="character"),
  make_option(c("-v", "--vcf"), type="character", default=NULL, 
              help="path of input vcf file", metavar="character"),
  make_option(c("-t", "--type"), type="character", default=NULL, 
              help="type of windows to use (coordinates or sites)", metavar="character"),
  make_option(c("-s", "--size"), type="numeric", default=NULL, 
              help="size of windows", metavar="numeric"),
  make_option(c("-d", "--phased"), type="logical", default=NULL, 
              help="whether the vcf file is phased", metavar="logical"),
  make_option(c("-n", "--nj"), type="logical", default=NULL, 
              help="whether to calculate neighbor-joining trees", metavar="logical"),
  make_option(c("-a", "--ali"), type="logical", default=NULL, 
              help="whether to output alignments in fasta format for each window", metavar="logical"),
  make_option(c("-i", "--incr"), type="numeric", default=0, 
              help="overlap between windows", metavar="numeric"),
  make_option(c("-r", "--regions"), type="character", default=NULL, 
              help="A bed file specifying regions in which to calculate trees", metavar="character")
)

opt_parser <- OptionParser(option_list = option_list)

opt <- parse_args(opt_parser)

PREF <- opt$prefix
VCF <- opt$vcf
TYPE <- opt$type
SIZE <- opt$size
PHASED <- opt$phased
NJ <- opt$nj
ALN <- opt$ali
INCR <- opt$incr
REG <- opt$regions

if(!is.null(REG)){
  tree.region(vcf=VCF, regions=REG, phased=PHASED, write.seq=ALN, nj=NJ, prefix=PREF)
} else if( TYPE == "s"){
  topo.windows.sites(vcf=VCF, size=SIZE, incr=INCR, prefix=PREF, phased=PHASED, write.seq=ALN, nj=NJ)
} else if( TYPE == "c" ){
  topo.windows.coord(vcf=VCF, size=SIZE, incr=INCR, prefix=PREF, phased=PHASED, write.seq=ALN, nj=NJ)
}
