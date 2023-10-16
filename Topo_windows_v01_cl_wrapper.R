getwd()

arguments <- commandArgs(trailingOnly = TRUE)

PREF <- arguments[1]
print(PREF)
VCF <- arguments[2]
print(VCF)
TYPE <- arguments[3]
print(TYPE)
SIZE <- as.numeric(arguments[4])
print(SIZE)
INCR <- as.numeric(arguments[5])
print(INCR)
PHASED <- as.logical(arguments[6])
print(PHASED)
NJ <- as.logical(arguments[7])
print(NJ)
ALN <- as.logical(arguments[8])
print(ALN)


if( TYPE == "s"){
  topo.windows.sites(vcf=VCF, size=SIZE, incr=INCR, phased=PHASED, write.seq = ALN, nj = NJ)
} else if( TYPE == "c" ){
  topo.windows.coord(VCF, SIZE, PHASED, write.seq = ALN, nj = NJ)
}
