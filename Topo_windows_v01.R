##### A script to calculate trees in sliding windows from a vcf file. Called from Topo_windows.sh

##### L. Rancilhac
# v. 0.1 - 15/09/2023

library(vcfR)
library(ape)

topo.windows.sites <- function(vcf, size, phased, PREF, write.seq = T, nj = T){

  vcf <- read.vcfR(vcf)
  if(phased == T){ dnabin <- vcfR2DNAbin(vcf) }
  else if(phased == F){ dnabin <- vcfR2DNAbin(vcf, extract.haps = F, consensus = T) }
  

  stats <- matrix(ncol=7)
  colnames(stats) <- c("CHR","CHR.START", "CHR.END", "CHR.SIZE", "NSITES", "PROP.MISS", "TREE")
  START <- 1

  seq.dir <- paste(PREF, "_sequences", sep = "")
  print(seq.dir)
  if(write.seq == T & !file.exists(seq.dir)){ dir.create(seq.dir) }

  nj.file.name <- paste(PREF, "_NJ_trees.trees", sep="")

  while (START < length(dnabin[1,])) {

    # define the end of the window
    END <- START + (size - 1)
    if(END > length(dnabin[1,])) { END <- length(dnabin[1,]) }

    # extract the data for the defined window
    curr.window <- dnabin[,START:END]

    # write sequence in fasta format if specified
    if(write.seq == T){
	      NAME <- paste(seq.dir, "/", colnames(curr.window)[1],"-",colnames(curr.window)[length(curr.window[1,])],"-",length(curr.window[1,]),"pos.fasta",sep="")
      write.FASTA(curr.window, NAME)
    }

    # calculate tree
    if(nj == T){
      curr.dist <- dist.dna(curr.window, model="JC69", pairwise.deletion = T)
      curr.tree <- try(nj(curr.dist), silent = T)
      if(class(curr.tree) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
        TREE <- "NA"}
      else{write.tree(curr.tree, nj.file.name, append=T)
        TREE <- "YES"}
    }
    else{ TREE <- "NA" }

    # collect window information
    CHR <- strsplit(colnames(curr.window)[1], "_")[[1]][1]
    CHR.START <- strsplit(colnames(curr.window)[1], "_")[[1]][2]
    CHR.END <- strsplit(colnames(curr.window)[length(curr.window[1,])], "_")[[1]][2]
    WIN.SIZE <- as.numeric(CHR.END) - as.numeric(CHR.START)
    NSITES <- length(curr.window[1,])
    PROP.MISS <- round(length(which(curr.window == "f0"))/length(curr.window), digits=4)

    # write window information
    stats <- rbind(stats, c(CHR, CHR.START, CHR.END, WIN.SIZE, NSITES, PROP.MISS, TREE))

    # increment to next window
    START <- START + size
  }

  # write information table
  STAT.NAME <- paste(PREF,"_",size,"_sites_windows_stats.tsv", sep="")
  stats <- stats[-1,]
  write.table(as.data.frame(stats), STAT.NAME, quote = F, row.names = F, col.names = T, sep="\t")

}

topo.windows.coord <- function(vcf, size, phased, PREF, write.seq = T, nj = T){

  vcf <- read.vcfR(vcf)

  stats <- matrix(ncol=7)
  colnames(stats) <- c("CHR","CHR.START", "CHR.END", "CHR.SIZE", "NSITES", "PROP.MISS", "TREE")

  NWIN <- 1
  START <- 1

  seq.dir <- paste(PREF, "_sequences", sep = "")
  if(write.seq == T & !file.exists(seq.dir)){ dir.create(seq.dir) }

  nj.file.name <- paste(PREF, "_NJ_trees.trees", sep="")

  LAST.POS <- max(as.numeric(vcf@fix[,2]))

  while (START < LAST.POS) {

    paste("processing window ", NWIN, sep="")
    # define the end of the window
    END <- START + (size - 1)

    if(END > LAST.POS) { END <- LAST.POS }

    # extract the data for the defined window

    curr.vcf <- vcf[which(as.numeric(vcf@fix[,2]) >= START & as.numeric(vcf@fix[,2]) <= END), ]
    if(phased == T){ curr.window <- vcfR2DNAbin(curr.vcf, verbose=F) }
    else if(phased == F){ curr.window <- vcfR2DNAbin(curr.vcf, verbose=F, extract.haps = F, consensus = T) }

    if(dim(curr.vcf)[1] == 0){
      CHR <- vcf@fix[1,1]
      CHR.START <- START
      CHR.END <- END
      WIN.SIZE <- as.numeric(CHR.END) - as.numeric(CHR.START)
      NSITES <- length(curr.window[1,])
      PROP.MISS <- round(length(which(curr.window == "f0"))/length(curr.window), digits=4)

      # write window information
      stats <- rbind(stats, c(CHR, CHR.START, CHR.END, WIN.SIZE, NSITES, PROP.MISS, TREE))

      # increment to next window
      START <- START + size
      NWIN <- NWIN+1
    }
    else{

    # write sequence in fasta format if specified
    if(write.seq == T){
      NAME <- paste(seq.dir, "/", colnames(curr.window)[1],"-",colnames(curr.window)[length(curr.window[1,])],"-",size,"pos.fasta",sep="")
      write.FASTA(curr.window, NAME)
    }

    # calculate tree if specified
    if(nj == T){
      curr.dist <- dist.dna(curr.window, model="JC69", pairwise.deletion = T)
      curr.tree <- try(nj(curr.dist), silent = T)
      if(class(curr.tree) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
        TREE <- "NA"}
      else{write.tree(curr.tree, nj.file.name, append=T)
        TREE <- "YES"}
    }
    else{ TREE <- "NA" }

    # collect window information
    CHR <- strsplit(colnames(curr.window)[1], "_")[[1]][1]
    CHR.START <- START
    CHR.END <- END
    WIN.SIZE <- as.numeric(CHR.END) - as.numeric(CHR.START)
    NSITES <- length(curr.window[1,])
    PROP.MISS <- round(length(which(curr.window == "f0"))/length(curr.window), digits=4)

    # write window information
    stats <- rbind(stats, c(CHR, CHR.START, CHR.END, WIN.SIZE, NSITES, PROP.MISS, TREE))

    # increment to next window
    START <- START + size
    NWIN <- NWIN+1
    }
  }

  # write information table
  STAT.NAME <- paste(PREF,"_",size,"_coordinates_windows_stats.tsv", sep="")
  stats <- stats[-1,]
  write.table(as.data.frame(stats), STAT.NAME, quote = F, row.names = F, col.names = T, sep="\t")

}
