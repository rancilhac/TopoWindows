##### A script to calculate trees in sliding windows from a vcf file. Called from Topo_windows.sh

##### L. Rancilhac
# v. 0.2 - 06/11/2023

library(vcfR)
library(ape)

# Main functions

topo.windows.sites <- function(vcf, size, incr=0, phased, prefix, write.seq = T, nj = T, dna.dist){
  if(incr == 0){ incr <- size }
  vcf <- read.vcfR(vcf)
  if(phased == T){ dnabin <- vcfR2DNAbin(vcf) }
  else if(phased == F){ dnabin <- vcfR2DNAbin(vcf, extract.haps = F, consensus = T) }
  
  stats <- matrix(ncol=9)
  colnames(stats) <- c("CHR","CHR.START", "CHR.END", "CHR.SIZE", "NSITES", "PROP.MISS", "PROP.PIS", "TREE", "NTIPS")
  CHR <- getFIX(vcf)[,1][1]
  START <- 1

  seq.dir <- paste(prefix, "_sequences", sep = "")
  if(write.seq == T & !file.exists(seq.dir)){ dir.create(seq.dir) }

  nj.file.name <- paste(prefix, "_NJ_trees.trees", sep="")

  while (START < length(dnabin[1,])) {

    # define the end of the window
    END <- START + (size - 1)
    if(END > length(dnabin[1,])) { END <- length(dnabin[1,]) }

    # extract the data for the defined window
    curr.window <- dnabin[,START:END]
    
    CHR.START <- strsplit(colnames(curr.window)[1], "_")[[1]][length(strsplit(colnames(curr.window)[1], "_")[[1]])]
    CHR.END <- strsplit(colnames(curr.window)[length(curr.window[1,])], "_")[[1]][length(strsplit(colnames(curr.window)[length(curr.window[1,])], "_")[[1]])]

    # write sequence in fasta format if specified
    if(write.seq == T){
	      NAME <- paste(seq.dir, "/", CHR,"-", CHR.START, "-", CHR.END,".fasta",sep="")
      write.FASTA(curr.window, NAME)
    }
    
    # calculate tree
    if(nj == T){
      curr.window <- as.character(curr.window)
      miss.seq <- which(apply(curr.window, 1, missing.seq) == 0)
      if(length(miss.seq) > 0){
      curr.window <- curr.window[-miss.seq, ]
      cat(paste("\nthe following sequences were removed from window ", CHR.START, "-", CHR.END, " because they contained only Ns\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
      cat(paste(rownames(curr.window)[miss.seq], "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
      }
      curr.window <- as.DNAbin(curr.window)
      curr.dist <- dist.dna(curr.window, model=dna.dist, pairwise.deletion = T)
      curr.dist[curr.dist == Inf] <- NA
      cat(paste("\nWarning: ", length(which(is.na(curr.dist)))/length(curr.dist),"% of the distances could not be calculated in window", CHR.START, "-", CHR.END, ". A NJ tree is still calculated. This is likely because some sequences are too divergent. You may want to use the \"raw\" distance model or use larger windows.\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
      curr.tree <- try(njs(curr.dist), silent = T)
      if(class(curr.tree) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
        TREE <- "NA"
        cat(paste("\nthe tree could no be calculated for window ", CHR.START, "-", CHR.END, "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
      else if(class(try(write.tree(curr.tree), silent=T)) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
        TREE <- "NA"
        cat(paste("\nthe tree could no be written for window ", CHR.START, "-", CHR.END, "; this may be because some sequences have too much missing data or are too distant from each other\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
      else{write.tree(curr.tree, nj.file.name, append=T)
        TREE <- "YES"
	NTIPS <- Ntip(curr.tree)}
    }
    else if(nj == F){ TREE <- NA }

    # collect window information
    WIN.SIZE <- as.numeric(CHR.END) - as.numeric(CHR.START)
    NSITES <- length(curr.window[1,])
    PROP.PIS <- prop.pis(as.character(curr.window))
    PROP.MISS <- round(length(which(curr.window == "f0"))/length(curr.window), digits=4)

    # write window information
    stats <- rbind(stats, c(CHR, CHR.START, CHR.END, WIN.SIZE, NSITES, PROP.MISS, PROP.PIS, TREE, NTIPS))

    # increment to next window
    print(START)
    START <- START + incr
  }

  # write information table
  STAT.NAME <- paste(prefix,"_windows_stats.tsv", sep="")
  stats <- stats[-1,]
  write.table(as.data.frame(stats), STAT.NAME, quote = F, row.names = F, col.names = T, sep="\t")

}

topo.windows.coord <- function(vcf, size, incr=0, phased, prefix, write.seq = T, nj = T, dna.dist){
  
  if(incr == 0){ incr <- size }
  
  vcf <- read.vcfR(vcf)

  stats <- matrix(ncol=9)
  colnames(stats) <- c("CHR","CHR.START", "CHR.END", "CHR.SIZE", "NSITES", "PROP.MISS", "PROP.PIS", "TREE", "NTIPS")

  NWIN <- 1
  START <- 1

  seq.dir <- paste(prefix, "_sequences", sep = "")
  if(write.seq == T & !file.exists(seq.dir)){ dir.create(seq.dir) }

  nj.file.name <- paste(prefix, "_NJ_trees.trees", sep="")

  LAST.POS <- max(as.numeric(vcf@fix[,2]))

  while (START < LAST.POS) {

    paste("processing window ", NWIN, sep="")
    # define the end of the window
    END <- START + (size - 1)
    
    print(paste(START, "-", END, sep=""))

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
      PROP.PIS <- 0
      PROP.MISS <- round(length(which(curr.window == "f0"))/length(curr.window), digits=4)

      # write window information
      stats <- rbind(stats, c(CHR, CHR.START, CHR.END, WIN.SIZE, NSITES, PROP.MISS, PROP.PIS, TREE))

      # increment to next window
      START <- START + size
      NWIN <- NWIN+1
    }
    else{
    
    CHR <- strsplit(colnames(curr.window)[1], "_")[[1]][1]
    CHR.START <- START
    CHR.END <- END
    # write sequence in fasta format if specified
    if(write.seq == T){
      NAME <- paste(seq.dir, "/", CHR,"-",CHR.START, "-", CHR.END,".fasta",sep="")
      write.FASTA(curr.window, NAME)
    }

    # calculate tree if specified
    if(nj == T){
      curr.window <- as.character(curr.window)
      miss.seq <- which(apply(curr.window, 1, missing.seq) == 0)
      if(length(miss.seq) > 0){
        curr.window <- curr.window[-miss.seq, ]
        cat(paste("\nthe following sequences were removed from window ", CHR.START, "-", CHR.END, " because they contained only Ns\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
        cat(paste(rownames(curr.window)[miss.seq], "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
      }
      curr.window <- as.DNAbin(curr.window)
      curr.dist <- dist.dna(curr.window, model=dna.dist, pairwise.deletion = T)
      curr.dist[curr.dist == Inf] <- NA
      cat(paste("\nWarning: ", length(which(is.na(curr.dist)))/length(curr.dist),"% of the distances could not be calculated in window", CHR.START, "-", CHR.END, ". A NJ tree is still calculated. This is likely because some sequences are too divergent. You may want to use the \"raw\" distance model or use larger windows.\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
      curr.tree <- try(njs(curr.dist), silent = T)
      if(class(curr.tree) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
        TREE <- "NA"
        cat(paste("\nthe tree could no be calculated for window ", CHR.START, "-", CHR.END, "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
      else if(class(try(write.tree(curr.tree), silent=T)) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
        TREE <- "NA"
        cat(paste("\nthe tree could no be written for window ", CHR.START, "-", CHR.END, "; this may be because some sequences have too much missing data or are too distant from each other\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
      else{write.tree(curr.tree, nj.file.name, append=T)
        TREE <- "YES"
	NTIPS <- Ntip(curr.tree)}
    }
    else if(nj == F){ TREE <- NA }

    # collect window information
    WIN.SIZE <- as.numeric(CHR.END) - as.numeric(CHR.START)
    NSITES <- length(curr.window[1,])
    PROP.PIS <- prop.pis(as.character(curr.window))
    PROP.MISS <- round(length(which(curr.window == "f0"))/length(curr.window), digits=4)

    # write window information
    stats <- rbind(stats, c(CHR, CHR.START, CHR.END, WIN.SIZE, NSITES, PROP.MISS, PROP.PIS, TREE, NTIPS))

    # increment to next window
    START <- START + incr
    NWIN <- NWIN+1
    }
  }

  # write information table
  STAT.NAME <- paste(prefix,"_",size,"_coordinates_windows_stats.tsv", sep="")
  stats <- stats[-1,]
  write.table(as.data.frame(stats), STAT.NAME, quote = F, row.names = F, col.names = T, sep="\t")

}

tree.region <- function(vcf, regions, phased, write.seq = T, nj = T, prefix, dna.dist){
  
  vcf <- read.vcfR(vcf)
  
  regions <- read.table(regions)
  
  for(i in 1:nrow(regions)){
  chr <- regions[i,1]
  curr.vcf <- vcf[vcf@fix[,1] == chr, ]
  
  start <- regions[i,2]
  end <- regions[i,3]
  
  nj.file.name <- paste(chr, "_", start, "-", end, "_NJ_tree.tre", sep="")
  
 if(start > max(as.numeric(curr.vcf@fix[,2]))) { print("error: starting position after the end of vcf") }
 if(end > max(as.numeric(curr.vcf@fix[,2]))) { end <- max(as.numeric(curr.vcf@fix[,2])) }
 curr.vcf <- curr.vcf[which(as.numeric(curr.vcf@fix[,2]) >= start & as.numeric(curr.vcf@fix[,2]) <= end), ]
 if(phased == T){ curr.window <- vcfR2DNAbin(curr.vcf, verbose=F, asterisk_as_del = F) }
 else if(phased == F){ curr.window <- vcfR2DNAbin(curr.vcf, verbose=F, extract.haps = F, consensus = T, asterisk_as_del = F) }
    
 if(dim(curr.vcf)[1] == 0){ print("error: vcf is empty for selected region") }
  else{
      # write sequence in fasta format if specified
      if(write.seq == T){
        NAME <- paste("./", prefix, chr,"-", start, "-", end, ".fasta",sep="")
        write.FASTA(curr.window, NAME)
      }
      
      # calculate tree if specified
      if(nj == T){
        curr.window <- as.character(curr.window)
        miss.seq <- which(apply(curr.window, 1, missing.seq) == 0)
        if(length(miss.seq) > 0){
          curr.window <- curr.window[-miss.seq, ]
          cat(paste("\nthe following sequences were removed from window ", start, "-", end, " because they contained only Ns\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
          cat(paste(rownames(curr.window)[miss.seq], "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
        }
        curr.window <- as.DNAbin(curr.window)
        curr.dist <- dist.dna(curr.window, model=dna.dist, pairwise.deletion = T)
        curr.dist[curr.dist == Inf] <- NA
	if(length(which(is.na(curr.dist)))/length(curr.dist)){
	cat(paste("\nWarning: ", length(which(is.na(curr.dist)))/length(curr.dist),"% of the distances could not be calculated in window", CHR.START, "-", CHR.END, ". A NJ tree is still calculated. This is likely because some sequences are too divergent. You may want to use the \"raw\" distance model or use larger windows.\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
        }
	curr.tree <- try(nj(curr.dist), silent = T)
        if(class(curr.tree) == "try-error"){ cat(paste("\nthe tree could no be calculated for region ", chr, ":", start, "-", end, "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
        else if(class(try(write.tree(curr.tree), silent=T)) == "try-error"){ cat(paste("\nthe tree could no be written for window ", chr, ":", start, "-", end, "; this may be because some sequences have too much missing data or are too distant from each other\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
        else{write.tree(curr.tree, nj.file.name, append=T)}
      }
  }
  }
  }

# Helper functions

prop.pis <- function(x){
  pars.inf <- function(x) {
    x <- table(x)
    x <- x[x > 1]
    n <- c("-", "n", "b", "h", "d", 
           "v", "k", "s", "r", "w", 
           "y")
    if (length(x[!names(x) %in% n]) > 1) 
      x <- TRUE
    else FALSE
  }
  nbchar <- dim(x)[2]
  out <- apply(x, 2, pars.inf)
  out <- round(length(out[out])/nbchar, digits=3)
  return(out)
  
}

missing.seq <- function(x){
  length(which(names(table(x)) != "n"))
  }
