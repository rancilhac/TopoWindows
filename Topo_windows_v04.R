##### Functions to calculate trees in sliding windows from a vcf file.

##### L. Rancilhac
# v. 0.4 - 10/2024

library(vcfR)
library(ape)
library(phangorn)

# Main functions

### addition of ML phylogeny (package phangorn)
# tree: three options: N (none), NJ or ML

topo.windows.sites <- function(vcf, size, incr=0, phased, prefix, write.seq = T, tree = "NJ", dna.model="JC",
                               missing.thresh=0.7, force = T){
  
  nj.file.name <- paste(prefix, ".trees", sep="")
  if(file.exists(nj.file.name) & force == F){ stop("Execution halted: output file already exists. Please delete existing file or use force=T")}
  else if (file.exists(nj.file.name) & force == T){ file.remove(nj.file.name) }
  
  if(!tree %in% c("N", "NJ", "ML")){ stop(paste("Execution halted: argument tree = ", tree, " not recognized, please use one of N, NJ and ML", sep="")) }

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
    if(tree != "N"){
      curr.window <- as.character(curr.window)
      rem.seq <- which(apply(curr.window, 1, prop.missing) > missing.thresh)
      if(length(rem.seq) > 0){
        curr.window <- curr.window[-rem.seq, ]
        cat(paste("\nthe following sequences were removed from window ", CHR.START, "-", CHR.END, " because they contained more than ", missing.thresh,"% of missing sites\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
        cat(paste(names(rem.seq), "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
      }
      curr.window <- as.DNAbin(curr.window)
      curr.dist <- dist.dna(curr.window, model=dna.model, pairwise.deletion = T)
      curr.dist[curr.dist == Inf] <- NA
      if(tree == "ML"){
        curr.tree <- try(ml.phylo(curr.window, curr.dist, dna.model), silent = T)
        if(class(curr.tree$tree) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
          TREE <- "NA"
          cat(paste("\nThe tree could no be calculated for window ", CHR.START, "-", CHR.END, "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
        else if(class(try(write.tree(curr.tree$tree), silent=T)) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
          TREE <- "NA"
          cat(paste("\nThe tree could no be written for window ", CHR.START, "-", CHR.END, "; this may be because some sequences have too much missing data or are too distant from each other\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
        else{write.tree(curr.tree$tree, nj.file.name, append=T)
          TREE <- "YES"
          NTIPS <- Ntip(curr.tree$tree)}
      }
      else if(tree == "NJ"){
        curr.tree <- try(bionjs(curr.dist), silent = T)
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
      }
    else { TREE <- NA }
    
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

topo.windows.coord <- function(vcf, size, incr=0, phased, prefix, write.seq = T, tree = "NJ", dna.model="JC",
                               missing.thresh=0.7, force=T){
  nj.file.name <- paste(prefix, ".trees", sep="")
  if(file.exists(nj.file.name) & force == F){ stop("Execution halted: output file already exists. Please delete existing file or use force=T")}
  else if (file.exists(nj.file.name) & force == T){ file.remove(nj.file.name) }
  
  if(!tree %in% c("N", "NJ", "ML")){ stop(paste("Execution halted: argument tree = ", tree, " not recognized, please use one of N, NJ and ML", sep="")) }
  
  if(incr == 0){ incr <- size }
  
  vcf <- read.vcfR(vcf)
  
  stats <- matrix(ncol=9)
  colnames(stats) <- c("CHR","CHR.START", "CHR.END", "CHR.SIZE", "NSITES", "PROP.MISS", "PROP.PIS", "TREE", "NTIPS")
  
  NWIN <- 1
  START <- 1
  
  seq.dir <- paste(prefix, "_sequences", sep = "")
  if(write.seq == T & !file.exists(seq.dir)){ dir.create(seq.dir) }
  
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
      if(tree != "N"){
        curr.window <- as.character(curr.window)
        rem.seq <- which(apply(curr.window, 1, prop.missing) > missing.thresh)
        if(length(rem.seq) > 0){
          curr.window <- curr.window[-rem.seq, ]
          cat(paste("\nthe following sequences were removed from window ", CHR.START, "-", CHR.END, " because they contained more than ", missing.thresh,"% of missing sites\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
          cat(paste(names(rem.seq), "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
        }
        curr.window <- as.DNAbin(curr.window)
        curr.dist <- dist.dna(curr.window, model=dna.model, pairwise.deletion = T)
        curr.dist[curr.dist == Inf] <- NA
        if(tree == "ML"){
          curr.tree <- try(ml.phylo(curr.window, curr.dist, dna.model), silent = T)
          if(class(curr.tree$tree) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
            TREE <- "NA"
            cat(paste("\nThe tree could no be calculated for window ", CHR.START, "-", CHR.END, "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
          else if(class(try(write.tree(curr.tree$tree), silent=T)) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
            TREE <- "NA"
            cat(paste("\nThe tree could no be written for window ", CHR.START, "-", CHR.END, "; this may be because some sequences have too much missing data or are too distant from each other\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
          else{write.tree(curr.tree$tree, nj.file.name, append=T)
            TREE <- "YES"
            NTIPS <- Ntip(curr.tree$tree)}
        }
        else if(tree == "NJ"){
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
      }
      else { TREE <- NA }
      
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

tree.region <- function(vcf, regions, phased, write.seq = T, tree = "NJ", dna.model = "JC",
                        missing.thresh=0.7, force = T){
  
  if(!tree %in% c("N", "NJ", "ML")){ stop(paste("Execution halted: argument tree = ", tree, " not recognized, please use one of N, NJ and ML", sep="")) }
  
  vcf <- read.vcfR(vcf)
  
  regions <- read.table(regions)
  
  for(i in 1:nrow(regions)){
    chr <- regions[i,1]
    curr.vcf <- vcf[vcf@fix[,1] == chr, ]
    
    start <- regions[i,2]
    end <- regions[i,3]
    
    nj.file.name <- paste(chr, "_", start, "-", end, ".tre", sep="")
    if(file.exists(nj.file.name) & force == F){ stop("Execution halted: output file already exists. Please delete existing file or use force=T")}
    else if (file.exists(nj.file.name) & force == T){ file.remove(nj.file.name) }
    
    if(start > max(as.numeric(curr.vcf@fix[,2]))) { stop("Execution halted: specified starting position is after the end of the vcf") }
    if(end > max(as.numeric(curr.vcf@fix[,2]))) { end <- max(as.numeric(curr.vcf@fix[,2])) }
    curr.vcf <- curr.vcf[which(as.numeric(curr.vcf@fix[,2]) >= start & as.numeric(curr.vcf@fix[,2]) <= end), ]
    if(phased == T){ curr.window <- vcfR2DNAbin(curr.vcf, verbose=F, asterisk_as_del = F) }
    else if(phased == F){ curr.window <- vcfR2DNAbin(curr.vcf, verbose=F, extract.haps = F, consensus = T, asterisk_as_del = F) }
    
    if(dim(curr.vcf)[1] == 0){ stop(paste("Execution halted: vcf is empty for chromosome ", chr, sep="")) }
    else{
      # write sequence in fasta format if specified
      if(write.seq == T){
        NAME <- paste("./", chr,"-", start, "-", end, ".fasta",sep="")
        write.FASTA(curr.window, NAME)
      }
      
      # calculate tree if specified
      if(tree != "N"){
        curr.window <- as.character(curr.window)
        rem.seq <- which(apply(curr.window, 1, prop.missing) > missing.thresh)
        if(length(rem.seq) > 0){
          curr.window <- curr.window[-rem.seq, ]
          cat(paste("\nthe following sequences were removed from window ", CHR.START, "-", CHR.END, " because they contained more than ", missing.thresh*100,"% of missing sites\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
          cat(paste(names(rem.seq), "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)
        }
        curr.window <- as.DNAbin(curr.window)
        curr.dist <- dist.dna(curr.window, model=dna.model, pairwise.deletion = T)
        curr.dist[curr.dist == Inf] <- NA
        if(tree == "ML"){
          curr.tree <- try(ml.phylo(curr.window, curr.dist, dna.model), silent = T)
          if(class(curr.tree$tree) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
            TREE <- "NA"
            cat(paste("\nThe tree could no be calculated for window ", CHR.START, "-", CHR.END, "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
          else if(class(try(write.tree(curr.tree$tree), silent=T)) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
            TREE <- "NA"
            cat(paste("\nThe tree could no be written for window ",chr,":", CHR.START, "-", CHR.END, "; this may be because some sequences have too much missing data or are too distant from each other\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
          else{write.tree(curr.tree$tree, nj.file.name, append=T)
            TREE <- "YES"
            NTIPS <- Ntip(curr.tree$tree)}
        }
        else if(tree == "NJ"){
          curr.tree <- try(bionjs(curr.dist), silent = T)
          if(class(curr.tree) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
            TREE <- "NA"
            cat(paste("\nThe tree could no be calculated for window ", chr, ":", CHR.START, "-", CHR.END, "\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
          else if(class(try(write.tree(curr.tree), silent=T)) == "try-error"){ cat("NA\n", file=nj.file.name, append=T)
            TREE <- "NA"
            cat(paste("\nThe tree could no be written for window ", chr, ":", CHR.START, "-", CHR.END, "; this may be because some sequences have too much missing data or are too distant from each other\n", sep=""), file=paste(prefix,".log",sep=""), append=T)}
          else{write.tree(curr.tree, nj.file.name, append=T)
            TREE <- "YES"
            NTIPS <- Ntip(curr.tree)}
        }
      }
      else { TREE <- NA }
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

prop.missing <- function(seq){
  return(length(which(seq == "n"))/length(seq))
}

ml.phylo <- function(dnabin, dist, model){
  fit_initial <- pml(njs(dist), data = as.phyDat(dnabin))
  fit_ml <- optim.pml(fit_initial, model = model)
}
