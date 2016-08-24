###################### BRASS circos plots ###########################
########################### Version 1.0 #############################

## Maximilian Stammnitz, mrs72
## Last Update - 24/08/2016 ##
##############################

library(stringr)
library(circlize)
library(GenomicRanges)


## 1. Load SVs (as .rearr file)
###############################

## Set of multi-output: 4 Tumours, 4 Hosts (already removed read-IDs with awk prior to loading, saving memory)
SVs <- read.table("All_samples_multi_small.rearr",sep = '\t', header=F)
SVs <- SVs[,c(1,3,4,5,7,8,2,6,9:16)] # might have to re-sort tumour and host columns (as specified in .rearr-header)

colnames(SVs) <- c("CHR-1","START-1","END-1","CHR-2","START-2","END-2","STRAND-1","STRAND-2", 
                   "T1","T2","T3","T4","H1","H2","H3","H4")


## 2. Annotate Genes
####################

annotate.genes <- function(x){
  
  # 1. Load file with gene annotations
  genes <- read.table('Genes.txt', header=T, sep="\t")
  
  # 2. Check First breakpoint
  BP1.Ranges <- GRanges(seqnames = Rle(x[,1]),
                        ranges = IRanges(start = as.numeric(x[,2]),
                                         end = as.numeric(x[,3])))
  Genes.Ranges <- GRanges(seqnames = Rle(genes[,1]),
                          ranges = IRanges(start = genes[,2],
                                           end = genes[,3]))
  OL <- findOverlaps(BP1.Ranges, Genes.Ranges)
  OL <- as.matrix(OL)
  colnames(OL) <- c("Coordinates", "Genes")
  OL[,2] <- paste(as.character(genes[OL[,2],'GENE']))
  
  # 3. Add hits to end
  add <- rep('',nrow(x))
  add[as.numeric(OL[,1])] <- as.character(OL[,2])
  x <- cbind(x, add)
  colnames(x)[ncol(x)] <- 'BP1'
  
  # 4. Check First breakpoint
  BP2.Ranges <- GRanges(seqnames = Rle(x[,4]),
                        ranges = IRanges(start = as.numeric(x[,5]),
                                         end = as.numeric(x[,6])))
  OL <- findOverlaps(BP2.Ranges, Genes.Ranges)
  OL <- as.matrix(OL)
  colnames(OL) <- c("Coordinates", "Genes")
  OL[,2] <- paste(as.character(genes[OL[,2],'GENE']))
  
  # 5. Add hits to end
  add <- rep('',nrow(x))
  add[as.numeric(OL[,1])] <- as.character(OL[,2])
  x <- cbind(x, add)
  colnames(x)[ncol(x)] <- 'BP2'
  
  # 6. Output
  return(x)
}
SVs <- annotate.genes(SVs)


## 3. Raw Filter
################

# 1. Ends of contigs
contigends <- function(x){
  
  devil.contig.gaps <- read.table("Contigends.txt", header=T)

  # 1. Left Breakpoint
  Input.left.Ranges <- GRanges(seqnames = Rle(x[,1]),
                               ranges = IRanges(start = as.integer(x[,2]),
                                                end = as.integer(x[,3])))
  
  Contigend.Ranges <- GRanges(seqnames = Rle(devil.contig.gaps[,1]),
                              ranges = IRanges(start = as.integer(devil.contig.gaps[,2]),
                                               end = as.integer(devil.contig.gaps[,3])))
  
  # 2. Match with positions in Input
  Overlaps.left <- findOverlaps(Input.left.Ranges, Contigend.Ranges)
  Overlaps.left <- as.matrix(Overlaps.left)
  colnames(Overlaps.left) <- c("Sets", "Contigends")
  
  # 3. Right Breakpoint
  Input.right.Ranges <- GRanges(seqnames = Rle(x[,4]),
                                ranges = IRanges(start = as.integer(x[,5]),
                                                 end = as.integer(x[,6])))
  
  # 4. Match with positions in Input
  Overlaps.right <- findOverlaps(Input.right.Ranges, Contigend.Ranges)
  Overlaps.right <- as.matrix(Overlaps.right)
  colnames(Overlaps.right) <- c("Sets", "Contigends")
  
  # 5. Remove the ones which match
  cat("\n Total removed: ", round(length(unique(c(Overlaps.left[,1],Overlaps.right[,1])))/nrow(x),4)*100, "%")
  x <- x[-unique(c(Overlaps.left[,1],Overlaps.right[,1])),]
  
  # 6. Output
  return(x)
}
SVs.filt <- contigends(SVs)

# 2. Simple repeats
simple.repeats <- function(x){
  
  devil.repeats <- read.table("Simplerepeats.txt",header=T)
  require(GenomicRanges)
  
  # 1. Left Breakpoint
  Input.left.Ranges <- GRanges(seqnames = Rle(x[,1]),
                               ranges = IRanges(start = as.integer(x[,2]),
                                                end = as.integer(x[,3])))
  
  Repeat.Ranges <- GRanges(seqnames = Rle(devil.repeats[,1]),
                           ranges = IRanges(start = as.integer(devil.repeats[,2]-5),
                                            end = as.integer(devil.repeats[,3])+5))
  
  # 2. Match with positions in Input
  Overlaps.left <- findOverlaps(Input.left.Ranges, Repeat.Ranges)
  Overlaps.left <- as.matrix(Overlaps.left)
  colnames(Overlaps.left) <- c("Sets", "Repeats")
  
  # 3. Right Breakpoint
  Input.right.Ranges <- GRanges(seqnames = Rle(x[,4]),
                                ranges = IRanges(start = as.integer(x[,5]),
                                                 end = as.integer(x[,6])))
  
  # 4. Match with positions in Input
  Overlaps.right <- findOverlaps(Input.right.Ranges, Repeat.Ranges)
  Overlaps.right <- as.matrix(Overlaps.right)
  colnames(Overlaps.right) <- c("Sets", "Repeats")
  
  # 5. Remove the ones which match
  cat("\n Total removed: ", round(length(unique(c(Overlaps.left[,1],Overlaps.right[,1])))/nrow(x),4)*100, "%")
  x <- x[-unique(c(Overlaps.left[,1],Overlaps.right[,1])),]
  
  # 6. Output
  return(x)
}
SVs.filt <- simple.repeats(SVs.filt)

# 3. Filter out breakpoints that self-overlap
self.overlap <- function(x){
  
  # 1. Make Granges objects of left and right
  require(GenomicRanges)
  Input.left.Ranges <- GRanges(seqnames = Rle(x[,1]),
                               ranges = IRanges(start = as.integer(x[,2]),
                                                end = as.integer(x[,3])))
  
  Input.right.Ranges <- GRanges(seqnames = Rle(x[,4]),
                                ranges = IRanges(start = as.integer(x[,5]),
                                                 end = as.integer(x[,6])))
  
  # 2. Pairwise overlap-testing
  x <- x[Input.left.Ranges %outside% Input.right.Ranges,]
  
  # 3. Output
  return(x)
}
SVs.filt <- self.overlap(SVs.filt)


## 4. Format
############

# 1. ChrX renaming
for (i in c(1,4)){
  SVs.filt[,i] <- sub('Chrx', 'ChrX', SVs.filt[,i])
}

# 2. Translate coordinates
contig.to.genomic <- function(x){
  
  # 1. Load contig-genomic position file
  Contiglengthfile <- "ContigPositions.txt"
  contig.pos <- read.table(Contiglengthfile, header=T)
  contig.pos[,1] <- sub('Chrx','ChrX',contig.pos[,1])
  contig.pos <- contig.pos[-grep('ChrU', contig.pos[,1]),]
  
  # 2. BREAKPOINT 1
  if(length(grep("_supercontig_000000000|MT", as.character(x[,1])))!=0){
    contig <- as.character(x[-grep("_supercontig_000000000|MT", as.character(x[,1])),1])
    add <- as.numeric(contig.pos[match(contig,contig.pos[,1])-1,3])
    x[-grep("_supercontig_000000000|MT", as.character(x[,1])),2] <- 
      as.numeric(x[-grep("_supercontig_000000000|MT", as.character(x[,1])),2])+add
    x[-grep("_supercontig_000000000|MT", as.character(x[,1])),3] <- 
      as.numeric(x[-grep("_supercontig_000000000|MT", as.character(x[,1])),3])+add
  }else{
    contig <- as.character(x[,1])
    add <- as.numeric(contig.pos[match(contig,contig.pos[,1])-1,3])
    x[,2] <- as.numeric(x[,2])+add
    x[,3] <- as.numeric(x[,3])+add
  }
  
  # 3. BREAKPOINT 2
  if(length(grep("_supercontig_000000000|MT", as.character(x[,4])))!=0){
    contig <- as.character(x[-grep("_supercontig_000000000|MT", as.character(x[,4])),4])
    add <- as.numeric(contig.pos[match(contig,contig.pos[,1])-1,3])
    x[-grep("_supercontig_000000000|MT", as.character(x[,4])),5] <- 
      as.numeric(x[-grep("_supercontig_000000000|MT", as.character(x[,4])),5])+add
    x[-grep("_supercontig_000000000|MT", as.character(x[,4])),6] <- 
      as.numeric(x[-grep("_supercontig_000000000|MT", as.character(x[,4])),6])+add
  }else{
    contig <- as.character(x[,4])
    add <- as.numeric(contig.pos[match(contig,contig.pos[,1])-1,3])
    x[,5] <- as.numeric(x[,5])+add
    x[,6] <- as.numeric(x[,6])+add
  }
  
  # 4. Chromosome renaming
  if(length(grep('MT', x[,1]))!=0){
    x[-grep("MT", as.character(x[,1])),1] <- str_split_fixed(x[-grep("MT", as.character(x[,1])),1],'_',2)[,1]
  }else{
    x[,1] <- str_split_fixed(x[,1],'_',2)[,1]
  }
  if(length(grep('MT', x[,4]))!=0){
    x[-grep("MT", as.character(x[,4])),4] <- str_split_fixed(x[-grep("MT", as.character(x[,4])),4],'_',2)[,1]
  }else{
    x[,4] <- str_split_fixed(x[,4],'_',2)[,1]
  }
  
  # 5. Output
  return(x)
}
SVs.filt <- contig.to.genomic(SVs.filt)


## 5. Subsetting
################

## Germline removal

# a. Less than or 3 reads from all hosts
SVs.filt.nogerm <- SVs.filt[rowSums(SVs.filt[,13:16])<=3,]

# b. Less than or 10 reads in total
SVs.filt.nogerm <- SVs.filt.nogerm[-which(rowSums(SVs.filt.nogerm[,9:16])<=10),]

## additional rules for subsetting tumour-shared and tumour-unique variants, e.g...

## Potentially somatic subset of shared T1/T2
## 10 or more reads in both T1 and T2, less than or 3 reads in T3 & T4
T1.T2.shared.only <- SVs.filt.nogerm[which(SVs.filt.nogerm[,'T1']>=10 & SVs.filt.nogerm[,'T2']>=10 & SVs.filt.nogerm[,'T3']<=3 & SVs.filt.nogerm[,'T4']<=3),]

## Somatic subset of T1-unique
## 10 or more reads in T1 but less than or 3 reads in T2, T3 & T4
T86.unique <- SVs.filt.nogerm[which(SVs.filt.nogerm[,'T1']>=10 & SVs.filt.nogerm[,'T2']<=3 & SVs.filt.nogerm[,'T3']<=3 & SVs.filt.nogerm[,'T4']<=3),]


## 6. Circos-Plots
##################

# Set image-path
setwd("XYZ/My_Plots/Circos/")

circos.shared <- function(x, title, sample1, sample2){
  
  # 1. Define borders
  chromosome.ranges <- matrix(0, ncol=2, nrow=8)
  chromosome.ranges[,1] <- 1
  chromosome.ranges[,2] <- c(680437123,732629474,635140686,480844255,297391021,262034688,85730105,16627)
  rownames(chromosome.ranges) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "ChrX", "MT")
  
  total_1 <- length(which(as.character(x[,'BP1'])!=''))
  total_2 <- length(which(as.character(x[,'BP2'])!=''))
  
  # 2. Start plotting
  png(paste0(title, '.png'), width = 2000, height = 1900)
  color.ramp.up <- max(c(x[,sample1]+x[,sample2])/2)   # scaled to 1.0
  color.ramp.down <- min(c(x[,sample1]+x[,sample2])/2) # scaled to 0.4
  delta <- c(1-0.4)/c(color.ramp.up-color.ramp.down)
  color.ramp <- 0.4+c(c(c(x[,sample1]+x[,sample2])/2)-color.ramp.down)*delta
  
  circos.par("track.height" = 0.2, cell.padding = c(0, 0, 0, 0))
  circos.initialize(factors = rownames(chromosome.ranges), 
                    xlim=chromosome.ranges)
  circos.trackPlotRegion(ylim = c(0, 1), 
                         panel.fun = function(x, y) {print(get.cell.meta.data("xlim"))}, 
                         track.height=0.02, bg.col=c(2:9), bg.border=c(2:9),
                         track.index=1)
  
  # Add gene labels of first breakpoint
  if(total_1>1){
    for (i in 1:total_1){
      chr <- x[which(as.character(x[,'BP1'])!=''),1][i]
      circos.axis(sector.index = chr, 
                  major.at = x[which(as.character(x[,'BP1'])!=''),2][i], 
                  labels = as.character(x[which(as.character(x[,'BP1'])!=''),'BP1'])[i],
                  direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                  labels.away.percentage=1/1.2, minor.ticks = 4)
    }    
  }else if(total_1==1){
    chr <- x[which(as.character(x[,'BP1'])!=''),1]
    circos.axis(sector.index = chr, 
                major.at = x[which(as.character(x[,'BP1'])!=''),2], 
                labels = as.character(x[which(as.character(x[,'BP1'])!=''),'BP1']),
                direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                labels.away.percentage=1/1.2, minor.ticks = 4)    
  }
  
  # Add gene labels of second breakpoint
  if(total_2>1){
    for (i in 1:total_2){
      chr <- x[which(as.character(x[,'BP2'])!=''),4][i]
      circos.axis(sector.index = chr, 
                  major.at = x[which(as.character(x[,'BP2'])!=''),5][i],
                  labels = as.character(x[which(as.character(x[,'BP2'])!=''),'BP2'])[i],
                  direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                  labels.away.percentage=1/1.2, minor.ticks = 4)
    }
  }else if(total_2==1){
    chr <- x[which(as.character(x[,'BP2'])!=''),4]
    circos.axis(sector.index = chr,
                major.at = x[which(as.character(x[,'BP2'])!=''),5],
                labels = as.character(x[which(as.character(x[,'BP2'])!=''),'BP2']),
                direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                labels.away.percentage=1/1.2, minor.ticks = 4)
  }
  
  # Add chromosome names
  for (i in 1:length(rownames(chromosome.ranges))){
    circos.axis(h='top',sector.index = rownames(chromosome.ranges)[i], 
                major.at = chromosome.ranges[i,2]/2, 
                labels = rownames(chromosome.ranges)[i],
                direction = "outside", major.tick.percentage = 1, labels.cex=3,
                labels.away.percentage=1/1.2, minor.ticks = 4)
  }
  
  circos.genomicLink(region1=x[,1:3],region2=x[,4:6], 
                     col=rgb(0,0,1,color.ramp), lwd=5,
                     rou=0.9)
  circos.clear()
  dev.off()
  
}
circos.shared(T1.T2.shared.only, 'T1.T2_shared', 'T1', 'T2', 'Y') # example

circos.unique <- function(x, sample_unique, y, sample_shared1, sample_shared2, title){
  
  # 1. Define borders
  chromosome.ranges <- matrix(0, ncol=2, nrow=8)
  chromosome.ranges[,1] <- 1
  chromosome.ranges[,2] <- c(680437123,732629474,635140686,480844255,297391021,262034688,85730105,16627)
  rownames(chromosome.ranges) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "ChrX", "MT")
  
  total_1.un <- length(which(as.character(x[,'BP1'])!=''))
  total_2.un <- length(which(as.character(x[,'BP2'])!=''))
  total_1.sh <- length(which(as.character(y[,'BP1'])!=''))
  total_2.sh <- length(which(as.character(y[,'BP2'])!=''))
  
  # 2. Start plotting
  png(paste0(title, '.png'), width = 2000, height = 1900)
  
  # shared 
  color.ramp.up <- max(c(y[,sample_shared1]+y[,sample_shared2])/2)   # scaled to 1.0
  color.ramp.down <- min(c(y[,sample_shared1]+y[,sample_shared2])/2) # scaled to 0.4
  delta <- c(1-0.4)/c(color.ramp.up-color.ramp.down)
  color.ramp.y <- 0.4+c(c(c(y[,sample_shared1]+y[,sample_shared2])/2)-color.ramp.down)*delta
  
  # unique 
  color.ramp.up <- max(x[,sample_unique])   # scaled to 1.0
  color.ramp.down <- min(x[,sample_unique]) # scaled to 0.4
  delta <- c(1-0.4)/c(color.ramp.up-color.ramp.down)
  color.ramp.x <- 0.4+c(x[,sample_unique]-color.ramp.down)*delta
  
  circos.par("track.height" = 0.2, cell.padding = c(0, 0, 0, 0))
  circos.initialize(factors = rownames(chromosome.ranges), 
                    xlim=chromosome.ranges)
  circos.trackPlotRegion(ylim = c(0, 1), 
                         panel.fun = function(x, y) {print(get.cell.meta.data("xlim"))}, 
                         track.height=0.02, bg.col=c(2:9), bg.border=c(2:9),
                         track.index=1)
  
  
  # Add gene labels of first breakpoint (shared)
  if(total_1.sh>1){
    for (i in 1:total_1.sh){
      chr <- y[which(as.character(y[,'BP1'])!=''),1][i]
      circos.axis(sector.index = chr, 
                  major.at = y[which(as.character(y[,'BP1'])!=''),2][i], 
                  labels = as.character(y[which(as.character(y[,'BP1'])!=''),'BP1'])[i],
                  direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                  labels.away.percentage=1/1.2, minor.ticks = 4)
    }    
  }else if(total_1.sh==1){
    chr <- y[which(as.character(y[,'BP1'])!=''),1]
    circos.axis(sector.index = chr, 
                major.at = y[which(as.character(y[,'BP1'])!=''),2], 
                labels = as.character(y[which(as.character(y[,'BP1'])!=''),'BP1']),
                direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                labels.away.percentage=1/1.2, minor.ticks = 4)    
  }
  
  # Add gene labels of second breakpoint (shared)
  if(total_2.sh>1){
    for (i in 1:total_2.sh){
      chr <- y[which(as.character(y[,'BP2'])!=''),4][i]
      circos.axis(sector.index = chr, 
                  major.at = y[which(as.character(y[,'BP2'])!=''),5][i],
                  labels = as.character(y[which(as.character(y[,'BP2'])!=''),'BP2'])[i],
                  direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                  labels.away.percentage=1/1.2, minor.ticks = 4)
    }
  }else if(total_2.sh==1){
    chr <- y[which(as.character(y[,'BP2'])!=''),4]
    circos.axis(sector.index = chr,
                major.at = y[which(as.character(y[,'BP2'])!=''),5],
                labels = as.character(y[which(as.character(y[,'BP2'])!=''),'BP2']),
                direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                labels.away.percentage=1/1.2, minor.ticks = 4)
  }
  
  # Add gene labels of first breakpoint (unique)
  if(total_1.un>1){
    for (i in 1:total_1.un){
      chr <- x[which(as.character(x[,'BP1'])!=''),1][i]
      circos.axis(sector.index = chr, 
                  major.at = x[which(as.character(x[,'BP1'])!=''),2][i], 
                  labels = as.character(x[which(as.character(x[,'BP1'])!=''),'BP1'])[i],
                  direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                  labels.away.percentage=1/1.2, minor.ticks = 4)
    }    
  }else if(total_1.un==1){
    chr <- x[which(as.character(x[,'BP1'])!=''),1]
    circos.axis(sector.index = chr, 
                major.at = x[which(as.character(x[,'BP1'])!=''),2], 
                labels = as.character(x[which(as.character(x[,'BP1'])!=''),'BP1']),
                direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                labels.away.percentage=1/1.2, minor.ticks = 4)    
  }
  
  # Add gene labels of second breakpoint  (unique)
  if(total_2.un>1){
    for (i in 1:total_2.un){
      chr <- x[which(as.character(x[,'BP2'])!=''),4][i]
      circos.axis(sector.index = chr, 
                  major.at = x[which(as.character(x[,'BP2'])!=''),5][i],
                  labels = as.character(x[which(as.character(x[,'BP2'])!=''),'BP2'])[i],
                  direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                  labels.away.percentage=1/1.2, minor.ticks = 4)
    }
  }else if(total_2.un==1){
    chr <- x[which(as.character(x[,'BP2'])!=''),4]
    circos.axis(sector.index = chr,
                major.at = x[which(as.character(x[,'BP2'])!=''),5],
                labels = as.character(x[which(as.character(x[,'BP2'])!=''),'BP2']),
                direction = "inside", major.tick.percentage = 1, labels.cex=1.5,
                labels.away.percentage=1/1.2, minor.ticks = 4)
  }
  
  
  # Add chromosome names
  for (i in 1:length(rownames(chromosome.ranges))){
    circos.axis(h='top',sector.index = rownames(chromosome.ranges)[i], 
                major.at = chromosome.ranges[i,2]/2, 
                labels = rownames(chromosome.ranges)[i],
                direction = "outside", major.tick.percentage = 1, labels.cex=3,
                labels.away.percentage=1/1.2, minor.ticks = 4)
  }
  
  # Add shared links
  circos.genomicLink(region1=y[,1:3],region2=y[,4:6], 
                     col=rgb(0,0,1,color.ramp.y), lwd=5,
                     rou=0.9)
  
  # Add unique links
  circos.genomicLink(region1=x[,1:3],region2=x[,4:6], 
                     col=rgb(1,0,0,color.ramp.x), lwd=5,
                     rou=0.9)
  
  circos.clear()
  dev.off()
  
}
circos.unique(T1.unique, 'T1', T1.T2.shared.only, 'T1', 'T2', 'T1-unique', 'Y') # example