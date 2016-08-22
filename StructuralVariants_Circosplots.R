###################### BRASS circos plots ###########################
########################### Version 1.0 #############################

## mrs72
## Last Update - 22/08/2016 ##
##############################

library(circlize)


## 1. Load SVs (as .bedpe file)
###############################

## Test (example with BRASS-output of just a single-BAM input run)
T_SVs <- read.table("/Users/ms37/Desktop/Labwork/2016/08_August/05_08_BRASS/08_08_Stringent_Hostfiltering (1 read)/86T_unique_SVs.bedpe",sep = '\t', header=F)
colnames(T_SVs) <- c("CHR-1",	"START-1",	"END-1", 
                     "CHR-2",	"START-2",	"END-2", 
                     "STRAND 1", "STRAND 2", "READS")


## 2. Filter
############

# e.g Ends of contigs, ends of supercontigs and full supercontig filter (devil-spec.)
# e.g. Simple repeats (devil spec.)
# e.g. Read-threshold: > 15 supports
thresh <- 15
T_SVs.filt <- T_SVs[which(T_SVs[,"READS"]>thresh),]


## 3. Formatting
################

# e.g. Translate contig-positions to chromosomal positions (devil-spec.)


## 4. Circos
############

# Set image-path
setwd("XYZ/My_Plots/Circos/")

circos.image <- function(x, title){
  
  # 1. Define borders
  chromosome.ranges <- matrix(0, ncol=2, nrow=8)
  chromosome.ranges[,1] <- 1 # define all chromosomal start positions (= 1)
  chromosome.ranges[,2] <- c(680437123,732629474,635140686,480844255,297391021,262034688,85730105,16627) # define all chromosomal lengths
  rownames(chromosome.ranges) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "ChrX", "MT")


  # 2. Plot
  png(paste0(title, '.png'), width = 2000, height = 1900)
  color.ramps <- 1-c(1/x[,grep('READS', colnames(x))[1]]) # darker colour for SVs with more read supports (could still think of a nicer function...)
  circos.par("track.height" = 0.2, cell.padding = c(0, 0, 0, 0))
  
  # Initialise
  circos.initialize(factors = rownames(chromosome.ranges), 
                    xlim=chromosome.ranges)
  
  # Build circos-'frame'
  circos.trackPlotRegion(ylim = c(0, 1), 
                         panel.fun = function(x, y) {print(get.cell.meta.data("xlim"))}, 
                         track.height=0.02, bg.col=c(2:c(nrow(chromosome.ranges)+1)), 
                         bg.border=c(2:c(nrow(chromosome.ranges)+1)),
                         track.index=1)
  
  # Add chromosome names
  for (i in 1:nrow(chromosome.ranges)){
    circos.axis(h='top',sector.index = rownames(chromosome.ranges)[i], 
                major.at = chromosome.ranges[i,2]/2, 
                labels = rownames(chromosome.ranges)[i],
                direction = "outside", major.tick.percentage = 1, labels.cex=3,
                labels.away.percentage=1/1.2, minor.ticks = 4)
  }

  # Add links
  circos.genomicLink(region1=x[,1:3],region2=x[,4:6], 
                     col=rgb(1,0,0,color.ramps*0.3),
                     lwd=5, rou=0.9)
  circos.clear()
  dev.off()
}
circos.image(T_SVs.filt, 'T_SVs_filtered')