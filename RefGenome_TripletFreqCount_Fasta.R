############ Processing Devil 7.1 Genome Fasta File #################
########################### Version 1.0 #############################

## Maximilian Stammnitz, mrs72
## Last Update - 08/12/2015 ##
##############################

library(Biostrings)

# Import Fasta-Genome
Devil_7.1 <- readDNAStringSet("/Users/ms37/Desktop/Data/References/Reference_7.1.fa", format="fasta", use.names=TRUE)

# Adjust Contig Regions (+- 999 BP) # other filters removed #
Devil_7.1 <- padAndClip(Devil_7.1, IRanges(1000, as.numeric(Fasta.regions[,3])), Lpadding.letter=".", Rpadding.letter=".")
                              
# Calculate
Devil_7.1.trinucleotides <- trinucleotideFrequency(Devil_7.1)
Devil_7.1.trinucleotides <- apply(Devil_7.1.trinucleotides, 2, sum)

# Plot
pdf("09_12_Devil7.1_Trinucleotides64.pdf", height=5, width=11)
barplot(Devil_7.1.trinucleotides, col="black",
        cex.names=0.2, 
        main="Trinucleotide Frequencies across Devil 7.1 (specified high-Quality contigs and regions)")
dev.off()

# Convert 64-tripletts into 32 'one-strand' tripletts
equivalent.triplett <- matrix(c('ACA', 'TGT', ### First 16/32
                                'ACC', 'TGG',
                                'ACG', 'TGC',
                                'ACT', 'TGA',
                                'CCA', 'GGT',
                                'CCC', 'GGG',
                                'CCG', 'GGC',
                                'CCT', 'GGA',
                                'GCA', 'CGT',
                                'GCC', 'CGG',
                                'GCG', 'CGC',
                                'GCT', 'CGA',
                                'TCA', 'AGT',
                                'TCC', 'AGG',
                                'TCG', 'AGC',
                                'TCT', 'AGA',
                                'ATA', 'TAT', ### Second 16/32
                                'ATC', 'TAG',
                                'ATG', 'TAC',
                                'ATT', 'TAA',
                                'CTA', 'GAT',
                                'CTC', 'GAG',
                                'CTG', 'GAC',
                                'CTT', 'GAA',
                                'GTA', 'CAT',
                                'GTC', 'CAG',
                                'GTG', 'CAC',
                                'GTT', 'CAA',
                                'TTA', 'AAT',
                                'TTC', 'AAG',
                                'TTG', 'AAC',
                                'TTT', 'AAA'),nrow=32, ncol=2, byrow=TRUE)
equivalent.triplett <- cbind(equivalent.triplett,c(1:nrow(equivalent.triplett)),c(1:nrow(equivalent.triplett)),c(1:nrow(equivalent.triplett)))
colnames(equivalent.triplett) <- c("Strand1", "Strand2", "Tripletts1", "Tripletts2", "Sum")
equivalent.triplett[,3] <- Devil_7.1.trinucleotides[match(equivalent.triplett[,1], names(Devil_7.1.trinucleotides))]
equivalent.triplett[,4] <- Devil_7.1.trinucleotides[match(equivalent.triplett[,2], names(Devil_7.1.trinucleotides))]
equivalent.triplett[,5] <- as.numeric(equivalent.triplett[,3])+as.numeric(equivalent.triplett[,4])

# Save
write.table(equivalent.triplett, "/Users/ms37/Desktop/Data/Info-Files/Devil_7.1Reference_Tripletts.txt", col.names=T, quote=F, row.names=F)