args = commandArgs(trailingOnly=TRUE)

# args[1] must be the snp file
# args[2] is the genome size in mb
# args[3] is the centimorgan by mb

snps <- read.table(args[1], stringsAsFactors = FALSE)
n.marker <- dim(snps)[1]

## Map file
# Marker names
marker1 <- "M"
marker2 <- 1:n.marker
marker2 <- sprintf("%03d", marker2)
marker <-paste0(marker1,marker2)

# Chromossome and position
tot = as.numeric(args[2])*as.numeric(args[3])

# The markers will be equally distribuited. There will be one marker each
by = (tot)/(n.marker-1)

pos <- seq(from = 0, to = tot, by = by)
chr <- rep("C1",length(pos))

map_file <- data.frame(marker=marker, chromosome=chr, position= pos)
write.table(map_file, file = paste0("mapfile.map"), quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

# Using sn1ps simulated by pirs genome described in ref_oryza.snp.lst file. Looking markers at the position:

P1_1 <- snps$V4
P1_2 <- snps$V4
P2_1 <- snps$V5
P2_2 <- snps$V5

chr <- snps$V1
pos <- snps$V2

founder_file <- data.frame(marker=marker, P1_1 , P1_2, P2_1, P2_2)

write.table(founder_file, file = paste0("founderfile.gen"), quote=FALSE, col.names = TRUE, row.names = FALSE, sep = "\t" )

## Parameters file

parameter <- paste0("PLOIDY = 2
                    MAPFUNCTION = HALDANE
                    MISSING = NA
                    CHROMFILE = inb.chrom
                    POPTYPE = F2
                    POPSIZE = 150
                    MAPFILE = mapfile.map
                    FOUNDERFILE = founderfile.gen
                    OUTPUT = sim_inb")

write.table(parameter, file = paste0("sim.par"), quote=FALSE, col.names = FALSE, row.names = FALSE, sep = "\t" )

chrom <- data.frame("chromosome"= "C1", "length"= tot, "centromere"=tot/2, "prefPairing"= 0.0, "quadrivalents"=0.0)

write.table(chrom, file= "inb.chrom", quote = F, col.names = T, row.names = F, sep= "\t")
