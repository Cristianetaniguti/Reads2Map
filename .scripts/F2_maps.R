
# Path
path <- "~/github/errors_workflow/cromwell-executions/F2/7f3c41c1-2ddf-4b08-b42f-39fa5016f04e/"

setwd(path)
cMByMb <- 4.63


# Packages

library(ggplot2)
library(reshape2)

## Gatk

# Map distances

maps_gatk <- read.table("call-all_maps/shard-0/execution/gatk_map_df.txt", header =T)

poscM <- (as.numeric(as.character(maps[,2]))/1000000)*cMByMb
poscM.norm <- poscM-poscM[1]

maps_gatk <- cbind(maps_gatk, poscM, poscM.norm)

diff.dis <- sqrt((maps_gatk$poscM.norm - maps_gatk$rf)^2)
diff.dis.mean <- mean(sqrt((maps_gatk$poscM.norm - maps_gatk$rf)^2))
diff.dis.median <- median(sqrt((maps_gatk$poscM.norm - maps_gatk$rf)^2))

diff.tot.dis <- sqrt((maps_gatk$poscM.norm[length(maps_gatk$poscM.norm)] - maps_gatk$rf[length(maps_gatk$rf)])^2)

mapfile <- read.table("call-pedsim_files/execution/mapfile.map", header=T)

tot_mks <- read.table("call-pedsim_files/execution/tot_mks.txt")

mapfile.pos <- cbind(mapfile, tot_mks[,2])

coverage <- maps_gatk$pos[length(maps_gatk$pos)]*100/mapfile.pos[,4][length(mapfile.pos[,4])]

## Depths

alt_depth <- read.table("call-aval_vcf/execution/gatk_alt_depth.txt", header = T)


ref_depth <- read.table("call-aval_vcf/execution/gatk_ref_depth.txt", header = T)

                                        # Pensar em como vai representar as profundidades

