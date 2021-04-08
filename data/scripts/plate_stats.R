get_plate_stats <- function(){

library(dplyr)
library(ggplot2)
library(platetools)
library(viridis)
args = commandArgs(trailingOnly=TRUE)

source=args[1]
output=args[2]

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}



df_counts <- read.table(paste0(source, "BC.stats.txt"), sep = "\t", header=F)
df_quality <- read.table(paste0(source, "collapsed.stats.txt"), fill = TRUE, sep="\t", header=T)

#df_counts processing

colnames(df_counts) <- c("barcode", "count", "wells")

df_counts[["ID"]] <- gsub("^results.*(P[0-9]{2}[A-Za-z]{1,10}[0-9]{0,2})\\..*fasta","\\1", df_counts$wells)
df_counts[["plate_id"]] <- gsub("P([0-9]{1,2}).*", "plate_\\1", gsub("P0", "P", df_counts[["ID"]]))
df_counts[["well_id"]] <- gsub("P[0-9]{1,2}", "", df_counts[["ID"]])

df_counts[["ID"]] <- paste0(df_counts[["plate_id"]], "_", df_counts[["well_id"]])

df_counts_unmatched <- df_counts[grep("unmatched" ,df_counts[["well_id"]]), c("count", "ID")]

complete_wells <- rep(num_to_well(1:96), 12)
complete_plate_id <- rep(sprintf("plate_%s",1:12), each = 96)
ID <- paste0(complete_plate_id, "_", complete_wells)

df_counts_plates <- data.frame(complete_plate_id=complete_plate_id, complete_wells=complete_wells, ID = ID )

final_plates <- left_join(df_counts_plates, df_counts, by = "ID")
final_plates[["complete_plate_id"]] <- factor(final_plates[["complete_plate_id"]], levels= sprintf("plate_%s",1:12))

rainbow_colours <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan","#7FFF7F",  "#FF7F00", "red", "#7F0000"))
b <- c(10, 1000, 2000, 3000, 4000, 5000)


pdf(file=paste0(output, ".demultiplexed.pdf"), width = 10, height = 20)
p <- raw_grid(data = final_plates[["count"]],
       well = final_plates[["complete_wells"]],
       plate_id = final_plates[["complete_plate_id"]],) +
    #scale_fill_distiller(type = "div") +
    ggtitle("Sequences demultiplexed per well") +
    scale_fill_gradientn(colors = rainbow_colours(5), breaks = b, labels = format(b))
print(p)
dev.off()

#df_quality processing

colnames(df_quality) <- gsub("^.*P[0]{0,1}([0-9]{1,2})([A-Z][0-9]{1,2})","plate_\\1_\\2",colnames(df_quality))
df_quality_pct <- apply(df_quality, MARGIN = 2, FUN = function(x){sum(x[1:3], na.rm = T )/sum(x,  na.rm = T )})
wells=gsub("^.*_","", names(df_quality_pct))
plate_id=gsub("_[A-Z].*","", names(df_quality_pct))

df_quality_pct <- data.frame(score=df_quality_pct,
                              plate_id=plate_id,
                              wells=wells,
                              ID=paste0(plate_id, "_", wells))

final_plates_qc <- left_join(df_counts_plates, df_quality_pct, by = "ID")
#final_plates_qc[["score"]][is.na(final_plates_qc[["score"]])] <- 0
final_plates_qc[["complete_plate_id"]] <- factor(final_plates_qc[["complete_plate_id"]], levels= sprintf("plate_%s",1:12))

b2 <- c(0.1,0.3,0.5)
colours <- colorRampPalette(c("red", "yellow", "#7FFF7F"))

pdf(file=paste0(output, ".quality.pdf"), width = 10, height = 20)
p2 <- raw_grid(data = final_plates_qc[["score"]],
          well = final_plates_qc[["complete_wells"]],
          plate_id = final_plates_qc[["complete_plate_id"]],
          plate = 96) +
       #scale_fill_distiller(type = "div") +
       ggtitle("Collapsed quality score") +
       scale_fill_gradientn(colors = colours(3), breaks = b2, labels = format(b2))
print(p2)
dev.off()


write.table(x = final_plates, file = paste0(output, ".df_counts.tsv"), quote = F, row.names = F)
write.table(x = final_plates_qc, file = paste0(output, ".df_quality.tsv"), quote = F, row.names = F)
write.table(x = df_counts_unmatched, file = paste0(output, ".df_unmatched.tsv"), quote = F, row.names = F)
}

get_plate_stats()
