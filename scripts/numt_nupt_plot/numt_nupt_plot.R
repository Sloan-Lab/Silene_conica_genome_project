library(tidyverse)
library(hash)
library(viridis)
library(scales)


plot.new()

tick_vert_dist = 0.02
chrom_line_weight=0.5
tick_weight=0.2
horizontal_compression_factor = 4
max_chrom_len = 101606156 * horizontal_compression_factor
chromosomes = read.table("Chromosomes.txt", header=TRUE)
plastid_hits= read.table("Silene_conica_blastSummary_e-6.plastid.hsps.hits.no_rRNA.txt", header=TRUE)
plastid_hits = plastid_hits[order(plastid_hits$PercentID),]
mito_hits = read.table("Silene_conica_blastSummary_e-6.mito.hsps.hits.no_rRNA.txt", header=TRUE)
mito_hits = mito_hits[order(mito_hits$PercentID),]
plastid_sw = read.table("plastid.sw.no_rRNA.txt", header=TRUE)
mito_sw = read.table("mito.sw.no_rRNA.txt", header=TRUE)

min_hit_len = 300
evalue_thresh = 1e-6

ChromNameHash = hash(chromosomes$SeqName, chromosomes$Chromosome)
ChromLengthHash = hash(chromosomes$SeqName, chromosomes$Length)
ChromHeightHash = hash(chromosomes$SeqName, chromosomes$HeightOffset)
ChromMitoPercentHash = hash(chromosomes$SeqName, chromosomes$mito)
ChromPlastidPercentHash = hash(chromosomes$SeqName, chromosomes$plastid)



palette(viridis(50))

for (i in 1:dim(mito_hits)[1]){
  if(mito_hits$HitLength[i] >= min_hit_len & mito_hits$E.value[i] <= evalue_thresh){
    segments (x0=mito_hits$HitStart[i]/max_chrom_len, y0=ChromHeightHash[[mito_hits$SeqName[i]]],x1=mito_hits$HitStart[i]/max_chrom_len,y1=ChromHeightHash[[mito_hits$SeqName[i]]]+tick_vert_dist, col=mito_hits$PercentID[i], lwd=tick_weight)
  }
}

for (i in 1:dim(plastid_hits)[1]){
  if(plastid_hits$HitLength[i] >= min_hit_len & plastid_hits$E.value[i] <= evalue_thresh){
    segments (x0=plastid_hits$HitStart[i]/max_chrom_len, y0=ChromHeightHash[[plastid_hits$SeqName[i]]],x1=plastid_hits$HitStart[i]/max_chrom_len,y1=ChromHeightHash[[plastid_hits$SeqName[i]]]-tick_vert_dist, col=plastid_hits$PercentID[i], lwd=tick_weight)
  }
}


for (i in 1:dim(chromosomes)[1]){
  segments(0, ChromHeightHash[[chromosomes$SeqName[i]]], ChromLengthHash[[chromosomes$SeqName[i]]] / max_chrom_len, ChromHeightHash[[chromosomes$SeqName[i]]], lwd=chrom_line_weight)
  mito_sw_filt = mito_sw %>% filter(Chromosome == chromosomes$SeqName[i])
  lines(x=(mito_sw_filt$WindowPosition)/max_chrom_len, y = mito_sw_filt$Coverage + ChromHeightHash[[chromosomes$SeqName[i]]] + 0.002, col=alpha(col='firebrick4',0.8), lwd=0.3)
  plastid_sw_filt = plastid_sw %>% filter(Chromosome == chromosomes$SeqName[i])
  lines(x=(plastid_sw_filt$WindowPosition)/max_chrom_len, y = -1*plastid_sw_filt$Coverage + ChromHeightHash[[chromosomes$SeqName[i]]] - 0.002, col=alpha(col='firebrick4',0.8), lwd=0.3)
  text(x=-0.03, y=ChromHeightHash[[chromosomes$SeqName[i]]], adj=0, cex=0.25, font=2, ChromNameHash[[chromosomes$SeqName[i]]])
  text(x=ChromLengthHash[[chromosomes$SeqName[i]]] / max_chrom_len + 0.005, y=ChromHeightHash[[chromosomes$SeqName[i]]] + 0.008, adj=0, cex=0.15, paste("Mito:    ", ChromMitoPercentHash[[chromosomes$SeqName[i]]]))
  text(x=ChromLengthHash[[chromosomes$SeqName[i]]] / max_chrom_len + 0.005, y=ChromHeightHash[[chromosomes$SeqName[i]]] - 0.008, adj=0, cex=0.15, paste("Plastid:", ChromPlastidPercentHash[[chromosomes$SeqName[i]]]))
}

legend_x = 0.315
legend_y_bottom = 0.3
scale_height = 0.05
scale_width=0.01
legend_y_top = 0.725
legend_y_top2 = legend_y_top
segments (x0=legend_x, y0=legend_y_bottom, x1=legend_x, y1=legend_y_bottom + scale_height, col=alpha(col='firebrick4',0.8), lend=2, lwd=0.5)
text (x=legend_x + 0.005, y= legend_y_bottom + scale_height/2, adj=0, cex=0.2, "5%")

for (i in 100:51){
  segments(x0=legend_x, y0=legend_y_top, x1=legend_x+scale_width, y1=legend_y_top, col=i, lend=2, lwd=8*tick_weight)
  legend_y_top = legend_y_top - 0.004
}

text(x=legend_x, y=legend_y_bottom + scale_height + 0.04, adj = 0, cex=0.2, font=2, "Shared")
text(x=legend_x, y=legend_y_bottom + scale_height + 0.02, adj = 0, cex=0.2, font=2, "DNA (%)")
text(x=legend_x, y=legend_y_top2 + 0.03, adj = 0, cex=0.2, font=2, "Identity (%)")
text(x=legend_x + 0.015, y=legend_y_top2, adj = 0, cex=0.2, 100)
gradient_length = 0.2
text(x=legend_x + 0.015, y=legend_y_top2 - gradient_length, adj = 0, cex=0.2, 50)
text(x=legend_x + 0.015, y=legend_y_top2 - gradient_length/2, adj = 0, cex=0.2, 75)


map_scale_offset = 0.05
segments(0, legend_y_bottom - map_scale_offset, 1/horizontal_compression_factor, legend_y_bottom - map_scale_offset, lwd=chrom_line_weight, col="darkgray")

scale_tick_height = 0.01
tick_pos=0
tick_label=0
tick_inc=20
while (tick_pos < 1/horizontal_compression_factor){
  segments(tick_pos, legend_y_bottom - map_scale_offset, tick_pos, legend_y_bottom - map_scale_offset - scale_tick_height, col="darkgray", lwd=chrom_line_weight)
  text (tick_pos, legend_y_bottom - map_scale_offset - scale_tick_height - 0.01, cex=0.2, paste(tick_label, "Mb"))
  tick_pos = tick_pos + tick_inc*1e6/max_chrom_len
  tick_label = tick_label + tick_inc
}



