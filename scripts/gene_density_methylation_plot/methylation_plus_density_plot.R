library(tidyverse)
library(hash)
library(viridis)
library(scales)



plot.new()

tick_vert_dist = 0.02
chrom_line_weight=0.5
tick_weight=0.2
horizontal_compression_factor = 3
max_chrom_len = 101369167 * horizontal_compression_factor
chromosomes = read.table("Chromosomes.txt", header=TRUE)
sw = read.table("methylation.txt", header=TRUE, fill=TRUE)
sw2 = read.table("genes_density.txt", header=TRUE)

min_hit_len = 300
evalue_thresh = 1e-6
methylation_compression_factor = 3000
gene_density_compression_factor = 7


ChromNameHash = hash(chromosomes$SeqName, chromosomes$Chromosome)
ChromLengthHash = hash(chromosomes$SeqName, chromosomes$Length)
ChromHeightHash = hash(chromosomes$SeqName, chromosomes$HeightOffset)
ChromMitoPercentHash = hash(chromosomes$SeqName, chromosomes$mito)
ChromPlastidPercentHash = hash(chromosomes$SeqName, chromosomes$plastid)



palette(viridis(50))



for (i in 1:dim(chromosomes)[1]){
  segments(0, ChromHeightHash[[chromosomes$SeqName[i]]], ChromLengthHash[[chromosomes$SeqName[i]]] / max_chrom_len, ChromHeightHash[[chromosomes$SeqName[i]]], lwd=chrom_line_weight)
  sw_filt = sw %>% filter(Chromosome == chromosomes$SeqName[i])
  lines(x=(sw_filt$WindowPosition)/max_chrom_len, y = sw_filt$MethProp / methylation_compression_factor + ChromHeightHash[[chromosomes$SeqName[i]]] + 0.002, col=alpha(col='goldenrod2',0.8), lwd=0.3)
  sw2_filt = sw2 %>% filter(Chromosome == chromosomes$SeqName[i])
  lines(x=(sw2_filt$WindowPosition)/max_chrom_len, y = sw2_filt$Genes_per_kb / gene_density_compression_factor + ChromHeightHash[[chromosomes$SeqName[i]]] + 0.002, col=alpha(col='dodgerblue3',0.8), lwd=0.3)
  text(x=-0.03, y=ChromHeightHash[[chromosomes$SeqName[i]]], adj=0, cex=0.25, font=2, ChromNameHash[[chromosomes$SeqName[i]]])
}

legend_x = 0.315
legend_y_bottom = 0.4
scale_height = 0.05 / gene_density_compression_factor
scale_width=0.01
legend_y_top = 0.725
legend_y_top2 = legend_y_top
segments (x0=legend_x, y0=legend_y_bottom, x1=legend_x, y1=legend_y_bottom + scale_height, col=alpha(col='dodgerblue3',0.8), lend=2, lwd=0.5)
text (x=legend_x + 0.005, y= legend_y_bottom + scale_height/2, adj=0, cex=0.2, "0.05 genes per kb")

legend_x = 0.315
legend_y_bottom = 0.3
scale_height = 50 / methylation_compression_factor
scale_width=0.01
legend_y_top = 0.725
legend_y_top2 = legend_y_top
segments (x0=legend_x, y0=legend_y_bottom, x1=legend_x, y1=legend_y_bottom + scale_height, col=alpha(col='goldenrod2',0.8), lend=2, lwd=0.5)
text (x=legend_x + 0.005, y= legend_y_bottom + scale_height/2, adj=0, cex=0.2, "0.5 CpG Methylation")




for (i in 100:51){
  legend_y_top = legend_y_top - 0.004
}

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



