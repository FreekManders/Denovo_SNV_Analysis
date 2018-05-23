library(reshape2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(optparse)


####________________Parse command line arguments_________________________####
option_list = list(
  make_option(c("-f", "--FILELIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", 
              help="The filelist", metavar="character"),
  make_option(c("-o", "--OUTPUT_PATH"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", 
              help="The output path", metavar="character"),
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

outpdf = paste0(opt$OUTPUT_PATH, "/WGS_QC/Output/WGS_QC.pdf")


####________________Read filelist and check for duplicates________________####
filelist = read.table(opt$FILELIST, header = T, stringsAsFactors = F, sep = "\t")
if (sum(duplicated(filelist$Sample)) >= 1){
  stop("The file list contains duplicate samples")
}

####_________________Read the WGS metrics and flagstat summary files._____####
#The wgs metrics and flagstat files are combined into one WGS and one flagstat file.
homeadress = ""
folders = unique(filelist$Folder)
WGS = data.frame()
flagstat = data.frame()
for (folder in folders){
  fnameWGS = paste(homeadress ,folder, "/QCStats/WGSMetrics_summary.txt", sep = "")
  WGStable = read.table(fnameWGS, header = T, stringsAsFactors = F)
  WGS = rbind(WGS, WGStable)
  
  fnameflag = paste(homeadress ,folder, "/QCStats/flagstat_summary.txt", sep = "")
  flagtable = read.table(fnameflag, header = F, stringsAsFactors = F, sep = "\t")
  flagtable = t(flagtable)
  colnames(flagtable) = gsub(" ", ".", flagtable[1,])
  flagtable = as.data.frame(flagtable[-1,], stringsAsFactors = F)
  flagstat = rbind(flagstat, flagtable)
}

flagstat$Total.number.of.reads = as.numeric(flagstat$Total.number.of.reads)
flagstat$Percentage.reads.mapped = as.numeric(gsub("%", "", flagstat$Percentage.reads.mapped))/100

####____Join the WGS and flagstat files and extract info from the sample names._____####
intermed_table = full_join(WGS, flagstat, by = "sample")
intermed_table$sampleshort = gsub("-DNA-1_dedup", "", intermed_table$sample)
intermed_table$family = gsub("-.*$", "", intermed_table$sample) %>% factor()

fammember = lapply(intermed_table$sample, function(x){
  sampleno = strsplit(x, "-")[[1]][2];
  if (substring(sampleno, 1, nchar("1")) == "1"){
    return("Father")
  }
  else if (substring(sampleno, 1, nchar("2")) == "2"){
    return("Mother")
  }
  else if (substring(sampleno, 1, nchar("3")) == "3" | substring(sampleno, 1, nchar("4")) == "4" ){
    return("Child")
  }
  else{
    return("NA")
  }
})
intermed_table$fammember = unlist(fammember)


####______Join the intermediate table with the filelist._______####
qctable = left_join(filelist, intermed_table, by = c("Sample" = "sampleshort"))


####_____Make different versions of the table for the different figures.____####
qctable$Runname = as.factor(qctable$Run)
levels(qctable$Runname) = 1:length(levels(qctable$Runname))
qctable$Runname = paste("Run: ", qctable$Runname, sep = "")
qctable$mapped.reads = qctable$Total.number.of.reads * qctable$Percentage.reads.mapped
qctable$unmapped.reads = qctable$Total.number.of.reads - qctable$mapped.reads

qctable2 = qctable %>% select(Sample, family, Runname, unmapped.reads, mapped.reads) %>% melt(id = c("Sample", "family", "Runname"), variable.name = "Legend", value.name = "nrreads")
qctable2$Legend = gsub("\\.", " ", qctable2$Legend)

cols_qctable3 = grep("PCT_.*X$", colnames(qctable))
qctable3 = qctable %>% select(Sample, cols_qctable3) %>% melt(id = c("Sample"), variable.name = "Xaxis", value.name = "Coverage")
qctable3$Xaxis = gsub("PCT_", "", qctable3$Xaxis) %>% gsub("X", "", .) %>% as.numeric()
qctable3$Coverage = qctable3$Coverage * 100

qctable4 = qctable
cols_qctable4 = grep("^PCT_EXC_(?!.*TOTAL)", colnames(qctable), perl = T)
#qctable4[,cols_qctable4] = qctable4[,cols_qctable4] * qctable4$Total.number.of.reads
qctable4 = qctable4 %>% select(Sample, Runname, cols_qctable4) %>% melt(id = c("Sample", "Runname"), variable.name = "Legend", value.name = "nrexcluded")
qctable4$Legend = gsub("PCT_EXC_", "", qctable4$Legend)
qctable4$nrexcluded = qctable4$nrexcluded * 100


####________Create figures and export them in a pdf file.________________####
maxyfig1 = 1.05 * max(qctable$MEAN_COVERAGE + qctable$SD_COVERAGE)
maxyfig2 = 1.05 * max(qctable2$nrreads)
maxyfig4 = 1.5 * max(qctable4$nrexcluded)

nrcolors = length(unique(qctable$family))
set.seed(001)
colors = sample(colorRampPalette(brewer.pal(8, "Accent"))(nrcolors))
MeanCovFig = ggplot(qctable, aes(x = Sample, y = MEAN_COVERAGE, fill = family)) + 
  geom_bar(stat = "identity") + 
  geom_errorbar(aes(x = Sample, ymin = MEAN_COVERAGE - SD_COVERAGE, ymax = MEAN_COVERAGE + SD_COVERAGE)) + 
  geom_hline(aes(yintercept = 30)) + 
  facet_grid(. ~ Runname, scales = "free", space = "free") + 
  theme_bw() + 
  labs(x = "Sample", y = "Mean coverage +/- SD") + 
  coord_cartesian(ylim = c(0, maxyfig1), expand = F) + 
  theme(axis.text.x = element_text(angle = 80, size = 6, margin = margin(t = 20))) + 
  guides(fill=guide_legend("Family")) +
  scale_fill_manual(values = colors)

NrReadsFig = ggplot(qctable2, aes(x = Sample, y = nrreads, fill = Legend)) + 
  facet_grid(. ~ Runname, scales = "free", space = "free") + 
  geom_bar(stat = "identity") + 
  theme_bw() + 
  labs(y = "Number of reads", x = "Sample") + 
  coord_cartesian(ylim = c(0, maxyfig2), expand = F) + 
  theme(axis.text.x = element_text(angle = 80, size = 6, margin = margin(t = 20)))

nrcolors = length(unique(qctable3$Sample))
set.seed(001)
colors = sample(colorRampPalette(brewer.pal(8, "Accent"))(nrcolors))
MinReadCovFig = ggplot(qctable3, aes(x = Xaxis, y=Coverage, color = Sample)) + 
  geom_line() + 
  labs(x = "Minimal read coverage", y = "% of genome") + 
  theme_bw() + 
  coord_cartesian(ylim = c(0,100), xlim = c(0,100), expand = F) + 
  theme(legend.text=element_text(size=6)) + 
  guides(color=guide_legend(title="Sample")) +
  scale_color_manual(values = colors)

ReadExclusionFig = ggplot(qctable4, aes(x = Sample, y = nrexcluded, fill = Legend)) + 
  geom_bar(stat = "identity") + 
  facet_grid(. ~ Runname, scales = "free", space = "free") + 
  theme_bw() + 
  labs(x = "Sample", y = "% reads excluded") + 
  coord_cartesian(ylim = c(0, maxyfig4), expand = F) + 
  theme(axis.text.x = element_text(angle = 80, size = 6, margin = margin(t = 20))) + 
  guides(fill=guide_legend("Exclusion filter"))

pdf(outpdf, paper = "a4")
MeanCovFig
NrReadsFig
MinReadCovFig
ReadExclusionFig
dev.off()
print("done")
