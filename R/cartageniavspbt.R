library(ggplot2)
library(dplyr)
library(VennDiagram)
library(gridExtra)
library(VariantAnnotation)
library(reshape2)
library(ggpubr)
library(optparse)

####___________________Parse command line arguments__________________________________________####
option_list = list(
  make_option(c("-f", "--FILELIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", 
              help="The filelist", metavar="character"),
  make_option(c("-o", "--OUTPUT_PATH"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", 
              help="The output path", metavar="character"),
  make_option(c("--CARTA_DATA"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Characteristics_denovo/Myfiltervscartagenia/Cartagenia_Exon_SNVs.txt", 
              help="The data from cartagenia", metavar="character"),
  make_option(c("--CARTA_FAMS"), type="character", default="1,2,3,4,6,7,8,9,10,11,15,16,17,18,19,20,21,22,23,24,25,26,27,29,30", 
              help="The families analyzed by cartagenia", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


####____________________Read in filelist______________________________________________________####
filelist = read.table(opt$FILELIST, header = T, stringsAsFactors = F, sep = "\t")
if (sum(duplicated(filelist$Sample)) >= 1){
  stop("The file list contains duplicate samples")
}
filelistchild = filelist %>% filter((grepl("Child|child", Name)))

folder = paste0(opt$OUTPUT_PATH, "/Callableregion_NrMendelViols/TrueDenovo/SNVs/")
filelistchild$vcf_files = paste(folder, filelistchild$Sample, "_Truedenovo.vcf", sep = "")

filelistchild$Runname = as.factor(filelistchild$Run)
levels(filelistchild$Runname) = 1:length(levels(filelistchild$Runname))
filelistchild$Runname = paste("Run: ", filelistchild$Runname, sep = "")


####____________________Read in cartagenia data and create plots______________________________####
carta = read.csv(opt$CARTA_DATA, sep = "\t", stringsAsFactors = F)
cartafams = as.integer(base::strsplit(opt$CARTA_FAMS, ",")[[1]])
cartafams = paste0("MP", sprintf("%02d", cartafams))

tpfig = ggplot(carta, aes(x = IGVcheck)) + 
  geom_bar(fill = "darkred") + 
  labs(x = "True or False positive", y = "Nr. of denovo SNVs", title = "True and false positives from cartagenia") + 
  coord_cartesian(ylim = c(0, 60), expand = F) + 
  theme_bw() +
  theme(panel.grid.major.x = element_blank(), text = element_text(size = 20))


cartatrue = carta %>% filter(IGVcheck == "True")

nrsyn = grep("^synonymous", cartatrue$Effect..codingEffect.) %>% length(.)
nrmis = grep("nonsynonymous", cartatrue$Effect..codingEffect.) %>% length(.)
nrstopgain = grep("stopgain", cartatrue$Effect..codingEffect.) %>% length(.)
mutationtypes = data.frame("Mutationtype" = c("Synonymous", "Nonsynonymous", "stopgain"), "Counts" = c(nrsyn, nrmis, nrstopgain))
maxy = max(mutationtypes$Counts) * 1.1
muttypefig = ggplot(data = mutationtypes, aes(x = Mutationtype, y = Counts, fill = Mutationtype)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  coord_cartesian(ylim = c(0, maxy), expand = F) +
  labs(x = "Type of mutation", y = "Nr. of denovo SNVs", title = "Type of exonic mutation in Cartagenia") +
  theme(panel.grid.major.x = element_blank(), text = element_text(size = 20))


nrsamples = length(cartafams)
nrtrue = nrow(cartatrue)
avgtrue = round(nrtrue / nrsamples, 3)

tpbyfam = cartatrue %>% group_by(family) %>% summarise("Nrmuts" = n())
cartafamsdf = data.frame("family" = cartafams, stringsAsFactors = F)
tpbyfam = left_join(cartafamsdf, tpbyfam, by = "family")
tpbyfam = tpbyfam %>% mutate_all(funs(replace(., is.na(.), 0)))

mutsbyfamfig = ggplot(tpbyfam, aes(x = family, y = Nrmuts)) + 
  geom_bar(stat = "identity", fill = "darkblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 80, size = 8, margin = margin(t = 12))) +
  coord_cartesian(ylim = c(0, 5), expand = F) +
  labs(x = "Family", y = "True positives", title = "True positives per sample in cartagenia") +
  annotate("text", x = 9, y = 4.5, label = paste("Average number of true positives per sample: ", avgtrue, sep = "")) +
  theme(panel.grid.major.x = element_blank(), text = element_text(size = 20))

pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/Myfiltervscartagenia/Cartagenia_exon.pdf"))
tpfig
muttypefig
mutsbyfamfig
dev.off()


####____________________Read in true denovo exonic snvs_______________________________________####
vcfs = NULL
for (i in 1:nrow(filelistchild)){
  vcfname = filelistchild$vcf_files[i]
  if (!file.exists(vcfname)){
    next
  }
  vcf = readVcf(vcfname, "hg19")
  effects = sapply(info(vcf)$ANN, function(x) {gsub("^.\\|", "", x) %>% gsub("\\|.*", "", .) %>% unique(.) %>% grep("synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant", ., value = T) %>% paste(., collapse = ";")})
  chrom = seqnames(vcf) %>% as.vector()
  start = ranges(vcf) %>% as.data.frame(.)
  start = start$start
  df = data.frame("chrom" = chrom, "start" = start, "Effects" = effects, "Sample" = filelistchild$Sample[i], "Run" = filelistchild$Runname[i], stringsAsFactors = F)
  vcfs = rbind(vcfs, df)
}
vcfs$family = gsub("-.*", "", vcfs$Sample)
vcfs$Effects = sapply(vcfs$Effects, function(x) {strsplit(x, "&|;") %>% .[[1]] %>% unique(.) %>% grep("synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant", ., value = T) %>% paste(., collapse = ";")})
vcfs$exonic = vcfs$Effects != ""
vcfs$Exoniclabel = ifelse(vcfs$exonic, "Exonic", "nonExonic")

exonic = vcfs %>% filter(exonic)
exonic2 = exonic %>% filter(family %in% cartafams) %>% dplyr::select(chrom, start, Effects, Sample, family)

write.table(exonic, paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/Myfiltervscartagenia/All_Exonic_SNVs.txt"), quote = F, row.names = F, sep = "\t")

####____________________Compare exonic snvs with cartagenia___________________________________####
#My exonic trudenovos have been filtered in igv. This removed only 4 sites.
nrshared = inner_join(exonic2, carta, by = c("family" = "family", "chrom" = "Chromosome", "start" = "Start")) %>% nrow()
nrsharedtrue = inner_join(exonic2, cartatrue, by = c("family" = "family", "chrom" = "Chromosome", "start" = "Start")) %>% nrow()

pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/Myfiltervscartagenia/CartageniasharedTPs.pdf"))
vennshared = draw.pairwise.venn(nrow(exonic2), nrow(carta), nrshared, category = c("WGS custom filter","Alissa Interpret"),inverted = F, ext.text = T, fill = c("darkred","skyblue"), alpha = c(0.8,0.95), cat.pos = c(180,180), cat.cex = c(2,2), cex = c(3,3,3), ind = F)
grid.arrange(gTree(children=vennshared), top="Shared denovo SNVs with all cartagenia snvs")
vennsharetrue = draw.pairwise.venn(nrow(exonic2), nrow(cartatrue), nrsharedtrue, category = c("WGS custom filter","Alissa Interpret"),inverted = F, ext.text = T, fill = c("darkred","skyblue"), alpha = c(0.8,0.95), cat.pos = c(0,0), cat.cex = c(2,2), cex = c(3,3,3), ind = F)
grid.arrange(gTree(children=vennsharetrue), top="Shared denovo SNVs with Alissa Interpret true positives")
dev.off()

missedbyme = anti_join(cartatrue, exonic2, by = c("family" = "family", "Chromosome" = "chrom", "Start" = "start"))
write.table(missedbyme, paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/Myfiltervscartagenia/TPsmissedbyme.txt"), quote = F, row.names = F, sep = "\t")
missedbycarta = anti_join(exonic2, cartatrue, by = c("family" = "family", "chrom" = "Chromosome", "start" = "Start"))
write.table(missedbycarta, paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/Myfiltervscartagenia/TPsmissedbycarta.txt"), quote = F, row.names = F, sep = "\t")

####____________________Characteristics of the true denovo exonic snvs________________________####
nrsyn = grep("synonymous", exonic$Effects) %>% length(.)
nrmis = grep("missense", exonic$Effects) %>% length(.)
nrstopgain = grep("stop_gained", exonic$Effects) %>% length(.)
nrsplice = grep("splice", exonic$Effects) %>% length(.)
mutationtypes = data.frame("Mutationtype" = c("Synonymous", "Missense", "Stop gain", "Splice site"), "Counts" = c(nrsyn, nrmis, nrstopgain, nrsplice))
maxy = max(mutationtypes$Counts) * 1.1
muttypefig = ggplot(data = mutationtypes, aes(x = Mutationtype, y = Counts, fill = Mutationtype)) + 
  geom_bar(stat = "identity") +
  theme_bw() +
  coord_cartesian(ylim = c(0, maxy), expand = F) +
  labs(x = "Type of mutation", y = "Nr of de novo SNVs", title = "Type of exonic mutation in WGS") +
  theme(panel.grid.major.x = element_blank(), text = element_text(size = 20))



nrsamples = nrow(filelistchild)
nrtrue = nrow(exonic)
avgtrue = round(nrtrue / nrsamples, 3)

mutsbys = exonic %>% group_by(Sample) %>% summarise("ExonicMuts" = n())
mutsbys = left_join(filelistchild, mutsbys, by = "Sample")
mutsbys = mutsbys %>% mutate_all(funs(replace(., is.na(.), 0)))

mutsbysfig = ggplot(mutsbys, aes(x = Sample, y = ExonicMuts)) + 
  geom_bar(stat = "identity", fill = "darkblue") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 80, size = 8, margin = margin(t = 20))) +
  coord_cartesian(ylim = c(0, 4), expand = F) +
  labs(x = "Family", y = "True positives", title = "True positives per sample in WGS") +
  annotate("text", x = 9, y = 3.75, label = paste("Average number of true positives per sample: ", avgtrue, sep = "")) +
  theme(panel.grid.major.x = element_blank(), text = element_text(size = 20))


pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/Myfiltervscartagenia/WGS_exon.pdf"))
muttypefig
mutsbysfig
dev.off()
