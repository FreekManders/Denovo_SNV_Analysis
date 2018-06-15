library(plyr)
library(dplyr)
library(ggplot2)
library(VariantAnnotation)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(cba)
library(circlize)
library(ggpubr)
library(VennDiagram)
library(gridExtra)
library(optparse)

####___________________Parse command line arguments__________________________________________####
option_list = list(
  make_option(c("-f", "--FILELIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", 
              help="The filelist", metavar="character"),
  make_option(c("-o", "--OUTPUT_PATH"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", 
              help="The output path", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

####_______________Read in the filelist_____________________________________________####
filelist = read.table(opt$FILELIST, header = T, stringsAsFactors = F, sep = "\t")
if (sum(duplicated(filelist$Sample)) >= 1){
  stop("The file list contains duplicate samples")
}
filelist = filelist %>% filter(grepl("Child|child", Name))
folder = paste0(opt$OUTPUT_PATH, "/Callableregion_NrMendelViols/TrueDenovo/SNVs/")
filelist$vcf_files = paste(folder, filelist$Sample, "_Truedenovo.vcf", sep = "")

filelist$Runname = as.factor(filelist$Run)
levels(filelist$Runname) = 1:length(levels(filelist$Runname))
filelist$Runname = paste("Run: ", filelist$Runname, sep = "")

####_______________Read in vcfs_____________________________________________________####
vcfs = NULL
for (i in 1:nrow(filelist)){
	vcfname = filelist$vcf_files[i]
  if (!file.exists(vcfname)){
    next
  }
  vcf = readVcf(vcfname, "hg19")
  effects = sapply(info(vcf)$ANN, function(x) {gsub("^.\\|", "", x) %>% gsub("\\|.*", "", .) %>% unique(.) %>% grep("synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant", ., value = T) %>% paste(., collapse = ";")})
  DANN = as.numeric(info(vcf)$DANN)
  regbuild = info(vcf)$Ensemblregbuild
  Fetal_Brain_Female = info(vcf)$Fetal_Brain_Female
  Fetal_Brain_Male = info(vcf)$Fetal_Brain_Male
  Brain_Hippocampus_Middle = info(vcf)$Brain_Hippocampus_Middle
  Brain_Germinal_Matrix = info(vcf)$Brain_Germinal_Matrix
  CD14_Primary_Cells = info(vcf)$CD14_Primary_Cells
  phastcons = info(vcf)$phastCons46way
  gnomAD_AC = as.numeric(info(vcf)$gnomAD_AC)
  GoNLv5_AC = as.numeric(info(vcf)$GoNLv5_AC)
  PON_COUNT = as.numeric(info(vcf)$PON_COUNT)
  df = data.frame("Effects" = effects, "regbuild" = regbuild, "Fetal_Brain_Female" = Fetal_Brain_Female, "Fetal_Brain_Male" = Fetal_Brain_Male, "Brain_Hippocampus_Middle" = Brain_Hippocampus_Middle, "Brain_Germinal_Matrix" = Brain_Germinal_Matrix, "CD14_Primary_Cells" = CD14_Primary_Cells, "DANN" = DANN, "phastcons" = phastcons, "PON_COUNT" = PON_COUNT, "GoNLv5_AC" = GoNLv5_AC, "gnomAD_AC" = gnomAD_AC, "Sample" = filelist$Sample[i], "Run" = filelist$Runname[i], stringsAsFactors = F)
  vcfs = rbind(vcfs, df)
}
vcfs$vsplit = as.factor(vcfs$Run)
nrlevels = nlevels(vcfs$vsplit)
levels(vcfs$vsplit) = c(rep(1, ceiling(nrlevels/2)), rep(2, floor(nrlevels/2)))
vcfs$regbuild[is.na(vcfs$regbuild)] = "Not regulatory"
chromhmms = c("Fetal_Brain_Female", "Fetal_Brain_Male", "Brain_Hippocampus_Middle", "Brain_Germinal_Matrix", "CD14_Primary_Cells")
vcfs[,chromhmms] = apply(vcfs[,chromhmms], 2, function(x) {mapvalues(x, from = c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "9_Het", "10_TssBiv", "11_BivFlnk", "12_EnhBiv", "13_ReprPC", "14_ReprPCWk", "15_Quies"), to = c("Active TSS", "Flanking Active TSS", "Transcr. at gene 5' and 3'", "Strong transcription", "Weak transcription", "Genic enhancers", "Enhancers", "ZNF genes & repeats", "Heterochromatin", "Bivalent/Poised TSS", "Flanking Bivalent TSS/Enh", "Bivalent Enhancer", "Repressed PolyComb", "Weak Repressed PolyComb", "Quiescent/Low"))})
vcfs$Effects = sapply(vcfs$Effects, function(x) {strsplit(x, "&|;") %>% .[[1]] %>% unique(.) %>% grep("synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant", ., value = T) %>% paste(., collapse = ";")})
vcfs$exonic = vcfs$Effects != ""
vcfs$Exoniclabel = ifelse(vcfs$exonic, "Exonic", "nonExonic")
vcfs$majorrows = apply(vcfs[,chromhmms], 1, function(x) {sum(duplicated(x)) >= 2})
vcfs$majorchromhmm = apply(vcfs[,chromhmms], 1, function(x) {names(sort(table(x), decreasing = T))[1]})
vcfs$gnomAD_AC[is.na(vcfs$gnomAD_AC)] = 0
vcfs$GoNLv5_AC[is.na(vcfs$GoNLv5_AC)] = 0
vcfs$PON_COUNT[is.na(vcfs$PON_COUNT)] = 0

chrstatecolors = c("Red", "Orange Red", "LimeGreen", "Green", "DarkGreen", "GreenYellow", "Yellow", "Medium Aquamarine", "PaleTurquoise", "IndianRed", "DarkSalmon", "DarkKhaki", "#C0C0C0", "Gainsboro", "Black")
names(chrstatecolors) = c("Active TSS", "Flanking Active TSS", "Transcr. at gene 5' and 3'", "Strong transcription", "Weak transcription", "Genic enhancers", "Enhancers", "ZNF genes & repeats", "Heterochromatin", "Bivalent/Poised TSS", "Flanking Bivalent TSS/Enh", "Bivalent Enhancer", "Repressed PolyComb", "Weak Repressed PolyComb", "Quiescent/Low")



####_______________Create figures for DANN scores___________________________________####

danndistrFig = ggplot(data = vcfs, aes(y = DANN, x = "")) + 
  geom_boxplot() +
  geom_dotplot(stackdir = "center", binaxis = "y", position = "dodge", binwidth = 0.025, binpositions = "all", dotsize = 0.4, stackratio = 0.8) +
  labs(x = "All samples", title = "DANN scores of all samples") +
  theme_bw() + 
  theme(text = element_text(size=20))

danndistrRunFig = list()
for (i in 1:2){
vcfs2 = vcfs %>% filter(vsplit == i)
dotplot = ggplot(data = vcfs2, aes(x = Sample, y = DANN)) + 
  geom_boxplot() +
  geom_dotplot(aes(fill = Sample), stackdir = "center", binaxis = "y", position = "dodge", binwidth = 0.025, binpositions = "all", dotsize = 0.4, stackratio = 0.75) + 
  facet_grid(. ~ Run, scales = "free", space = "free") +
  labs(title = "DANN scores per sample") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, size = 10, margin = margin(t = 20)), text = element_text(size=20))
  danndistrRunFig[[i]] = dotplot
}

dannbyregionFig = ggplot(data = vcfs, aes(y = DANN, x = Exoniclabel)) + 
  geom_boxplot() +
  geom_dotplot(stackdir = "center", binaxis = "y", position = "dodge", binwidth = 0.025, binpositions = "all", dotsize = 0.4, stackratio = 0.8) +
  labs(x = "Genomic region", title = "DANN scores of all samples by genomic region") +
  stat_compare_means(label.x.npc = 0.5, label.y = 1.05, size = 10) + 
  theme_bw() +
  theme(text = element_text(size=20))

exonic = vcfs %>% dplyr::filter(exonic)
dannbyeffectFig = ggplot(data = exonic, aes(y = DANN, x = "", fill = Effects)) + 
  geom_dotplot(stackgroups = T, stackdir = "center", binaxis = "y", binwidth = 0.025, binpositions = "all", dotsize = 0.4, stackratio = 0.8) +
  labs(x = "", title = "DANN scores of all exonic samples") +
  theme_bw() +
  theme(axis.ticks.x = element_blank()) +
  theme(text = element_text(size=20))

###Dann by regbuild
reglevels = levels(as.factor(vcfs[,"regbuild"]))
pvals = c()
for (reg in reglevels){
  if (reg != "Not regulatory"){
    regdanns = vcfs %>% filter(regbuild == reg) %>% dplyr::select(DANN)
    nonregdanns = vcfs %>% filter(regbuild == "Not regulatory") %>% dplyr::select(DANN)
    test = wilcox.test(regdanns[,,drop = T], nonregdanns[,,drop = T], alternative = "two.sided", exact = F)
    pvals = c(pvals, test[[3]])
    }
}
padj = p.adjust(pvals, method = "fdr")
padj = sprintf("fdr: %1.3f", padj)
notpval = grep("Not regulatory", reglevels)
padj = c(padj[1:(notpval-1)], "", padj[notpval:length(padj)])

ann_text = data.frame(regbuild = reglevels, DANN = 1, lab = "Test")
danndistrregbuildFig = ggplot(data = vcfs, aes(x = regbuild, y = DANN)) + 
  geom_boxplot() +
  geom_dotplot(aes(fill = regbuild), stackdir = "center", binaxis = "y", position = "dodge", binwidth = 0.025, binpositions = "all", dotsize = 0.4, stackratio = 0.75) + 
  facet_grid(. ~ regbuild, scales = "free", space = "free") +
  labs(title = "DANN scores per Ensembl regbuild annotation", x = "") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, size = 12, margin = margin(t = 30), vjust = 0.75), text = element_text(size = 20)) +
  geom_text(data = ann_text, label = padj)

###Dann by chromhmm sites occuring in more than 3 celltypes.
vcfs3 = vcfs %>% filter(majorrows == T)
reglevels = levels(as.factor(vcfs3[,"majorchromhmm"]))
pvals = c()
notpvals = c()
for (reg in reglevels){
  if (reg != "Quiescent/Low"){
    regdanns = vcfs3 %>% filter(majorchromhmm == reg) %>% dplyr::select(DANN)
    if (nrow(regdanns) > 1){
      nonregdanns = vcfs3 %>% filter(majorchromhmm == "Quiescent/Low") %>% dplyr::select(DANN)
      test = wilcox.test(regdanns[,,drop = T], nonregdanns[,,drop = T], alternative = "two.sided", exact = F)
      pvals = c(pvals, test[[3]])
    }
    else{
      notpval = grep(reg, reglevels)
      notpvals = c(notpvals, notpval)
    }
  }
}
padj = p.adjust(pvals, method = "fdr")
padj = sprintf("fdr: %1.3f", padj)
notpvals = c(notpvals, grep("Quiescent/Low", reglevels))
for (notpval in notpvals){
  padj = c(padj[1:(notpval-1)], "", padj[notpval:length(padj)])
}

ann_text = data.frame(majorchromhmm = reglevels, DANN = 1.1, lab = "Test")
danndistrChromhmmFig = ggplot(data = vcfs3, aes(x = majorchromhmm, y = DANN)) + 
  geom_boxplot() +
  geom_dotplot(aes(fill = majorchromhmm), stackdir = "center", binaxis = "y", position = "dodge", binwidth = 0.025, binpositions = "all", dotsize = 0.4, stackratio = 0.75) + 
  facet_grid(. ~ majorchromhmm, scales = "free", space = "free") +
  labs(title = "DANN scores per ChromHMM segmentation", x = "") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, size = 12, margin = margin(t = 30), vjust = 0.73), text = element_text(size = 20)) +
  geom_text(data = ann_text, label = padj)

pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/DANNscores.pdf"), width = 15, height = 10)
danndistrFig
danndistrRunFig[[1]]
danndistrRunFig[[2]]
dannbyregionFig
dannbyeffectFig
danndistrregbuildFig
danndistrChromhmmFig
dev.off()

###Create quantile table
cutoffs = seq(0, 1, 0.1)
quant = quantile(vcfs$DANN, probs = cutoffs)
quantnames = names(quant)
quant = c("Allsamples", quant)
quants = vcfs %>% group_by(Sample) %>% do(data.frame(t(quantile(.$DANN, probs = cutoffs)))) %>% as.data.frame()
quants = rbind(quant, quants)
colnames(quants) = c("Sample", quantnames)
write.table(quants, paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/DANNscoresquantiles.txt"), sep = "\t", quote = F, row.names = F)
write.table(vcfs[,c("DANN")], paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/alldanns.txt"), quote = F, row.names = F, col.names = F)

####_______________Regulatory features and chromhmm_________________________________####
vcfs2 = vcfs %>% dplyr::select("Sample", "Run", "regbuild", chromhmms) %>% melt(id = c("Sample", "Run", "regbuild"), value.name = "Chromatin_state", variable.name = "Celltype_or_tissue") %>% group_by(Celltype_or_tissue, Chromatin_state) %>% dplyr::summarise(count = n())
mutsbystate = ggplot(data = vcfs2, aes(x = Celltype_or_tissue, fill = Chromatin_state, y = count)) +
  geom_bar(stat = "identity", position = "dodge") +
  coord_cartesian(ylim = c(1, 1600), expand = F) + 
  labs(x = "Celtype or tissue", y = "Nr. of denovo SNVS", title = "Nr of denovo SNVs in different chromatin states") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, size = 12, margin = margin(t = 15), vjust = 0.75), text = element_text(size = 20)) + 
  scale_fill_manual(name = "Chromatin state", values = chrstatecolors) +
  scale_y_log10()

vcfs3 = vcfs %>% dplyr::select("Sample", "Run", "regbuild", chromhmms) %>% melt(id = c("Sample", "Run", "regbuild"), value.name = "Chromatin_state", variable.name = "Celltype_or_tissue") %>% group_by(regbuild, Celltype_or_tissue, Chromatin_state) %>% dplyr::summarise(count = n())
mutsbystateandbuild = ggplot(data = vcfs3, aes(x = Celltype_or_tissue, fill = Chromatin_state, y = count)) +
  geom_bar(stat = "identity", position = "dodge") +
  facet_grid(. ~ regbuild, scales = "free", space = "free") +
  coord_cartesian(ylim = c(1, 1600), expand = F) + 
  labs(x = "Celtype or tissue", y = "Nr. of denovo SNVs", title = "Nr of denovo SNVs in different chromatin states") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 20, size = 8, margin = margin(t = 15))) + 
  scale_fill_manual(name = "Chromatin state", values = chrstatecolors) +
  scale_y_log10()


quiescent = vcfs[rowSums(vcfs[,chromhmms] == "Quiescent/Low") == 5, ]
sum(quiescent$regbuild == "Not regulatory" & quiescent$exonic == F)

#Cooccurence of regbuild with Chromhmm heatmaps
makeheats = function(cooccur, title){
cooccur[is.na(cooccur)] = 0
row.names(cooccur) = cooccur$regbuild
cooccur = dplyr::select(cooccur, -regbuild)
cooccur2 = log10(cooccur)

sort.cor.col = cor(cooccur) #calculate Pearson correlations
sort.cor.dist.col = as.dist(1-sort.cor.col) #create distance table based on correlations
sort.cor.dist.hclust.col = hclust(sort.cor.dist.col, method="complete")
sort.cor = cor(t(cooccur)) #calculate Pearson correlations
sort.cor.dist = as.dist(1-sort.cor) #create distance table based on correlations
sort.cor.dist.hclust = hclust(sort.cor.dist, method="complete")
sort.cor.dist.hclust.ddg = as.dendrogram(sort.cor.dist.hclust)
sort.cor.dist.hclust.col.ddg = as.dendrogram(sort.cor.dist.hclust.col)

opt <- order.optimal(sort.cor.dist, sort.cor.dist.hclust$merge)
optimal <- sort.cor.dist.hclust
optimal$merge <- opt$merge
optimal$order <- opt$order
opt.col <- order.optimal(sort.cor.dist.col, sort.cor.dist.hclust.col$merge)
optimal.col <- sort.cor.dist.hclust.col
optimal.col$merge <- opt.col$merge
optimal.col$order <- opt.col$order
optimal.ddg = as.dendrogram(optimal)
optimal.col.ddg = as.dendrogram(optimal.col)

regbuildvschromhmm = Heatmap(cooccur2, col = colorRamp2(seq(0,4,0.4),rev(brewer.pal(11, "RdBu"))), column_dend_reorder = F,
              cluster_columns = as.dendrogram(optimal.col.ddg), cluster_rows = as.dendrogram(optimal.ddg), na_col = "grey",
              row_names_side = "left", column_names_side = "top", show_row_names = T, show_column_names = T, column_title = paste(title, "\nEpigenomics Roadmap", sep = ""), row_title = "Ensembl regulatory build",
              heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal", legend_width = unit(4, "cm"), title_position = "topcenter"),
              name = "Log10(Regbuild vs ChromHMM)", row_names_gp = gpar(fontsize = 18), column_names_gp = gpar(fontsize = 18))

cooccur3 = cooccur / rowSums(cooccur)
regbuildvschromhmmrel = Heatmap(cooccur3, col = colorRamp2(seq(0,0.7,0.07),rev(brewer.pal(11, "RdBu"))), column_dend_reorder = F,
              cluster_columns = as.dendrogram(optimal.col.ddg), cluster_rows = as.dendrogram(optimal.ddg), na_col = "grey",
              row_names_side = "left", column_names_side = "top", show_row_names = T, show_column_names = T, column_title = paste(title, " relative\nEpigenomics Roadmap", sep = ""), row_title = "Ensembl regulatory build",
              heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal", legend_width = unit(4, "cm"), title_position = "topcenter"),
              name = "Fraction of regbuild", row_names_gp = gpar(fontsize = 18), column_names_gp = gpar(fontsize = 18))

cooccur4 = apply(cooccur, 2, function(x) x/sum(x))/rowSums(cooccur)*sum(cooccur) #Divide by row and col sums and then multiply by total counts to find enrichment.
cooccur4 = log2(cooccur4)
cooccur4[is.infinite(cooccur4)] = NA
regbuildvschromhmmenrich = Heatmap(cooccur4, col = colorRamp2(seq(-6,6,1.2),rev(brewer.pal(11, "RdBu"))), column_dend_reorder = F,
                                cluster_columns = as.dendrogram(optimal.col.ddg), cluster_rows = as.dendrogram(optimal.ddg), na_col = "grey",
                                row_names_side = "left", column_names_side = "top", show_row_names = T, show_column_names = T, column_title = paste(title, " relative\nEpigenomics Roadmap", sep = ""), row_title = "Ensembl regulatory build",
                                heatmap_legend_param = list(color_bar = "continuous", legend_direction = "horizontal", legend_width = unit(4, "cm"), title_position = "topcenter", at = c(-6, -3, 0, 3, 6)),
                                name = "Enrichment (log2)", row_names_gp = gpar(fontsize = 18), column_names_gp = gpar(fontsize = 18))

return(c(regbuildvschromhmm, regbuildvschromhmmrel, regbuildvschromhmmenrich))
}

cooccur = vcfs %>% dplyr::select("Sample", "Run", "regbuild", chromhmms) %>% melt(id = c("Sample", "Run", "regbuild"), value.name = "Chromatin_state", variable.name = "Celltype_or_tissue") %>% group_by(regbuild, Chromatin_state) %>% dplyr::summarise(count = n())
cooccur = dcast(cooccur, regbuild~Chromatin_state, value.var = "count")
write.table(cooccur, paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/Regbuild_ChromHMM_Cooccurence.txt"), sep = "\t", quote = F, row.names = F)
figscooccur = makeheats(cooccur, "Regbuild vs ChromHMM")

cooccurmajor = vcfs %>% filter(majorrows == T) %>% dplyr::select(regbuild, majorchromhmm) %>% group_by(regbuild, majorchromhmm) %>% dplyr::summarise(count = n())
cooccurmajor = dcast(cooccurmajor, regbuild~majorchromhmm, value.var = "count")
write.table(cooccurmajor, paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/Regbuild_MajorChromHMM_Cooccurence.txt"), sep = "\t", quote = F, row.names = F)
figscooccurmajor = makeheats(cooccurmajor, "Regbuild vs ChromHMM occuring in at least 3 sites")


pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/regbuild_chromhmm.pdf"), width = 15, height = 10)
mutsbystate
mutsbystateandbuild
figscooccur[[1]]
figscooccur[[2]]
figscooccur[[3]]
figscooccurmajor[[1]]
figscooccurmajor[[2]]
figscooccurmajor[[3]]
dev.off()

####_______________Phastcons________________________________________________________####

nrNA = sum(is.na(vcfs$phastcons))
nrnotNA = sum(!is.na(vcfs$phastcons))
navsnotna = data.frame("category" = c("No score (no alignment)", "Score present"), "counts" = c(nrNA, nrnotNA))
navsnotnaFig = ggplot(data = navsnotna, aes(x = category, y = counts)) +
  geom_bar(stat = "identity", fill = "darkblue") +
  labs(x = "", y = "Nr. of denovo SNVs", title = "Phastscore present") +
  theme_bw() + 
  theme(text = element_text(size=20))

vcfs4 = vcfs %>% dplyr::filter(!is.na(phastcons)) %>% mutate(phastcons = as.numeric(phastcons))

phastdistrFig = ggplot(data = vcfs4, aes(y = phastcons, x = "")) + 
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(stackdir = "center", binaxis = "y", position = "dodge", binwidth = 0.005, binpositions = "all", dotsize = 0.5, stackratio = 0.8) +
  labs(x = "All samples", y = "Phastcons 46-way scores", title = "Phastcons 46-way scores of all samples") +
  theme_bw() + 
  theme(text = element_text(size=20))

phastdistrRunFig = list()
for (i in 1:2){
  vcfs5 = vcfs4 %>% filter(vsplit == i)
  dotplot = ggplot(data = vcfs5, aes(x = Sample, y = phastcons)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_dotplot(aes(fill = Sample), stackdir = "center", binaxis = "y", position = "dodge", binwidth = 0.01, binpositions = "all", dotsize = 0.6, stackratio = 0.75) + 
    facet_grid(. ~ Run, scales = "free", space = "free") +
    labs(title = "Phastcons 46-ways scores per sample", y = "Phastcons 46-way scores") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 75, size = 10, margin = margin(t = 20)), text = element_text(size=20))
  phastdistrRunFig[[i]] = dotplot
}

DannPhastCor = round(cor(vcfs4$DANN, vcfs4$phastcons), 3)
phastvsdannFig = ggplot(data = vcfs4, aes(x = phastcons, y = DANN)) +
  geom_point() +
  geom_smooth(method = "lm") +
  labs(title = "Phastcons 46-ways vs DANN", x = "Phastcons 46-way scores", y = "DANN") + 
  theme_bw() +
  theme(text = element_text(size=20)) +
  annotate("text", x = 0.04, y = 1.1, size = 8, label = paste0("Cor: ", DannPhastCor))

phastbyregionFig = ggplot(data = vcfs4, aes(y = phastcons, x = Exoniclabel)) + 
  geom_boxplot() +
  labs(x = "Genomic region", title = "Phastcons 46-way scores of all samples by genomic region") +
  stat_compare_means(label.x.npc = 0.5, label.y = 1.05, size = 10) + 
  theme_bw() +
  theme(text = element_text(size=20))

exonic = vcfs4 %>% dplyr::filter(exonic)
phastbyeffectFig = ggplot(data = exonic, aes(y = phastcons, x = "", fill = Effects)) + 
  geom_dotplot(stackgroups = T, stackdir = "center", binaxis = "y", binwidth = 0.025, binpositions = "all", dotsize = 0.4, stackratio = 0.8) +
  labs(x = "", title = "Phastcons 46-way scores of all exonic samples") +
  theme_bw() +
  theme(axis.ticks.x = element_blank(), text = element_text(size=20))


###Phastcons by regbuild
reglevels = levels(as.factor(vcfs4[,"regbuild"]))
pvals = c()
for (reg in reglevels){
  if (reg != "Not regulatory"){
    regdanns = vcfs4 %>% filter(regbuild == reg) %>% dplyr::select(phastcons)
    nonregdanns = vcfs4 %>% filter(regbuild == "Not regulatory") %>% dplyr::select(phastcons)
    test = wilcox.test(regdanns[,,drop = T], nonregdanns[,,drop = T], alternative = "two.sided", exact = F)
    pvals = c(pvals, test[[3]])
  }
}
padj = p.adjust(pvals, method = "fdr")
padj = sprintf("fdr: %1.3f", padj)
notpval = grep("Not regulatory", reglevels)
padj = c(padj[1:(notpval-1)], "", padj[notpval:length(padj)])

ann_text = data.frame(regbuild = reglevels, phastcons = 1, lab = "Test")
phastdistrregbuildFig = ggplot(data = vcfs4, aes(x = regbuild, y = phastcons)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(aes(fill = regbuild), stackdir = "center", binaxis = "y", position = "dodge", binwidth = 0.025, binpositions = "all", dotsize = 0.4, stackratio = 0.75) + 
  facet_grid(. ~ regbuild, scales = "free", space = "free") +
  labs(title = "Phastcons 46-way scores per Ensembl regbuild annotation", x = "", y = "Phastcons 46-way scores") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, size = 12, margin = margin(t = 30), vjust = 0.75), text = element_text(size = 20)) +
  geom_text(data = ann_text, label = padj)

###Phastcons by chromhmm sites occuring in more than 3 celltypes.
vcfs6 = vcfs4 %>% filter(majorrows == T)
reglevels = levels(as.factor(vcfs6[,"majorchromhmm"]))
pvals = c()
notpvals = c()
for (reg in reglevels){
  if (reg != "Quiescent/Low"){
    regdanns = vcfs6 %>% filter(majorchromhmm == reg) %>% dplyr::select(phastcons)
    if (nrow(regdanns) > 1){
      nonregdanns = vcfs6 %>% filter(majorchromhmm == "Quiescent/Low") %>% dplyr::select(phastcons)
      test = wilcox.test(regdanns[,,drop = T], nonregdanns[,,drop = T], alternative = "two.sided", exact = F)
      pvals = c(pvals, test[[3]])
    }
    else{
      notpval = grep(reg, reglevels)
      notpvals = c(notpvals, notpval)
    }
  }
}
padj = p.adjust(pvals, method = "fdr")
padj = sprintf("fdr: %1.3f", padj)
notpvals = c(notpvals, grep("Quiescent/Low", reglevels))
for (notpval in notpvals){
  padj = c(padj[1:(notpval-1)], "", padj[notpval:length(padj)])
}

ann_text = data.frame(majorchromhmm = reglevels, phastcons = 1.1, lab = "Test")
phastdistrChromhmmFig = ggplot(data = vcfs6, aes(x = majorchromhmm, y = phastcons)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_dotplot(aes(fill = majorchromhmm), stackdir = "center", binaxis = "y", position = "dodge", binwidth = 0.025, binpositions = "all", dotsize = 0.4, stackratio = 0.75) + 
  facet_grid(. ~ majorchromhmm, scales = "free", space = "free") +
  labs(title = "DANN scores per ChromHMM segmentation", x = "") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle = 75, size = 12, margin = margin(t = 30), vjust = 0.73), text = element_text(size = 20)) +
  geom_text(data = ann_text, label = padj)


pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/Phastconscores.pdf"), width = 20, height = 10)
navsnotnaFig
phastdistrFig
phastdistrRunFig[[1]]
phastdistrRunFig[[2]]
phastvsdannFig
phastbyregionFig
phastbyeffectFig
phastdistrregbuildFig
phastdistrChromhmmFig
dev.off()


###Create quantile table
cutoffs = seq(0, 1, 0.1)
quant = quantile(vcfs4$phastcons, probs = cutoffs)
quantnames = names(quant)
quant = c("Allsamples", quant)
quants = vcfs4 %>% group_by(Sample) %>% do(data.frame(t(quantile(.$phastcons, probs = cutoffs)))) %>% as.data.frame()
quants = rbind(quant, quants)
colnames(quants) = c("Sample", quantnames)
write.table(quants, paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/Phastcons_46way_scoresquantiles.txt"), sep = "\t", quote = F, row.names = F)


####_______________Occurence denovo SNVs in databases(gnomAD, GoNL and HMF PON)_____####

###gnomAD + GoNL exonic vs non exonic
vcfs$pathogenic[vcfs$GoNLv5_AC == 0 & vcfs$gnomAD_AC == 0 & vcfs$PON_COUNT == 0] = "New variant"
vcfs$pathogenic[vcfs$GoNLv5_AC == 1 | vcfs$gnomAD_AC == 1 | vcfs$PON_COUNT == 1] = "Only occurs once in database"
vcfs$pathogenic[vcfs$GoNLv5_AC > 1 | vcfs$gnomAD_AC > 1 | vcfs$PON_COUNT > 1] = "Known variant (filtered out)"

mutsindb = vcfs %>% dplyr::select(Exoniclabel, pathogenic) %>% group_by(Exoniclabel, pathogenic) %>% summarise(counts = n()) %>% mutate(groupsize = sum(counts), countsbygroup = counts / groupsize)
dbsbyregionsAbsFig = ggplot(data = mutsindb, aes(y = counts, x = Exoniclabel, fill = pathogenic)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
  labs(y = "Nr of denovo SNVs", x = "Genomic region", title = "gnomAD, GoNL and HMF-PON occurences") + 
  coord_cartesian(ylim = c(0, 2000), expand = F) +
  theme(axis.text.x = element_text(angle = 80, size = 10, margin = margin(t = 25))) +
  theme(panel.grid.major.x = element_blank())


dbsbyregionsRelFig = ggplot(data = mutsindb, aes(y = countsbygroup, x = Exoniclabel, fill = pathogenic)) + 
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() + 
  labs(y = "Percentage of denovo SNVs", x = "Genomic region", title = "gnomAD, GoNL and HMF-PON occurences") + 
  scale_y_continuous(labels = scales::percent) +
  coord_cartesian(ylim = c(0, 1), expand = F) +
  theme(axis.text.x = element_text(angle = 80, size = 10, margin = margin(t = 25)), panel.grid.major.x = element_blank())

pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/gnomAD_GoNL_exonsvsnonexon.pdf"))
dbsbyregionsAbsFig
dbsbyregionsRelFig
dev.off()


###gnomAD + GoNL AC distribution
nrGoNL = sum(vcfs$GoNLv5_AC != 0)
nrgnomAD = sum(vcfs$gnomAD_AC != 0)
nrPON = sum(vcfs$PON_COUNT != 0)
nrGoNLgnomAD = sum(vcfs$GoNLv5_AC != 0 & vcfs$gnomAD_AC != 0)
nrGoNLPON = sum(vcfs$GoNLv5_AC != 0 & vcfs$PON_COUNT != 0)
nrgnomADPON = sum(vcfs$gnomAD_AC != 0 & vcfs$PON_COUNT != 0)
nrall = sum(vcfs$GoNLv5_AC != 0 & vcfs$gnomAD_AC != 0 & vcfs$PON_COUNT != 0)

vennshared = draw.triple.venn(nrGoNL, nrgnomAD, nrPON, n12 = nrGoNLgnomAD, n13 = nrGoNLPON, n23 = nrgnomADPON, n123 = nrall, category = c("GoNL","gnomAD", "HMF-PON"),inverted = F, ext.text = T, fill = c("green","mediumorchid","skyblue"), cat.cex = c(2,2,2), cex = c(rep(3,7)), ind = F)

acCounts = vcfs %>% dplyr::rename(GoNLv5 = GoNLv5_AC, gnomAD = gnomAD_AC, HMF_PON = PON_COUNT) %>% dplyr::select(GoNLv5, gnomAD, HMF_PON, Sample) %>% melt(id = "Sample") %>% filter(value != 0)
databasesdistFig = ggplot(data = acCounts, aes(y = value, x = variable)) + 
  geom_boxplot() +
  geom_dotplot(stackdir = "center", binaxis = "y", position = "dodge", binwidth = 0.04, binpositions = "all", dotsize = 1, stackratio = 0.8) +
  scale_y_log10() +
  labs(y = "Allele count", x = "Database", title = "gnomAD, GoNL and HMF-PON allele counts in the denovo SNVs") +
  theme_bw()

pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/gnomAD_GoNL_distribution"))
grid.arrange(gTree(children=vennshared), top="Denovo SNVs in gnomAD, GoNL and HMF-PON")
databasesdistFig
dev.off()

