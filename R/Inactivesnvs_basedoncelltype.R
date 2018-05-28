###This script looks at the number of snvs, that are inactive and outside the exons in different cell types.
###inactivity is based on chromhmm in different celltypes and on the ensembl regulatory build.

library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(VariantAnnotation)
library(optparse)


###Parse command line arguments
option_list = list(
  make_option(c("-f", "--FILELIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", 
              help="The filelist", metavar="character"),
  make_option(c("-o", "--OUTPUT_PATH"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", 
              help="The output path", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###Read in the filelist
filelist = read.table(opt$FILELIST, header = T, stringsAsFactors = F, sep = "\t")
if (sum(duplicated(filelist$Sample)) >= 1){
  stop("The file list contains duplicate samples")
}
filelist = filelist %>% filter(grepl("Child|child", Name))
folder = paste0(opt$OUTPUT_PATH, "/Possibly_pathogenic/SNVs/Filteredonfreq/")
filelist$vcf_files_pathogenic = paste(folder, filelist$Sample, "_freqfilter.vcf", sep = "")

###Determine the different states and celltypes
activestates = c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "10_TssBiv", "11_BivFlnk", "12_EnhBiv")
inactivestates = c("9_Het", "13_ReprPC", "14_ReprPCWk", "15_Quies")
chromhmms = c("Fetal_Brain_Female", "Fetal_Brain_Male", "Brain_Hippocampus_Middle", "Brain_Germinal_Matrix", "CD14_Primary_Cells")
chromhmmsbrainonly = c("Fetal_Brain_Female", "Fetal_Brain_Male", "Brain_Hippocampus_Middle", "Brain_Germinal_Matrix")
chromhmmsfetalonly = c("Fetal_Brain_Female", "Fetal_Brain_Male")

###Read in the vcfs
infovcfstate = DataFrame()
for (i in 1:nrow(filelist)){
  vcfname = filelist$vcf_files_pathogenic[i]
  sample = filelist$Sample[i]
  if (file.exists(vcfname)){
    vcf = readVcf(vcfname, "hg19")
    infovcfstate_sample = info(vcf)[,c(chromhmms, "Ensemblregbuild", "ANN")]
    infovcfstate = rbind(infovcfstate, infovcfstate_sample)
  }
}
infovcfstate$ANN = sapply(infovcfstate$ANN, function(x) {gsub("^.\\|", "", x) %>% gsub("\\|.*", "", .) %>% unique(.) %>% grep("synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant", ., value = T) %>% paste(., collapse = ";")})
infovcfstate$ANN= sapply(infovcfstate$ANN, function(x) {strsplit(x, "&|;") %>% .[[1]] %>% unique(.) %>% grep("synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant", ., value = T) %>% paste(., collapse = ";")})
infovcfstate$exonic = infovcfstate$ANN != ""
infovcfstate$Exoniclabel = ifelse(infovcfstate$exonic, "Exonic", "nonExonic")

###Determine how many snvs are inactive
inactive_bool = apply(infovcfstate[,chromhmms], 2, function(x) x %in% inactivestates)
inactives = data.frame(stringsAsFactors = F)
for (i in 1:5){
  all = sum(rowSums(inactive_bool) >= i & is.na(infovcfstate$Ensemblregbuild) & infovcfstate$exonic == F)
  brain = sum(rowSums(inactive_bool[,chromhmmsbrainonly]) >= i & is.na(infovcfstate$Ensemblregbuild) & infovcfstate$exonic == F)
  fetal = sum(rowSums(inactive_bool[,chromhmmsfetalonly]) >= i & is.na(infovcfstate$Ensemblregbuild) & infovcfstate$exonic == F)
  inactives_i = data.frame("nrcellsinactive" = i, "celltype" = c("all", "brain", "fetal"), "sum" = c(all, brain, fetal))
  inactives = rbind(inactives, inactives_i)
}

###Create figure showing how many snvs are inactive
inactivesnvsfig = ggplot(data = inactives, aes(x = nrcellsinactive, y = sum, fill = celltype)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(y = "Nr. inactive snvs", x = "Nr. of celltypes in which the snv is inactive", title = "Number of snvs that are inactive based on different celltypes") +
  theme_bw() +
  coord_cartesian(ylim = c(0,1500), expand = F)

pdf(paste0(opt$OUTPUT_PATH, "/Possibly_pathogenic/inactivesnvs.pdf"))
inactivesnvsfig
dev.off()
