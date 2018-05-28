###This script removes sites from vcf files that are in inactive non exonic regions,
### based on fetal brain chromhmm and the ensembl regulatory build. The vcfs need to already contain this annotation.
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(VariantAnnotation)
library(optparse)

###Parse command line arguments
option_list = list(
  make_option(c("-f", "--FILELIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", 
              help="The filelist", metavar="character"),
  make_option(c("--OVERWRITE"), type="character", default="true", 
              help="Whether or not to overwrite existing data files", metavar="character"),
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
folder = paste0(opt$OUTPUT_PATH, "/Possibly_pathogenic/SNVs/")
filelist$vcf_files_pathogenic = paste(folder, "Filteredonfreq/", filelist$Sample, "_freqfilter.vcf", sep = "")
filelist$vcf_files_pathogenic_out = paste(folder, "freq_chromfunction/", filelist$Sample, "_freq_chromfunction_filter.vcf", sep = "")


###Look at the number of snvs, that are inactive in different cell types.
activestates = c("1_TssA", "2_TssAFlnk", "3_TxFlnk", "4_Tx", "5_TxWk", "6_EnhG", "7_Enh", "8_ZNF/Rpts", "10_TssBiv", "11_BivFlnk", "12_EnhBiv")
inactivestates = c("9_Het", "13_ReprPC", "14_ReprPCWk", "15_Quies")
chromhmmsfetalonly = c("Fetal_Brain_Female", "Fetal_Brain_Male")

for (i in 1:nrow(filelist)){
  vcfname = filelist$vcf_files_pathogenic[i]
  sample = filelist$Sample[i]
  vcfname_out = filelist$vcf_files_pathogenic_out[i]
  if (file.exists(vcfname) & (!file.exists(vcfname_out) | opt$OVERWRITE == "true")){
    vcf = readVcf(vcfname, "hg19", row.names = F)
    infovcfstate = info(vcf)[,c(chromhmmsfetalonly, "Ensemblregbuild", "ANN")]
    infovcfstate$ANN = sapply(infovcfstate$ANN, function(x) {gsub("^.\\|", "", x) %>% gsub("\\|.*", "", .) %>% unique(.) %>% grep("synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant", ., value = T) %>% paste(., collapse = ";")})
    infovcfstate$ANN= sapply(infovcfstate$ANN, function(x) {strsplit(x, "&|;") %>% .[[1]] %>% unique(.) %>% grep("synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant", ., value = T) %>% paste(., collapse = ";")})
    infovcfstate$exonic = infovcfstate$ANN != ""
    inactive_bool = apply(infovcfstate[,chromhmmsfetalonly], 2, function(x) x %in% inactivestates)
    inactiverows = rowSums(inactive_bool[,chromhmmsfetalonly]) == 2 & is.na(infovcfstate$Ensemblregbuild) & infovcfstate$exonic == F
    vcf_out = vcf[!inactiverows]
    writeVcf(vcf_out, vcfname_out)
    print(paste("Filtered sample: ", sample, sep = ""))
  }
}
