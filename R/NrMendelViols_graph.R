library(ggplot2)
library(reshape2)
library(dplyr)
library(optparse)

####________Parse command line arguments_____________________________________________________________####
option_list = list(
  make_option(c("-f", "--FILELIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", 
              help="The filelist", metavar="character"),
  make_option(c("-o", "--OUTPUT_PATH"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", 
              help="The output path", metavar="character"),
  make_option(c("--FAMS3KIDS"), type="character", default="MP14", 
              help="Families with 3 kids", metavar="character"),
  make_option(c("--COMPAREPRIOR"), type="character", default="true", 
              help="Whether or not to compare the results with a different PBT prior used in a different run.", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);



filelist_fname = opt$FILELIST
output_path = opt$OUTPUT_PATH
famMultChild = strsplit(opt$FAMS3KIDS, ",")[[1]]
famMultChild = c("MP14")
comparePrior = opt$COMPAREPRIOR
folder = paste0(output_path, "/Callableregion_NrMendelViols")
pathogenicfolder = paste0(output_path, "/Possibly_pathogenic/SNVs")

####____________Read the file list___________________________________________####

filelist = read.table(filelist_fname, header = T, stringsAsFactors = F, sep = "\t")
if (sum(duplicated(filelist$Sample)) >= 1){
  stop("The file list contains duplicate samples")
}
filelistchildren = filelist %>% filter(grepl("Child|child", Name))


####____________Create the filepaths_________________________________________####

vcfdirs = c()
for (i in 1:nrow(filelistchildren)){
  vcfdir = paste(filelistchildren[i,"Folder"], "/vcf_new_filter", sep = "")
  if (!dir.exists(vcfdir)){
    vcfdir = paste(filelistchildren[i,"Folder"], "/vcf", sep = "")
  }
  vcfdirs = c(vcfdirs, vcfdir)
}
filelistchildren$vcfdir = vcfdirs

MendelViol = paste(filelistchildren$vcfdir, "/", filelistchildren$Family, "/", filelistchildren$Family, ".MendelViol", sep = "")
CallableMendelViol = paste(folder, "/PBT/", filelistchildren$Family, ".MendelViol", sep = "")
CallableMendelVioldefaultprior = paste(folder, "/PBT/Defaultprior/", filelistchildren$Family, ".MendelViol", sep = "")
MendelViolsfullgt = paste(folder, "/PBT/", filelistchildren$Family, "_fullgenotypes.MendelViol", sep = "")
MendelViolsfullgtdefaultprior = paste(folder, "/PBT/Defaultprior/", filelistchildren$Family, "_fullgenotypes.MendelViol", sep = "")
MendelViolssnv = paste(folder, "/Denovo/SNVs/", filelistchildren$Sample, "_denovo_snv.recode.vcf", sep = "")
MendelViolssnvdefaultprior = paste(folder, "/Denovo/Defaultprior/SNVs/", filelistchildren$Sample, "_denovo_snv.recode.vcf", sep = "")
MendelViolssnvtrue = paste(folder, "/TrueDenovo/SNVs/", filelistchildren$Sample, "_Truedenovo.vcf", sep = "")
MendelViolssnvtruedefaultprior = paste(folder, "/TrueDenovo/SNVs/Defaultprior/", filelistchildren$Sample, "_Truedenovo.vcf", sep = "")
Filteredonfreq = paste(pathogenicfolder, "/Filteredonfreq/", filelistchildren$Sample, "_freqfilter.vcf", sep = "")
freq_chromfunction = paste(pathogenicfolder, "/freq_chromfunction/", filelistchildren$Sample, "_freq_chromfunction_filter.vcf", sep = "")
vcfbefore = paste(filelistchildren$vcfdir, "/", filelistchildren$Family, "/", filelistchildren$Family, ".filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.vcf", sep = "")
vcfafter = paste(folder, "/CallableVCF/", filelistchildren$Sample, "_filtered_callable_annotated.recode.vcf", sep = "")

#Create specific filepaths for families with multiple children
for (fam in famMultChild){
  famrows = grep(fam, filelistchildren$Family)
  famdir = paste0(folder, "/", fam)
  nrchildren = length(grep("filtered_variants_snpEff_snpSift_Cosmicv76_GoNLv5.recode.vcf" ,list.files(famdir)))
  childSeq = 1:nrchildren
  CallableMendelViol[famrows] = paste0(famdir, "/", fam, "_child", childSeq ,".MendelViol")
  CallableMendelVioldefaultprior[famrows] = paste0(famdir, "/Defaultprior/", fam, "_child", childSeq ,".MendelViol")
  MendelViolsfullgt[famrows] = paste0(famdir, "/", fam, "_child", childSeq, "_fullgenotypes.MendelViol", sep = "")
  MendelViolsfullgtdefaultprior[famrows] = paste0(famdir, "/Defaultprior/", fam, "_child", childSeq ,"_fullgenotypes.MendelViol")
  vcfafter[famrows] = paste0(famdir, "/", fam, "_child", childSeq, "_filtered_callable_annotated.recode.vcf")
}
counttable = data.frame(stringsAsFactors = F)

####____________Get number of sites at different points in the pipeline______####
iapviols = c()
for (i in 1:nrow(filelistchildren)){
  if (file.exists(MendelViol[i])){
    command = paste("wc -l <", MendelViol[i])
    iapviol = system(command, intern = T)
  }
  else{iapviol = 0}
  iapviols = c(iapviols, iapviol)
}
iapviols = as.numeric(iapviols)


callableviols = c()
for (i in 1:nrow(filelistchildren)){
  if (file.exists(CallableMendelViol[i])){
    command = paste("wc -l <", CallableMendelViol[i])
    callableviol = system(command, intern = T)
  }
  else{callableviol = 0}
  callableviols = c(callableviols, callableviol)
}
callableviols = as.numeric(callableviols)


fullgtviols = c()
for (i in 1:nrow(filelistchildren)){
  if (file.exists(MendelViolsfullgt[i])){
    command = paste("wc -l <", MendelViolsfullgt[i])
    fullgtviol = system(command, intern = T)
  }
  else{fullgtviol = 0}
  fullgtviols = c(fullgtviols, fullgtviol)
}
fullgtviols = as.numeric(fullgtviols)


snvviols = c()
for (i in 1:nrow(filelistchildren)){
  if (file.exists(MendelViolssnv[i])){
    command = paste("grep -v '^#'", MendelViolssnv[i], "| wc -l")
    snvviol = system(command, intern = T)
  }
  else{snvviol = 0}
  snvviols = c(snvviols, snvviol)
}
snvviols = as.numeric(snvviols)


snvviolstrue = c()
for (i in 1:nrow(filelistchildren)){
  if (file.exists(MendelViolssnvtrue[i])){
    command = paste("grep -v '^#'", MendelViolssnvtrue[i], "| wc -l")
    snvvioltrue = system(command, intern = T)
  }
  else{snvvioltrue = 0}
  snvviolstrue = c(snvviolstrue, snvvioltrue)
}
snvviolstrue = as.numeric(snvviolstrue)


pathos_freqf = c()
for (i in 1:nrow(filelistchildren)){
  if (file.exists(Filteredonfreq[i])){
    command = paste("grep -v '^#'", Filteredonfreq[i], "| wc -l")
    patho_freq_chromfunc_f = system(command, intern = T)
  }
  else{patho_freq_chromfunc_f = 0}
  pathos_freqf = c(pathos_freqf, patho_freq_chromfunc_f)
}
pathos_freqf = as.numeric(pathos_freqf)

pathos_freq_chromfunc_f = c()
for (i in 1:nrow(filelistchildren)){
  if (file.exists(freq_chromfunction[i])){
    command = paste("grep -v '^#'", freq_chromfunction[i], "| wc -l")
    patho_freq_chromfunc_f = system(command, intern = T)
  }
  else{patho_freq_chromfunc_f = 0}
  pathos_freq_chromfunc_f = c(pathos_freq_chromfunc_f, patho_freq_chromfunc_f)
}
pathos_freq_chromfunc_f = as.numeric(pathos_freq_chromfunc_f)



#This part is slow
sitesbefore = c()
for (i in 1:nrow(filelistchildren)){
  vcfbeforesample = paste(vcfbefore[i], ".gz", sep = "")
  if (file.exists(vcfbefore[i])){
    command = paste("grep -v '^#'", vcfbefore[i], "| wc -l")
    sitebefore = system(command, intern = T)
  }
  else if (file.exists(vcfbeforesample)){
    command = paste("zgrep -v '^#'", vcfbeforesample, "| wc -l")
    sitebefore = system(command, intern = T)
  }
  else{sitebefore = 0}
  sitesbefore = c(sitesbefore, sitebefore)
}
sitesbefore = as.numeric(sitesbefore)

sitesafter = c()
for (i in 1:nrow(filelistchildren)){
  if (file.exists(vcfafter[i])){
    command = paste("grep -v '^#'", vcfafter[i], "| wc -l")
    siteafter = system(command, intern = T)
  }
  else{siteafter = 0}
  sitesafter = c(sitesafter, siteafter)
}
sitesafter = as.numeric(sitesafter)

pathosfiltered = snvviolstrue - pathos_freqf
pathosfiltered = pathosfiltered[pathosfiltered != 0]
avgpathosfiltered = mean(pathosfiltered)

counttable = cbind(filelistchildren, iapviols, callableviols, fullgtviols, snvviols, snvviolstrue, pathos_freqf, pathos_freq_chromfunc_f, sitesbefore, sitesafter)

####____________Compare to the default mutation prior________________________####
if (comparePrior == "true"){
  callableviolsdefaultprior = c()
  for (i in 1:nrow(filelistchildren)){
    if (file.exists(CallableMendelVioldefaultprior[i])){
      command = paste("wc -l <", CallableMendelVioldefaultprior[i])
      callablevioldefaultprior = system(command, intern = T)
    }
    else{callablevioldefaultprior = 0}
    callableviolsdefaultprior = c(callableviolsdefaultprior, callablevioldefaultprior)
  }
  callableviolsdefaultprior = as.numeric(callableviolsdefaultprior)
  
  fullgtviolsdefaultprior = c()
  for (i in 1:nrow(filelistchildren)){
    if (file.exists(MendelViolsfullgtdefaultprior[i])){
      command = paste("wc -l <", MendelViolsfullgtdefaultprior[i])
      fullgtvioldefaultprior = system(command, intern = T)
    }
    else{fullgtvioldefaultprior = 0}
    fullgtviolsdefaultprior = c(fullgtviolsdefaultprior, fullgtvioldefaultprior)
  }
  fullgtviolsdefaultprior = as.numeric(fullgtviolsdefaultprior)
  
  snvviolsdefaultprior = c()
  for (i in 1:nrow(filelistchildren)){
    if (file.exists(MendelViolssnvdefaultprior[i])){
      command = paste("grep -v '^#'", MendelViolssnvdefaultprior[i], "| wc -l")
      snvvioldefaultprior = system(command, intern = T)
    }
    else{snvvioldefaultprior = 0}
    snvviolsdefaultprior = c(snvviolsdefaultprior, snvvioldefaultprior)
  }
  snvviolsdefaultprior = as.numeric(snvviolsdefaultprior)
  
  snvviolstruedefaultprior = c()
  for (i in 1:nrow(filelistchildren)){
    if (file.exists(MendelViolssnvtruedefaultprior[i])){
      command = paste("grep -v '^#'", MendelViolssnvtruedefaultprior[i], "| wc -l")
      snvvioltruedefaultprior = system(command, intern = T)
    }
    else{snvvioltruedefaultprior = 0}
    snvviolstruedefaultprior = c(snvviolstruedefaultprior, snvvioltruedefaultprior)
  }
  snvviolstruedefaultprior = as.numeric(snvviolstruedefaultprior)
  
  counttable = cbind(counttable, callableviolsdefaultprior, fullgtviolsdefaultprior, snvviolsdefaultprior, snvviolstruedefaultprior)
}

####____________Write the results in a table_________________________________####
time = format(Sys.time(), "%H:%M_%d-%m-%y")
write.table(counttable, paste0(output_path, "/NrMendelViols/NrMendelViolstable_", time, ".txt"), row.names = F, quote = F, sep = "\t")

####____________Visualize the results________________________________________####
NrMendelViols = counttable
NrMendelViols$Run = as.factor(NrMendelViols$Folder)
levels(NrMendelViols$Run) = paste("Run: ", 1:length(levels(NrMendelViols$Run)), sep = "")

NrMendelViols2 = NrMendelViols %>% select(Sample, Run, iapviols, callableviols, fullgtviols, snvviols, snvviolstrue) %>% melt(id = c("Sample", "Run"))
levels(NrMendelViols2$variable) = c("IAP", "Callable regions", "Full GT", "SNV", "True denovo SNV")
NrMendelViols2$value[NrMendelViols2$value == 0] = 1
overview_fig = ggplot(NrMendelViols2, aes(y = value, x = Sample, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(. ~ Run, scales = "free", space = "free") + 
  scale_y_log10() + 
  theme_bw() + 
  labs(y = "Mendelian violations", x = "Sample", title = "Mendelian violations") + 
  coord_cartesian(ylim = c(10, 100000), expand = F) + 
  theme(axis.text.x = element_text(angle = 80, size = 10, margin = margin(t = 25)), legend.position = "bottom", panel.grid.major.x = element_blank(), text = element_text(size=20)) + 
  scale_fill_discrete(name = "Regions")

NrMendelViols2_5 = NrMendelViols %>% select(Sample, Run, snvviolstrue, pathos_freqf, pathos_freq_chromfunc_f) %>% melt(id = c("Sample", "Run"))
levels(NrMendelViols2_5$variable) = c("True denovo SNV", "Filtered on gnomAD, PONN and GoNL", "Filtered out inactive sites")
maxy2_5 = 1.1 * max(NrMendelViols2_5$value)
prioritize_fig = ggplot(NrMendelViols2_5, aes(y = value, x = Sample, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(. ~ Run, scales = "free", space = "free") + 
  theme_bw() + 
  labs(y = "Mendelian violations", x = "Sample", title = "Mendelian violations") + 
  coord_cartesian(ylim = c(0, maxy2_5), expand = F) + 
  theme(axis.text.x = element_text(angle = 80, size = 10, margin = margin(t = 25)), legend.position = "bottom", panel.grid.major.x = element_blank(), text = element_text(size=20)) + 
  scale_fill_discrete(name = "Regions")

NrMendelViols3 = NrMendelViols %>% select(Sample, Run, sitesbefore, sitesafter) %>% melt(id = c("Sample", "Run"))
levels(NrMendelViols3$variable) = c("All regions", "Callable regions")
maxy3 = 1.1 * max(NrMendelViols3$value)
callsites_fig = ggplot(NrMendelViols3, aes(y = value, x = Sample, fill = variable)) + 
  geom_bar(stat = "identity", position = "dodge") + 
  facet_grid(. ~ Run, scales = "free", space = "free") + 
  theme_bw() + 
  labs(y = "Sites", x = "Sample", title = "Number of sites in unfiltered/filtered vcfs") + 
  coord_cartesian(ylim = c(0, maxy3), expand = F) + 
  theme(axis.text.x = element_text(angle = 80, size = 10, margin = margin(t = 25)), panel.grid.major.x = element_blank(), text = element_text(size=20)) + 
  scale_fill_discrete(name = "Regions")

if (comparePrior == "true"){
  NrMendelViols4 = NrMendelViols %>% select(Run, Sample, snvviols, snvviolsdefaultprior) %>% melt(id = c("Sample", "Run"))
  levels(NrMendelViols4$variable) = c("1.0E-04", "1.0E-08")
  maxy4 = 1.1 * max(NrMendelViols4$value)
  prior_fig = ggplot(NrMendelViols4, aes(y = value, x = Sample, fill = variable)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    facet_grid(. ~ Run, scales = "free", space = "free") + 
    theme_bw() + 
    labs(y = "Nr of denovo SNVs", x = "Sample", title = "Number of denovo snvs with different priors") + 
    coord_cartesian(ylim = c(0, maxy4), expand = F) + 
    theme(axis.text.x = element_text(angle = 80, size = 10, margin = margin(t = 25)), panel.grid.major.x = element_blank(), text = element_text(size=20)) + 
    scale_fill_discrete(name = "Prior")
  
  NrMendelViols5 = NrMendelViols %>% select(Run, Sample, snvviolstrue, snvviolstruedefaultprior) %>% melt(id = c("Sample", "Run"))
  levels(NrMendelViols5$variable) = c("1.0E-04", "1.0E-08")
  maxy5 = 1.1 * max(NrMendelViols5$value)
  prior_f_fig = ggplot(NrMendelViols5, aes(y = value, x = Sample, fill = variable)) + 
    geom_bar(stat = "identity", position = "dodge") + 
    facet_grid(. ~ Run, scales = "free", space = "free") + 
    theme_bw() + 
    labs(y = "Nr of denovo SNVs", x = "Sample", title = "Number of denovo snvs with different priors after filtering") + 
    coord_cartesian(ylim = c(0, maxy5), expand = F) + 
    theme(axis.text.x = element_text(angle = 80, size = 10, margin = margin(t = 25)), panel.grid.major.x = element_blank(), text = element_text(size=20)) + 
    scale_fill_discrete(name = "Prior")
}

pdf(paste0(output_path, "/NrMendelViols/NrMendelViols_", time, ".pdf"), width = 15)
overview_fig
prioritize_fig
callsites_fig
prior_fig
prior_f_fig
dev.off()
