library(dplyr)
library(ggplot2)
library(scales)
library(reshape2)
library(ggrepel)
library(Gviz)
library(BSgenome.Hsapiens.UCSC.hg19)
library(VariantAnnotation)
library(karyoploteR)
library(optparse)


####___________________Parse command line arguments____________________________________________####
option_list = list(
  make_option(c("-f", "--FILELIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", 
              help="The filelist", metavar="character"),
  make_option(c("-o", "--OUTPUT_PATH"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", 
              help="The output path", metavar="character"),
  make_option(c("--AGE_PARENTS"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/age_parents.txt", 
              help="A file containing the ages of the parents", metavar="character"),
  make_option(c("--BREAKPOINTS"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Complex_svs/SV/Data/Breakpoints.txt", 
              help="A file containing breakpoints", metavar="character")               
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


####____________________Read filelist___________________________________________________________####
filelist = read.table(opt$FILELIST, header = T, stringsAsFactors = F, sep = "\t")
if (sum(duplicated(filelist$Sample)) >= 1){
  stop("The file list contains duplicate samples")
}
filelist = filelist %>% filter((grepl("Child|child", Name)))
#filelistchild$effectfile = paste(folder, filelistchild$Sample, "_effect.txt", sep = "")
folder = paste0(opt$OUTPUT_PATH, "/Callableregion_NrMendelViols/TrueDenovo/SNVs/")
filelist$vcf_files = paste(folder, filelist$Sample, "_Truedenovo.vcf", sep = "")

####____________________known parent of origin._________________________________________________####
###Create the snv table
snvtable = data.frame(stringsAsFactors = F)
for (i in 1:nrow(filelist)){
  vcfname = filelist$vcf_files[i]
  sample = filelist$Sample[i]
  if (!file.exists(vcfname)){
    next
  }
  vcf = readVcf(vcfname, "hg19")
  chroms = as.vector(seqnames(rowRanges(vcf)))
  starts = start(rowRanges(vcf))
  parents = info(vcf)$PARENT
  if (isEmpty(parents)){
    parents = "not phased"
    print(paste0("Sample ", sample, " has not yet been phased."))
  }
  df_sample = data.frame("Sample" = sample, "chrom" = chroms, "pos" = starts, "Parent" = parents, stringsAsFactors = F)
  snvtable = rbind(snvtable, df_sample)
}
snvtableori = snvtable %>% filter(Parent != "not phased")

###Create figures for the denovo snvs with known parents of origin.
snvtableori1a = snvtableori %>% group_by(Sample, Parent) %>% summarise(count = n()) %>% mutate(percentage = count / sum(count))
originrichfig = ggplot(data = snvtableori1a, aes(x = Parent, fill = Sample)) + 
  geom_bar(aes(y = percentage), position = "dodge", stat = "identity") + 
  scale_y_continuous(labels = percent) + 
  scale_fill_hue(l=45) +
  labs(title = "Parental origin of denovo SNVs", x = "", y = "Percentage of denovo SNVs") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  theme(text = element_text(size = 20))



notphasedorevidence = snvtableori$Parent == "Notenoughreadsandsnps" | snvtableori$Parent == "Notphased"
snvtableori1b = snvtableori
snvtableori1b[notphasedorevidence, "Parent"] = "Not phased properly"
originfig = ggplot(data = snvtableori1b, aes(x = Parent)) + 
  geom_bar(aes(y = (..count..)/sum(..count..)), fill = "darkred") + 
  scale_y_continuous(labels = percent) + 
  labs(title = "Parental origin of denovo SNVs", x = "", y = "Percentage of denovo SNVs") +
  theme_bw() +
  theme(panel.grid.major.x = element_blank()) +
  theme(text = element_text(size = 20))



snvtableori2 = snvtableori %>% group_by(Sample) %>% summarise(Fathercount = sum(Parent == "Father"), Mothercount = sum(Parent == "Mother")) %>% mutate(Ratio = Fathercount / Mothercount)
ttest = t.test(snvtableori2$Fathercount, snvtableori2$Mothercount, paired = T)
pval = round(ttest[[3]], 4)
avgratio = round(mean(snvtableori2$Ratio), 3)
parentsratiofig =  ggplot(data = snvtableori2, aes(y = Ratio, x = "")) +
  geom_boxplot(width = 0.5) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.2) +
  labs(title = "Parental origin of denovo SNVs", x = "", y = "Father / Mother") +
  theme_bw() +
  annotate("text", x = 0.6, y = 20, size = 5, label = paste("Mean: ", avgratio, sep = "")) +
  theme(panel.grid.major.x = element_blank(), text = element_text(size = 20))



snvtableori3 = snvtableori2 %>% dplyr::select(Father = Fathercount, Mother = Mothercount, Sample) %>% melt(id = "Sample", variable.name = "Parent", value.name = "count")
perparentrichfig = ggplot(data = snvtableori3, aes(x = Parent, y = count)) +
  geom_boxplot(width = 0.5) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.6) +
  geom_text_repel(aes(label = Sample), box.padding = 0.5, size = 3) +
  labs(title = "Parental origin of denovo SNVs", y = "Nr. of denovo SNVs") + 
  theme_bw() +
  annotate("text", x = 0.6, y = 24, size = 5, label = paste("p: ", pval, sep = "")) +
  theme(text = element_text(size = 20))

perparentfig = ggplot(data = snvtableori3, aes(x = Parent, y = count)) +
  geom_boxplot(width = 0.5) +
  geom_dotplot(binaxis = "y", stackdir = "center", binwidth = 0.6) +
  labs(title = "Parental origin of denovo SNVs", y = "Nr. of denovo SNVs") + 
  theme_bw() +
  annotate("text", x = 0.5, y = 33, size = 5, label = paste("p: ", pval, sep = "")) +
  theme(text = element_text(size = 20))

pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/Parentalorigin.pdf"), paper = "a4")
originrichfig
originfig
parentsratiofig
perparentrichfig
perparentfig
dev.off()


####____________________The number of snvs per sample combined with the age of the parents______####
###Combine data
snvcounts = snvtable %>% group_by(Sample) %>% summarise(count = n())
ageparents = read.table(opt$AGE_PARENTS, sep = "\t", stringsAsFactors = F, header = T)
ageparents[,c("age_father", "age_mother")] = apply(ageparents[,c("age_father", "age_mother")], 2, function(x) as.numeric(gsub(",", ".", x)))
snvcounts = inner_join(snvcounts, ageparents, by = c("Sample" = "ID"))
corparentsage = cor(snvcounts$age_father, snvcounts$age_mother)

###Plot the number of snvs per Sample versus the age of the parents
mutsbyage = function(parent){
  ylims = c(min(snvcounts$count), max(snvcounts$count))
  xlims = c(min(snvcounts[,c("age_father", "age_mother")]), max(snvcounts[,c("age_father", "age_mother")]))
  parentcounts = snvcounts %>% dplyr::select(Sample, count, age = parent)
  lmmodel = summary(lm(count ~ age, data = parentcounts))
  pval = round(lmmodel$coefficients[8], 3)
  a = round(lmmodel$coefficients[2], 3)
  b = round(lmmodel$coefficients[1], 3)
  stderrora = round(lmmodel$coefficients[4], 3)
  stderrorb = round(lmmodel$coefficients[3], 3)
  rsq = round(lmmodel$r.squared, 3)
  
  if (parent == "age_father"){
    mytitle = "Father"
  }
  else if (parent == "age_mother"){
    mytitle = "Mother"
  }
  text.x = xlims[1] + (0.00 * (xlims[2] - xlims[1]))
  text.y = ylims[2] - (0.07 * (ylims[2] - ylims[1]))
  fig = ggplot(parentcounts, aes(x = age, y = count)) +
    geom_point(color = "red") + 
    geom_smooth(method = "lm") +
    coord_cartesian(ylim = ylims, xlim = xlims) +
    labs(title = mytitle, x = "Age", y = "Nr. Denovo SNVs") +
    theme_bw() +
    annotate("text", x = text.x, y = text.y, label = paste0("p: ", pval, "\nR^2: ", rsq, "\ny = ", a, " +/- ", stderrora, " x", " + ", b, " +/- ", stderrorb), size = 5, hjust = 0) +
    theme(text = element_text(size=20))
  return(fig)
}

mutsbyfatheragefig = mutsbyage("age_father")
mutsbymotheragefig = mutsbyage("age_mother")


###plot the number of phased snvs per sample vs the age of the parents
snvcounts2 = snvtable %>% filter(Parent == "Mother" | Parent == "Father") %>% group_by(Sample, Parent) %>% summarise(count = n())
snvcounts2 = inner_join(snvcounts2, ageparents, by = c("Sample" = "ID"))
snvcounts2$age = ifelse(snvcounts2$Parent == "Father", snvcounts2$age_father, snvcounts2$age_mother)

stats = data.frame(stringsAsFactors = F) #calculate pvalue and y~x formula.
for (parent in c("Father", "Mother")){
  snvcountsparent = snvcounts2[snvcounts2$Parent == parent,]
  lmmodel = summary(lm(count ~ age, data = snvcountsparent))
  pval = round(lmmodel$coefficients[8], 3)
  a = round(lmmodel$coefficients[2], 3)
  b = round(lmmodel$coefficients[1], 3)
  stderrora = round(lmmodel$coefficients[4], 3)
  stderrorb = round(lmmodel$coefficients[3], 3)
  rsq = round(lmmodel$r.squared, 3)
  df = data.frame("parent" = parent, "pval" = pval, "a" = a, "stderrora" = stderrora, "b" = b, "stderrorb" = stderrorb, "rsq" = rsq)
  stats = rbind(stats, df)
}
xlims = c(min(snvcounts2$age), max(snvcounts2$age))
ylims = c(min(snvcounts2$count) - 5, max(snvcounts2$count) + 5)
text.x = xlims[1] + (0.45 * (xlims[2] - xlims[1]))
text.y = ylims[2] - (0.04 * (ylims[2] - ylims[1]))

phasedmutsbyagefig = ggplot(snvcounts2, aes(x = age, y = count, fill = Parent, color = Parent)) +
  geom_point() + 
  geom_smooth(method = "lm") +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  scale_fill_manual(values = c("Father" = "grey", "Mother" = "grey")) +
  scale_color_discrete(l = 40) +
  theme_bw() +
  labs(y = "Nr. of phased denovo snvs", x = "age", title = "Effect of age parents on denovo snvs") +
  annotate("text", x = text.x, y = text.y, size = 4, label = paste0("Father:  p: ", stats$pval[1], "  R^2 = ", stats$rsq[1], "   y = ", stats$a[1], " +/- ", stats$stderrora[1], " x + ", stats$b[1], " +/- ", stats$stderrorb[1], "\nMother:  p: ", stats$pval[2], "  R^2 = ", stats$rsq[2], "   y = ", stats$a[2], " +/- ", stats$stderrora[2], " x + ", stats$b[2], " +/- ", stats$stderrorb[2])) +
  theme(text = element_text(size = 20))

pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/mutsbyage.pdf"), width = 8)
mutsbyfatheragefig
mutsbymotheragefig
phasedmutsbyagefig
dev.off()


####____________________Location of the denovo snvs_____________________________________________####
breakpoints <- read.delim(opt$BREAKPOINTS, header = T)

plotsnvs = function(snvtablesample, breakpointssample, sample){
  pp = getDefaultPlotParams(plot.type=2)
  pp$ideogramheight = 15
  pp$data1height = 100
  pp$data2height = 100
  pp$data1outmargin = 25
  pp$data2outmargin = 25
  
  title = paste("SV breakpoints per mb cohort and denovo snvs from sample: ", sample, sep = "")
  kp = plotKaryotype(main = title, plot.type = 2, plot.params = pp, labels.plotter = NULL)
  kp = kpDataBackground(kp, data.panel=1)
  kp = kpDataBackground(kp, data.panel=2)
  kp = kpAddChromosomeNames(kp, offset(1.8))
  kp = kpAddLabels(kp, labels = "SVs(1mb)", cex = 0.6, label.margin = 0.025, data.panel = 1)
  kp = kpAddLabels(kp, labels = "SNVs(1kb)", cex = 0.6, label.margin = 0.025, data.panel = 2)
  kp = kpAxis(kp, ymin = 0, ymax = 10, numticks = 2, data.panel = 1, cex = 0.5, side = 1)
  kp = kpAxis(kp, ymin = 0, ymax = 10, numticks = 2, data.panel = 1, cex = 0.5, side = 2)
  kp = kpAxis(kp, ymin = 0, ymax = 6, numticks = 2, data.panel = 2, cex = 0.5, side = 1)
  kp = kpAxis(kp, ymin = 0, ymax = 6, numticks = 2, data.panel = 2, cex = 0.5, side = 2)
  
  if (nrow(breakpointssample) >= 1){
    breakpoints_g = GRanges(seqnames = paste("chr", breakpointssample$Chr, sep = ""), IRanges(start = breakpointssample$Breakpoint, end = breakpointssample$Breakpoint+1))
    kp = kpPlotDensity(kp, data = breakpoints_g, border = "black", col = "black", window.size = 1e6, ymax = 10)
  }
  snvs_g = makeGRangesFromDataFrame(data.frame(chr = paste("chr", snvtablesample$chrom, sep = ""), start = snvtablesample$pos, end = snvtablesample$pos + 1))
  kp = kpPlotDensity(kp, data = snvs_g, window.size = 1e3, data.panel = 2, ymax = 6, border = "black", col = "black")
  return(kp)
}

pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/vcfannotations/Plotsnvs_location.pdf"), width = 16, height = 10)
for (i in 1:nrow(filelist)){
  sample = filelist$Sample[i]
  snvtablesample = snvtable %>% filter(Sample == sample)
  breakpointssample = breakpoints %>% filter(Patient == sample)
  if (nrow(snvtablesample) >= 1){
	  plotsnvs(snvtablesample, breakpointssample, sample)
	  print("Finished a sample")
  }
}
kp = plotsnvs(snvtable, breakpoints, "allsamples")
dev.off()
print("Done plotting the snv locations")
