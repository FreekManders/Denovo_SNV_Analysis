###This script reads in vcfs with denovo snvs (that have been filtered on frequency).
###It then measures how many genes (ensembl annotated) are near these snvs.
###Different distances around the snvs are measured.


library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(GenomicRanges)
library(VariantAnnotation)
library(karyoploteR)
library(optparse)

###Parse command line arguments
option_list = list(
  make_option(c("-f", "--FILELIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", 
              help="The filelist", metavar="character"),
  make_option(c("-o", "--OUTPUT_PATH"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", 
              help="The output path", metavar="character"),
  make_option(c("--GENE_LIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Complex_svs/Common_data/Genes/Genelist_HGNC_v4.txt", 
              help="A list of known genes, with annotations like PLi", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###Read in the filelist
filelist = read.table(opt$FILELIST, header = T, stringsAsFactors = F, sep = "\t")
if (sum(duplicated(filelist$Sample)) >= 1){
  stop("The file list contains duplicate samples")
}
filelist = filelist %>% filter(grepl("Child|child", Name))
folder = paste0(opt$OUTPUT_PATH, "/Possibly_pathogenic/SNVs/freq_chromfunction/")
filelist$vcf_files_pathogenic = paste(folder, filelist$Sample, "_freq_chromfunction_filter.vcf", sep = "")


genes = read.delim(opt$GENE_LIST)
genes$strand = mapvalues(genes$strand, from = c("1", "-1"), to = c("+", "-"))
genesgr = makeGRangesFromDataFrame(genes, keep.extra.columns = T, seqnames.field = "chromosome_name", start.field = "start_position", end.field = "end_position")
genesgr = genesgr[seqnames(genesgr) %in% 1:22]

###Determine overlap between snvs and genes
nr_overlaps = function(extensions, log_scale, extension_names){
  
  ###Read in vcfs and determine number of nearby genes.
  allvcfgrmatched = GRanges()
  nrgenesdf = data.frame(stringsAsFactors = F)
  snvsdf = data.frame(stringsAsFactors = F)
  for (i in 1:nrow(filelist)){
    vcfname = filelist$vcf_files_pathogenic[i]
    sample = filelist$Sample[i]
    if (file.exists(vcfname)){
      vcf = readVcf(vcfname, "hg19")
      vcfgr = rowRanges(vcf)
      for (extend in extensions){
        vcfgr2 = vcfgr + extend
        mtch = as.matrix(findOverlaps(vcfgr2, genesgr))
        nrsnvs = length(vcfgr2)
        
        #Nr genes per family (divided by nr. snvs)
        matchedgenesgr = genesgr[mtch[,2]]
        nrgenes = length(unique(mcols(matchedgenesgr)$ensembl_gene_id))
        nrgenesbysnv = nrgenes / nrsnvs
        sampledf = data.frame("sample" = sample, "extension" = extend, "nrgenes" = nrgenes, "nrgenesbysnv" = nrgenesbysnv, stringsAsFactors = F)
        nrgenesdf = rbind(nrgenesdf, sampledf)
        
        #counts per snv
        snvseq = seq(1, nrsnvs)
        counts = as.data.frame(table(mtch[,1]), stringsAsFactors = F)
        if (nrow(counts) != nrsnvs){
          notcounted = (snvseq)[!(snvseq %in% counts$Var1)]
          noncounts = data.frame("Var1" = notcounted, "Freq" = 0)
          counts = rbind(counts, noncounts)
        }
        counts$sample = sample
        counts$extension = extend
        snvsdf = rbind(snvsdf, counts)
        }
      }
  }
  
  ###Create figures showing the number of the genes that are near snvs.
  nrcolors = length(unique(nrgenesdf$sample))
  set.seed(001)
  colors = sample(colorRampPalette(brewer.pal(8, "Set1"))(nrcolors))
  
  nrgenesfig = ggplot(data = nrgenesdf, aes(x = extension, y = nrgenes, color = sample)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = colors) +
    labs(x = "Size around snv (bp)", y = "Unique nearby genes", title = "Genes near denovo SNVs") +
    theme_bw()
  if (log_scale == T){
    nrgenesfig = nrgenesfig + scale_x_log10()
    }
  
  nrgenesbysnvfig = ggplot(data = nrgenesdf, aes(x = extension, y = nrgenesbysnv, color = sample)) +
    geom_point() +
    geom_line() +
    scale_color_manual(values = colors) +
    labs(x = "Size around snv (bp)", y = "Mean unique nearby genes per snv", title = "Genes near denovo SNVs divided by the number of snvs") +
    theme_bw()
  if (log_scale == T){
    nrgenesbysnvfig = nrgenesbysnvfig + scale_x_log10()
  }
  
  if (log_scale == T){
    snvsdf$extensiondisplay = as.factor(log10(snvsdf$extension))
    labfill = "Max distance (log10)"
  } else {
    snvsdf$extensiondisplay = as.factor(snvsdf$extension)
    labfill = "Max distance"
  }
  nrgenespersnvfig = ggplot(data = snvsdf, aes(x = Freq, fill = extensiondisplay)) +
    facet_grid(extensiondisplay ~ .) +
    geom_histogram(binwidth = 1) +
    coord_cartesian(xlim = c(0,100), ylim = c(0,500), expand = F) +
    theme_bw() +
    labs(fill = labfill, x = "Nr. of nearby genes", y = "Nr. of snvs", title = "Number of snvs with the number of nearby genes using different distances") +
    scale_fill_hue(l=40) +
    theme(panel.grid.minor.x = element_blank())
  
  pdfname = paste0(opt$OUTPUT_PATH, "/Possibly_pathogenic/SNVs/Near_gene/snv_gene_overlap_", extension_names, ".pdf")
  pdf(pdfname)
  print(nrgenesfig)
  print(nrgenesbysnvfig)
  print(nrgenespersnvfig)
  dev.off()
}


extensions_log = 10^seq(1, 6, by = 0.5)
extensions_100kb = 100000*(1:10)

nr_overlaps(extensions = extensions_log, log_scale = T, extension_names = "log10_steps")
nr_overlaps(extensions = extensions_100kb, log_scale = F, extension_names = "100kb_steps")


###Look at snvs that aren't matched to any genes.
nonmatched = GRanges()
for (i in 1:nrow(filelist)){
  vcfname = filelist$vcf_files_pathogenic[i]
  sample = filelist$Sample[i]
  if (file.exists(vcfname)){
    vcf = readVcf(vcfname, "hg19")
    vcfgr = rowRanges(vcf)
    extend = 2e+05
    vcfgr2 = vcfgr + extend
    mtch = as.matrix(findOverlaps(vcfgr2, genesgr))
    vcfgrnotmatched = vcfgr[-unique(mtch[,1])]
    if (length(vcfgrnotmatched) != 0){
      mcols(vcfgrnotmatched)$sample = sample
      infovcf = info(vcf)[-unique(mtch[,1]),c("DANN"), drop = F]
      mcols(vcfgrnotmatched) = cbind(mcols(vcfgrnotmatched), infovcf)
      nonmatched = c(nonmatched, vcfgrnotmatched)
    }
  }
}

seqlevelsStyle(nonmatched) <- "UCSC"
pdf(paste0(opt$OUTPUT_PATH, "/Possibly_pathogenic/SNVs/Near_gene/nonmatchedsnvs.pdf"))
kp = plotKaryotype()
kpPlotDensity(kp, data = nonmatched, window.size = 1e3, ymax = 3)
dev.off()
