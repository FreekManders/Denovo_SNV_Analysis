library("BSgenome.Hsapiens.UCSC.hg19", character.only = T)
library(MutationalPatterns)
library(reshape2)
library(gridExtra)
library(dplyr)
library(ggplot2)
library(biomaRt)
library(rtracklayer)
library(optparse)
ref_genome = "BSgenome.Hsapiens.UCSC.hg19"


###Parse command line arguments
option_list = list(
  make_option(c("-f", "--FILELIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", 
              help="The filelist", metavar="character"),
  make_option(c("-o", "--OUTPUT_PATH"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", 
              help="The output path", metavar="character"),
  make_option(c("--CANCER_SIGNATURES"), type="character", default="http://cancer.sanger.ac.uk/cancergenome/assets/signatures_probabilities.txt", 
              help="An url to a file containing the cancer signatures", metavar="character")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

chroms = paste0("chr", 1:22)

###Run mutational patterns function
runmutpatterns = function(vcfs, factors, name){
  
  #Create mutation spectra
  type_occurences = mut_type_occurrences(vcfs, ref_genome)
  pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/MutationalPatterns/", name, "Mutation_spectra.pdf"))
  for (i in 1:ncol(factors)){
    specfigbyfactor = plot_spectrum(type_occurences, CT = T, by = factors[,i])
    print(specfigbyfactor)
  }
  dev.off()
  
  #Create 96 mutation spectra
  mut_mat = mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
  pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/MutationalPatterns/", name, "96Mutation_spectra.pdf"))
  for (i in 1:ncol(factors)){
    mut_mat_byfactor = mut_mat %>% t(.) %>% as.data.frame(.) %>% mutate(myfactor = factors[,i]) %>% group_by(myfactor) %>% summarise_all(funs(sum)) %>% dplyr::select(-myfactor) %>% t(.)
    colnames(mut_mat_byfactor) = levels(factors[,i])
    spec96figbyfactor = (plot_96_profile(mut_mat_byfactor, condensed = T))
    print(spec96figbyfactor)
  }
  dev.off()
  
  #Format the cosmic signatures based on the mutational spectra
  new_order = match(row.names(mut_mat), cancer_signatures_db$Somatic.Mutation.Type)
  cancer_signatures = cancer_signatures_db[as.vector(new_order),]
  row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
  cancer_signatures = as.matrix(cancer_signatures[,4:33])
  hclust_cosmic = cluster_signatures(cancer_signatures, method = "average")
  cosmic_order = colnames(cancer_signatures)[hclust_cosmic$order]
  
  
  plot_cosmics = function(mut_mat, combined){
    #Look at the similarity between mutational spectra and cosmic signatures
    cos_sim_samples_signatures = cos_sim_matrix(mut_mat, cancer_signatures)
    svssigfig = plot_cosine_heatmap(cos_sim_samples_signatures, col_order = cosmic_order, cluster_rows = T)
  
    #Look at contribution of cosmic signatures to mutational spectra
    fit_res = fit_to_signatures(mut_mat, cancer_signatures)
    select = which(rowSums(fit_res$contribution) > 10)
    if (combined == ""){
      maxy = max(colSums(fit_res$contribution[select,])) * 1.1
    } else {
      maxy = sum(fit_res$contribution[select,]) * 1.1
    }
    
    sigcontrifig1 = plot_contribution(fit_res$contribution[select,,drop = F], cancer_signatures[,select], coord_flip = F, mode = "absolute") +
      coord_cartesian(ylim = c(0, maxy), expand = F) +
      theme(axis.text.x = element_text(angle = 85, size = 10, margin = margin(t = 26)))
    sigcontrifig2 = plot_contribution_heatmap(fit_res$contribution, cluster_samples = T, method = "complete")
    
    #Look at similarity between the reconstructed spectra (based on the cosmic signatures) and the original spectra
    cos_sim_ori_rec = cos_sim_matrix(mut_mat, fit_res$reconstructed)
    cos_sim_ori_rec = as.data.frame(diag(cos_sim_ori_rec))
    colnames(cos_sim_ori_rec) = "cos_sim"
    cos_sim_ori_rec$sample = row.names(cos_sim_ori_rec)
    orivsrecfig = ggplot(cos_sim_ori_rec, aes(y = cos_sim, x = sample)) +
      geom_bar(stat = "identity", fill = "skyblue3") +
      coord_cartesian(ylim = c(0.6, 1), expand = F) +
      labs(y = "Cosine similarity\n original VS reconstructed", x = "") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 85, size = 10, margin = margin(t = 25))) +
      geom_hline(aes(yintercept = 0.95))
    
    
    pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/MutationalPatterns/", name, combined, "Cosmic_signatures.pdf"))
    print(svssigfig)
    print(sigcontrifig1)
    print(sigcontrifig2)
    print(orivsrecfig)
    dev.off()
  }
  plot_cosmics(mut_mat = mut_mat, combined = "")
  mut_mat_combined = data.frame("All samples" = rowSums(mut_mat))
  plot_cosmics(mut_mat = mut_mat_combined, combined = "combined_")
  
  #Rainfall plots for each sample
  chromosomes = seqnames(get(ref_genome))[1:22]
  pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/MutationalPatterns/", name, "rainfall_plots.pdf"))
  for (i in 1:length(vcfs)){
    if (sum(duplicated(as.character(seqnames(vcfs[[i]])))) > 0){
      rainfig = plot_rainfall(vcfs[[i]], title = names(vcfs[i]), chromosomes = chromosomes, cex = 2, ylim = 1e+09) +
        coord_cartesian(ylim = c(1, 1e+09))
      print(rainfig)
    }
  }
  dev.off()
  
  # #Distribution regulatory regions. First read in bedfiles and put them in a list.
  # #Bedfiles that are needed multiple times will only be read in once, because this is slow.
  # callable_list = vector("list", length(bed_files))
  # duplibeds = duplicated(bed_files)
  # for (i in 1:length(bed_files)){
  #   bed_file = bed_files[i]
  #   if (duplibeds[i] == "FALSE"){
  #     bed = import(bed_file)
  #     seqlevelsStyle(bed) = "UCSC"
  #     callable_list[[i]] = bed
  #   }
  #   else if (duplibeds[i] == "TRUE"){
  #     importedbed = grep(bed_file, bed_files)[1]
  #     callable_list[[i]] = callable_list[[importedbed]]
  #   }
  # }
  # print("Succeeded in reading in the bedfiles")
  
  #Read bedfiles and calculate genomic distributions
  distr = data.frame(stringsAsFactors = F)
  for (i in 1:length(vcfs)){
    vcf_list = vcfs[i]
    
    bed = import(bed_files[i])
    seqlevelsStyle(bed) = "UCSC"
    bed = bed[as.vector(seqnames(bed)) %in% chroms]
    seqlevels(bed) = chroms
    bed_list = list(bed)
    
    distr_sample = genomic_distribution(vcf_list, bed_list, regregions_list)
    distr = rbind(distr, distr_sample)
    print("Calculated the genomic distribution for one sample")
  }
  print("Finished calculating the genomic distributions")
  
  ##Calculate the actual genomic distribution. This is very slow.
  #distr = genomic_distribution(vcfs, callable_list, regregions_list)
  #print("Calculated the genomic distribution")
  
  #Create the genomic distribution figures
  pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/MutationalPatterns/", name, "Distribution.pdf"))
  for (i in 1:ncol(factors)){
    distr_test = enrichment_depletion_test(distr, by = factors[,i])
    enrdeplfig = plot_enrichment_depletion(distr_test)
    print(enrdeplfig)
  }
  dev.off()
  
  #Create the genomic distribution figures seperately for all samples
  samples = unique(distr$sample) 
  pdf(paste0(opt$OUTPUT_PATH, "/Characteristics_denovo/MutationalPatterns/", name, "Sample_distribution.pdf"))
  out = sapply(samples, function(mysample){
    distrsample = filter(distr, sample == mysample)
    distr_test = enrichment_depletion_test(distrsample)
    enrdeplfig = plot_enrichment_depletion(distr_test)
    print(enrdeplfig)
  })
  dev.off()
}

###Read in the filelist
filelist = read.table(opt$FILELIST, header = T, stringsAsFactors = F, sep = "\t")
if (sum(duplicated(filelist$Sample)) >= 1){
  stop("The file list contains duplicate samples")
}
filelist = filelist %>% filter(grepl("Child|child", Name))
folder = paste0(opt$OUTPUT_PATH, "/Callableregion_NrMendelViols/TrueDenovo/SNVs/")
vcf_files = paste(folder, filelist$Sample, "_Truedenovo.vcf", sep = "")
existing = file.exists(vcf_files)
filelist = filelist[existing,]
vcf_files = vcf_files[existing]

###Create bedfile paths
bedfolder = paste0(opt$OUTPUT_PATH, "/CallableLoci/Output/")
bed_files = paste(bedfolder, filelist$Family, "_intersect_CallableLoci.bed", sep = "")
famsmultikids = unique(filelist$Family[duplicated(filelist$Family)])
for (fam in famsmultikids){
  famrows = grep(fam, bed_files)
  bed_files[famrows] = paste0(bedfolder, fam, "_child", 1:length(famrows), "_intersect_CallableLoci.bed")
}

###Read in the cosmic signatures. This is done via an url for a remote file
cancer_signatures_db = read.table(opt$CANCER_SIGNATURES, sep = "\t", header = TRUE)

###Read in the regulatory regions
regulatory = useEnsembl(biomart = "regulation", dataset = "hsapiens_regulatory_feature", GRCh = 37)
regsites = getBM(attributes = c("chromosome_name", "chromosome_start",  "chromosome_end", "feature_type_name"), mart = regulatory)
features = unique(regsites$feature_type_name)
regregions_list = lapply(features, function(feature){
  regsitesfeature = filter(regsites, feature_type_name == feature)
  Grangeobj = reduce(GRanges(regsitesfeature$chromosome_name, IRanges(regsitesfeature$chromosome_start, regsitesfeature$chromosome_end)))
  return(Grangeobj)
})
regregions_list = GRangesList(regregions_list)
names(regregions_list) = features
seqlevelsStyle(regregions_list) = "UCSC"
for (i in 1:length(features)){
  sites_correct_chrom = as.vector(seqnames(regregions_list[[i]])) %in% chroms
  regregions_list[[i]] = regregions_list[[i]][sites_correct_chrom]
  seqlevels(regregions_list[[i]]) = chroms
}
  
###Read in the vcfs for all mutations and run mutational patterns.
print("Reading in vcfs for all mutations")
vcfs = read_vcfs_as_granges(vcf_files, filelist$Sample, ref_genome)
factors = data.frame(filelist$Sample, filelist$Sex, rep("All samples", length(vcf_files)))

print("Running mutational patterns for all mutations.")
runmutpatterns(vcfs = vcfs, factors = factors, name = "")
print("Finished mutational patterns for all samples")

###Read in the vcfs for phased mutations and run mutational patterns.
vcf_filesfather = paste(folder, "Parentaloriginmuts/", filelist$Sample, "_Truedenovofather.vcf", sep = "")
vcf_filesmother = paste(folder, "Parentaloriginmuts/", filelist$Sample, "_Truedenovomother.vcf", sep = "")
vcf_files = c(vcf_filesfather, vcf_filesmother)
phasedsamples = file.exists(vcf_files)

vcf_files = vcf_files[phasedsamples]
Parent = c(rep("Father", length(vcf_filesfather)), rep("Mother", length(vcf_filesmother)))[phasedsamples]
samplename = rep(filelist$Sample, 2)[phasedsamples]
sex = rep(filelist$Sex, 2)[phasedsamples]
bed_files = rep(bed_files, 2)[phasedsamples]
samplename_parent = paste(samplename, Parent, sep = "_")
sex_parent = paste(sex, Parent, sep = "_")
phasedfactors = data.frame(rep("All samples", length(vcf_files)), Parent, samplename, samplename_parent, sex, sex_parent)
print("Reading in vcfs for phased mutations")
phasedvcfs = read_vcfs_as_granges(vcf_files, samplename_parent, ref_genome)

print("Running mutational patterns for phased mutations")
runmutpatterns(vcfs = phasedvcfs, factors = phasedfactors, name = "phased_")

