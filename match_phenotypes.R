library(dplyr)
library(GenomicRanges)
library(VariantAnnotation)
library(gdata)
library(ggplot2)
library(scales)
library(reshape2)
library(optparse)
library(RColorBrewer)

####________Parse command line arguments_____________________________________________________________####
option_list = list(
  make_option(c("-f", "--FILELIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", 
              help="The filelist", metavar="character"),
  make_option(c("--OVERWRITE"), type="character", default="True", 
              help="Whether or not to overwrite existing data files", metavar="character"),
  make_option(c("-o", "--OUTPUT_PATH"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", 
              help="The output path", metavar="character"),
  make_option(c("--HPO_PATIENTS_EXCELL"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Complex_svs/Patients/Phenotype_overview2_Freek.xlsx", 
              help="The excell file containing the hpo terms of the patients", metavar="character"),
  make_option(c("--EXTRA_PHENOTYPE"), type="character", default="HP:0000006", 
              help="HPO term added to all patients", metavar="character"),
  make_option(c("--PLOT"), type="character", default="True", 
              help="Whether or not to plot", metavar="character"),
  make_option(c("--GENES2PHENO"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Complex_svs/Common_data/Genes/HPO/20171211_ALL_SOURCES_ALL_FREQUENCIES_genes_to_phenotype.txt", 
              help="A genes to pheno type file from the HPO consortium", metavar="character"),
  make_option(c("--OVERLAP_GENE_SNV_DIST"), type="integer", default=500000, 
              help="area around the snv which will be overlapped with genes."),
  make_option(c("--PHENOMATCH"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Complex_svs/Common_scripts/Phenomatch/bin/phenomatch.jar", 
              help="Location of phenomatch", metavar="character"),
  make_option(c("--KNOWNGENES_ENTREZ"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Complex_svs/Common_scripts/Phenomatch/data/knownGene.txt.entrez_id.tab.unique", 
              help="File containing known entrez genes", metavar="character"),
  make_option(c("--HPO_OBO"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Complex_svs/Common_data/Genes/HPO/20171211_HPO.obo", 
              help="HPO OBO file", metavar="character")
  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

####________Get the phenotypes from the patients from excell into a flat file________________________####
read_phenotype_from_excell = function(fname, sheetname){
  
  #Read excel file. This can be either a xls or a xlsx file.
  phenotypes = read.xls(fname, sheet = sheetname, stringsAsFactors = F)
  
  #The first row contains columns with subdivisions of the "Abnormality of the nervous system" column.
  #These subdivisions are added to the column names and the first row is removed. Any empty columns are also removed.
  nervous_sub_cols = phenotypes[1,] != ""
  nervous_sub_cols[is.na(nervous_sub_cols)] = F
  nervous_sub_colnames = paste("Nervous_system:",phenotypes[1,nervous_sub_cols], sep = "")
  phenotypes = phenotypes[-1,]
  colnames(phenotypes)[nervous_sub_cols] = nervous_sub_colnames
  emptycols = grep("^X.[0-9]+$", colnames(phenotypes), perl = T)
  if (length(emptycols) != 0){
    phenotypes = phenotypes[,-emptycols]
  }
  return(phenotypes)
}

write_phenotypes = function(fname_excel, fname_filelist, output_path, overwrite, plot, extra_phenotype = ""){
  
  filelist = read.table(fname_filelist, header = T, stringsAsFactors = F, sep = "\t")
  if (sum(duplicated(filelist$Sample)) >= 1){
    stop("The file list contains duplicate samples")
  }
  filelist = filelist %>% filter(grepl("Child|child", Name))
  
  pheno_sv = read_phenotype_from_excell(fname_excel, "Phenotype_overview")
  pheno_xl = read_phenotype_from_excell(fname_excel, "Xlinked_Phenotype_overview")
  pheno_snv = read_phenotype_from_excell(fname_excel, "Other MP patients")
  
  #combine excell sheets and remove wrong samples
  phenotypes = rbind(pheno_sv, pheno_xl, pheno_snv)
  phenotypes = phenotypes %>% filter(grepl("^MP[0-9]{2}", Patient, perl = T) & !duplicated(Patient) & Patient != "MP28-3041") %>% dplyr::select(-Male, -Female)
  phenotypes = inner_join(phenotypes, filelist, by = c("Patient" = "Sample"))
  phenotypes = replace(phenotypes, is.na(phenotypes), "")
  nrpatients = nrow(phenotypes)
  if (sum(lengths(strsplit(phenotypes$HPO_TERM, ",")) == lengths(strsplit(phenotypes$HPO_ID, ","))) != nrpatients){
    stop("The length of the HPO ids and terms is not the same for all patients.")
  }
  if (startsWith(extra_phenotype, "HP:")){
  phenotypes$HPO_ID = paste(phenotypes$HPO_ID, extra_phenotype, sep = ",")
  }
  phenotypes_out = paste0(output_path, "/Phenotypes/hpo_patients.txt")
  if (!file.exists(phenotypes_out) | overwrite == "True"){
    write.table(phenotypes[,c("Patient", "DGAP", "HPO_ID")], phenotypes_out, row.names = F, quote = F, sep = "\t")
  }
   if (plot == "True"){
    #Count the occurences of the different phenotypes
    countcols = grep("Abnormality|Nervous_system|Neoplasm", colnames(phenotypes))
    counts = apply(phenotypes[,countcols], 2, function(x) {sum(x != "")})
    nrmales = sum(phenotypes$Sex == "M")
    nrfemales = sum(phenotypes$Sex == "F")
    counts = c("Male" = nrmales, "Female" = nrfemales, counts)
    counts = counts / nrpatients
    
    #Put the occurences in a df and create formatted names
    phenotype_overview = data.frame("type" = names(counts), "counts" = counts, row.names = NULL)
    phenotype_overview$category = ifelse(grepl("Nervous_system", phenotype_overview$type), "Abnormality\nof the\n nervous system", "Abnormality of ...")
    phenotype_overview$category = ifelse(grepl("Male|Female", phenotype_overview$type), "Sex", phenotype_overview$category)
    phenotype_overview$category = factor(phenotype_overview$category, levels = c("Sex", "Abnormality\nof the\n nervous system", "Abnormality of ..."))
    
    phenotype_overview$typename = gsub("Nervous_system:|Abnormality.of.", "", phenotype_overview$type)
    phenotype_overview$typename = gsub("the.", "", phenotype_overview$typename)
    phenotype_overview$typename = gsub("\\.", " ", phenotype_overview$typename)
    phenotype_overview = phenotype_overview %>% filter(counts > 0.05)
    
    phenotype_overview_fig = ggplot(data = phenotype_overview, aes(x = typename, y = counts, fill = category)) + 
      geom_bar(stat = "identity", col = "black", width = 0.8) + theme_minimal() + 
      theme(axis.text.x = element_text(angle = 90, hjust = 1), 
            axis.title = element_text(face = "bold", size = 12),
            panel.spacing = unit(0, "lines"), 
            strip.placement = "outside",
            strip.background = element_rect(colour = "white", fill = "lightgrey"),
            legend.position="none",
            plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
            panel.grid.major.x = element_line(colour="darkgrey")) +  
      coord_flip(ylim = c(0,1)) + 
      facet_grid(category ~ ., switch="y", scales = "free_y", space = "free") +
      scale_y_continuous(labels = percent, expand = c(0.01,0.01), breaks = c(0,0.2,0.4,0.6,0.8,1)) + 
      labs(x  = "Clinical phenotype", y = "Frequency in cohort", title = paste("Overview of phenotypes (n=", nrpatients, ")", sep = ""))
    
    pdf(paste0(output_path, "/Phenotypes/phenotype_overview.pdf"))
    print(phenotype_overview_fig)
    dev.off()
  }
}
####________Change genes2phenotype file, so phenomatcher will include entrez gene ids in its output._####
mod_genes2pheno_file = function(fname, output_path, overwrite){
	fname_out = paste0(output_path, "/Phenotypes/genes2phenotype.txt")
  con = file(fname, "r")
	headerline = readLines(con, n = 1)
	close(con)
	genes2pheno = read.delim(fname, header = F, skip = 1)
	genes2pheno$V2 = paste(genes2pheno$V1, genes2pheno$V2, sep = "__")
	
	if (!file.exists(fname_out) | overwrite == "True"){
	  cat(paste0(headerline, "\n"), file = fname_out)
	  write.table(genes2pheno, fname_out , row.names = F, quote = F, sep = "\t", append = T, col.names = F)
  }
}

####________Get the snv windows. phenomatch will overlap genes with these windows.___________________####
get_snv_windows = function(extension, fname_filelist, output_path){
  filelist = read.table(fname_filelist, header = T, stringsAsFactors = F, sep = "\t")
  if (sum(duplicated(filelist$Sample)) >= 1){
    stop("The file list contains duplicate samples")
  }
  filelist = filelist %>% filter(grepl("Child|child", Name))
  folder = paste0(output_path, "/Possibly_pathogenic/SNVs/freq_chromfunction/")
  filelist$vcf_files_pathogenic = paste(folder, filelist$Sample, "_freq_chromfunction_filter.vcf", sep = "")
  
  vcfsdf = data.frame(stringsAsFactors = F)
  for (i in 1:nrow(filelist)){
    vcfname = filelist$vcf_files_pathogenic[i]
    sample = filelist$Sample[i]
    if (!file.exists(vcfname)){
      next
    }
    vcf = readVcf(vcfname, "hg19")
    vcfgr = rowRanges(vcf)
    vcfgr = vcfgr + extension
    vcfgr = GenomicRanges::trim(vcfgr)
    vcfdf = as.data.frame(vcfgr)[,c("seqnames", "start", "end")]
    vcfdf$sample = sample
    vcfdf$nr = seq_along(vcfdf[,1])
    vcfsdf = rbind(vcfsdf, vcfdf)
  }
  return(vcfsdf)
}


####________Run phenomatch___________________________________________________________________________####
Run_Phenomatch <- function(phenomatch, knowngenes_entrez, output_path, snv, extension, HPO_obo_file, overwrite, plot){
  HPO_patients_file = paste0(output_path, "/Phenotypes/hpo_patients.txt")
  HPO_genes_phenotype = paste0(output_path, "/Phenotypes/genes2phenotype.txt")
  output_folder = paste0(output_path, "/Phenotypes/")
  HPO_patients <- read.delim(HPO_patients_file, stringsAsFactors = F)
  snv = inner_join(snv, HPO_patients, by = c("sample" = "Patient"))
  patients = unique(snv$sample)
  snv = snv %>% mutate(chr = paste0("chr", seqnames), ID = paste(sample, nr, sep = "_"), HPO_ID = gsub(",", ";", HPO_ID)) %>% dplyr::select(chr, start, end, ID, HPO_ID)
  
  phenomatch_overview_all = data.frame(stringsAsFactors = F)
  for(patient in patients){
    patientrows = grep(patient, snv$ID)
    snv_patient = snv[patientrows,]
    if(nrow(snv_patient) == 0){
      next
    }
    
    #split into single hpo terms so their individual effects can be seen.
    HPO_split <- unlist(strsplit(snv_patient$HPO_ID[1], split = ";"))
    snv_patient_singlehpos = data.frame(stringsAsFactors = F)
    for (hpo in HPO_split){
      snv_patient_singlehpo = snv_patient %>% dplyr::select(-HPO_ID) %>% mutate(HPO_ID = hpo)
      snv_patient_singlehpos = rbind(snv_patient_singlehpos, snv_patient_singlehpo)
    }
    snv_patient = rbind(snv_patient, snv_patient_singlehpos)
    
    phenomatch_input_file = paste0(output_folder, "Phenomatch_input/", patient, "_phenomatch_input", extension, ".txt")
    write.table(snv_patient, phenomatch_input_file, quote = F, sep = "\t", col.names = F, row.names = F)
    
    phenomatch_output_file <- paste0(output_folder, "Phenomatch_raw_output/", patient, "_phenomatch_output", sep = "")
    
    # Phenomatch automatically adds ".overlapped_genes.txt" to the output filename
    phenomatch_file <- paste(phenomatch_output_file, ".overlapped_genes.txt", sep = "")
    if (file.exists(phenomatch_file) & overwrite != "True"){
      next
    }
    
    # This is the command to run Phenomatch.
    phenomatch_command <- paste("java -jar ", 
                                phenomatch," -i ", 
                                phenomatch_input_file, " -g ", 
                                knowngenes_entrez, " -O ", 
                                HPO_obo_file, " -a ", 
                                HPO_genes_phenotype, " -o ", 
                                phenomatch_output_file, sep = "")
    print(phenomatch_command)
    system(phenomatch_command)
    
    
    # Read the Phenomatch output
    phenomatch_patient <- read.delim(phenomatch_file)
    
    # Obtain the phenomatch score for the full phenotype (all HPO terms together) by removing the single HPO terms; except when here is only one HPO-term 
    if(length(HPO_split) > 1){
      phenomatch_patient_full <- phenomatch_patient[!as.vector(phenomatch_patient$phenotypes) %in% HPO_split,]
    }else{
      phenomatch_patient_full <- phenomatch_patient
    }
    
    phenomatch_overview <- data.frame(gene_symbol = phenomatch_patient_full$gene_symbol, phenoMatchScore = phenomatch_patient_full$phenoMatchScore, stringsAsFactors = F)
    phenomatch_overview <- phenomatch_overview[!duplicated(phenomatch_overview$gene_symbol),]
    phenomatch_overview <- phenomatch_overview[which(phenomatch_overview$phenoMatchScore > 0),]
    
    # Add the phenomatch scores for single HPO terms as columns to the phenomatch_overview
    for(hpo_term in HPO_split){
      phenomatch_single_HPO <- phenomatch_patient[phenomatch_patient$phenotypes == hpo_term,]
      
      phenomatch_single_HPO_filter <- phenomatch_single_HPO[,c("gene_symbol","phenoMatchScore")]
      phenomatch_single_HPO_filter <- phenomatch_single_HPO_filter[!duplicated(phenomatch_single_HPO_filter$gene_symbol),]
      phenomatch_single_HPO_filter <- phenomatch_single_HPO_filter[which(phenomatch_single_HPO_filter$phenoMatchScore > 0),]
      names(phenomatch_single_HPO_filter)[2] <- hpo_term
      phenomatch_overview <- merge(phenomatch_overview, phenomatch_single_HPO_filter, all.x = T)
    }
    
    # Translate the HPO ids to HPO terms:
    HPO_raw <- read.delim(HPO_genes_phenotype, skip = 1, header = F)
    names(HPO_raw) <- c("entrezgene", "Entrez_gene_name", "HPO_Term", "HPO_Term_ID")
    HPO_raw$entrezgene <- factor(HPO_raw$entrezgene)
    
    all_hpos <- HPO_raw[,c("HPO_Term", "HPO_Term_ID")]
    all_hpos <- all_hpos[!duplicated(all_hpos$HPO_Term_ID),]
    
    hpos <- data.frame(HPO_Term_ID = HPO_split)
    hpos_patient <- merge(hpos, all_hpos, by = "HPO_Term_ID")
    
    hpo <- c(names(phenomatch_overview))
    hpo[1] <- "HPO_ID"
    phenomatch_overview$gene_symbol <- as.vector(phenomatch_overview$gene_symbol)
    phenomatch_overview <- rbind(phenomatch_overview, hpo)
    
    for(names in names(phenomatch_overview)){
      hpo_id <- hpos_patient[hpos_patient$HPO_Term_ID == names,]
      if(nrow(hpo_id) > 0){
        names(phenomatch_overview)[which(names(phenomatch_overview) == names)] <- as.vector(hpo_id$HPO_Term)
      }
    }
    
    phenomatch_overview <- phenomatch_overview[which(phenomatch_overview$gene_symbol != "HPO_ID"),]
    phenomatch_overview$phenoMatchScore <- as.numeric(phenomatch_overview$phenoMatchScore)
    phenomatch_overview$gene_symbol <- factor(phenomatch_overview$gene_symbol, levels = phenomatch_overview$gene_symbol[order(phenomatch_overview$phenoMatchScore)])
    phenomatch_overview$sample = patient
    
    write.table(file = paste(output_folder, "Phenomatch_output/", patient, "_phenomatch_out.txt", sep = ""), x = phenomatch_overview, row.names = F, sep = "\t", quote = F)
    
    
    #Perform plotting
    if(plot == "True"){
      phenomatch_overview_melted = phenomatch_overview %>%  dplyr::select(-sample) %>% melt(id.vars = "gene_symbol")
      phenomatch_overview_melted$value <- as.numeric(phenomatch_overview_melted$value)
      phenomatch_overview_melted$Score <- phenomatch_overview_melted$value
      phenomatch_overview_melted$Score <- ifelse(phenomatch_overview_melted$Score > 30, 30, phenomatch_overview_melted$Score)
      
      phenomatch_overview_melted$value <- round(phenomatch_overview_melted$value, 1)
      
      names(phenomatch_overview_melted) <- c("Gene", "Phenotype", "phenoMatchScore", "Score")
      phenomatch_overview_melted$Gene = gsub(".*__", "", phenomatch_overview_melted$Gene)
      output_height <- nrow(phenomatch_overview) / 3 + 2
      
      output_width <- ncol(phenomatch_overview) / 2 + 2
      
      data_table <- ggplot(phenomatch_overview_melted, aes(x = Phenotype, y = Gene, fill = Score), label = phenoMatchScore) + geom_tile(lwd = 0.8, colour  = "black") +
        geom_text(size = 3.5, label = phenomatch_overview_melted$phenoMatchScore) + 
        scale_fill_gradient(low ="white", high = "red", limits=c(0,30)) +
        scale_x_discrete(position = "top") + 
        theme_minimal() + 
        theme(axis.text.x = element_text(angle = 60, hjust = 0), plot.title = element_text(hjust = 0.5)) + ggtitle(label = patient)
      
      
      pdf(file = paste(output_folder, "Phenomatch_output/", patient, "_overview.pdf", sep = ""), height = output_height, width = output_width)
      print(data_table)
      dev.off()
    }
    phenomatch_overview_score = phenomatch_overview[,c("gene_symbol", "phenoMatchScore", "sample")]
    phenomatch_overview_all = rbind(phenomatch_overview_all, phenomatch_overview_score)
  }
  write.table(phenomatch_overview_all, paste0(output_folder, "Phenomatch_output/allsamples.txt"), sep ="\t", row.names = F, quote = F)
}

####________Run functions____________________________________________________________________________####

write_phenotypes(output_path = opt$OUTPUT_PATH, extra_phenotype = opt$EXTRA_PHENOTYPE, fname_filelist = opt$FILELIST, fname_excel = opt$HPO_PATIENTS_EXCELL, overwrite = opt$OVERWRITE, plot = opt$PLOT)
mod_genes2pheno_file(fname = opt$GENES2PHENO, output_path = opt$OUTPUT_PATH, overwrite = opt$OVERWRITE)
vcfsdf = get_snv_windows(extension = opt$OVERLAP_GENE_SNV_DIST, fname_filelist = opt$FILELIST, output_path = opt$OUTPUT_PATH)
Run_Phenomatch(snv = vcfsdf, extension = opt$OVERLAP_GENE_SNV_DIST, output_path = opt$OUTPUT_PATH, HPO_obo_file = opt$HPO_OBO, phenomatch = opt$PHENOMATCH, knowngenes_entrez = opt$KNOWNGENES_ENTREZ, overwrite = opt$OVERWRITE, plot = opt$PLOT)

