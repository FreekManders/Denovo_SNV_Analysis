###This script overlaps the snvs with all genes that lie within a certain distance.
###The effect of different distances on the number of overlapping genes, can be visualized with the snv_gene_overlap_count.R script.

library(plyr)
library(dplyr)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(biomaRt)
library(GenomicRanges)
library(VariantAnnotation)
library(karyoploteR)
library(gdata)
library(scales)
library(stringr)
library(VennDiagram)
library(gridExtra)
library(optparse)

####________Parse command line arguments_____________________________________________________________####
option_list = list(
  make_option(c("-f", "--FILELIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/FilelistFreek.txt", 
              help="The filelist", metavar="character"),
  make_option(c("-o", "--OUTPUT_PATH"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV", 
              help="The output path", metavar="character"),
  make_option(c("--PLOT"), type="character", default="true", 
              help="Whether or not to plot", metavar="character"),
  make_option(c("--OVERLAP_GENE_SNV_DIST"), type="integer", default=500000, 
              help="area around the snv which will be overlapped with genes."),
  make_option(c("--PCHIC_CELLTYPES"), type="character", default="All", 
              help="Celltypes to use for PCHiC", metavar="character"),
  make_option(c("--GENE_LIST"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Complex_svs/Common_data/Genes/Genelist_HGNC_v4.txt", 
              help="Gene list", metavar="character"),
  make_option(c("--RNA_FOLDER1"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Complex_svs/RNA/Data/DE/", 
              help="Folder with rna files", metavar="character"),
  make_option(c("--RNA_FOLDER2"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/RNA/MP01_06/DE/", 
              help="Second folder with rna files. This is for samples sequenced in a different way.", metavar="character"),
  make_option(c("--PCHIC_FOLDER"), type="character", default="/hpc/cog_bioinf/cuppen/project_data/Complex_svs/Common_data/PCHiC/", 
              help="Folder with PCHiC data.", metavar="character")
  
  ); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


folder = opt$OUTPUT_PATH
setwd(folder)
filelist_f = opt$FILELIST
geneList = opt$GENE_LIST
folder_rna = opt$RNA_FOLDER1
folder_rna_2 = opt$RNA_FOLDER2
pchic_folder = opt$PCHIC_FOLDER

extension = opt$OVERLAP_GENE_SNV_DIST
plot = opt$PLOT
cell_typesPCHiC = opt$PCHIC_CELLTYPES


####_______Read in the filelist._______________________________________________________________####
filelist = read.table(filelist_f, header = T, stringsAsFactors = F, sep = "\t")
if (sum(duplicated(filelist$Sample)) >= 1){
  stop("The file list contains duplicate samples")
}
filelist = filelist %>% filter(grepl("Child|child", Name))
filelist$vcf_files_pathogenic = paste0("Possibly_pathogenic/SNVs/freq_chromfunction/", filelist$Sample, "_freq_chromfunction_filter.vcf")


####_______Read in gene annotations file.______________________________________________________####
genes = read.delim(geneList)
genes$strand = mapvalues(genes$strand, from = c("1", "-1"), to = c("+", "-"))
genesgr = makeGRangesFromDataFrame(genes, keep.extra.columns = T, seqnames.field = "chromosome_name", start.field = "start_position", end.field = "end_position")
genesgr = genesgr[seqnames(genesgr) %in% 1:22]

####_______Combine snvs from vcf with genes. results in snv*gene granges table.________________####
allvcfgrmatched = GRanges()
for (i in 1:nrow(filelist)){
  vcfname = filelist$vcf_files_pathogenic[i]
  sample = filelist$Sample[i]
  if (!file.exists(vcfname)){
    next
  }
  vcf = readVcf(vcfname, "hg19")
  vcfgr = rowRanges(vcf)
  vcfgr2 = vcfgr + extension
  mtch = as.matrix(GenomicRanges::findOverlaps(vcfgr2, genesgr))
  
  vcfgrmatched = vcfgr[mtch[,1]]
  genesgrmatched = genesgr[mtch[,2]]
  vcfgrmatched$start_gene = start(genesgrmatched)
  vcfgrmatched$end_gene = end(genesgrmatched)
  vcfgrmatched$strand_gene = as.vector(strand(genesgrmatched))
  vcfgrmatched$sample = sample
  mcols(vcfgrmatched) = cbind(mcols(vcfgrmatched), mcols(genesgrmatched))
  infovcf = info(vcf)[mtch[,1], c("GoNLv5_AC", "GoNLv5_AF", "GoNLv5_AN", "ANN", "DANN", "gnomAD_AC", "gnomAD_AF", "gnomAD_AN", "PON_COUNT", "Brain_Hippocampus_Middle", "Fetal_Brain_Female", "CD14_Primary_Cells", "Brain_Germinal_Matrix", "Fetal_Brain_Male", "Ensemblregbuild", "phastCons46way")]
  mcols(vcfgrmatched) = cbind(mcols(vcfgrmatched), infovcf)
  allvcfgrmatched = c(allvcfgrmatched, vcfgrmatched)
}

####_______Add columns to the snv_gene table.__________________________________________________####
allvcfgrmatched$snv_gene_id = seq_len(length(allvcfgrmatched))

#Change the ALT column from a DNAStrinSetList to a DNAStrinSet. This ensures it can be displayed normally in a data.frame
allvcfgrmatched$ALT = unstrsplit(allvcfgrmatched$ALT)

#Add column showing inheritance pattern of gene
allvcfgrmatched$dominant = ifelse(grepl("HP:0000006", allvcfgrmatched$HPO_Term_IDs), "Autosomal dominant", "Not autosomal dominant")

#Add columns based on the snpeff annotation
allvcfgrmatched$snpsift_gene = sapply(allvcfgrmatched$ANN, function(x) {annvector = strsplit(x[1], "\\|"); annvector[[1]][5]})
allvcfgrmatched$same_gene = allvcfgrmatched$snpsift_gene == allvcfgrmatched$ensembl_gene_id
allvcfgrmatched$same_gene[is.na(allvcfgrmatched$same_gene)] = F

effects = sapply(allvcfgrmatched$ANN, function(x) {gsub("^.\\|", "", x) %>% gsub("\\|.*", "", .) %>% unique(.) %>% grep("synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant", ., value = T) %>% paste(., collapse = ";")})
allvcfgrmatched$Effects_exonic = sapply(effects, function(x) {strsplit(x, "&|;") %>% .[[1]] %>% unique(.) %>% grep("synonymous_variant|missense_variant|stop_gained|stop_lost|stop_gained|start_lost|stop_retained_variant|splice_acceptor_variant|splice_donor_variant|splice_region_variant|protein_altering_variant|incomplete_terminal_codon_variant|coding_sequence_variant|NMD_transcript_variant", ., value = T) %>% paste(., collapse = ";")})
allvcfgrmatched$exonic = ifelse(allvcfgrmatched$Effects_exonic != "", "Exonic", "nonExonic")

allvcfgrmatched$snpsift_impact = sapply(allvcfgrmatched$ANN, function(x) {annvector = strsplit(x[1], "\\|"); annvector[[1]][3]})

#Add columns showing distance
allvcfgrmatched$distancedownstream = allvcfgrmatched$start_gene - start(allvcfgrmatched)
allvcfgrmatched$distanceupstream = start(allvcfgrmatched) - allvcfgrmatched$end_gene
allvcfgrmatched$distance = ifelse(allvcfgrmatched$distancedownstream >= 0, allvcfgrmatched$distancedownstream, allvcfgrmatched$distanceupstream)
allvcfgrmatched$distance[allvcfgrmatched$distance < 0] = 0 #snv is located between start and end of gene.
allvcfgrmatched$distancetostart = ifelse(allvcfgrmatched$strand_gene == "+", abs(allvcfgrmatched$start_gene - start(allvcfgrmatched)), abs(allvcfgrmatched$end_gene - start(allvcfgrmatched)))
allvcfgrmatched$distancetotss = abs(allvcfgrmatched$transcript_start - start(allvcfgrmatched))

#Add columns showing whether the snv is in an enhancer or promoter according to chromhmm
allvcfgrmatched$Fetal_Brain_Female_enhprom = allvcfgrmatched$Fetal_Brain_Female %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv")
allvcfgrmatched$Fetal_Brain_Male_enhprom = allvcfgrmatched$Fetal_Brain_Male %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv")

#Short notation for the chromhmm
allvcfgrmatched$Fetal_Brain_Female_short = gsub(".*_", "", allvcfgrmatched$Fetal_Brain_Female)
allvcfgrmatched$Fetal_Brain_Male_short = gsub(".*_", "", allvcfgrmatched$Fetal_Brain_Male)
allvcfgrmatched$Brain_Hippocampus_Middle_short = gsub(".*_", "", allvcfgrmatched$Brain_Hippocampus_Middle)
allvcfgrmatched$Brain_Germinal_Matrix_short = gsub(".*_", "", allvcfgrmatched$Brain_Germinal_Matrix)
allvcfgrmatched$CD14_Primary_Cells_short = gsub(".*_", "", allvcfgrmatched$CD14_Primary_Cells)

#change NAs to 0s so they look better in a table
allvcfgrmatched$GoNLv5_AC[is.na(allvcfgrmatched$GoNLv5_AC)] = 0
allvcfgrmatched$gnomAD_AC[is.na(allvcfgrmatched$gnomAD_AC)] = 0
allvcfgrmatched$PON_COUNT[is.na(allvcfgrmatched$PON_COUNT)] = 0

if (plot == "true"){
  ###Determine number of hpo terms per gene
  nrhpoterms = allvcfgrmatched$Number_HPO_Terms_Gene %>% as.data.frame(.)
  nrhpopergene_fig = ggplot(data = nrhpoterms, aes(x = .)) +
    geom_histogram(fill = "coral4", binwidth = 1) +
    scale_y_continuous(trans = "log1p", breaks = c(1,10,100,1000,10000), expand = c(0,0)) +
    labs(x = "nr hpo terms", y = "nr of genes (log10)", title = "Number of hpo terms per gene") +
    theme_bw()
  
  pdf("Possibly_pathogenic/SNVs/Near_gene/hpopergene_allsnvgenecombis.pdf")
  print(nrhpopergene_fig)
  dev.off()
}

####_______Add how high the DANN scores are, compared to the entire denovo snv set.____________####
allDanns = read.delim("Characteristics_denovo/vcfannotations/alldanns.txt", header = F)
dannquantsFunc = ecdf(allDanns$V1)
dannquants = dannquantsFunc(allvcfgrmatched$DANN)
allvcfgrmatched$dannPercentTop = percent((1 - dannquants))


####_______Count nr of genes between snv and matched gene._____________________________________####
#Upstream genes, downstream genes and genes around an snv are grouped seperately.
snv_gene = as.data.frame(allvcfgrmatched, row.names = seq_along(allvcfgrmatched))
snv_gene$snv_id = paste(snv_gene$sample, snv_gene$seqnames, snv_gene$start, sep = "_")

snv_gene$updown = ifelse(snv_gene$distancedownstream > 0, "down", "around_snv")
snv_gene$updown[snv_gene$distanceupstream > 0] = "up"
snv_gene = snv_gene %>% dplyr::group_by(snv_id, updown) %>% dplyr::arrange(distance) %>% dplyr::mutate(genes_between = row_number())
snv_gene$genes_between = snv_gene$genes_between - 1

#Count how many genes are around each snv. And add them the the nr of genes between snv and matched gene
nraround = snv_gene %>% dplyr::ungroup() %>% dplyr::filter(updown == "around_snv") %>% dplyr::group_by(snv_id) %>% dplyr::count()
snv_gene = left_join(snv_gene, nraround, by = "snv_id")
snv_gene$n[is.na(snv_gene$n)] = 0
snv_gene$genes_between = snv_gene$genes_between + snv_gene$n
snv_gene$genes_between[snv_gene$updown == "around_snv"] = 0
snv_gene = snv_gene %>% dplyr::select(-distanceupstream, -distancedownstream, -n)

####_______Read the phenomatch scores.________________________________________________________#### 
###uses functions from sourced file. Overlapping is done seperately by phenomatcher###
phenos = read.table(paste0(folder, "/Phenotypes/Phenomatch_output/allsamples.txt"), header = T, sep = "\t", stringsAsFactors = F)
phenos[,c("entrez_gene_id", "entrez_gene_symbol")] = str_split_fixed(phenos$gene_symbol, "__", 2)
phenos$entrez_gene_id = as.integer(phenos$entrez_gene_id)
phenos = phenos %>% dplyr::select(-gene_symbol)

####_______Merge the genes linked to the snvs with the genes from phenomatcher.________________####
####This merge works first on gene symbol. Then on entrez id.###
snv_gene = as.data.frame(snv_gene)
snv_gene$hgnc_symbol = as.character(snv_gene$hgnc_symbol)
snv_gene_hpo = inner_join(snv_gene, phenos, by = c("sample" = "sample", "hgnc_symbol" = "entrez_gene_symbol"))

snv_not_merged = anti_join(snv_gene, phenos, by = c("sample" = "sample", "hgnc_symbol" = "entrez_gene_symbol"))
phenos_not_merged = anti_join(phenos, snv_gene, by = c("sample" = "sample", "entrez_gene_symbol" = "hgnc_symbol"))
snv_gene_hpo2 = inner_join(snv_not_merged, phenos_not_merged, by = c("sample" = "sample", "entrezgene" = "entrez_gene_id")) #merge remaining genes on entrez id

snv_not_merged2 = anti_join(snv_not_merged, phenos_not_merged, by = c("sample" = "sample", "entrezgene" = "entrez_gene_id")) #Keep genes that aren't in phenomatch output.
snv_not_merged2$phenoMatchScore = 0

snv_gene_hpo = dplyr::select(snv_gene_hpo, -entrez_gene_id)
snv_gene_hpo2 = dplyr::select(snv_gene_hpo2, -entrez_gene_symbol)
snv_gene_hpo = rbind(snv_gene_hpo, snv_gene_hpo2, snv_not_merged2) #Bind them all together.
snv_gene_hpo[snv_gene_hpo$phenoMatchScore == 0, "dominant"] = "Not in HPO"


snv_gene_hpo$ANN = gsub("\n", "", snv_gene_hpo$ANN)
write.table(snv_gene_hpo, "Possibly_pathogenic/SNVs/Near_gene/snv_gene.txt", quote = F, sep = "\t", row.names = F)

###Check for exonic genes, that aren't in de snv_gene table. This is done by comparing with the vcfs from cartageniavspbt.R
# exonic_matchedgene = snv_gene %>% dplyr::select(sample, seqnames, start, exonic) %>% filter(exonic == "Exonic") %>% unique(.)
# exonic = vcfs %>% dplyr::select(Sample, chrom, start, Exoniclabel) %>% filter(Exoniclabel == "Exonic")
# anti_join(exonic, exonic_matchedgene, by = c("Sample" = "sample", "seqnames" = "chrom", "start"))
# exonic_matchedgene_hpo = snv_gene_hpo %>% dplyr::select(sample, seqnames, start, exonic) %>% filter(exonic == "Exonic") %>% unique(.)
# anti_join(exonic_matchedgene, exonic_matchedgene_hpo, by = c("sample", "seqnames", "start"))


if (plot == "true"){
###Plot characteristics of snv_gene combinations with hpo phenomatch scores
  snv_gene_counts = data.frame("group" = c("all", "with phenomatch score"), "counts" = c(nrow(snv_gene), nrow(snv_gene_hpo[snv_gene_hpo$phenoMatchScore != 0,])))
  snv_gene_counts_fig = ggplot(data = snv_gene_counts, aes(x = group, y = counts, fill = group)) +
    geom_bar(stat = "identity") +
    coord_cartesian(ylim = c(0, 1.1 * max(snv_gene_counts$counts))) +
    labs(x = "", y = "nr. of snv gene combinations", title = "snv gene combinations") +
    theme_bw() +
    scale_fill_discrete(l = 30)
  
  phenomatch_fig = ggplot(data = snv_gene_hpo, aes(y = phenoMatchScore, x = "")) + 
    geom_boxplot(outlier.shape = NA) +
    geom_dotplot(stackdir = "center", binaxis = "y", position = "dodge", binwidth = 1, binpositions = "all", dotsize = 0.8, stackratio = 0.8) +
    labs(x = "All samples", title = "Phenomatch scores of all combinations") +
    theme_bw()
  
  pdf("Possibly_pathogenic/SNVs/Near_gene/snv_gene_characteristics.pdf")
  print(snv_gene_counts_fig)
  print(phenomatch_fig)
  dev.off()
}

if (plot == "true"){
  ###Plot inheritance of the genes.
  fname_genes2pheno = "Phenotypes/genes2phenotype.txt"
  genes2pheno = read.delim(fname_genes2pheno, header = F, skip = 1)
  genes2pheno = group_by(genes2pheno, V2) %>% summarise(hpo_ids = paste(V4, collapse = ";"))
  nrboth_all = sum(grepl("HP:0000006", genes2pheno$hpo_ids) & grepl("HP:0000007", genes2pheno$hpo_ids))
  nrdominant_all = sum(grepl("HP:0000006", genes2pheno$hpo_ids))
  nrrecessive_all = sum(grepl("HP:0000007", genes2pheno$hpo_ids))
  inh_genes2pheno = data.frame("Inheritance" = c("Autosomal recessive", "Autosomal dominant", "Both"), "Nrgenes" = c(nrrecessive_all, nrdominant_all, nrboth_all), "genes" = "all")
  inh_genes2pheno$Nrgenes = inh_genes2pheno$Nrgenes / nrow(genes2pheno)
  
  bothgenes_r = grepl("HP:0000006", snv_gene_hpo$HPO_Term_IDs) & grepl("HP:0000007", snv_gene_hpo$HPO_Term_IDs)
  nrboth = length(unique(snv_gene_hpo[bothgenes_r,"ensembl_gene_id"]))
  dominantgenes_r = grepl("HP:0000006", snv_gene_hpo$HPO_Term_IDs)
  nrdominant = length(unique(snv_gene_hpo[dominantgenes_r,"ensembl_gene_id"]))
  recessivegenes_r = grepl("HP:0000007", snv_gene_hpo$HPO_Term_IDs)
  nrrecessive = length(unique(snv_gene_hpo[recessivegenes_r,"ensembl_gene_id"]))
  
  inh_snv_gene_hpo = data.frame("Inheritance" = c("Autosomal recessive", "Autosomal dominant", "Both"), "Nrgenes" = c(nrrecessive, nrdominant, nrboth), "genes" = "linked to snv")
  inh_snv_gene_hpo$Nrgenes = inh_snv_gene_hpo$Nrgenes / length(unique(snv_gene_hpo[snv_gene_hpo$phenoMatchScore != 0, "ensembl_gene_id"]))
  
  inh_genes = rbind(inh_genes2pheno, inh_snv_gene_hpo)
  inh_genes_fig = ggplot(data = inh_genes, aes(x = Inheritance, y = Nrgenes, fill = genes)) +
    geom_bar(stat = "identity", position = "dodge") +
    coord_cartesian(ylim = c(0, 1.1 * max(inh_genes$Nrgenes)), expand = F) +
    scale_y_continuous(labels = percent) +
    scale_x_discrete(limits = c("Autosomal recessive", "Autosomal dominant", "Both")) +
    theme_bw() +
    scale_fill_discrete(l = 30) +
    labs(x = "Mode of inheritance", y = "Nr of genes", title = "Inheritance of genes (with hpo terms)")
  
  inh_HI_fig = ggplot(data = snv_gene_hpo, aes(x = dominant, y = HI)) +
    geom_boxplot(outlier.shape = NA) +
    geom_dotplot(stackdir = "center", binaxis = "y", position = "dodge", binwidth = 1, binpositions = "all", dotsize = 0.8, stackratio = 0.8) +
    labs(x = "HPO inheritance", title = "HI versus inheritance HPO") +
    theme_bw()
  inh_pLI_fig = ggplot(data = snv_gene_hpo, aes(x = dominant, y = pLI)) +
    geom_boxplot(outlier.shape = NA) +
    geom_dotplot(stackdir = "center", binaxis = "y", position = "dodge", binwidth = 0.005, binpositions = "all", dotsize = 0.8, stackratio = 0.8) +
    labs(x = "HPO inheritance", title = "pLI versus inheritance HPO") +
    theme_bw()
  
  pdf("Phenotypes/inheritance.pdf", width = 15)
  print(inh_genes_fig)
  print(inh_HI_fig)
  print(inh_pLI_fig)
  dev.off()
}


####_______Match with RNA._____________________________________________________________________####
filelist$fname_rna = paste0(folder_rna, "MP_", filelist$Sample, "_DE_results.bed")
filelist$fname_rna_2 = paste0(folder_rna_2, "MP_", filelist$Sample, "_DE_results.bed")
rna = data.frame(stringsAsFactors = F)
for (i in 1:nrow(filelist)){
  fname_rna = filelist$fname_rna[i]
  sample = filelist$Sample[i]
  if (!file.exists(fname_rna)){#check if rna file exists. If it doesn't exist check if it exists in folder 2.
    fname_rna = filelist$fname_rna_2[i]
    if (!file.exists(fname_rna)){
      next
    }
  }
  rna_sample = read.delim(fname_rna)
  rna_sample$sample = sample
  rna_sample$ensembl_gene_id = rownames(rna_sample)
  if (length(grep("shrunk", colnames(rna_sample))) == 0){
    rna_sample[,c("log2FoldChange_shrunk", "log2lfcSE_shrunk")] = NA
  }
  rna = rbind(rna, rna_sample)
}
rownames(rna) = seq_along(rownames(rna))

snv_gene_hpo$ensembl_gene_id = as.character(snv_gene_hpo$ensembl_gene_id)
snv_gene_table = left_join(snv_gene_hpo, rna, by = c("sample" = "sample", "ensembl_gene_id" = "ensembl_gene_id"))

write.table(snv_gene_table, "Possibly_pathogenic/SNVs/Near_gene/snv_gene_rna.txt", quote = F, row.names = F, sep = "\t")


####_______Match with PCHIC.___________________________________________________________________####
genes_transcribed = genesgr[!is.na(genesgr$transcription_start_site)]
genes_tss = GRanges(seqnames = seqnames(genes_transcribed), IRanges( start = genes_transcribed$transcription_start_site, end = genes_transcribed$transcription_start_site +1), ensembl_gene_id = genes_transcribed$ensembl_gene_id)

snv_gene_table_pos = GRanges(seqnames = snv_gene_table$seqnames,
                             IRanges(start = snv_gene_table$start, end = snv_gene_table$end))


###Either supply celltypes manually or use all celltypes in the folder
if (cell_typesPCHiC == "All"){
	cell_typesPCHiC = list.files(pchic_folder, ".txt") %>% gsub("PCHiC_", "", .) %>% gsub(".txt", "", .)
} else {
	cell_typesPCHiC = strsplit(cell_typesPCHiC, ",")[[1]]
}

###Read the pchic data
pchic = data.frame(stringsAsFactors = F)
for (cell_type in cell_typesPCHiC){
  pchic_fname = paste0(pchic_folder, "PCHiC_", cell_type, ".txt")
  if (!file.exists(pchic_fname)){
    print(paste("The PCHiC file:", pchic_fname, "does not exist."))
    next
  }
  pchic_sample = read.delim(pchic_fname, stringsAsFactors = F)
  pchic = rbind(pchic, pchic_sample) 
}
pchic$Cell_type = paste0(pchic$Cell_type, "_celltype_PCHiC")
pchic$bait_chr = gsub(pattern = "chr", replacement = "", x = pchic$bait_chr)
pchic$PIR_chr = gsub(pattern = "chr", replacement = "", x = pchic$PIR_chr)

###Overlap the baits with the tss of the genes, to add the genes associated with the baits to the pchic data.frame.
pchic_bait = GRanges(seqnames = pchic$bait_chr,
                         IRanges(start = pchic$bait_start, 
                                 end = pchic$bait_end))
bait_tss_overlap = findOverlaps(pchic_bait, genes_tss)
pchic = pchic %>% dplyr::slice(queryHits(bait_tss_overlap)) %>% dplyr::mutate(ensembl_gene_id = genes_tss$ensembl_gene_id[subjectHits(bait_tss_overlap)])

###Overlap the captures with the snvs in the snv_gene_table. The snv_gene_id is added to the pchic data.frame.
pchic_capture = GRanges(seqnames = pchic$PIR_chr,
                        IRanges(start = pchic$PIR_start, end = pchic$PIR_end))
capture_snvgene_overlap = findOverlaps(pchic_capture, snv_gene_table_pos)
pchic = pchic %>% dplyr::slice(queryHits(capture_snvgene_overlap)) %>% dplyr::mutate(snv_gene_id = snv_gene_table$snv_gene_id[subjectHits(capture_snvgene_overlap)])

###Merge the pchic score to the snv_gene_table based on the snv_gene_id and the gene name columns.
pchic = pchic %>% dplyr::select(snv_gene_id, ensembl_gene_id, Cell_type, Score) %>% dcast(., snv_gene_id + ensembl_gene_id ~ Cell_type, value.var = "Score")
pchic$ensembl_gene_id = as.character(pchic$ensembl_gene_id)
snv_gene_table = left_join(snv_gene_table, pchic, by = c("snv_gene_id", "ensembl_gene_id"))

###Calculate in how many cells the interactions are present.
cellsmissed = snv_gene_table %>% dplyr::select(contains("_celltype_PCHiC")) %>% is.na(.) %>% rowSums(.)
nrcells = grep("_celltype_PCHiC", colnames(snv_gene_table)) %>% length(.)
nrinteractions = nrcells - cellsmissed
snv_gene_table$fraction_cells_PCHiC = nrinteractions / nrcells
snv_gene_table$fraction_cells_PCHiC_makeup = paste0(nrinteractions, "/", nrcells)

write.table(snv_gene_table, "Possibly_pathogenic/SNVs/Near_gene/snv_gene_rna_pchic.txt", quote = F, row.names = F, sep = "\t")

if (plot == "true"){
  
  fraction_cells_PCHiC = as.data.frame(table(snv_gene_table$fraction_cells_PCHiC))
  pchic_fig = ggplot(data = fraction_cells_PCHiC, aes(x = Var1, y = Freq)) +
    geom_bar(stat = "identity", fill = "skyblue4") +
    coord_cartesian(ylim = c(1, 1.1 * max(fraction_cells_PCHiC$Freq)), expand = F) +
    scale_y_log10() +
    labs(x = paste0("Interaction found in nr of cell types (", nrcells, " cell types)"), y = "Nr. of snv gene combinations", title = "SNV-Gene interactions according th Promoter Capture HiC") +
    scale_x_discrete(labels = seq(0, nrcells)) +
    theme_bw() +
    theme(panel.grid.major.x = element_blank())
  
  pchic_regbuild_fig = ggplot(data = snv_gene_table, aes(x = fraction_cells_PCHiC, color = Ensemblregbuild)) +
    geom_density() +
    coord_cartesian(expand = F) +
    labs(y = "Density", x = paste0("Interaction found in fraction of cells (", nrcells, " cell types)"), title = "SNV-Gene interactions according to Promoter Capture HiC") +
    theme_bw()
  
  pchic_distance_fig = ggplot(data = snv_gene_table, aes(x = distancetostart, y = fraction_cells_PCHiC)) +
    geom_point() +
    geom_smooth(method = "gam") +
    labs(x = "distance between snv and start of gene", y = "Interaction found in fraction of cell types", title = "Effect of distance between snv and gene on the Promoter Capture HiC interactions") +
    theme_bw()
  
  pdf("Possibly_pathogenic/SNVs/Near_gene/pchic.pdf")
  print(pchic_fig)
  print(pchic_regbuild_fig)
  print(pchic_distance_fig)
  dev.off()
}

####_______Classify snvs as possibly pathogenic________________________________________________####
#snv_gene_table = read.delim("/hpc/cog_bioinf/cuppen/project_data/Freek_SNV/Possibly_pathogenic/SNVs/Near_gene/snv_gene_rna_pchic.txt", stringsAsFactors = F)

#Remove less interesting columns, so the table becomes clearer for viewing..
#snv_gene_table_view = snv_gene_table %>%  dplyr::select(seqnames, start, REF, ALT, sample, hgnc_symbol, pLI, RVIS, HI, DDG2P, GoNLv5_AC, gnomAD_AC, PON_COUNT, DANN, dannPercentTop, phastCons46way, snpsift_impact, distance, distancetostart, distancetotss, genes_between, Ensemblregbuild, Effects_exonic, exonic, dominant, same_gene, phenoMatchScore, baseMean, log2FoldChange, log2FoldChange_shrunk, pvalue, padj, fraction_cells_PCHiC, fraction_cells_PCHiC_makeup, Fetal_Brain_Female, Fetal_Brain_Female_enhprom, Fetal_Brain_Male, Fetal_Brain_Male_enhprom, Brain_Hippocampus_Middle, Brain_Germinal_Matrix, CD14_Primary_Cells, snv_gene_id)

###For non exonic variants
#Genetic
#genetic = snv_gene_table_view %>% dplyr::filter(exonic == "nonExonic" & ((fraction_cells_PCHiC > 0.2) | (distancetotss < 100000 & (!is.na(Ensemblregbuild) | Fetal_Brain_Female %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv") | Fetal_Brain_Male %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv")))))

#Phenotypic
#phenotypic = snv_gene_table_view %>% dplyr::filter(exonic == "nonExonic" & (as.numeric(phenoMatchScore) >= 10 | !is.na(DDG2P)) & dominant == "Autosomal dominant")

#Possibly pathogenic (genetic + phenotypic)
possPathogenic = snv_gene_table %>% dplyr::filter((as.numeric(phenoMatchScore) >= 10 | !is.na(DDG2P)) & dominant == "Autosomal dominant" & exonic == "nonExonic" & ((fraction_cells_PCHiC > 0.2) | (distancetotss < 100000 & (!is.na(Ensemblregbuild) | Fetal_Brain_Female %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv") | Fetal_Brain_Male %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv")))))
possPathogenic$snv_gene_id = as.factor(possPathogenic$snv_gene_id)
possPathogenic$DANN = as.numeric(possPathogenic$DANN)
possPathogenic$phastCons46way = as.numeric(possPathogenic$phastCons46way)
possPathogenic$phenoAsId = as.factor(possPathogenic$phenoMatchScore + runif(nrow(possPathogenic), 0, 0.0000000001))
  
#Perform local p-adjustment. Look at pvalues for genes linked with a snv, that passed deseq2 independent filtering.
possPathogenic$pvalue2 = possPathogenic$pvalue
possPathogenic$pvalue2[is.na(possPathogenic$padj)] = NA
#possPathogenic = possPathogenic %>% group_by(sample) %>% mutate(padj_local = p.adjust(pvalue2, method = "fdr")) %>% dplyr::ungroup()
possPathogenic$padj_local = p.adjust(possPathogenic$pvalue2, method = "fdr")
possPathogenic$padj_local_nona = ifelse(is.na(possPathogenic$padj_local), 1, possPathogenic$padj_local)

write.table(possPathogenic, "Possibly_pathogenic/SNVs/Near_gene/possibly_pathogenic.txt", quote = F, row.names = F, sep = "\t")
possPathogenicTableReport = possPathogenic %>%  dplyr::select(seqnames, start, REF, ALT, sample, hgnc_symbol, pLI, RVIS, HI, DDG2P, GoNLv5_AC, gnomAD_AC, PON_COUNT, DANN, phastCons46way, snpsift_impact, distance, distancetostart, distancetotss, genes_between, Ensemblregbuild, same_gene, phenoMatchScore, baseMean, log2FoldChange, log2FoldChange_shrunk, pvalue, padj_local, fraction_cells_PCHiC_makeup, Fetal_Brain_Female, Fetal_Brain_Male, Brain_Hippocampus_Middle_short, Brain_Germinal_Matrix_short, CD14_Primary_Cells_short)
write.table(possPathogenic, "Possibly_pathogenic/SNVs/Near_gene/possibly_pathogenic_interestingcollumns.txt", quote = F, row.names = F, sep = "\t")

###For exonic variants
possPathogenicExonic = snv_gene_table %>% dplyr::filter(exonic == "Exonic" & same_gene & as.numeric(phenoMatchScore) > 0 & Effects_exonic != "synonymous_variant")

#Perform local p-adjustment. Look at pvalues for genes linked with a snv, that passed deseq2 independent filtering.
possPathogenicExonic$pvalue2 = possPathogenicExonic$pvalue
possPathogenicExonic$pvalue2[is.na(possPathogenicExonic$padj)] = NA
possPathogenicExonic$padj_local = p.adjust(possPathogenicExonic$pvalue2, method = "fdr")
possPathogenicExonic$padj_local_nona = ifelse(is.na(possPathogenicExonic$padj_local), 1, possPathogenicExonic$padj_local)

write.table(possPathogenicExonic, "Possibly_pathogenic/SNVs/Near_gene/possibly_pathogenic_Exonic.txt", quote = F, row.names = F, sep = "\t")
possPathogenicExonicTableReport = possPathogenicExonic %>%  dplyr::select(seqnames, start, REF, ALT, sample, hgnc_symbol, pLI, RVIS, HI, DDG2P, GoNLv5_AC, gnomAD_AC, PON_COUNT, Effects_exonic, DANN, phastCons46way, snpsift_impact, distance, distancetostart, distancetotss, genes_between, Ensemblregbuild, same_gene, phenoMatchScore, baseMean, log2FoldChange, log2FoldChange_shrunk, pvalue, padj_local, fraction_cells_PCHiC_makeup, Fetal_Brain_Female_short, Fetal_Brain_Male_short, Brain_Hippocampus_Middle_short, Brain_Germinal_Matrix_short, CD14_Primary_Cells_short)
write.table(possPathogenicExonicTableReport, "Possibly_pathogenic/SNVs/Near_gene/possibly_pathogenic_Exoncic_interestingcollumns.txt", quote = F, row.names = F, sep = "\t")

if (plot == "true"){
#Venn diagrams showing on which criteria non exonic snv-gene combinations were classified as pathogenic.
  nrPCHiCandDist = possPathogenic %>% dplyr::filter((fraction_cells_PCHiC > 0.2) & (distancetotss < 100000 & (!is.na(Ensemblregbuild) | Fetal_Brain_Female %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv") | Fetal_Brain_Male %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv")))) %>% nrow(.)
  nrPCHiC = possPathogenic %>% dplyr::filter((fraction_cells_PCHiC > 0.2)) %>% nrow(.)
  nrDist= possPathogenic %>% dplyr::filter((distancetotss < 100000 & (!is.na(Ensemblregbuild) | Fetal_Brain_Female %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv") | Fetal_Brain_Male %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv")))) %>% nrow(.)
  vennPCHiCorDist = draw.pairwise.venn(nrPCHiC, nrDist, nrPCHiCandDist, category = c("PCHiC","Distance + regulatory region"),inverted = F, ext.text = F, fill = c("darkred","skyblue"), alpha = c(0.8,0.95), cat.pos = c(180,180), cat.cex = c(2,2), cex = c(3,3,3), ind = F)
  
  nrDistRegbuildandChromhmm = possPathogenic %>% dplyr::filter(distancetotss < 100000 & !is.na(Ensemblregbuild) & (Fetal_Brain_Female %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv") | Fetal_Brain_Male %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv"))) %>% nrow(.)
  nrDistRegbuild = possPathogenic %>% dplyr::filter(distancetotss < 100000 & !is.na(Ensemblregbuild)) %>% nrow(.)
  nrDistChromhmm = possPathogenic %>% dplyr::filter(distancetotss < 100000 & (Fetal_Brain_Female %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv") | Fetal_Brain_Male %in% c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv"))) %>% nrow(.)
  vennRegbuildorChromhmm = draw.pairwise.venn(nrDistRegbuild, nrDistChromhmm, nrDistRegbuildandChromhmm, category = c("Ensembl regulatory build", "ChromHMM: Fetal Brain Male/Female"),inverted = F, ext.text = F, fill = c("dodgerblue","lightblue"), alpha = c(0.8,0.95), cat.pos = c(0,0), cat.dist = c(0.07,0.07), cat.cex = c(2,2), cex = c(3,3,3), ind = F)
  
  nrPhenoandDDG = possPathogenic %>% dplyr::filter(as.numeric(phenoMatchScore) >= 10 & !is.na(DDG2P)) %>% nrow(.)
  nrPheno = possPathogenic %>% dplyr::filter(as.numeric(phenoMatchScore) >= 10) %>% nrow(.)
  nrDDG = possPathogenic %>% dplyr::filter(!is.na(DDG2P)) %>% nrow(.)
  vennPhenoorDDG = draw.pairwise.venn(nrPheno, nrDDG, nrPhenoandDDG, category = c("PhenoMatch score", "DDG2P"),inverted = F, ext.text = F, fill = c("darkred","skyblue"), alpha = c(0.8,0.95), cat.pos = c(0,0), cat.dist = c(0.05,0.05), cat.cex = c(2,2), cex = c(3,3,3), ind = F)
  
  pdf("Possibly_pathogenic/SNVs/Near_gene/possPathogencClassifiedCriteria.pdf")
  grid.arrange(gTree(children = vennPCHiCorDist), top="Nr of snv-gene combinations classified as possibly pathogenic according to:")
  grid.arrange(gTree(children = vennRegbuildorChromhmm), top="Regulatory region:")
  grid.arrange(gTree(children = vennPhenoorDDG), top="Nr of snv-gene combinations classified as possibly pathogenic according to:")
  dev.off()
}

 
####_______Create final report figure, for the non exonic possibly pathogenic snvs_____________####
#Facet 1: patient info
if (plot == "true"){
  Patient = data.frame(id = possPathogenic$phenoAsId,
                        Variable = "Patient ID",
                        Facet = "Patient info",
                        Value = 0,
                        Label = possPathogenic$sample)
  
  #Facet 2: gene info
  Genesymbol = data.frame(id = possPathogenic$phenoAsId,
                     Variable = "Gene symbol",
                     Facet = "Gene info",
                     Value = 0,
                     Label = possPathogenic$hgnc_symbol)
  
  Distance = data.frame(id = possPathogenic$phenoAsId,
                         Variable = "Distance to TSS\n(Kb)",
                         Facet = "Gene info",
                         Value = ifelse(possPathogenic$distancetotss < 1e5, (1e5 - possPathogenic$distancetotss) / 1e5 * 100 , 0),
                         Label = as.character(round(possPathogenic$distancetotss / 1e3, 1)))
  
  Genesbetween = data.frame(id = possPathogenic$phenoAsId,
                        Variable = "Nr genes between",
                        Facet = "Gene info",
                        Value = ifelse(possPathogenic$genes_between < 10 , (10 - possPathogenic$genes_between) / 10 * 100, 0),
                        Label = as.character(possPathogenic$genes_between))
  
  #Facet 3: SNV info
  Ensemblregbuild = data.frame(id = possPathogenic$phenoAsId,
                                Variable = "Ensembl regulatory build",
                                Facet = "SNV region",
                                Value = ifelse(!is.na(possPathogenic$Ensemblregbuild), 50, 0),
                                Label = possPathogenic$Ensemblregbuild)
  
  Fetal_Brain_Female = data.frame(id = possPathogenic$phenoAsId,
                                   Variable = "Fetal brain female",
                                   Facet = "SNV region",
                                   Value = ifelse(possPathogenic$Fetal_Brain_Female_enhprom, 50, 0),
                                   Label = possPathogenic$Fetal_Brain_Female_short)
  
  Fetal_Brain_Male = data.frame(id = possPathogenic$phenoAsId,
                                 Variable = "Fetal brain male",
                                 Facet = "SNV region",
                                 Value = ifelse(possPathogenic$Fetal_Brain_Male_enhprom, 50, 0),
                                 Label = possPathogenic$Fetal_Brain_Male_short)
  
  DannTop = data.frame(id = possPathogenic$phenoAsId,
                        Variable = "Dann",
                        Facet = "SNV region",
                        Value = possPathogenic$DANN * 100,
                        Label = as.character(round(possPathogenic$DANN, 2)))
  
  phastcons46way = data.frame(id = possPathogenic$phenoAsId,
                               Variable = "phastCons46way",
                               Facet = "SNV region",
                               Value = possPathogenic$phastCons46way * 100,
                               Label = as.character(round(possPathogenic$phastCons46way, 2)))
  
  #Facet 4: Phenotype info
  Phenomatch_Score = data.frame(id = possPathogenic$phenoAsId,
                                 Variable = "Phenomatch",
                                 Facet = "Phenotype",
                                 Value = ifelse(possPathogenic$phenoMatchScore > 30, 100, possPathogenic$phenoMatchScore/30*100),
                                 Label = as.character(round(possPathogenic$phenoMatchScore, 1)))
  
  possPathogenic$DDG2P = factor(possPathogenic$DDG2P)
  levels(possPathogenic$DDG2P) = list("Yes" = "confirmed", "Prob" = "probable", "Pos"= "possible", "No" = NA)
  DDG2P = data.frame(id = possPathogenic$phenoAsId,
                      Variable = "DDG2P",
                      Facet = "Phenotype",
                      Value = ifelse(possPathogenic$DDG2P %in% c("Yes", "Prob"), 75, 0),
                      Label = as.character(possPathogenic$DDG2P))
  
  pLI = data.frame(id = possPathogenic$phenoAsId,
                    Variable = "pLI",
                    Facet = "Phenotype",
                    Value = ifelse(possPathogenic$pLI > 0.5, (possPathogenic$pLI-0.5)*200 ,0),
                    Label = as.character(round(possPathogenic$pLI, 2)))
  
  RVIS = data.frame(id = possPathogenic$phenoAsId,
                     Variable = "RVIS",
                     Facet = "Phenotype",
                     Value = ifelse(possPathogenic$RVIS < 50, 100 - possPathogenic$RVIS*2,0),
                     Label = as.character(round(possPathogenic$RVIS, 1)))

  HI = data.frame(id = possPathogenic$phenoAsId,
                   Variable = "HI",
                   Facet = "Phenotype",
                   Value = ifelse(possPathogenic$HI < 50, 100 - possPathogenic$HI*2, 0),
                   Label = as.character(round(possPathogenic$HI, 1)))
  
  
  #Facet 5: RNA
  RNAfoldchange = data.frame(id = possPathogenic$phenoAsId,
                      Variable = "RNA fold change\n(log2)",
                      Facet = "RNA",
                      Value = ifelse(possPathogenic$baseMean > 100,
                                     ifelse(possPathogenic$log2FoldChange < 0, 
                                            possPathogenic$log2FoldChange / 2 * -100,
                                            0),
                                              NA),
                    Label = ifelse(possPathogenic$baseMean > 100,
                                   ifelse(possPathogenic$padj_local_nona < 0.1, 
                                          paste(as.character(round(possPathogenic$log2FoldChange, 1)), "*",sep= ""),
                                          as.character(round(possPathogenic$log2FoldChange, 1))),
                                   "NE"))
  
  #Facet 6: PCHIC
  PCHiCInteractions = data.frame(id = possPathogenic$phenoAsId,
                        Variable = "PCHIC",
                        Facet = "Interaction",
                        Value = possPathogenic$fraction_cells_PCHiC * 100,
                        Label = possPathogenic$fraction_cells_PCHiC_makeup)
  #Combine all variables.
  plotData = rbind(Patient, Genesymbol, Distance, Genesbetween, Ensemblregbuild, Fetal_Brain_Female, Fetal_Brain_Male, DannTop, phastcons46way, Phenomatch_Score, DDG2P, pLI, RVIS, HI, RNAfoldchange, PCHiCInteractions)
  
  #plotting parameters
  bottommargin = 0.5
  nrmuts = nrow(possPathogenic)
  width = c(rep(1.2, nrmuts), rep(1, 3 * nrmuts), rep(4.4, nrmuts), rep(1.2, 2 * nrmuts), rep(1, 9 * nrmuts))
  hjust = c(rep(0.52, nrmuts), rep(0.5, 3 * nrmuts), rep(0.85, nrmuts), rep(0.6, nrmuts), rep(0.55, nrmuts), rep(0.5, 9 * nrmuts))
  
  ReportFig = ggplot(data = plotData, aes(x = Variable, y = id, fill = Value, width = width)) + 
    geom_tile(col = "black") + 
    geom_text(aes(label = Label, hjust = hjust), size = 3, na.rm = T) + 
    facet_grid(~ Facet, scales = "free", space = "free") + 
    scale_fill_gradient2(low = "darkblue", mid = "white", high = rgb(178,34,34, maxColorValue = 255), na.value = "gray", guide = F)+
    scale_x_discrete(position = "top") +
    theme_minimal() +
    theme(axis.ticks.y = element_blank(), 
          axis.title.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.text.x = element_text(angle = 60, hjust = 0, size = 12), plot.title = element_text(hjust = 0.5),
          strip.placement = "outside", strip.background = element_rect(fill = "lightgray", colour = "white"), strip.text = element_text(face = "bold", size = 10), 
          plot.margin = unit(c(0.5,0.5,bottommargin,0.5), "cm"))
  
  pdf("Possibly_pathogenic/SNVs/Near_gene/reportFig.pdf", width = 12, height = 10)
  print(ReportFig)
  dev.off()
}
####_______Effect of adding autosomal dominant to the hpo terms of the patients._______________####
#Run earlier part of script with and without dominant inheritance added to phenotypes of patients.
# no_inh= snv_gene_hpo_noinh %>% dplyr::select(sample, start, ensembl_gene_id, phenoMatchScore)
# with_inh= snv_gene_hpo_withinh %>% dplyr::select(sample, start, ensembl_gene_id, phenoMatchScore, HPO_Term_IDs)
# change_phenomatch = inner_join(no_inh, with_inh, by = c("sample", "start", "ensembl_gene_id"))
# change_phenomatch$change = change_phenomatch$phenoMatchScore.y - change_phenomatch$phenoMatchScore.x
# change_phenomatch$Dominant = ifelse(grepl("HP:0000006", change_phenomatch$HPO_Term_IDs), "Autosomal dominant", "not Autosomal dominant")
# 
# effect_inh_phenomatch_fig = ggplot(data = change_phenomatch, aes(y = change, x = "", fill = Dominant)) + 
#   geom_dotplot(stackdir = "center", binaxis = "y", stackgroups = T, binwidth = 0.05, binpositions = "all", dotsize = 0.2, stackratio = 0.8) +
#   labs(x = "", title = "Change in phenoMatchScore after using inheritance", y = "Change in phenoMatchScore") +
#   theme_bw()
# 
# pdf("Phenotypes/effect_inheritance_phenomatch.pdf", width = 20)
# effect_inh_phenomatch_fig
# dev.off()
