# Data Parsing
lapply(c("data.table","devtools","taigr","VariantAnnotation","readr","stringr"),
       library, character.only = T)

#read in clinvar variants 
clinvar <- fread("/Users/jaredcollins/Downloads/variant_summary.txt")
clinvar_snv <- clinvar[clinvar$Type == "single nucleotide variant",]

# Read in breast cancer GWAS data
GWAS <- fread("/Users/jaredcollins/Downloads/EFO_0000305_associations_export.tsv")

# read in all COSMIC noncoding variant info
COSMIC_noncoding <- fread("/Users/jaredcollins/Downloads/CellLinesProject_NonCodingVariants_Tsv_v100_GRCh38/CellLinesProject_NonCodingVariants_v100_GRCh38.tsv")

#subset noncoding variants for MCF7 variants
COSMIC_noncoding_sub <- COSMIC_noncoding[,c(3,8,9,10,11,12,16)]
COSMIC_noncoding_sub <- COSMIC_noncoding_sub[COSMIC_noncoding_sub$COSMIC_SAMPLE_ID %in% "COSS905946",]
#save the subset
write.table(COSMIC_noncoding_sub, file = "/Users/jaredcollins/Dropbox (Partners HealthCare)/Broad/broadhacks24/COSMIC_noncoding_MCF7.tsv", row.names=FALSE, sep="\t")

#read back in the subset to check it works 
MCF7_noncoding <- fread("/Users/jaredcollins/Dropbox (Partners HealthCare)/Broad/broadhacks24/COSMIC_noncoding_MCF7.tsv")

#read in and format the risk alleles overlapping with Hi-C data
riskAlleles <- fread("/Users/jaredcollins/Downloads/riskAlleles.txt",sep="}")
riskAlleles = as.data.frame(t(as.data.frame(riskAlleles)))
names(riskAlleles)<-"col1"
rm(risk_redone)
risk_redone <- data.frame(do.call('rbind', strsplit(riskAlleles$col1,',',fixed=TRUE)))

#screen GWAS-HiC overlap hits in MCF7 VCF and ClinVar
hits<-NULL
for(ID in COSMIC_noncoding_sub$GENOME_START){
  hits<-append(hits,grep(ID,risk_redone$X2))
}

hits_clinvar<-NULL
for(rsID in risk_redone$X2){
  rsID<-str_sub(rsID,start=16,end=-4)
  hits_clinvar<-append(hits_clinvar,grep(rsID,clinvar$`RS# (dbSNP)`))
}

