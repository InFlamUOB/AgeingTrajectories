#(Internal)Run in my R desktop in ukbb_ageonset-master project with folder Jan24FINAL_20230124_1752 as key results. Previous
#versions contained other GTEx analysis (ProximityAndEQTL.R)



# Need this file from 2_SummaryStats.R (ProximityWith8.tsv or Supplementary 4 table) - found here /Users/lxb732/Desktop/Multimorbidity_Paper_FINAL/Jan2023/ProximityWith8.tsv
#Also FUMA files generated using software (https://fuma.ctglab.nl/tutorial) , can be found appended too ( Healthy_Mix, Healthy_ONLY, Mix_Mix, Unhealthy_Mix, Unhealthy_ONLY)

library(tidyverse)
library(data.table)

mainDir <- "/Users/lxb732/Downloads/ukbb_ageonset-master"
NameRun <- 'Jan25'
subDir <- paste0(sub('\\..*', '', NameRun), format(Sys.time(), '_%Y%m%d_%H%M'))
dir.create(file.path(mainDir, subDir))
IncludeFigHere <- file.path(mainDir, subDir)
#IncludeFigHere <- "/Users/lxb732/Downloads/ukbb_ageonset-master/Jan18_20230120_1541"
#IncludeFigHere <- "/Users/lxb732/Downloads/ukbb_ageonset-master/Jan25_20230125_1333"
IncludeFigHere <- "/Users/lxb732/Downloads/ukbb_ageonset-master/Jan24FINAL_20230124_1752"

set.seed(132)

############### SNP associations description

#sig_res <-  read_tsv('/Users/lxb732/Desktop/Jan2023/ProximityWith8.tsv')

sig_res <-  read_tsv("/Users/lxb732/Desktop/Multimorbidity_Paper_FINAL/Jan2023/ProximityWith8.tsv")
# so associations of variants are unique by rsID but also by specific allele change - therefore to specify number of variants 
#should say unique SNPs and alleles!

SNPsPerCategory <- sig_res %>%
  mutate(Common= paste0(ID,"_", ALLELE0,"_", ALLELE1))

length(unique(SNPsPerCategory$Common)) #128 
length(unique(SNPsPerCategory$ID)) #128 
#one to one - no different alleless

#see <- SNPsPerCategory %>%
#  filter(ID == "rs174549") # one with two alleles

ForVenn <- data.frame(table( SNPsPerCategory$Common,  SNPsPerCategory$Variable)) %>%
  filter(Var1 != "") %>%
  mutate_if(is.factor, as.character) %>%
  mutate(Freq2 = ifelse(Freq > 0 , Var1, 0)) %>%
  dplyr::select(-Freq) %>%
  spread(Var2, Freq2) 

#see2 <- ForVenn %>%
#  filter(Var1 == see$Common)

ForVenn[ForVenn==0] <- NA

Venn <- list( PRP = na.omit(ForVenn$Healthy),  ToHealthy = na.omit(ForVenn$HealthyCross), NRN = na.omit(ForVenn$Unhealthy),  ToUnhealthy = na.omit(ForVenn$UnhealthyCross))

#see2 <- intersect(Venn$Healthy, see$Common)


library(ggVennDiagram)

VennPlot <- ggVennDiagram(Venn, label = "count", label_size = 4, edge_size = 5) +  
  scale_fill_gradient(low="blue",high = "red") + 
  theme_void() + 
  theme(legend.position = "none")

GenesPerCat <- process_region_data(Venn(Venn))

HealthySNPs <- SNPsPerCategory %>% filter(Common %in% GenesPerCat$item[[1]]) %>% select(ID) %>% unique()
UnhealthySNPs <- SNPsPerCategory %>% filter(Common %in% GenesPerCat$item[[2]]) %>% select(ID) %>% unique()
MixSNPs <- SNPsPerCategory %>% filter(Common %in% GenesPerCat$item[[3]]) %>% select(ID) %>% unique()

write.table(HealthySNPs$ID, paste0(IncludeFigHere, "/Healthy_SNPs.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(UnhealthySNPs$ID, paste0(IncludeFigHere, "/Unhealthy_SNPs.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(MixSNPs$ID, paste0(IncludeFigHere, "/Mix_SNPs.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

fwrite(GenesPerCat, paste0(IncludeFigHere, "/GenesPerCatInTableSNPs.txt"))
#GenesPerCat <- fread(paste0(IncludeFigHere, "/GenesPerCatInTable.txt"))

library(phenoscanner)

#GenesPerCat$item[[2]] #Unhealthy
#incredible!!!
Unhealthy <- phenoscanner(
  #snpquery =  c("rs174541", "rs174544", "rs174548", "rs174549", "rs174555", "rs174556",  "rs174557", "rs174560",     "rs174561", "rs174576",  "rs4246215", "rs751229834"), 
  snpquery = UnhealthySNPs$ID,
  # genequery = ,
  regionquery = NULL,
  catalogue = "eQTL", #instead of GWAS get eQTL info an dGWAS literature from database: https://www.ncbi.nlm.nih.gov/gap/phegeni
  pvalue = 1e-05,
  proxies = "None",
  r2 = 0.8,
  build = 37
)

SNPsUnhealthy <- Unhealthy$snps
ResultsUnhealthy <- Unhealthy$results

UnderstandBetas1 <- sig_res %>%
  filter(ID %in% UnhealthySNPs$ID  ) %>% 
  dplyr::rename(rsid = ID) %>%
  inner_join(SNPsUnhealthy %>% select(rsid, consequence, hgnc))

fwrite(UnderstandBetas1, paste0(IncludeFigHere, "/UnderstandBetasUnhealthy.txt"))

########### Healthy

Healthy1 <- phenoscanner(
  snpquery =   HealthySNPs$ID[1:80][grepl("rs", HealthySNPs$ID[1:80])], # genequery = ,
  regionquery = NULL,
  catalogue = "eQTL",
  pvalue = 1e-05,
  proxies = "None",
  r2 = 0.8,
  build = 37
)

Healthy2 <- phenoscanner(
  snpquery =   HealthySNPs$ID[81:108][grepl("rs", HealthySNPs$ID[81:108])], # genequery = ,
  regionquery = NULL,
  catalogue = "eQTL",
  pvalue = 1e-05,
  proxies = "None",
  r2 = 0.8,
  build = 37
)

SNPsHealthy <-  rbind(Healthy1$snps, Healthy2$snps)
ResultsHealthy <- rbind(Healthy1$results, Healthy2$results)

UnderstandBetas3 <- sig_res %>%
  filter(ID %in% HealthySNPs$ID)%>% 
  dplyr::rename(rsid = ID) %>%
  inner_join(SNPsHealthy %>% select(rsid, consequence, hgnc)) 

#checked in dbSNP 

UnderstandBetas3[75, "hgnc"] <- "RCCD1"
UnderstandBetas3[75, "consequence"] <- "intron"

write.table(HealthySNPs$ID, quote = FALSE, row.names = FALSE, col.names = FALSE,  paste0(IncludeFigHere, "/FUMA_Healthy.txt"))

fwrite(UnderstandBetas1, paste0(IncludeFigHere, "/UnderstandBetasUnhealthy.txt"))


####################### Mix
#GenesPerCat$item[[3]] #UnhealthyAndHealthy
#incredible!!!
MixUnAndHealthy <- phenoscanner(
  snpquery =  MixSNPs$ID,
  # genequery = ,
  regionquery = NULL,
  catalogue = "eQTL",
  pvalue = 1e-05,
  proxies = "None",
  r2 = 0.8,
  build = 37
)

SNPsMixUnAndHealthy <- MixUnAndHealthy$snps
ResultsMixUnAndHealthy <- MixUnAndHealthy$results

UnderstandBetas2 <- sig_res %>%
  filter(ID %in% MixSNPs$ID,
  ) %>% 
  dplyr::rename(rsid = ID) %>%
  inner_join(SNPsMixUnAndHealthy %>% select(rsid, consequence, hgnc))

fwrite(UnderstandBetas2, paste0(IncludeFigHere, "/UnderstandBetasHealthyAndUnhealthy.txt"))


##############


pdf(paste0(IncludeFigHere, "/VennSNPs_Names.pdf"), 10,5)
print(VennPlot)
dev.off()

################## GWAS Catalog from https://www.ncbi.nlm.nih.gov/gap/phegeni

MixGWAS <- fread("/Users/lxb732/Downloads/Mix_SNPS_GWAS.tab")
UnhealthyGWAS <- fread("/Users/lxb732/Downloads/Unhealthy_SNPS_GWAS.tab")
HealthyGWAS <- fread("/Users/lxb732/Downloads/Healthy_SNPS_GWAS.tab")

GWASCatalog <- list(MixGWAS %>% mutate(Type = "Mix"), 
                    UnhealthyGWAS %>% mutate(Type = "Unhealthy"), 
                    HealthyGWAS %>% mutate(Type = "Healthy")) %>%
  bind_rows()

fwrite(GWASCatalog, "GWASCatalogResults.csv")

#Gene Annotation ################
################
################
################

#before did not have the eQTL results associated to mix so only proxy ( 3 genes) would appear
eQTL <- list(ResultsUnhealthy %>% mutate(Type = "Unhealthy"), 
             ResultsHealthy %>% mutate(Type = "Healthy"),
             ResultsMixUnAndHealthy %>% mutate(Type = "Unhealthy"), 
             ResultsMixUnAndHealthy %>% mutate(Type = "Healthy")) %>% #mix is included in both!
  bind_rows()

write.table(eQTL, paste0(IncludeFigHere, "/eQTL_info.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE)

#Proximity 

gwas2GRanges <- function(gwasRes, SNP = "RefSNP_id", start = "BP", chr = "CHR", cols2retain = c('ALLELE0', 'ALLELE1'), genome='hg19' ){
  library(tidyverse)
  library(GenomicRanges)
  gwasRes <- gwasRes %>%
    dplyr::rename(RefSNP_id = SNP,
                  start = start,
                  chr = chr) %>%
    dplyr::mutate(end=`start`)
  gwasRes <- gwasRes %>%
    dplyr::select( RefSNP_id, chr, start, end, cols2retain)
  snpinf <- makeGRangesFromDataFrame(gwasRes)
  if (substr(seqlevels(snpinf),1,3)!='chr') {
    seqlevels(snpinf)=paste('chr',seqlevels(snpinf),sep='')
  }
  genome(snpinf)=genome
  for( colnm in c('RefSNP_id', cols2retain)){
    values(snpinf)[[colnm]]=gwasRes[[colnm]]
  }
  return(snpinf)
}

martx=biomaRt::useMart('ensembl','hsapiens_gene_ensembl')
genemap <- biomaRt::getBM(attributes = c('entrezgene_id','hgnc_symbol','ensembl_gene_id','description'),
                          mart = martx) %>%
  dplyr::rename(entrezgene = entrezgene_id)


library(dplyr)
#detach(package:plyr)

# Summary Stats GWAS ################
################
################
################


gwasResOrig <- read_tsv('/Users/lxb732/Desktop/Multimorbidity_Paper_FINAL/Jan2023/ProximityWith8.tsv') %>% #random summary stats?
  dplyr::rename(    "CHR" = "CHROM") %>%
  dplyr::rename(      "BP" = "GENPOS") %>%
  #  dplyr::rename(     "P" = "p") %>%
  #  dplyr::rename(  "A1" = "ALLELE1") %>%
  dplyr::rename( "SNP" = "ID")
#sometimes opposite side? with dplyr this is the correct one


allvarFinal <- list()
proxyFinal <- list()

for (i in unique(gwasResOrig$Variable)){
  
  print(i)
  gwasRes <- gwasResOrig %>%
    filter(Variable == i)
  
  #gwasRes <- read_tsv('data/processed/ukbb/gwas/bolt/a1092.imp.stats')
  gwas_as_GR <- gwas2GRanges(gwasRes, SNP = 'SNP',start = 'BP',chr = 'CHR',genome = 'hg19')
  
  library(VariantAnnotation)
  library(GenomicRanges)
  
  txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  varInfo <- unique(granges(gwas_as_GR))
  allvar <- locateVariants(varInfo, txdb,
                           AllVariants(promoter=PromoterVariants(upstream=1000,
                                                                 downstream=1000)))
  
  allvar <- allvar[which(!is.na(allvar$GENEID)),]
  overs <- findOverlaps(gwas_as_GR,allvar)
  for(colnm in colnames(mcols(gwas_as_GR))){
    mydat <-(mcols(gwas_as_GR)[[colnm]])[queryHits(overs)]
    mcols(allvar)[[colnm]] = '*'
    mcols(allvar)[[colnm]][subjectHits(overs)] <- mydat
  }
  geneids <- unique(allvar$GENEID)
  geneids <- geneids[complete.cases(geneids)]
  genemap <- genemap %>% unique() %>%
    mutate(entrezgene=as.character(entrezgene))
  allvar <- as.tibble(allvar) %>% #spliceSite intron fiveUTR threeUTR coding intergenic promoter
    filter(LOCATION!='intergenic') %>%  # no intergenic either way
    mutate(entrezgene=GENEID)%>%
    left_join(genemap) %>%
    mutate(Type = i)
  
  allvarFinal[[i]] <- allvar
  
  #two stages, all variants and then map to entrez genes! Docuement both 
  
  #rm(list=setdiff(ls(),c('allvar','gwas_as_GR','gwasRes','gwas2Granges','genemap')))
  #saveRDS(allvar,paste0("./allvar_",i,".rds"))
  proxyres <- allvar %>%
    dplyr::rename(SNP = RefSNP_id,
                  CHR = seqnames,
                  BP = start,
                  proxy_entrez = entrezgene,
                  proxy_hgnc = hgnc_symbol,
                  proxy_ensembl = ensembl_gene_id,
                  proxy_type = LOCATION,
                  Ref = ALLELE1,
                  Alt = ALLELE0) %>%
    dplyr::mutate(CHR = gsub('chr','',CHR)) %>%
    dplyr::select( SNP, CHR, BP, Ref, Alt, proxy_entrez, proxy_ensembl, proxy_hgnc, proxy_type, Type) %>%
    unique()
  print('proxy mapping is done')
  # system('mkdir ./data/processed/genomicAnalysis')
  #saveRDS(proxyres, file=paste0("./proxyRes_",i,".rds"))
  print('proxy file is saved')
  
  proxyFinal[[i]] <- proxyres
  #rm(allvar)
  #rm(proxyres)
  #rm(gwas_as_GR)
  
}

#proxy analysis for each category

Proxy <- proxyFinal %>%
  bind_rows()

#fwrite(Proxy, paste0(IncludeFigHere, "/Proxy8.csv"))
Proxy <- fread(paste0("/Users/lxb732/Downloads/ukbb_ageonset-master/Jan24FINAL_20230124_1736/Proxy8.csv"))

Vars <- allvarFinal %>%
  bind_rows()

fwrite(Vars, paste0(IncludeFigHere, "/Vars8.csv"))

#################### Merge eQTL and Proximity


#eQTL <-  fread(paste0(IncludeFigHere, "/eQTL_info.txt"))
#Proxy <- fread(paste0(IncludeFigHere, "/Proxy8.csv"))

MergeAll <- list( #Info %>% dplyr::select(SNP = Uploaded_variation, Variable, gene) %>% mutate(Type = "VEP"), 
  eQTL = eQTL %>% dplyr::select(SNP = rsid,  gene = exp_gene, Variable = Type, ALLELE0 = a2, ALLELE1 = a1) %>% mutate(Type = "eQTL"), #changes a2 and a1 to coin cide with proxy 
  Proximity = Proxy %>% dplyr::select(SNP, gene = proxy_hgnc, Variable = Type, ALLELE0 = Ref, ALLELE1 = Alt) %>% mutate(Type = "Proxy") ) %>%
  bind_rows() %>%
  unique() %>%
  filter(gene != "-") %>%
  separate_rows(gene, sep = "[;]") %>%
  mutate(Common = paste0(SNP,"_", ALLELE0, "_", ALLELE1)) %>% #other way round? 
  mutate(gene = toupper(gene)) %>%  #C11orf10
  mutate(gene = gsub("GPR44", "PTGDR2", gene)) %>% #according to gene cards - as not recognized by FUMA
  mutate(gene = gsub("C11ORF10", "TMEM258", gene)) %>%
  mutate(gene = gsub("NECTIN2", "PVRL2", gene))


fwrite(MergeAll, paste0(IncludeFigHere, "/MergeAll8.csv"))
MergeAll <- fread(paste0(IncludeFigHere, "/MergeAll8.csv"))
#MergeAll <- fread("/Users/lxb732/Downloads/ukbb_ageonset-master/Jan18_20230120_1541/MergeAll8.csv")
############ understand SNP to gene ratio


#real unique associations - 128 associations from GWAS
length(unique(MergeAll$SNP)) #125  - lost 3 because no mapping to genes? 
setdiff( SNPsPerCategory$ID, MergeAll$SNP) #9:45413234_GGT_G, rs157592, 15:91508579_AG_A
length(unique(MergeAll$Common)) #166 - ALLELE1 or ALLELE0 issue? but in dbSNP just one allele
#setdiff(MergeAll$Common, SNPsPerCategory$Common) 

SNPsPerGeneType <- MergeAll %>% 
  filter(gene != "") %>%
  dplyr::select(gene, Type) %>% 
  unique() %>%
  group_by(gene) %>%
  add_tally()

length(unique(SNPsPerGeneType$gene)) #100

GenesPaper <- SNPsPerGeneType %>% 
  filter(Type == "Proxy") 

write.table(GenesPaper$gene , paste0(IncludeFigHere, "/GenesPaper.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE, sep = ",")

PRC1 <- MergeAll %>%
  filter(gene == "PRC1")

length(unique(PRC1$SNP)) #78


SNPsPerGene <-  data.frame(table( MergeAll$gene,  MergeAll$Variable)) %>%
  filter(Var1 != "") %>%
  mutate_if(is.factor, as.character) %>%
  mutate(Freq2 = ifelse(Freq > 0 , Var1, 0)) %>%
  filter(Freq != 0) 

fwrite(SNPsPerGene, paste0(IncludeFigHere, "/SNPsPerGene_All.csv"))

ForVenn <- data.frame(table( MergeAll$gene,  MergeAll$Variable)) %>%
  filter(Var1 != "") %>%
  mutate_if(is.factor, as.character) %>%
  mutate(Freq2 = ifelse(Freq > 0 , Var1, 0)) %>%
  dplyr::select(-Freq) %>%
  spread(Var2, Freq2) 

ForVenn[ForVenn==0] <- NA

Venn <- list( PRP = na.omit(ForVenn$Healthy),  ToHealthy = na.omit(ForVenn$HealthyCross), NRN = na.omit(ForVenn$Unhealthy),  ToUnhealthy = na.omit(ForVenn$UnhealthyCross))


VennPlot <- ggVennDiagram(Venn, label = "count", label_size = 4, edge_size = 5) +  
  scale_fill_gradient(low="blue",high = "red") + 
  theme_void() + 
  theme(legend.position = "none")


GenesPerCat <- process_region_data(Venn(Venn))

fwrite(GenesPerCat, paste0(IncludeFigHere, "/GenesPerCatIn_All.csv"))

pdf(paste0(IncludeFigHere, "/Venn_All_ChangeNames.pdf"), 10,5)
print(VennPlot)
dev.off()

AllGenes <- c(as.character(Venn[["Healthy"]]), as.character(Venn[["ToHealthy"]]),  as.character(Venn[["Unhealthy"]]),  as.character(Venn[["ToUnhealthy"]])) %>% setdiff(.,"0")

write.table(GenesPerCat$item[[1]], paste0(IncludeFigHere, "/Healthy_Genes.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(GenesPerCat$item[[2]], paste0(IncludeFigHere, "/Unhealthy_Genes.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

write.table(GenesPerCat$item[[3]], paste0(IncludeFigHere, "/Mix_Genes.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

#MergeAll8 to ageing databases!

#####################
######################

#Done!

########################
########################

#map it with FUMA results. 

#- get genes.txt


#redo but with GO processes - marked differences immunological versus metabolic signature??!!

FUMAPlots <-  function(FUMA, Title){
  
  
  SupplementaryPlot <- FUMA %>%  #no Wikipathways, GO_mf, Immunologic_signatures
    #filter(Category %in% c( "KEGG", "TF_targets","Reactome", "GO_bp" )) %>%
    mutate(GeneSet = gsub("KEGG_", "", GeneSet)) %>%
    mutate(GeneSet = gsub("REACTOME_", "", GeneSet)) %>% 
    mutate(GeneSet = gsub("GO_", "", GeneSet)) %>% 
    mutate(GeneSet = gsub("_", " ", GeneSet)) %>%
    mutate(GeneSet = tolower(GeneSet)) %>%
    arrange(desc(-adjP)) %>% 
    group_by(Category) %>% 
    dplyr::slice(1:12)
  
  pdf(paste0(IncludeFigHere, "/Categories_", Title, ".pdf"), 10, 5)
  
  for (i in unique(SupplementaryPlot$Category)){
    
    SupplementaryPlotFil <-  SupplementaryPlot %>% 
      filter(Category == i)
    
    print(
      ggplot(SupplementaryPlotFil, aes( fct_reorder(GeneSet, -log10(adjP)), y = -log10(adjP),   size = N_overlap/N_genes)) + 
        geom_point() + 
        coord_flip() + 
        theme_bw() + 
        labs(y = "logAdjP", x = "GeneSet", title = paste0( i," from ", Title)) +  #or instead Fold Enrichment for logAdjP
        guides( color = FALSE) 
    )
  }
  
  dev.off()
  
  
  
  ######3 plot of enriched datasets!!
  ReducedFUMA <- FUMA %>%  #no Wikipathways, GO_mf, Immunologic_signatures
    filter(Category %in% c( "Canonical_Pathways", "TF_targets", "GO_bp" )) %>% #c( "Canonical_Pathways", "Wikipathways","Reactome", "GO_bp" )
    mutate(GeneSet = gsub("KEGG_", "KEGG: ", GeneSet)) %>%
    mutate(GeneSet = gsub("REACTOME_", "REACTOME: ", GeneSet)) %>% 
    mutate(GeneSet = gsub("GO_", "", GeneSet)) %>% 
    mutate(GeneSet = gsub("BIOCARTA_", "BIOCARTA: ", GeneSet)) %>% 
    mutate(GeneSet = gsub("_", " ", GeneSet)) %>%
    mutate(GeneSet = tolower(GeneSet)) %>%
    arrange(desc(-adjP)) %>% 
    group_by(Category) %>% 
    dplyr::slice(1:20)
  
  Plot1 <-  ggplot(ReducedFUMA, aes( fct_reorder(GeneSet, -log10(ReducedFUMA$adjP)), y = -log10(ReducedFUMA$adjP),
                                     size = N_overlap/N_genes)) + #color = Category
    geom_point() + 
    coord_flip() + 
    facet_wrap(.~Category, scales = "free") + 
    theme_bw() + 
    labs(y = "logAdjP", x = "GeneSet") +  #or instead Fold Enrichment for logAdjP
    guides( color = FALSE)  #size = FALSE
  
  #plt.tight_layout(Plot1)
  #pdf("Plot1.pdf", 20, 10)
  #print(Plot1)
  #dev.off()
  
  #ReducedFUMA2 <- FUMA %>% 
  #  filter(Category %in% c( "GWAScatalog" )) #Wikipathways
  #
  #
  #Plot2 <-  ggplot(ReducedFUMA2 %>% arrange(adjP) %>% slice_head(n = 12), 
  #                 aes( fct_reorder(GeneSet, -log10(adjP)), y = -log10(adjP), 
  #                      size = N_overlap/N_genes), color = Category) + #color = N_overlap/N_genes)
  #  geom_point() + 
  #  coord_flip()+
  #  labs(title = "GWAS Catalog", y = "logAdjP", x = "GeneSet" #or instead Fold Enrichment for logAdjP
  #  ) +
  #  theme_bw() 
  
  library(patchwork)
  
  
  
  #pdf(paste0(IncludeFigHere, "/PathwaysOnly_",Title,".pdf"), 35, 12)
  #print( Plot1)
  #dev.off()
  
  return(Plot1) #+ Plot2 +  plot_layout(widths = c(3, 1))
  
}





# Mix
Mix <- fread(paste0(IncludeFigHere, "/Mix_Mix/GS.txt"))
MixPlot  <- FUMAPlots(Mix,"MixI") #+  labs(title = "GWAS Catalog")

# HealthyWithMix

FUMAHeMix <- fread(paste0(IncludeFigHere, "/Healthy_Mix/GS.txt"))
FUMAHeMixPlot  <- FUMAPlots(FUMAHeMix,"FUMAHeMixI") #+  labs(title = "GWAS Catalog")

# UnhealthyWithMix

FUMAUnHeMix <- fread(paste0(IncludeFigHere, "/Unhealthy_Mix/GS.txt"))
FUMAUnHeMixPlot  <- FUMAPlots(FUMAUnHeMix,"FUMAUnHeMixI") #+  labs(title = "GWAS Catalog")

#HealthyOnly 

Healthy <- fread(paste0(IncludeFigHere, "/Healthy_ONLY/GS.txt"))
HealthyPlot  <- FUMAPlots(Healthy,"Healthy_OnlyI") #+  labs(title = "GWAS Catalog")


#UnhealthyOnly 

UnHealthy <- fread(paste0(IncludeFigHere, "/Unhealthy_ONLY/GS.txt"))
UnhealthyPlot  <- FUMAPlots(UnHealthy,"Unhealthy_OnlyI") #+  labs(title = "GWAS Catalog")

################## Plots


pdf(paste0(IncludeFigHere, "/MixNOI.pdf"), 22, 15)
print( HealthyPlot/UnhealthyPlot +  plot_layout(guides = 'collect'))
dev.off()

pdf(paste0(IncludeFigHere, "/MixJUstI.pdf"), 22, 15)
print( MixPlot +  plot_layout(guides = 'collect'))
dev.off()

pdf(paste0(IncludeFigHere, "/MixYESI.pdf"), 22, 15)

print( FUMAHeMixPlot/FUMAUnHeMixPlot +  plot_layout(guides = 'collect') )

dev.off()

pdf(paste0(IncludeFigHere, "/MixFinal1Thesis.pdf"), 22, 4)

print( FUMAHeMixPlot)

dev.off()

pdf(paste0(IncludeFigHere, "/MixFinal2Thesis.pdf"), 22, 4)

print(FUMAUnHeMixPlot )

dev.off()

#Info on TF 



AllFUMAResults <- list(Mix %>% mutate(Type = "Mix_Only"), 
                       FUMAHeMix %>% mutate(Type = "HealthyAndMix"),
                       FUMAUnHeMix %>% mutate(Type = "UnhealthyAndMix"), 
                       Healthy %>% mutate(Type = "Healthy_Only"), 
                       UnHealthy %>% mutate(Type = "Unhealthy_Only")) %>% 
  bind_rows()

fwrite(AllFUMAResults, "AllFUMAResults.csv")

TF_targets <- AllFUMAResults %>% 
  filter(Category == "TF_targets") #%>% 
#filter(Type == "UnhealthyAndMix")

HRH <- AllFUMAResults %>% 
  filter(Type == "HealthyAndMix") %>%
  arrange(adjP)

HRH$genes
############

see <- FUMAUn %>% 
  filter(Category %in% c( "KEGG", "TF_targets","Reactome", "Wikipathways" )) %>%
  select(c("Category", "GeneSet", "adjP", "genes")) %>%
  separate_rows(., genes, convert = TRUE, sep = "[:]") %>%
  rename(gene = genes) %>%
  filter(gene != "")

############### claim on different betas and categories

sig_res <-  read_tsv('/Users/lxb732/Desktop/Multimorbidity_Paper_FINAL/Jan2023/ProximityWith8.tsv')
MergeAll <- fread(paste0(IncludeFigHere, "/MergeAll8.csv"))

AgeGenesLDLR <- MergeAll %>%
  filter(gene %in% c( "LDLR")) %>%
  select(ID = SNP, gene, Type) %>% #VarGene = Variable
  left_join(., sig_res %>% 
              select(VarSumm = Variable, ID, BETA, SE, LOG10P))

unique(AgeGenesLDLR$ID) #11

AgeGenesTOMM40 <- MergeAll %>%
  filter(gene %in% c("TOMM40")) %>%
  select(ID = SNP, gene, Type) %>% #VarGene = Variable
  left_join(., sig_res %>% 
              select(VarSumm = Variable, ID, BETA, SE, LOG10P))

unique(AgeGenesTOMM40$ID) #28
unique(filter(AgeGenesTOMM40, Type == "eQTL")$ID) #28
unique(filter(AgeGenesTOMM40, Type == "Proxy")$ID) #10

########## Longevity  - add the beta stats myself

Zenin <- sig_res %>%
  filter(ID %in% c("rs4420638", "rs429358", "rs2075650")) 

APOE <- sig_res %>%
  filter(ID %in% filter(MergeAll,gene %in% c("APOE"))$SNP ) 

