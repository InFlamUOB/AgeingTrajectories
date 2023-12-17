#(Internal) Run as slurm job (RuthEthnicityScript3GWAS.R by RuthEthnicity2.sh)


#!/usr/bin/env Rscript

.libPaths( c( .libPaths(), "/rds/homes/l/lxb732/R/library/4.0.0/haswell") )



library(tidyverse)
library(data.table)
library(ggrepel)
library(qqman)
library(RSQLite)
options(bitmapType='cairo')

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)


fn=commandArgs(trailingOnly=TRUE)

print(fn)

IncludeFigHere <- fn

print(IncludeFigHere)

#IncludeFigHere <- "/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930"


filename_assoc_results=paste0(fn,".txt")
graphical_output=".png"
myYlim = "" #max -log P value for the yaxis - any out of bounds points will be brought down to within the plot. Default is top value plus 2
num_anno = 20 #number of annotations (below p = 1e-5)
title_tail = fn
man_title = paste0("Manhattan Plot for UK biobank cohort for ",title_tail)
qq_title = paste0("QQ Plot for UK biobank cohort for ",title_tail)
file_name_out=paste0(fn,graphical_output) #ggsave for manhattan plot so will change format according to string - not true for QQ plot this is hard coded png below


#setwd(paste0(IncludeFigHere, "/results/"))
file_list <-  list.files(path = paste0(IncludeFigHere, "/results/"), pattern = paste0("*",".regenie")) #taken fn from here?

print(file_list[1])

ldf <- lapply(paste0(IncludeFigHere, "/results/",file_list), fread) #got warnings?

names(ldf) <- file_list %>%
  str_replace_all(., "ukb_step2_BT","") %>% #ukb_step2_Laura before
  str_replace_all(., ".regenie","") %>%
  str_replace_all(., "_Healthy","&Healthy") %>% #change clusters here
  str_replace_all(., "_Unhealthy","&Unhealthy")

ldf2 <- ldf %>%
  bind_rows(.id = "Names") %>%
  separate(Names, c("Cont", "Variable"), sep = "&") %>%
  #dplyr::select(-Eliminate) %>%
  mutate( p = 10^(-LOG10P) )

ldf3 <- ldf2 %>%
  filter(A1FREQ > 0.01 & INFO > 0.3)

s <- filter(ldf3, Variable == "HealthyCross" ) %>%
  unique()

s1 <- s %>% 
  filter(LOG10P > 5) 

s2 <- s %>% 
  filter(LOG10P > 8) 

saveRDS(ldf3, paste0(IncludeFigHere,"/ldf3.RDS")) #maybe better the fread and fwirte? 

#ldf3 <- readRDS("/rds/projects/g/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930/ldf3.RDS")
manPlot <- list()



for (j in levels(as.factor(ldf3$Variable))) {  #"Healthy"        "HealthyCross"   "Unhealthy"      "UnhealthyCross"
  
  print(paste0(".................... ", j))
  ldf4 <- ldf3 %>%
    filter(Variable == j) %>%
    select(-c(Cont, Variable))
  
  #for each 9,946,350 associations!
  
  #CHROM GENPOS          ID ALLELE0 ALLELE1   A1FREQ     INFO TEST        BETA        SE     CHISQ   LOG10P          p
  #1:     1  10177 rs367896724      AC       A 0.599062 0.467098  ADD  0.00510208 0.0110620 0.2127290 0.190685 0.64463666
  #2:     1  10352 rs201106462      TA       T 0.605565 0.447098  ADD  0.00182131 0.0113843 0.0255952 0.059039 0.87289298
  #3:     1  11008 rs575272151       G       C 0.914028 0.492883  ADD -0.00711804 0.0188852 0.1420630 0.151048 0.70623949
  #4:     1  11012 rs544419019       G       C 0.914028 0.492883  ADD -0.00711804 0.0188852 0.1420630 0.151048 0.70623949
  #5:     1  13110 rs540538026       A       G 0.942804 0.389850  ADD -0.05661640 0.0256201 4.8833900 1.566770 0.02711627
  #6:     1  13116  rs62635286       G       T 0.813275 0.407047  ADD -0.03633350 0.0150192 5.8522400 1.808070 0.01555715
  
  if(num_anno > 0){
    sig = ldf4 %>% 
      filter(LOG10P > 5) %>%
      group_by(ID) %>% #no variable
      slice_max(n=1, order_by = desc(LOG10P)) %>%
      group_by(CHROM) %>%
      arrange(desc(LOG10P)) %>%
      mutate(keep = NA) %>%
      ungroup()
    i = 1
    counter = 1
    var_list = c()
    while (i <= num_anno){
      snp = sig$ID[counter]
      c = sig$CHROM[sig$ID == snp]
      p_low = sig$GENPOS[sig$ID == snp] - 2000000
      p_up = sig$GENPOS[sig$ID == snp] + 2000000
      top_of_list = sig %>%
        filter(CHROM == c & GENPOS > p_low & GENPOS < p_up) %>%
        slice_max(n = 1, order_by = LOG10P)
      if (snp == top_of_list$ID){
        sig$keep[sig$ID == snp] = "yes"
        i = i+1
      }
      counter = counter + 1
    }
    
    
    sig = sig %>%
      filter(keep == "yes") %>%
      dplyr::select(ID)
    
    
    #to_anno = rbind(cassi1,cassi2,sig) %>% mutate(is_annotate="yes")
    to_anno = rbind(sig) %>% mutate(is_annotate="yes")
    
    chroms = ldf4$CHROM[ldf4$ID %in% sig$ID]
    snp_int = paste(sig$ID, collapse = '","')
    snp_int = gsub('^','"',snp_int)
    snp_int = gsub('$','"',snp_int)
    
    mydb <- dbConnect(RSQLite::SQLite(), "/rds/projects/k/karwatha-rate-af/datacentre/GWAS/annotation/imputed/vep_imputed.sqlite") #Dom maintenance!
    
    for (c in unique(chroms)){
      print(c)
      if (c == chroms[1]){
        an1 = dbGetQuery(mydb, paste0('SELECT "#Uploaded_variation", SYMBOL, Consequence, NEAREST FROM chr',c,' WHERE "#Uploaded_variation" IN (',snp_int,')'))
      } else {
        an2 = dbGetQuery(mydb, paste0('SELECT "#Uploaded_variation", SYMBOL, Consequence, NEAREST FROM chr',c,' WHERE "#Uploaded_variation" IN (',snp_int,')'))
        an1 = rbind(an1,an2)
      }
      print(dim(an1))
    }
    
    dbDisconnect(mydb)
    
    colnames(an1) = c("Uploaded_variation","gene","Consequence", "NEAREST")
    
    
    ukb_anno = an1
    #ukb_anno = fread("/castles/nr/projects/2017/gkoutosg-variant-prediction/UKBiobankPlink/genetic_variants_thin.tsv", sep ="\t") 
    
    ukb_anno = ukb_anno %>%
      mutate(gene = ifelse(Consequence == "intergenic_variant", paste0(NEAREST,"*"), NEAREST)) %>% 
      dplyr::select(ID = Uploaded_variation, gene) %>%
      group_by(ID) %>%
      slice_head(n=1) %>%
      ungroup
    anno = to_anno %>% left_join(ukb_anno, by = "ID")
  } else {
    anno = data.frame(ID = "rsX", gene = NA)
  }
  
  ldf4 = ldf4 %>% 
    left_join(anno, by = "ID") %>%
    mutate(is_annotate = ifelse(is.na(gene),"no","yes"))
  
  sig_res = ldf4 %>%
    dplyr::select(ID,gene,LOG10P) %>%
    filter(LOG10P>5)
  
  write_tsv(sig_res, paste0(IncludeFigHere, "/sig_res_",j,".tsv"))
  
  
  don = ldf4 %>%
    
    # Compute chromosome size
    group_by(CHROM) %>%
    summarise(chr_len=max(GENPOS)) %>%
    
    # Calculate cumulative position of each chromosome
    mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
    dplyr::select(-chr_len) %>%
    
    # Add this info to the initial dataset
    left_join(ldf4, ., by=c("CHROM"="CHROM")) %>%
    
    # Add a cumulative position of each SNP
    arrange(CHROM, GENPOS) %>%
    mutate( BPcum=GENPOS+tot)
  
  axisdf = don %>% group_by(CHROM) %>% summarize(center=( max(BPcum) + min(BPcum) )/ 2 )
  if(!myYlim==""){
    y_limit <- myYlim
  }else{
    y_limit <- max(ldf4$LOG10P[!is.na(ldf4$LOG10P)]) + 2
  }
  
  don <- don %>% mutate(
    out_of_bounds = (LOG10P - 1) > y_limit,
    LOG10P = ifelse(out_of_bounds, (y_limit -1), LOG10P)
  )
  
  saveRDS(don, paste0(IncludeFigHere, "/don_",j,".RDS"))
  
  
  manPlot[[j]] <- ggplot(don, aes(x=BPcum, y=LOG10P)) +
    
    # Show all points
    geom_point( aes(color=as.factor(CHROM)), alpha=0.75, size=0.001) +
    scale_color_manual(values = rep(c("aquamarine4", "slateblue"), 22 )) +
    
    # custom X axis:
    scale_x_continuous( label = axisdf$CHROM, breaks= axisdf$center ) +
    scale_y_continuous(expand = c(0, 0), limits=c(0,y_limit) ) +     # remove space between plot area and x axis
    geom_hline(yintercept=-log10(5E-8), linetype="dashed", color = "red") +
    # Add label using ggrepel to avoid overlapping
    geom_label_repel( data=subset(don, is_annotate=="yes"), aes(label=gene), size=1.5) +
    
    # Custom the theme:
    theme_bw() +
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(), 
      axis.text.x = element_text(color = "grey20", size = 8, angle = 90, hjust = .5, vjust = .5, face = "plain")
    ) +
    xlab("Chromosome") +
    ylab("-log p") 
  #ggtitle(man_title) +
  
  #+ facet_wrap(.~ Variable, scales = "free") #+
  #  geom_line(data = don[!is.na(don$grp),], aes(group = grp), size = 0.25, colour = "gold2")
  
  #pdf(paste0("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/Ruth_20210914_2347/ManhattanClust_",j,".pdf"))
  #print(manPlot)
  #dev.off()
  
  
  ggsave(paste0(IncludeFigHere,"/Manhattan_",j,".png"), plot = manPlot, device = "png")
  
  #png(paste0("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/Ruth_20210914_2347/Manhattan_",j,".png"), width=12, height=7, units="in",res=300)
  #print(manPlot)
  #dev.off()
  
  
  
}

save(manPlot, file = paste0(IncludeFigHere, "/manPlot.RData"))

library(patchwork)
#png(paste0(IncludeFigHere, "/MannPlotAll",  format(Sys.time(), '_%Y%m%d_%H%M'),".png"), width = 1200, height = 1800, units = "px")
Merge <- (manPlot[[1]] + manPlot[[2]]) / (manPlot[[3]]  + manPlot[[4]]) +
  plot_annotation(tag_levels = 'A')

ggsave(paste0(IncludeFigHere,"/Manhattan_",  format(Sys.time(), '_%Y%m%d_%H%M'),".png"), 
       plot = Merge, width = 13, height = 14, dpi = 600, units = "in", device = "png")



#gets significant IDS


#ToEnrich <- ldf3  %>% 
# distinct(gene, Variable, .keep_all = T)


sig = ldf3 %>% 
  filter(LOG10P > 8) %>% #normally 8! Had 5
 # group_by(ID) %>%
  #slice_max(n=1, order_by = desc(LOG10P)) %>% #noooooo - not just want to keep 1?? - no difference at all? 
 # group_by(CHROM) %>%
  arrange(desc(LOG10P)) %>%
  select(ID)


ToProximity <- ldf3 %>% 
  filter(LOG10P > 8) %>% #normally 8! Had 5
  # group_by(ID) %>%
  #slice_max(n=1, order_by = desc(LOG10P)) %>% #noooooo - not just want to keep 1?? - no difference at all? 
  # group_by(CHROM) %>%
  arrange(desc(LOG10P))

write_tsv(ToProximity, paste0(IncludeFigHere, "/ProximityWith8.tsv")) 
#ToProximity <- read_tsv(paste0(IncludeFigHere, "/Proximity.tsv")) 
#length(unique(ToProximity$ID))

#ldf3 <- ToProximity

chroms = ldf3$CHROM[ldf3$ID %in% sig$ID]
snp_int = paste(sig$ID, collapse = '","')
snp_int = gsub('^','"',snp_int)
snp_int = gsub('$','"',snp_int)

library(dplyr)
library(dbplyr)


#load annotation db and find annotations
mydb <- dbConnect(RSQLite::SQLite(), "/rds/projects/k/karwatha-rate-af/datacentre/GWAS/annotation/imputed/vep_imputed.sqlite")

src_dbi(mydb)
surveys <- tbl(mydb, "chr1")

for (c in unique(chroms)){
  print(c)
  if (c == chroms[1]){
    an1 = dbGetQuery(mydb, paste0('SELECT "#Uploaded_variation", SYMBOL, Consequence, NEAREST FROM chr',c,' WHERE "#Uploaded_variation" IN (',snp_int,')'))
  } else {
    an2 = dbGetQuery(mydb, paste0('SELECT "#Uploaded_variation", SYMBOL, Consequence, NEAREST FROM chr',c,' WHERE "#Uploaded_variation" IN (',snp_int,')'))
    an1 = rbind(an1,an2)
  }
  print(dim(an1))
}

dbDisconnect(mydb)

colnames(an1) = c("Uploaded_variation","gene","Consequence", "NEAREST")

#add annotations
ukb_anno = an1

#rs8111416 what about the genes associated? dbSNP gives none!

#just coding genes!!
#ukb_anno = ukb_anno %>%
#  mutate(gene = ifelse(Consequence == "intergenic_variant", paste0(NEAREST,"*"), NEAREST)) %>% 
#  select(ID = Uploaded_variation, gene) %>%
#  group_by(ID) %>%
#  slice_head(n=1) %>% #whyyy?
#  ungroup()

AllAnos = ukb_anno %>%
  mutate(gene = ifelse(Consequence == "intergenic_variant", paste0(NEAREST,"*"), NEAREST)) #%>% 
#select(ID = Uploaded_variation, gene) %>%
#group_by(ID) %>%
#slice_head(n=1) %>% #whyyy?
#ungroup

write_tsv(AllAnos, paste0(IncludeFigHere, "/AllAnnotationsVEPWith8.tsv")) #not created here!

sig_res = ldf3 %>% 
  filter(LOG10P>5) %>% #before 5
  left_join(ukb_anno, by = "ID")  %>%
  select(Variable,CHROM, ID,gene,GENPOS,LOG10P) #%>%
  #filter(LOG10P>5) same

#################3write_tsv(sig_res, paste0(IncludeFigHere, "/AllClust_all_annos.tsv")) #not created here!
write_tsv(sig_res, paste0(IncludeFigHere, "/AllClust_all_annos2.tsv")) #not created here!
AllClust_all_annos2 <- read_tsv(paste0(IncludeFigHere, "/AllClust_all_annos2.tsv"))


#sig_res <- read_tsv("AllClust_all_annos.tsv")

sig_res2 <-  sig_res %>% 
  distinct(gene, Variable, .keep_all = T)


sig_res2 %>%
  select(Variable, gene) %>%
  group_by(Variable) %>%
  count()

table(sig_res2$Variable, sig_res2$gene)
#venn diagram

#ForVenn <- data.frame(table( sig_res2$gene,  sig_res2$Variable)) %>%
#  mutate_if(is.factor, as.character) %>%
#  mutate(Freq2 = ifelse(Freq == 1, Var1, 0)) %>%
#  spread(Var2, Freq2) %>%
#  filter(Freq == 1)
#
#
##All <- list(Control= ForVenn$Control, Baseline = ForVenn$DataBase,  NAFLD = ForVenn$NAFLD,  #TransplantMaster = ForVenn$TransplantMaster, TransplantPostTx = ForVenn$TransplantPostTx, #TransplantPreTx = ForVenn$TransplantPreTx )
#
#All <- list( Clust1 = na.omit(ForVenn$Clust1),  Clust2 = na.omit(ForVenn$Clust2), Clust3 = na.omit(ForVenn$Clust3),  Clust4 = na.omit(ForVenn$Clust4),  Clust5 = na.omit(ForVenn$Clust5 ))
#
#saveRDS(All,paste0(IncludeFigHere, "/ForVenn.RDS"))


#VennDiagram.R file in Desktop ( cannot download functions here)

#VennDiagramList <- readRDS("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/#AllComorbiditiesIssue_20210930_2242/Ruth_20210914_2347/VennDiagramList.RDS")





#############################QQ plot 



fn=commandArgs(trailingOnly=TRUE)

filename_assoc_results=paste0(fn,".txt")
graphical_output=".png"
myYlim = "" #max -log P value for the yaxis - any out of bounds points will be brought down to within the plot. Default is top value plus 2
num_anno = 20 #number of annotations (below p = 1e-5)
title_tail = fn
man_title = paste0("Manhattan Plot for UK biobank cohort for ",title_tail)
qq_title = paste0("QQ Plot for UK biobank cohort for ",title_tail)
file_name_out=paste0(fn,graphical_output) #ggsave for manhattan plot so will change format according to string - not true for QQ plot this is hard coded png below

library(tidyverse)
library(data.table)
library(ggrepel)
library(qqman)
library(RSQLite)
options(bitmapType='cairo')

library(data.table)
library(dplyr)
library(stringr)
library(tidyr)


ldf3 <- readRDS(paste0(IncludeFigHere, "/ldf3.RDS") )
QQ <- list()

for (j in levels(as.factor(ldf3$Variable))) { 
  
  print(paste0(".................... ", j))
  ldf4 <- ldf3 %>%
    filter(Variable == j) %>%
    select(-c(Cont, Variable))
  
  QQ[[j]] <- qq(ldf4$p)
  
  #print(qq(ldf4$LOG10P, main = paste0("QQplot - ", j)))
  
  
  png(paste0(IncludeFigHere, "/QQplot_",j ,".png"))
  print(qq(ldf4$p))
  dev.off()
  
}

library(patchwork)

#save(QQ, file = paste0(IncludeFigHere, "/QQ.RData"))
#
#png(paste0(IncludeFigHere, "/QQplotAll.png"))
#(QQ[[1]] + QQ[[2]]) / (QQ[[3]]  + QQ[[4]])
#dev.off()

print("QQPLOT")







