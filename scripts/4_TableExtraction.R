#Analysis run to generate tables from (1) GWAS (2) Rest 
#IncludeFigHere <- "/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930"
IncludeFigHere <- "/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/Final_20231217_1813" #Final_20231216_1613

library(xlsx)
library(data.table)
library(tidyverse)
select <- dplyr::select

#Supp Table
Supp1 <- fread(paste0(IncludeFigHere, "/Appendix1.txt"))
Supp2A <- fread(paste0(IncludeFigHere, "/BiomarkersPerParticipantInfo.txt"))
Supp2B <- fread(paste0(IncludeFigHere, "/TotalNumOfReadings.txt"))
Supp2C <- fread(paste0(IncludeFigHere, "/ParticipantsWithBiomarkers.txt"))
Supp3 <- fread(paste0("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/TreeWAS/FinalICD10CodesAndMultiCategoriesTOREPORT.tsv"))
Supp4 <- read.table(paste0(IncludeFigHere, "/LegendCategoriesAppendix.txt")) 

#change titles

Extra1 <- fread(paste0(IncludeFigHere, "/RRTableFINAL.csv"))


names(Extra1) <- names(Extra1) %>%
  str_replace_all(., "Healthy_Low", "PRP") %>%
  str_replace_all(., "Healthy_Cross", "NBP") %>%
  str_replace_all(., "Unhealthy_Cross", "PBN") %>%
  str_replace_all(., "Unhealthy_High", "NRN") %>%
  str_replace_all(., "dis", "Disease") 

Extra2 <- fread(paste0(IncludeFigHere, "/AssocDiab.csv"))


names(Extra2) <- names(Extra2) %>%
  str_replace_all(., "Healthy_Low", "PRP") %>%
  str_replace_all(., "Healthy_Cross", "NBP") %>%
  str_replace_all(., "Unhealthy_Cross", "PBN") %>%
  str_replace_all(., "Unhealthy_High", "NRN") %>%
  str_replace_all(., "dis", "Disease") 

#Variable and title

Supp5 <- fread(paste0("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930/ProximityWith8.tsv")) %>% #fread(paste0(IncludeFigHere, "/ProximityWith8.tsv")) if run GWAS!!
  mutate(Variable = str_replace(Variable, "Healthy", "PRP")) %>%
  mutate(Variable = str_replace(Variable, "Unhealthy", "NRN"))

#Type and Title


Extra3 <- fread(paste0("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930/GWASCatalogResults.csv")) %>% #fread(paste0(IncludeFigHere, "/GWASCatalogResults.csv")) if run GWAS!!
  mutate(Type = str_replace(Type, "Healthy", "PRP")) %>%
  mutate(Type = str_replace(Type, "Unhealthy", "NRN"))

#Variable and title - betas

#careful with associations of healthy and unhealthy - did well taking Varaiable of dataset? Same results?



Extra4 <- fread(paste0("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930/MergeAll8.csv")) %>% #fread(paste0(IncludeFigHere, "/MergeAll8.csv")) if run GWAS!!  
  select(-c(Common, Variable)) %>%
  inner_join(., Supp5 %>% select(Variable, ID, BETA:p) %>% rename(SNP = ID)) %>%
  distinct() %>% 
  mutate(Variable = str_replace(Variable, "Healthy", "PRP")) %>%
  mutate(Variable = str_replace(Variable, "Unhealthy", "NRN"))

#-- in ProximityAndEQTLFinal3.R

#Type and Title


Extra5 <- fread(paste0("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930/AllFUMAResults.csv")) %>% #fread(paste0(IncludeFigHere, "/AllFUMAResults.csv")) %>% #- in ProximityAndEQTLFinal3.R -UPLOAD manually
  mutate(Type = str_replace(Type, "Healthy", "PRP")) %>%
  mutate(Type = str_replace(Type, "Unhealthy", "NRN"))

#Variable  title and betas 

Extra6 <- fread(paste0(IncludeFigHere, "/Longevity_Table2.tsv")) %>% 
  mutate(Variable = str_replace(Variable, "Healthy", "PRP")) %>%
  mutate(Variable = str_replace(Variable, "Unhealthy", "NRN"))

Supp6 <- fread(paste0(IncludeFigHere, "/Longevity_Table.tsv")) %>% 
  select(-Variable) %>%
  inner_join(., Supp5 %>% select(Variable, ID, BETA:p) %>% rename(SNP = ID)) %>%
  distinct() %>% 
  mutate(Variable = str_replace(Variable, "Healthy", "PRP")) %>%
  mutate(Variable = str_replace(Variable, "Unhealthy", "NRN"))



wb <- createWorkbook()

sheet <- createSheet(wb, "Supp Table 1 - Read codes")
addDataFrame(Supp1, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Supp Table 2A - Biomarkers Per Participant")
addDataFrame(Supp2A, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Supp Table 2B - Total Num Of Readings")
addDataFrame(Supp2B, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Supp Table 2C - Participants With Biomarkers")
addDataFrame(Supp2C, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Supp Table 3 - ICD10 and Categories")
addDataFrame(Supp3, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Supp Table 4 - Legend")
addDataFrame(Supp4, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Supp Table 5 - Summary Stats Significant")
addDataFrame(Supp5, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Supp Table 6 - Longevity SNPs")
addDataFrame(Supp6, sheet, row.names = FALSE, col.names = TRUE)


saveWorkbook(wb, paste0(IncludeFigHere, "/SupplementaryTable.xlsx"))




#Appendix/Github

#Appendix1 <- fread(paste0(IncludeFigHere, "/Appendix1.txt"))
#Appendix2 <- read.table(paste0(IncludeFigHere, "/Appendix2A.txt")) - info on remaining participants
#Appendix3 <- read.table(paste0(IncludeFigHere, "/LegendCategoriesAppendix.txt")) 


wb <- createWorkbook()
sheet <- createSheet(wb, "Extra 1 - Relative Risk")
addDataFrame(Extra1, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Extra 2 - Diabetes Association - Relative Risk")
addDataFrame(Extra2, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Extra 3 - GWAS Catalog rsIDs")
addDataFrame(Extra3, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Extra 4 - Gene annotations")
addDataFrame(Extra4, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Extra 5 - FUMA results")
addDataFrame(Extra5, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Extra 5 - All databases longevity")
addDataFrame(Extra6, sheet, row.names = FALSE, col.names = TRUE)

saveWorkbook(wb, paste0(IncludeFigHere, "/Extra.xlsx"))


#Table
Tab1 <-  fread(paste0(IncludeFigHere, "/TableLifestyeFactors.csv"))[,-1] %>%
  mutate(across(everything(), ~gsub("Healthy", "Positive", ., fixed = TRUE))) %>%
  mutate(across(everything(), ~gsub("Unhealthy", "Negative", ., fixed = TRUE))) %>%
  mutate(across(everything(), ~gsub("unhealthy", "negative", ., fixed = TRUE))) %>%
  mutate(across(everything(), ~gsub("healthy", "positive", ., fixed = TRUE)))

Tab2 <-  fread(paste0(IncludeFigHere, "/AllPrevalence4CatsToo.csv"))[,-1] 

names(Tab2) <-  names(Tab2) %>% 
  str_replace(., "Healthy_Cross", "NBP" ) %>%
  str_replace(., "Healthy_Low", "PRP" ) %>%
  str_replace(., "Unhealthy_Cross", "PBN" ) %>%
  str_replace(., "Unhealthy_High", "NRN" )



wb <- createWorkbook()
sheet <- createSheet(wb, "Table 1 - Lifestyle")
addDataFrame(Tab1, sheet, row.names = FALSE, col.names = TRUE)

sheet <- createSheet(wb, "Table 2 - Prevalence")
addDataFrame(Tab2, sheet, row.names = FALSE, col.names = TRUE)


saveWorkbook(wb, paste0(IncludeFigHere, "/Table.xlsx"))

