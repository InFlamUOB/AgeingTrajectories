setwd("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity")


#Main analysis conducted here ( run sbatch (.sh attached) and able to reproduce everything in a single go), with extras performed separately (GWAS results compilation found in 2_SumamryStats.R,
#GWAS results enrichment with Proxiimty, eQTL analysis and enrichment found in 3_GWASEnrichment.R and final tables compiled uniformly in 
#4_TableExtraction). Extra functions for missing value analysis found in 1_MainAnalysisFunction.R


#Unhealthy = Negative
#Healthy = Positive

#Figures
# - Biorender
#- Fig2: PlotLines.png y PlotPPT.pdf
# - Fig3: RRnet_log2_sig_log2RR1Size x 4 categories
# - Fig 4:  MixFinal1.pdf + MixFinal2.pdf (FUMA results #-- in ProximityAndEQTLFinal3.R)
# - Fig 5: AgeingDatabase8.pdf+ AgeingDatabase_Legend.pdf

#SuppFigures

#Fig1: MultimorbidityHistogram.pdf
#Fig2A: Manhattan_xxx.png x 4
#Fig2B: QQplot_xx.png x 4
#Fig3A: VennSNPs.pdf (Venn Diagram rsIDs) - in ProximityAndEQTLFinal3.R
#Fig3B: Venn_All.pdf (Venn Diagram gene annotation) - in ProximityAndEQTLFinal3.R


# More analysis found in ManuscriptDecRun.R (rules, multimrbidity, disgenet, TreeWAS, Knowledge Graph, jaccard, correlation, clustering - here paper analysis only).
#2187 line pointing at GWAS results from previous run (sig_res Supp4) - change to run in a single go!

################################################### ################################################### ################################################### 
################################################### Directory and packages ################################################### ################################################### 
################################################### ################################################### ################################################### 

mainDir <- "/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov"
NameRun <- "Final"
subDir <- paste0(sub('\\..*', '', NameRun), format(Sys.time(), '_%Y%m%d_%H%M'))
dir.create(file.path(mainDir, subDir))
IncludeFigHere <- file.path(mainDir, subDir)


print(IncludeFigHere)


library(data.table)
library(dplyr)
library(stringr)
library(tidyr)
library(lubridate)
library(ggplot2)
library(readr)
library(ggpubr)
library(dunn.test) ####
library(broom)
library(purrr)
library(ggsignif)
library(tibble)
library(RColorBrewer)

library(patchwork)

############

#library(ggnet)
library(GGally)
library(igraph)
#install.packages("intergraph")
library(intergraph)  #required to get ggnet2 working 
library(ggrepel)
library(pheatmap)
library(forcats)

#install.packages("GGally")
#library(GGally) #substitute ggnet

######

library(DescTools) ####
library(viridis)
library(forcats)
library(gridExtra)
library(correlation) ####
library(tibble)

#https://mran.microsoft.com/snapshot/2017-12-11/web/packages/ggCompNet/vignettes/examples-from-paper.html
#https://journal.r-project.org/archive/2017/RJ-2017-023/RJ-2017-023.pdf
#https://journal.r-project.org/archive/2017/RJ-2017-023/RJ-2017-023.pdf

library(finalfit)
library(gt)

library(mice)
library(tidylog)
library(scales)
library(readxl)

'%ni%' <- Negate('%in%') #check good - 4 patients out 


################################################### ################################################### ################################################### 
################################################### Filtering GP data and blood tests ################################################### ################################################### 
################################################### ################################################### ################################################### 


#GP_Clinical <- fread( '/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/data_files_ukb/hes/update_20210415/gp_clinical.txt')
#Then #1.GpClinicalWDescriptionFileCreationCopy18Jan.Rmd - aggregate gpclinical info from UKBiobank


WithDescription <- fread("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/GpClinicalWDescription.csv")

WithDescriptionFilter <- WithDescription %>%
  filter(str_detect(term_description, "albumin|alk. phos.|creatinine|glucose|CRP|C-reactive|C-reactive protein|lymphocyte|mean cell volume|MCV|Red blood|alkaline phosphatase|White blood|white count")) %>%
  filter(!str_detect(term_description, "WAVESENSE JAZZ|diabetes|test|Review|diseases|monitoring|meter|Microalbumin|asting|Urine|minute|RBC - Red blood cell count|Red blood cell (RBC) count|Red blood cell count|(RBC)")) %>%
  mutate(FinalValue = ifelse(value1 == "OPR003" | value1 == "OPR001" | value1 == "HPS002" | value1 == "OPR004" | value1 == "OPR002",value2, value1)) %>% # get common column with final values
  filter(FinalValue != "") %>% #lots of values with nothing!
  #filter(FinalValue != "0.000") %>% #lots of values with nothing! -
  #filter(FinalValue != "0.00") %>% #lots of values with nothing!-
  #filter(FinalValue != "0") %>% #lots of values with nothing!-
  filter(FinalValue != "^") %>% #lots of values with nothing!
  filter(!is.na(event_dt)) %>% #lots of values with nothing!
  group_by(term_description) %>%
  mutate(Count= n_distinct(eid)) %>%
  ungroup() %>%
  mutate(Age = time_length(ymd("2022-07-28") - dmy(event_dt) , unit = "years")) %>%
  filter(Age > 0) %>%
  filter(Age < 100) #delete strange dates
#5 minutes

SelectedTerms <- WithDescriptionFilter %>%
  select(Count, term_description) %>%
  unique() #more precise filtering direction

Filter2CountWithDescription <- WithDescriptionFilter%>%
  filter(Count > 10000 | term_description %in% c("Plasma random glucose level", "Total alkaline phosphatase", "Serum random glucose level", "Plasma creatinine level", "Serum albumin normal", "Serum creatinine normal", "Plasma albumin level", "Serum alk. phos. normal", "Red blood cell size", "Plasma alkaline phosphatase level", "Serum creatinine raised", "Serum creatinine low", "Blood glucose normal", "Leucocytosis -high white count", "Serum alk. phos. raised", "Leucopenia - low white count", "Serum albumin low", "MCV - normal", "MCV - low", "MCV - raised", "Serum creatinine abnormal", "Total white blood count"))

Filter2CountWithDescription$term_description <- as.factor(Filter2CountWithDescription$term_description)
Filter2CountWithDescription$event_dt <- dmy(Filter2CountWithDescription$event_dt)
Filter2CountWithDescription$FinalValue <- as.numeric(Filter2CountWithDescription$FinalValue) 
Filter2CountWithDescription$eid <- as.factor(Filter2CountWithDescription$eid)


#PhenoAge paper

# Albumin, mg/dL	42.7 ± 4.2
# Creatinine, µmol/L	75.2 ± 21.9
# Glucose, mmol/L	5.4 ± 1.5
# C-reactive protein, g/L	0.4 ± 1.0
# Lymphocyte percent, %	30.0 ± 9.6
# Mean cell volume, fL	90.0 ± 6.0
# Red cell distribution width, %	12.7 ± 1.1
# Alkaline phosphatase, U/L	69.4 ± 28.0
# White blood cell count, 1000cell/uL	6.8 ± 2.4


#Figures

Hist1 <- list()

for (i in levels(as.factor(Filter2CountWithDescription$Value))){
  
  data <- filter(Filter2CountWithDescription,Value == i )
  data$Value <- droplevels(as.factor(data$Value))
  
  
  Hist1[[i]] <-  ggplot(data, aes( x = event_dt, y = FinalValue, fill =  term_description, colour = term_description)) + 
    geom_point(alpha= 0.2)+
    theme_bw() +
    ggtitle(paste0(i)) + 
    labs(y = paste0(i)) +
    #theme(legend.position = "none") +
    theme(plot.title = element_text(size = 7),
          axis.title=element_text(size=7))
  
  
}


png(paste0(IncludeFigHere, "/MeanClusterBox2.png"),type = "cairo", width = 2250, height = 1580,
    units = "px",  res = NA)

wrap_plots(Hist1, ncol = 6,nrow = 8)

dev.off()

#Help decide units


Hist2 <- list()

for (i in levels(as.factor(Filter2CountWithDescription$term_description))){
  
  data <- filter(Filter2CountWithDescription,term_description == i )
  data$term_description <- droplevels(as.factor(data$term_description)) 
  data$value3 <-  data$value3 %>%
    str_replace_all(c("µmol/L" = "Greekmol/L"))
  
  
  Hist2[[i]] <-  ggplot(data, aes( x = event_dt, y = FinalValue, fill =  value3, colour = value3)) + 
    geom_point(alpha= 0.4)+
    theme_bw() +
    ggtitle(paste0(i)) + 
    labs(y = paste0(i)) +
    #theme(legend.position = "none") +
    theme(plot.title = element_text(size = 7),
          axis.title=element_text(size=7))
  
  
}

#options(bitmapType = 'cairo')
png(paste0(IncludeFigHere, "/Units4.png"),type = "cairo", width = 2050, height = 1580,
    units = "px",  res = NA)

wrap_plots(Hist1, ncol = 6,nrow = 8)

dev.off()

Units <- Filter2CountWithDescription %>% 
  group_by(term_description, value3,Value ) %>%
  summarise( mean = mean(FinalValue), sd = sd(FinalValue), n()) %>%
  mutate(AllTogether = paste0(term_description,"_", value3, "_", Value ))

write.table(Units, paste0(IncludeFigHere, "/Units4.txt"))
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5940111/ table 1 plasma glucose out - it is serum 

TakeOut <- filter(Units, is.na(sd))


CuratedData2 <- Filter2CountWithDescription %>%
  filter(!c(term_description == "Percentage lymphocytes" & value3 == "MEA000" & Value == "42b1.")) %>%
  filter(!c(term_description == "Plasma random glucose level" & value3 == "%" & Value == "44g0.")) %>%
  filter(!c(term_description == "Plasma random glucose level" & value3 == "bm" & Value == "44g0.")) %>%
  filter(!c(term_description == "Plasma random glucose level" & value3 == "mg/mmol" & Value == "44g0.")) %>%
  filter(!c(term_description == "Plasma random glucose level" & value3 == "ml/min" & Value == "44g0.")) %>%
  filter(!c(term_description == "Plasma random glucose level" & value3 == "nmol/l" & Value == "44g0.")) %>%
  filter(!c(term_description == "Percentage lymphocytes" & value3 == "MEA037" & Value == "42b1.")) %>%
  filter(!c(term_description == "Red blood cell distribution width" & value3 == "MEA000" & Value == "42Z7.")) %>%
  #filter(!c(term_description == "Serum albumin" & value3 == "" & Value == "44M4.")) %>%
  #filter(!c(term_description == "Serum creatinine" & value3 == "" & Value == "XE2q5")) %>%
  #filter(!c(term_description == "Serum creatinine level" & value3 == "" & Value == "XE2q5")) %>%
  filter(!c(term_description == "White blood count" & value3 == "MEA153" & Value == "42H..")) %>%
  filter(!c(term_description == "White blood count" & value3 == "/ul" & Value == "42H..")) %>%
  mutate(AllTogether = paste0(term_description,"_", value3, "_", Value )) %>%
  filter(AllTogether %ni% TakeOut$AllTogether)


ff <- CuratedData2 %>%
  mutate(CommonLabels = case_when(str_detect(term_description, regex("albumin", ignore_case = TRUE))  ~ "Albumin", 
                                  str_detect(term_description, regex("creatinine", ignore_case = TRUE))  ~ "Creatinine", 
                                  str_detect(term_description, regex("glucose", ignore_case = TRUE))  ~ "Glucose", 
                                  str_detect(term_description, regex("white", ignore_case = TRUE))  ~ "White blood count",  
                                  str_detect(term_description, regex("reactive|protein|CRP", ignore_case = TRUE))  ~ "Plasma C-reactive protein level",
                                  str_detect(term_description, regex("lymphocyte", ignore_case = TRUE))  ~ "Percentage lymphocytes",  
                                  str_detect(term_description, regex("mean cell volume|volume|MCV", ignore_case = TRUE))  ~ "Mean corpuscular volume (MCV)",  
                                  str_detect(term_description, regex("red cell|width|size", ignore_case = TRUE))  ~ "Red blood cell dist width",  
                                  str_detect(term_description, regex("alk|alkaline phosphatase", ignore_case = TRUE))  ~ "Alkaline Phosphatase",  
                                  TRUE  ~ "ELSE")
  )

#many with repeated measurements per date! Leave just 1!
#https://www.healthline.com/health/rheumatoid-arthritis-crp-levels


Compare <- ff %>%
  #filter(Value %in% unique(Filter2CountWithDescription$Value)) %>%
  select(data_provider, Read, Value, term_description, CommonLabels) %>%
  unique() %>%
  arrange((term_description))


write.table(paste0("Data provider: ", Compare$data_provider, ", Read code: ", Compare$Read, ", Value: ", Compare$Value, ", Term description: ", Compare$term_description, ", Label: ", Compare$CommonLabels ), paste0(IncludeFigHere, "/Appendix1.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)



#Outliers and final data


ff2 <- ff %>%
  select(eid, event_dt ,FinalValue, Age, CommonLabels) %>%
  unique() %>%
  group_by(CommonLabels ) %>%
  mutate(ThresholdUpper =  ifelse(CommonLabels %in% c("C-reactive protein"), 1000, ifelse(CommonLabels %in% c("Lymphocyte percentage"), 100, quantile(FinalValue, 0.995)[[1]])), 
         ThresholdLower =   quantile(FinalValue, 0.005)[[1]],  
         Delete = ifelse(FinalValue < ThresholdUpper & FinalValue > ThresholdLower, "In", "Out")
  ) %>%
  ungroup() 

Threshold <- ff2 %>%
  select(CommonLabels, ThresholdUpper, ThresholdLower) %>%
  unique()

write.table(Threshold, paste0(IncludeFigHere, "/Threshold.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)


HistA <- list()

for (i in levels(as.factor(ff2$CommonLabels))){
  
  data <- filter(ff2,CommonLabels == i )
  data$CommonLabels <- droplevels(as.factor(data$CommonLabels)) 
  
  ValueUpper <- filter(Threshold, CommonLabels == i)$ThresholdUpper
  ValueLower <- filter(Threshold, CommonLabels == i)$ThresholdLower
  print(paste0(ValueUpper,",", ValueLower))
  
  HistA[[i]] <-  ggplot(data, aes( x = event_dt, y = FinalValue, fill =  Delete, colour = Delete)) + 
    geom_point(alpha= 0.4)+
    theme_bw() +
    ggtitle(paste0(i)) + 
    labs(y = paste0(i)) +
    # geom_hline( aes(yintercept = ValueUpper), linetype='solid', size=1) + #not working?
    #geom_hline( aes(yintercept = ValueLower), linetype='solid', size=1) +
    #theme(legend.position = "none") +
    theme(plot.title = element_text(size = 7),
          axis.title=element_text(size=7))
  
  
}


ff3 <- ff2 %>%
  mutate(FinalValue = ifelse(CommonLabels == "C-reactive protein", FinalValue/10, FinalValue)) 


HistA <- list()

for (i in levels(as.factor(ff3$CommonLabels))){
  
  data <- filter(ff3,CommonLabels == i )
  data$CommonLabels <- droplevels(as.factor(data$CommonLabels)) 
  
  ValueUpper <- filter(Threshold, CommonLabels == i)$ThresholdUpper
  ValueLower <- filter(Threshold, CommonLabels == i)$ThresholdLower
  print(paste0(ValueUpper,",", ValueLower))
  
  HistA[[i]] <-  ggplot(data, aes( x = event_dt, y = FinalValue, fill =  Delete, colour = Delete)) + 
    geom_point(alpha= 0.4)+
    theme_bw() +
    ggtitle(paste0(i)) + 
    labs(y = paste0(i)) +
    # geom_hline( aes(yintercept = ValueUpper), linetype='solid', size=1) + #not working?
    #geom_hline( aes(yintercept = ValueLower), linetype='solid', size=1) +
    #theme(legend.position = "none") +
    theme(plot.title = element_text(size = 7),
          axis.title=element_text(size=7))
  
  
}


png(paste0(IncludeFigHere, "/Thresh2WithCRP.png"),type = "cairo", width = 1650, height = 1080,
    units = "px",  res = NA)

wrap_plots(HistA, ncol = 3,nrow = 3)

dev.off()



CuratedData3 <- ff3 %>%
  filter(Delete == "In")

fwrite(CuratedData3, paste0(IncludeFigHere, "/CuratedDataJuly.csv"))  


################################################### ################################################### ################################################### 
################################################### GP data summary, supplementary ################################################### ################################################### 
################################################### ################################################### ################################################### 

CuratedData <- fread(paste0(IncludeFigHere, "/CuratedDataJuly.csv")) 


write.table(paste0("All participants: ", length(unique(WithDescription$eid)), ", Participants without read codes: ", length(unique(CuratedData$eid))),  paste0(IncludeFigHere, "/Appendix2A.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)

length(unique(WithDescription$eid)) #230,030
length(unique(CuratedData$eid)) #206,218


#Table and Plots 

Table <- CuratedData %>%
  group_by(eid, CommonLabels) %>%
  summarise(n = n()) #%>%
#ungroup() %>%
#group_by(CommonLabels) %>%
#summarise(mean = mean(n), sd = sd(n), Nums = n())


TableA <- CuratedData %>%
  group_by(CommonLabels, eid) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  group_by(CommonLabels) %>%
  summarise( max = max(n), mean = round(mean(n),2), sd = round(sd(n),2),  NumberPeopleWithAtLeastOneValue = n()) #at least one of those values is the n() here.

#write.csv(TableA, "BiomarkersPerParticipantInfo.csv")
fwrite(TableA, paste0(IncludeFigHere, "/BiomarkersPerParticipantInfo.txt"))


Table2 <- CuratedData %>%
  group_by(CommonLabels) %>%
  summarise(TotalNumberOfReadings = n())

fwrite(Table2, paste0(IncludeFigHere, "/TotalNumOfReadings.txt"))

Table3 <- CuratedData %>%
  group_by(eid) %>%
  summarise(n = n())

Table4 <- Table3 %>%
  select(n) %>%
  ungroup() %>%
  group_by(n) %>%
  count()%>%
  ungroup() %>%
  mutate(Groups = case_when(n >= 1000 ~ "1000-1532", 
                            n < 1000 & n >= 500  ~ "500-1000",
                            n < 500 & n >= 300 ~ "300-500",
                            n < 300 & n >= 200 ~ "200-300",
                            n < 200 & n >= 100 ~ "100-200",
                            n < 100 & n >= 80 ~ "80-100",
                            n < 80 & n >= 60 ~ "60-80",
                            n < 60 & n >= 50 ~ "50-60",
                            n < 50 & n >= 40 ~ "40-50",
                            n < 40 & n >= 30 ~ "30-40",
                            n < 30 & n >= 20 ~ "20-30",
                            n < 20 & n >= 15 ~ "15-20",
                            n < 15 & n >= 10 ~ "10-15",
                            n < 10 & n >= 5 ~ "5-10",
                            n < 5 & n >= 0 ~ "0-5",
                            TRUE ~ "Nothing"
  )) %>%
  group_by(Groups) %>%
  mutate(Sum = sum(nn))

Table5 <- Table4 %>%
  select(Groups, Sum) %>%
  unique()

fwrite(Table5, paste0(IncludeFigHere, "/ParticipantsWithBiomarkers.txt"))



################################################### ################################################### ################################################### 
################################################### Add Biobank data to GP  ################################################### ################################################### 
################################################### ################################################### ################################################### 


IDS <- read_csv("/rds/projects/g/gkoutosg-variant-prediction/Laura2/oldlxb732/BiobankPhenoAge/DeleteIDsBiobankFeb2.csv", col_names = FALSE) # file they sent us with the IDs we had to remove. 
SPANKPhenoAgeHESS1 <- read.csv('/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/SPANKComorbidities.csv') #UKBiobank data biomarker info to create phenoage!


SPANKPhenoAgeHESS <- SPANKPhenoAgeHESS1 %>%
  mutate(phenoage1.age_recruitment = phenoage.age_recruitment + time_length(ymd(phenoage1.date_recruitment_center) - ymd(phenoage.date_recruitment_center), unit = "years")) %>% 
  mutate(phenoage2.age_recruitment = phenoage.age_recruitment + time_length(ymd(phenoage2.date_recruitment_center) - ymd(phenoage.date_recruitment_center), unit = "years"))

# have to get bmi from clinicla gp codes too!

SelectHESSAge <- SPANKPhenoAgeHESS %>%
  select(-c(contains("comorbidity"), "n_admissions","time_until_readmission" , contains("death"), phenoage.sex_male, phenoage.smoking_status_coding, phenoage.ethnic_background_coding, phenoage.weight, phenoage.genetic_sex_male:phenoage.imd_scotland, phenoage1.smoking_status_coding, phenoage1.ethnic_background_coding,phenoage1.weight, phenoage2.height:phenoage2.ethnic_background_coding, phenoage2.weight,age_disease))

SelectHESSAgeVars1 <- SelectHESSAge   %>%   
  select(c(eid,phenoage.age_recruitment, phenoage.date_recruitment_center:phenoage.glucose_mmolL)) %>%
  pivot_longer(., cols = phenoage.white_blood_cell_1000uL:phenoage.glucose_mmolL,  names_to = "PhenoAgeValues", values_to= "Values") %>%
  rename(c("age" = "phenoage.age_recruitment",  "date" = "phenoage.date_recruitment_center",  "bmi" = "phenoage.bmi"))

SelectHESSAgeVars2 <- SelectHESSAge   %>%   
  select(c(eid, phenoage1.age_recruitment,phenoage1.date_recruitment_center:phenoage1.glucose_mmolL)) %>%
  pivot_longer(., cols = phenoage1.white_blood_cell_1000uL:phenoage1.glucose_mmolL,  names_to = "PhenoAgeValues", values_to= "Values") %>%
  rename(c("age" = "phenoage1.age_recruitment",  "date" = "phenoage1.date_recruitment_center",  "bmi" = "phenoage1.bmi"))

SelectHESSAgeVars3 <- SelectHESSAge   %>%   
  select(c(eid, phenoage2.age_recruitment,phenoage2.date_recruitment_center:phenoage2.lymphocyte_perc)) %>% 
  pivot_longer(., cols = phenoage2.white_blood_cell_1000uL:phenoage2.lymphocyte_perc,  names_to = "PhenoAgeValues", values_to= "Values") %>%
  rename(c("age" = "phenoage2.age_recruitment",  "date" = "phenoage2.date_recruitment_center",  "bmi" = "phenoage2.bmi"))


PhenoAgeValues <- rbind(SelectHESSAgeVars1, SelectHESSAgeVars2)

PhenoAgeValues <- rbind(PhenoAgeValues, SelectHESSAgeVars3) %>%
  drop_na(Values) 

length(unique(PhenoAgeValues$eid)) #490942

FirstVisit <- SelectHESSAge   %>%   
  select(c(eid,phenoage1.date_recruitment_center)) %>%
  filter(phenoage1.date_recruitment_center != "")

SecondVisit <- SelectHESSAge   %>%   
  select(c(eid,phenoage2.date_recruitment_center)) %>%
  filter(phenoage2.date_recruitment_center != "")

Merge <- unique(c(FirstVisit$eid, SecondVisit$eid))


see <- PhenoAgeValues %>%
  filter(grepl("phenoage1.", PhenoAgeValues) | grepl("phenoage2.", PhenoAgeValues))

length(unique(see$eid)) #number of patients with at least a second visit 

PhenoAgeValues$PhenoAgeValues <-  PhenoAgeValues$PhenoAgeValues %>% 
  str_replace_all(c( "e1." = "e.", "e2." = "e.")) %>%
  str_replace_all(c( "phenoage.albumin_gL" = "Albumin",
                     "phenoage.alkaline_phosphatase_UL" = "Alkaline Phosphatase",
                     "phenoage.creat_umolL" = "Creatinine" , 
                     "phenoage.crp_mgL" = "Plasma C-reactive protein level", 
                     "phenoage.glucose_mmolL" = "Glucose", 
                     "phenoage.lymphocyte_perc" = "Percentage lymphocytes", 
                     "phenoage.mean_corpuscular_volume_fL" = "Mean corpuscular volume (MCV)", 
                     "phenoage.red_blood_cell_RDW_perc" = "Red blood cell dist width",
                     "phenoage.white_blood_cell_1000uL" = "White blood count"))



CuratedData <- fread(paste0(IncludeFigHere, "/CuratedDataJuly.csv")) 

CuratedDataSelected <- CuratedData %>%
  select(eid,event_dt,CommonLabels, FinalValue  ) %>%
  mutate(Type = "GP")

# do an age in column here


Age <- SPANKPhenoAgeHESS1 %>%
  rename(c("age" = "phenoage.age_recruitment",  "date" = "phenoage.date_recruitment_center")) %>%
  #select(date, age, eid) %>%
  select(date, age, eid) %>%
  group_by(eid) %>%
  slice(1) %>%
  ungroup()


CuratedDataSelected2 <- left_join(CuratedDataSelected, Age) %>% 
  dplyr::mutate(AgeGP2 = case_when(time_length(ymd(event_dt) - ymd(date), unit = "years") >= 0 ~ age + time_length(ymd(event_dt) - ymd(date), unit = "years"), TRUE ~ age - time_length(ymd(date) - ymd(event_dt), unit = "years") )) %>%
  select(-c(date, age)) %>%
  rename(c("date" = "event_dt", "PhenoAgeValues" = "CommonLabels", "Values" = "FinalValue", "age" = "AgeGP2")) 

#takes a bit

CuratedDataSelected2$age <- as.numeric(CuratedDataSelected2$age)
CuratedDataSelected2$date <- as.Date(CuratedDataSelected2$date)

PhenoAgeValues2 <- PhenoAgeValues %>%
  mutate(Type = "PhenoAge") %>%
  select(-bmi) 

PhenoAgeValues2$date <- as.Date(PhenoAgeValues2$date)

MergeData <- rbind(CuratedDataSelected2, PhenoAgeValues2)


fwrite(MergeData, paste0(IncludeFigHere, "/MergeData.csv"))  #############################################################################
#############################################################################
#############################################################################


MergeData <- fread(paste0(IncludeFigHere, "/MergeData.csv"))   #############################################################################
#############################################################################
#############################################################################

length(unique(filter(MergeData, Type == "PhenoAge")$eid)) #490942
length(unique(filter(MergeData, Type == "GP")$eid)) #206218

MergeDataWide <- MergeData %>% 
  pivot_wider(names_from = PhenoAgeValues,
              values_from = Values, values_fn = mean)

fwrite(MergeDataWide, paste0(IncludeFigHere, "/MergeDataWide.csv"))  

################################################### ################################################### ################################################### 
################################################### Extra missing data analysis ################################################### ################################################### 
################################################### ################################################### ################################################### 


MergeDataWide <- fread(paste0(IncludeFigHere, "/MergeDataWide.csv"))  



pdf(paste0(IncludeFigHere, "/Missing.pdf"), 50, 55) 
print(MergeDataWide %>%
        missing_pattern(explanatory = names(MergeDataWide))
)
dev.off()


source("1_MainAnalysisFunction.R")


pdf(paste0(IncludeFigHere, "/PlotMissingSwitch3.pdf"), 10,15)

PlotMissingMagda(MergeDataWide %>% 
                   filter(Type == "GP"))

dev.off()

#just 1 single one that has a complete case!

CompleteCases <- MergeDataWide %>%
  filter(Type == "GP") %>%
  drop_na()

length(unique(CompleteCases$eid))

################################################### ################################################### ################################################### 
################################################### Misisng value imputation window ################################################### ################################################### 
################################################### ################################################### ################################################### 

MergeDataWide$date <- ymd(MergeDataWide$date)

MergeDataWide2 <- MergeDataWide %>%
  filter(date < as.Date('2030-05-02')  & date > as.Date('1915-05-02')) %>%
  arrange((date)) 

MergeDataWideFill <- MergeDataWide2 %>%
  group_by(eid) %>%  #change to group_by(eid, Type) to see how much the analysis changes!
  fill('Mean corpuscular volume (MCV)':'Plasma C-reactive protein level', .direction = "updown") %>%
  ungroup()

#MissingwindowPlot ()


Window1  <- MergeDataWide2 %>%
  filter(eid %in% c("1000048", "1000055", "1000072", "1000080", "1000099")) %>%
  select(eid, Type, `Mean corpuscular volume (MCV)`:`Plasma C-reactive protein level`) 
  
Window1$eid <- factor(Window1$eid , levels = c("1000048", "1000055", "1000072", "1000080", "1000099"), 
                             labels = c("Participant 1", "Participant 2", "Participant 3", "Participant 4", "Participant 5"))

Window2  <- MergeDataWideFill %>%
  filter(eid %in% c("1000048", "1000055", "1000072", "1000080", "1000099")) %>%
  select(eid, Type, `Mean corpuscular volume (MCV)`:`Plasma C-reactive protein level`)

Window2$eid <- factor(Window2$eid , levels = c("1000048", "1000055", "1000072", "1000080", "1000099"), 
                      labels = c("Participant 1", "Participant 2", "Participant 3", "Participant 4", "Participant 5"))



fwrite(Window1, paste0(IncludeFigHere, "/Window1.csv"))
fwrite(Window2, paste0(IncludeFigHere, "/Window2.csv"))
#mixing GP and PhenoAge!


################################################### ################################################### ################################################### 
################################################### Calculate PhenoAge  ################################################### ################################################### 
################################################### ################################################### ################################################### 


# there are missing values if there has never been a record of a specific marker 
PhenoAgeMergeDataWideFill <- MergeDataWideFill %>%
  dplyr::mutate(funsPhen = -19.907 - 0.0336*Albumin + 0.0095*Creatinine + 0.1953*Glucose + 0.0954*log(`Plasma C-reactive protein level`) - 0.0120*`Percentage lymphocytes` + 0.0268*`Mean corpuscular volume (MCV)` + 0.3306*`Red blood cell dist width` +  0.00188*`Alkaline Phosphatase` + 0.0554*`White blood count` + 0.0804*age) %>%
  dplyr::mutate(PhenoAge = 141.50 + (log((-0.00553)*(((-1.51714)*exp(funsPhen))/0.0076927)))/(0.09165)) %>%
  drop_na() %>% #length(unique(PhenoAgeMergeDataWideFill$eid))
  filter(`Plasma C-reactive protein level` != 0) %>% # added this log of 0 = Inf
  group_by(eid) %>%
  mutate(Count = n_distinct(date)) %>%
  ungroup() %>%
  filter(Count > 1) # take away patients with just one reading (phenoage? Makes sense GP info only had these 200,000 patients!)


#see <- filter(PhenoAgeMergeDataWideFill, eid %in% c(filter(SlopesOnly, Situation2 == "Unknown")$eid)) #too strange ups and downs!!


length(unique(PhenoAgeMergeDataWideFill$eid))#199228

length(levels(as.factor(PhenoAgeMergeDataWideFill$eid)))
#432395 with 1 and all 198459 with more than 1 reading! old
length(levels(as.factor(PhenoAgeMergeDataWideFill$date)))
#9147 old


################################################### ################################################### ################################################### 
################################################### Calculate Slopes  ################################################### ################################################### 
################################################### ################################################### ################################################### 


FinalFiltered <- PhenoAgeMergeDataWideFill %>%
  # filter(PhenoAge < 200 & PhenoAge > -100) %>% #filtered by sudden changes filter later? - takes out Inf values so in 
  #filter(Count == 60) %>%
  #filter(eid %in% c("5078410","3880505","2938678" )) %>%
  mutate(size = ifelse(Type == "GP", 2, 3))

fwrite(FinalFiltered, paste0(IncludeFigHere, "/FinalFiltered.csv"))        #############################################################################


TableInfo <- FinalFiltered %>%
  select(eid, date, age, Count) %>%
  group_by(eid) %>%
  mutate(Length = max(age) - min(age)) %>%
  ungroup() %>%
  select(eid, Count, Length) %>%
  unique()

write.table(TableInfo,paste0(IncludeFigHere, "/TableInfoLengthTimeANdCountFinalFiltered.txt"))

#s <- fread(paste0(IncludeFigHere, "/TableInfoLengthTimeANdCountFinalFiltered.txt"))

FinalFiltered$eid <- as.factor(FinalFiltered$eid )
#Unhealthy <- length(unique(filter(FinalFiltered, age < PhenoAge)$eid))
#Unhealthy <- filter(FinalFiltered, age < PhenoAge)
#Healthy <- length(unique(filter(FinalFiltered, age > PhenoAge)$eid))

WithSlopes <- FinalFiltered %>%
  group_by(eid) %>%
  arrange((date)) %>% do(
    {
      r <- residuals(lm((PhenoAge) ~ age, data = .))
      data.frame(., r)
    }) %>%
  filter(abs(r) < 10) %>% do( #deleted those completely out of range. Maybe include a measure for most varaibel least variable? Residual info - changes with age? Number of samples. 
    {
      s <- (lm((PhenoAge) ~ age, data = .))
      a1 <- data.frame(., Intercept=s[["coefficients"]][1])
      a2 <- data.frame(a1, AgeCoeffs = s[["coefficients"]][2])
    }) 
#big



SlopesOnly <- WithSlopes %>%
  ungroup() %>%
  dplyr::mutate( Result = as.numeric(-Intercept/(AgeCoeffs -1))) %>%
  dplyr::mutate( Situation = case_when(  (Intercept > 0 &  AgeCoeffs >= 1) ~ "Unhealthy_High", 
                                         (Intercept > 0 &  AgeCoeffs < 1) ~ "Healthy_Cross",
                                         (Intercept < 0 &  AgeCoeffs > 1) ~ "Unhealthy_Cross",
                                         (Intercept < 0 &  AgeCoeffs <= 1 & AgeCoeffs > 0) ~ "Healthy_Low", 
                                         TRUE ~ "Unknown"        
  )) %>%
  dplyr::mutate( Situation2 = case_when( (Situation == "Unhealthy_Cross" & Result > 100) ~ "Healthy_Low", 
                                         (Situation == "Healthy_Cross" & Result > 100) ~ "Unhealthy_High", 
                                         TRUE ~ Situation))




fwrite(SlopesOnly, paste0(IncludeFigHere, "/SlopesAndClassifictaion3.csv"))        

################################################### ################################################### ################################################### 
################################################### Slopes Description  ################################################### ################################################### 
################################################### ################################################### ################################################### 


SlopesOnly <- fread(paste0(IncludeFigHere, "/SlopesAndClassifictaion3.csv")) 

SlopesOnly <- SlopesOnly %>% 
  mutate(Situation2 = case_when(Situation2 == "Healthy_Cross" ~ "NBP", 
                                Situation2 == "Healthy_Low" ~ "PRP", 
                                Situation2 == "Unhealthy_Cross" ~ "PBN", 
                                Situation2 == "Unhealthy_High" ~ "NRN"))
PlotLines <-  ggplot(filter(SlopesOnly, Situation != "Unknown")) +
  geom_abline(aes(
    intercept = Intercept, 
    slope = AgeCoeffs, 
    color=factor(Situation2))) + #factor(Situation)
  xlim(0, 100) +
  ylim(0, 100) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed", size = 2 ) +
  geom_abline(slope=1, intercept = 0, color = "yellow", size = 2, alpha = 0.4 ) +
  theme_bw() + facet_grid(~ Situation2, scales = "free") +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(size = 20),#15
        axis.text.x = element_text( size = 20 ),
        axis.text.y = element_text( size = 20 ),
        axis.title.x = element_text( size = 20 ),
        axis.title.y = element_text( size = 20 )) +
  scale_colour_manual(values = c("#F8766D","#00BA38", "#F564E3", "#00BFC4")) +
  xlab("Age") +
  ylab("PhenoAge")

png(paste0(IncludeFigHere, "/PlotLines2.png"), type = "cairo", width = 1050, height = 520,
    units = "px")

#width = 1250, height = 580,

print(PlotLines)

dev.off()



aa <- SlopesOnly %>% 
  select(eid, Situation2) %>% 
  unique() 

write.table(table(aa$Situation2), paste0(IncludeFigHere, "/TableNumbersSlopes.txt"))

see <- filter(SlopesOnly, eid %in% c(filter(SlopesOnly, Situation2 == "Unknown")$eid)) #12 participants with no fitted lines, weird data points. 


FinalFiltered <- fread(paste0(IncludeFigHere, "/FinalFiltered.csv")) #############################################################################
FinalFiltered$eid <- as.factor(FinalFiltered$eid )

Plots2 <- ggplot(filter(FinalFiltered, eid %in% c("3054375", "2038638","2461370","1263795","1163740")), aes(x = age, y = PhenoAge, group = eid, col = eid, size = size)) + 
  geom_line(size = 1, alpha=0.4) + 
  geom_point(alpha = 0.7) +
  theme_bw() +
  theme(axis.line = element_line(colour = "grey50", size = 1),
        panel.grid.minor = element_blank()) +
  #scale_x_date(date_breaks = "5 years", date_minor_breaks = "1 year") +
  #ggtitle("PhenoAge components") + 
  labs(y="PhenoAge", x = "Age") + 
  theme(axis.title.x = element_text( size = 15 ),
        axis.title.y = element_text( size = 15),
        axis.text.x = element_text( size = 15 ),
        axis.text.y = element_text( size = 15 )) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed", size = 2 ) +
  geom_abline(slope=1, intercept = 0, color = "yellow", size = 2, alpha = 0.4 ) 



kk <- FinalFiltered %>% filter(eid %in% c("3054375", "2038638","2461370","1263795","1163740")) %>% filter(Type == "PhenoAge")


#Plots3 <- ggplot(kk, aes(x = age, y = PhenoAge, group = eid, col = eid, size = size)) + 
#  geom_line(size = 1,alpha=0.4) + 
#  geom_point(alpha = 0.7) +
#  theme_bw() +
#  theme(axis.line = element_line(colour = "grey50", size = 1),
#        panel.grid.minor = element_blank()) +
#  labs(y="PhenoAge", x = "Age") + 
#  geom_abline(slope=1, intercept = 0, linetype = "dashed", size = 2 ) +
#  geom_abline(slope=1, intercept = 0, color = "yellow", size = 2, alpha = 0.4 ) 

Plots4 <- ggplot( filter(FinalFiltered, eid %in% c("3054375", "2038638","2461370")), aes(x = age, y = PhenoAge, group = eid, col = eid, size = size)) + #1263795, 1163740
  geom_line(size = 1) + 
  geom_point(alpha = 0.8) +
  theme_bw() +
  theme(axis.line = element_line(colour = "grey50", size = 1),
        panel.grid.minor = element_blank()) +
  labs(y="PhenoAge", x = "Age") + 
  #scale_color_discrete(name = "Patient Id") +
  #theme(legend.position = "none") + 
  scale_colour_manual(values = c("#00BA38", "#00BFC4", "#F564E3")) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed", size = 2 ) +
  geom_abline(slope=1, intercept = 0, color = "yellow", size = 2, alpha = 0.4 ) +
  #geom_abline(slope= Lines$AgeCoeffs[1], intercept = Lines$Intercept[1], color = "dark orange",alpha = 0.3 , size = 1.5) +00BFC4
  # geom_abline(slope= Lines$AgeCoeffs[2], intercept = Lines$Intercept[2], color = "#556b2f",alpha = 0.3, size = 1.5 ) +
  geom_abline(slope= unique(filter(SlopesOnly, eid == "3054375")$AgeCoeffs), intercept = unique(filter(SlopesOnly, eid == "3054375")$Intercept), color = "#F564E3",alpha = 0.3, size = 1.5 ) + #Lines
  geom_abline(slope= unique(filter(SlopesOnly, eid == "2038638")$AgeCoeffs), intercept = unique(filter(SlopesOnly, eid == "2038638")$Intercept),color = "#00BA38",alpha = 0.3 , size = 1.5) +
  geom_abline(slope= unique(filter(SlopesOnly, eid == "2461370")$AgeCoeffs), intercept = unique(filter(SlopesOnly, eid == "2461370")$Intercept), color = "#00BFC4",alpha = 0.3 , size = 1.5) + 
  theme(axis.title.x = element_text( size = 15 ),
        axis.title.y = element_text( size = 15),
        axis.text.x = element_text( size = 15 ),
        axis.text.y = element_text( size = 15 ))


#show_col(hue_pal()(6))


pdf(paste0(IncludeFigHere, "/PlotPPT2.pdf"), 14,5 )

print(Plots2 + Plots4) #Plots3

dev.off()


################################################### ################################################### ################################################### 
################################################### Add diseases  ################################################### ################################################### 
################################################### ################################################### ################################################### 

Slopes <- fread(paste0(IncludeFigHere, "/SlopesAndClassifictaion3.csv"))

#HESS <- fread("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/AllHESS.csv", na.strings = c("",NA)) #diagnosis!
#NEW HESS!! old?
IDS <- read_csv("/rds/projects/g/gkoutosg-variant-prediction/Laura2/oldlxb732/BiobankPhenoAge/DeleteIDsBiobankFeb2.csv", col_names = FALSE) # file they sent us with the IDs we had to remove. 

SelectSelf <- fread("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/CompareHESSandSelf/HESSCheckedAugust_SN.csv") %>% #lots of extra instances!
  filter(eid %ni% c(IDS$X1))


ToDate <- function(x){
  y <- ymd(as.character(x))
  return(y)
}



SlopeMerge <- Slopes %>% 
  group_by(eid) %>%
  mutate(MaximumAgeValueBlood = max(age)) %>% #death is eithe oldest reading age or death age
  ungroup() %>% #include max age value 
  rename(., TypePlot = Type) %>%
  rename(., Age = age) %>%
  #select(-c(Mean.corpuscular.volume..MCV.:funsPhen, date, r, Count)) %>% #rename Type for TypePlot - age, PhenoAge, size and Type why take this out?? - took it out for Gwen
  select(-c(Mean.corpuscular.volume..MCV.:funsPhen, date, PhenoAge, r, Count)) %>% #rename Type for TypePlot - age, PhenoAge, size and Type why take this out?? - took it out for Gwen
  select(-c(Age, size, TypePlot)) %>%
  unique() 

#all renames are new
# PhenoAge,size, Type  why if i take his it goes crazy 

SelectSelf2 <- SelectSelf 

Rules <- inner_join(SlopeMerge, SelectSelf2) %>% #before just SelectSelf
  group_by(eid) %>%
  dplyr::mutate(Multimorb = n_distinct(Disease)) %>%
  ungroup() %>%
  mutate(FactorDeath = ifelse(is.na(Death), "No", "Yes")) %>%
  group_by(eid) %>%
  mutate(MaximumAgeValueDisease = max(Age)) %>% #death is eithe oldest reading age or death age 
  mutate(MaxAgeMix = ifelse(FactorDeath == "No", max(MaximumAgeValueDisease,MaximumAgeValueBlood), Death)) %>%
  ungroup() %>%
  filter(Situation2 != "Unknown") %>% 
  select(-Type) 

Rules2 <- left_join(SlopeMerge, SelectSelf) %>%
  group_by(eid) %>%
  dplyr::mutate(Multimorb = n_distinct(Disease)) %>%
  ungroup() %>%
  mutate(FactorDeath = ifelse(is.na(Death), "No", "Yes")) %>%
  group_by(eid) %>%
  mutate(MaximumAgeValueDisease = max(Age)) %>% #death is eithe oldest reading age or death age
  mutate(MaxAgeMix = ifelse(FactorDeath == "No", max(MaximumAgeValueDisease,MaximumAgeValueBlood), Death)) %>%
  ungroup() %>%
  filter(Situation2 != "Unknown") %>% 
  select(-Type)

#fwrite(Rules2, "Rules2.csv")

fwrite(Rules2, paste0(IncludeFigHere, "/Rules2AllDataNoDiseasesIn.csv") )
#fwrite(Rules, "Rules.csv")

#length(unique(Rules$eid)) #143230
#table(Rules$Situation2)
table(filter(SlopeMerge, eid %in% unique(Rules$eid))$Situation2)
table(filter(SlopeMerge, eid %in% unique(Rules2$eid))$Situation2)

Out <- unique(setdiff(Rules2$eid, Rules$eid))
#HESS$eid <- as.character(HESS$eid)

#h <- filter(select(HESS, eid,Num_self_noncancer) , eid %in% as.character(Out))
#table(h$Num_self_noncancer) #55986

#what to do? - deleting 55000 participants with non chronic diseases by our standards. 22000 (none!) but other 25000?
#0     1     2     3     4     5     6     7     8     9    10    11    12    13    14    15 
#21665 17340  9569  4250  1786   732   345   148    69    31    17     5     7     4     1     1 

#Healthy_Cross     Healthy_Low Unhealthy_Cross  Unhealthy_High 
#249356           43259          237114           44289 

fwrite(Rules, paste0(IncludeFigHere, "/Rules.csv")) 
Rules <- fread(paste0(IncludeFigHere, "/Rules.csv"))
Rules3 <- fread(paste0(IncludeFigHere, "/Rules2AllDataNoDiseasesIn.csv"))
Rules2a <- fread("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930/Rules.csv")
Rules4a <- fread("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930/Rules2AllDataNoDiseasesIn.csv")
Slopes <- fread("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930/SlopesAndClassifictaion3.csv")
Slopes2 <- fread(paste0(IncludeFigHere, "/SlopesAndClassifictaion3.csv"))


MergePlotSlopDiseaseFilter <- Rules %>%
  mutate(DiseasePhenoAge = as.numeric(Age*AgeCoeffs + Intercept )) %>%
  #filter(eid %in% SituationEid$eid) %>%
  filter(eid %in% c("3054375", "2038638","2461370"))

MergePlotSlopDiseaseFilter$eid <- as.factor(MergePlotSlopDiseaseFilter$eid)


Plot5 <-  ggplot(MergePlotSlopDiseaseFilter, aes(x = Age, y = DiseasePhenoAge)) +
  geom_abline(aes(
    intercept = Intercept, 
    slope = AgeCoeffs, 
    color=factor(Situation2)), size = 1.5, alpha = 0.4) + #factor(Situation)
  geom_point(shape = 5, size = 2)  + 
  xlim(0, 100) +
  ylim(0, 100) +
  scale_colour_manual(values = c("#00BA38", "#F564E3","#00BFC4")) +
  geom_abline(slope=1, intercept = 0, linetype = "dashed", size = 2 ) +
  geom_abline(slope=1, intercept = 0, color = "yellow", size = 2, alpha = 0.4 ) +
  theme_bw() +
  theme(axis.line = element_line(colour = "grey50", size = 1),
        panel.grid.minor = element_blank()) +
  ylab("Proxy of Disease PhenoAge")
#+
#theme(legend.position = "none")



png(paste0(IncludeFigHere, "/PlotPPT2.png"), type = "cairo",width = 420, height = 320,
    units = "px",  res = NA)


print(Plot5)

dev.off()



################################################### ################################################### ################################################### 
################################################### Description Slopes and Extras  ################################################### ################################################### 
################################################### ################################################### ################################################### 




#library(gridExtra)
#library(tidyverse)

#aa <- fread("/rds/projects/g/gkoutosg-variant-prediction/ukbiobank/data_files_ukb/ukb46518.tab")

#CompareParents <- aa %>%
#  select(eid = "f.eid", FatherDeath = "f.1807.0.0", MotherDeath = "f.3526.0.0", Sleep = "f.1160.0.0", 
#         MET_Physical = "f.22033.0.0", PA = "f.884.0.0", Summed_Minutes = "f.22034.0.0",
#         MET_IPAQ = "f.22032.0.0", MET_Physical_Summed = "f.22040.0.0",
#         AgeFather = "f.2946.0.0", AgeMother = "f.1845.0.0") %>%
#  mutate(LasteRecordedAgeFather = fcoalesce(FatherDeath,AgeFather)) %>%
#  mutate(LasteRecordedAgeMother = fcoalesce(MotherDeath,AgeMother)) #%>%
##select(-c(FatherDeath,AgeFather,MotherDeath,AgeMother))
##combine? 
#
#AllBiomarkers <- fread("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/AllBiomarkersNames.csv")
##table(AllBiomarkers$Pregnant, useNA = "always")
#
#AllBiomarkers2 <- AllBiomarkers %>%
#  left_join(.,CompareParents )
#
#OutcomesSpank <- AllBiomarkers2 %>%
#  select(eid, names(CompareParents), genetic_sex_male, smoking_status_coding, body_fat_percentage,Townsend_deprivation_index_at_recruitment, Overall_health_rating, Had_major_operations,Ethnic_background, 'Body_mass_index_(BMI)', Basal_metabolic_rate) #ethnic_genetic_background ,Fluid_intelligence_score
#
#fwrite(OutcomesSpank, "/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/AllBiomarkersNames2.csv")
#

#fwrite(Rules2, paste0(IncludeFigHere, "/Rules2AllDataNoDiseasesIn.csv") )
AllBiomarkers <- fread("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/AllBiomarkersNames2.csv")
#check all instances and coding 
#sleep not interesting!!

#mothers death and mothers age the same 

Rules2 <- fread(paste0(IncludeFigHere, "/Rules2AllDataNoDiseasesIn.csv"))

FinalData <-  Rules2 %>%  #complet set, not just diseases
  #filter(Situation != "Unknown") %>% 
  select(-c(Disease,Age, Dates)) %>%
  unique() %>%
  select(eid, Situation2, FactorDeath, MaxAgeMix, Multimorb) %>%
  inner_join(AllBiomarkers) %>% #not all have this!
  mutate_each_(funs(factor(.)),c( "FactorDeath", "genetic_sex_male", "smoking_status_coding","Had_major_operations", "Ethnic_background")) %>%
  #select(-c(eid, genetic_sex_male, smoking_status_coding, FactorDeath, Had_major_operations, body_fat_percentage, Sleep, Basal_metabolic_rate, Percentage_Had_major_operations)) %>%
  mutate_if(is.numeric, round, digits = 2) %>%
  mutate_at(vars(starts_with("Overall")),funs(factor)) %>%
  mutate_at(vars(starts_with("PA")),funs(factor)) %>%
  rename("Categories" = "Situation2", 
         "Age" = "MaxAgeMix", 
         "Comorbidities" = "Multimorb", 
         # "Fathers Lifespan" = "FatherDeath", 
         # "Mothers Lifespan" = "MotherDeath",
         # "Physical MET" = "MET_Physical", 
         # "Physical PA" = "PA", 
         # "Townsend deprivation index" = "Townsend_deprivation_index_at_recruitment", 
         # "Health rating" = "Overall_health_rating", 
         "BMI" = "Body_mass_index_(BMI)", 
         # "`Males %`" = "genetic_sex_male", 
         "Mortality" = "FactorDeath")#, 
#"Never smoked" = "smoking_status_coding") 


explanatory = setdiff(names(FinalData), c("Categories", "eid",  "Ethnic_background","Basal_metabolic_rate", "Sleep", "Had_major_operations", "body_fat_percentage" ))
dependent = "Categories"

t <- FinalData %>%
  mutate(Overall_health_rating = recode_factor(Overall_health_rating, 
                                               "1" = "Excellent", 
                                               "2" = "Good", 
                                               "3" = "Fair", 
                                               "4" = "Poor", 
                                               "-1" = "Do not know", 
                                               "-3" = "Prefer not to answer")) %>%
  mutate(smoking_status_coding = recode_factor(smoking_status_coding, 
                                               "1" = "Previous", 
                                               "2" = "Current", 
                                               "0" = "Never", 
                                               "-3" = "Prefer not to answer")) %>%
  mutate(genetic_sex_male = recode_factor(genetic_sex_male, 
                                          "0" = "Female", 
                                          "1" = "Male" )) %>%
  mutate(PA = recode_factor(PA, 
                            "-1" = "Do not know", 
                            "-3" = "Prefer not to answer")) %>%
  mutate(MET_IPAQ = recode_factor(MET_IPAQ, 
                                  "0" = "low", 
                                  "1" = "moderate", 
                                  "2" = "high")) %>%
  summary_factorlist(dependent, explanatory, p = TRUE, na_include = TRUE,  p_cont_para = "aov", add_col_totals = TRUE)  %>%
  relocate("Lifestye factors" = "label", levels,
           "Healthy remaining healthy" = Healthy_Low, 
           "Unhealthy becoming healthy" = Healthy_Cross, 
           "Healthy becoming unhealthy" = Unhealthy_Cross,
           "Unhealthy remaining unhealthy" = Unhealthy_High, p) 
#check names?


write.csv(t,paste0(IncludeFigHere, "/TableLifestyeFactors.csv"))


################################################### ################################################### ################################################### 
################################################### Prevalence  ################################################### ################################################### 
################################################### ################################################### ################################################### 


Prevalence2 <- function(PhenoAgeBaseline){
  
  PrevalenceBaseline <- PhenoAgeBaseline %>%
    group_by(Disease) %>%
    add_tally() %>%
    select(Disease, n) %>%
    unique()
  
  
  Prevalence <- sapply(1:length(PrevalenceBaseline$n), function(i) { BinomCI(PrevalenceBaseline$n[i], length(levels(as.factor(PhenoAgeBaseline$eid))) ,
                                                                             conf.level = 0.95,
                                                                             method = "clopper-pearson")}
  ) %>% 
    t() %>%
    as.data.frame()
  
  PrevalenceBaseline <- PrevalenceBaseline %>%
    ungroup() %>%
    add_column( Prevalence = Prevalence$V1*100 )
  
  attr(PrevalenceBaseline$Prevalence, "ATT") <- NULL
  
  PrevalenceBaseline <- PrevalenceBaseline %>% 
    mutate(InOut = case_when(Prevalence <0 ~ "Out", TRUE ~ "In")) ##########no filter by prevalence
  
  PrevalenceBaseline$Disease <- stringr::str_to_sentence(PrevalenceBaseline$Disease)
  
  return(PrevalenceBaseline)
  
}


PrevalenceBaseline <- Prevalence2(Rules) %>%
  mutate(Prevalence = round(Prevalence,2)) %>%
  arrange(-Prevalence)

write.csv(PrevalenceBaseline, paste0(IncludeFigHere, "/PrevalenceAll.csv"))
#PrevalenceBaseline <- fread( paste0(IncludeFigHere, "/PrevalenceAll.csv"))

PlotPrev <-  PrevalenceBaseline %>% 
  filter(InOut == "In") %>%
  mutate(Disease = fct_reorder(Disease, Prevalence )) %>%
  ggplot( aes(x=Disease , y=Prevalence, fill=Prevalence)) + 
  geom_bar(stat="identity") + 
  theme_classic() +coord_flip() +
  scale_fill_viridis(name="Prevelance %") + xlab("") + ylab("") + ggtitle("Baseline")

aAll <- PrevalenceBaseline %>% 
  select(Disease, n, Prevalence) %>%
  arrange(-Prevalence)

aAll$Prevalence <- format(aAll$Prevalence, digits = 2)

myt <- ttheme_default(
  core = list(fg_params=list(hjust = 1, x=1),
              bg_params=list(fill=c("white", "pink"))),
  
  colhead = list(fg_params=list(col="white"),
                 bg_params=list(fill="black"))
)


pdf(paste0(IncludeFigHere, "/PrevalenceAll.pdf"), 18,19)
PlotPrev + gridExtra::tableGrob(aAll, theme=myt,rows = NULL)
dev.off()

PlotPrevAge <- list()
aSituation <- list()

for (i in c("Healthy_Low", "Unhealthy_High", "Healthy_Cross", "Unhealthy_Cross")){
  
  print(i)
  
  AddDate2Age <- Rules %>% 
    filter(Situation2 == i)
  
  PrevalenceAge <- Prevalence2(AddDate2Age) 
  
  #if want other filters - filter AddDate2 (e.g above below/60 - sex)
  
  Names <- data.frame("Healthy remaining healthy", 
                      "Unhealthy becoming healthy" ,
                      "Healthy becoming unhealthy",
                      "Unhealthy remaining unhealthy")
  
  names(Names) <- c("Healthy_Low", "Healthy_Cross", "Unhealthy_Cross","Unhealthy_High")
  
  PlotPrevAge[[i]] <-  PrevalenceAge %>% 
    filter(InOut == "In") %>%
    mutate(Disease = fct_reorder(Disease, Prevalence )) %>%
    ggplot( aes(x=Disease , y=Prevalence, fill=Prevalence)) + 
    geom_bar(stat="identity") + theme_classic() +
    coord_flip() +scale_fill_viridis(name="Prevelance %") + xlab("") + ylab("") + 
    ggtitle(paste0(Names[[i]]))
  
  a <- PrevalenceAge %>% 
    select(Disease, n, Prevalence) %>%
    arrange(-Prevalence)
  
  a$Prevalence <- format(a$Prevalence, digits = 2)
  
  aSituation[[i]] <- a
  
  
  myt <- ttheme_default(
    core = list(fg_params=list(hjust = 1, x=1),
                bg_params=list(fill=c("white", "pink"))),
    
    colhead = list(fg_params=list(col="white"),
                   bg_params=list(fill="black"))
  )
  
  
  
  pdf(paste0(IncludeFigHere, "/Prev_",i, ".pdf"),  18,19)
  print(PlotPrevAge[[i]] + gridExtra::tableGrob(a, theme=myt,rows = NULL))
  dev.off()
  
}


ASituation<- aSituation %>%
  bind_rows(.id="Names")

AFin <- rbind(ASituation,aAll %>% mutate(Names = "All"))

pdf(paste0(IncludeFigHere, "/Prevalence.pdf"),  18,19)

ggplot(AFin, aes(fill=Names, y=Disease, x=as.numeric(Prevalence))) + 
  geom_bar(position="dodge", stat="identity") + theme_bw()

dev.off()



pdf(paste0(IncludeFigHere, "/PrevalenceAllBars.pdf"),  18,19)
(PlotPrevAge[[1]] + PlotPrevAge[[2]]) / (PlotPrevAge[[3]] + PlotPrevAge[[4]]) +  
  plot_layout(guides = "collect")
dev.off()

AllPrevalence <- AFin %>% #rbind(AFin,PrevalenceBaseline %>% rename(Names = InOut)) %>%
  pivot_wider(names_from = Names, values_from = c(n, Prevalence), names_glue = "{Names}_{.value}") %>%
  select(order(colnames(.))) %>% 
  arrange(desc(All_n))

write.csv(AllPrevalence, paste0(IncludeFigHere, "/AllPrevalence4CatsToo.csv"))


################################################### ################################################### ################################################### 
################################################### Pairwise Disease Analysis  ################################################### ################################################### 
################################################### ################################################### ################################################### 


#Correlation

Rules <-  fread( paste0(IncludeFigHere, "/Rules.csv")) 


ParticipantsWithNoComorbidities <- Rules %>%
  group_by(eid) %>%
  count() %>%
  filter(n == 1)

Rules <- Rules %>%
  filter(!eid %in%  ParticipantsWithNoComorbidities$eid)


#Did not take out non comorbid participants! - count in number of diseases? 

LegendCategories <- read.csv("/rds/projects/g/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/LegendCategories.csv")%>%
  mutate(DiseaseNew = toupper(Disease))

LegendCategories$Categories <- LegendCategories$Categories  %>%
  str_replace_all(c("Neoplasm" = "Genetic")) 

LegendCategories$DiseaseNew <- LegendCategories$DiseaseNew  %>%
  str_replace_all(.,"\\.", ",") %>% 
  str_replace_all(.,"_", " ") %>%
  str_replace_all(., "NEIROTIC STRESS RELATED AND SOMATOFORM DISEASES", "NEUROTIC, STRESS-RELATED AND SOMATOFORM DISEASES")

write.table(LegendCategories %>% select(Categories, DiseaseNew), paste0(IncludeFigHere,"/LegendCategoriesAppendix.txt"))

#Plot 

DiseasesInTime <- Rules %>% 
  select(eid,Situation2,Disease, Age, Death ) %>%
  mutate(AgeCat = cut(Age, breaks=seq(0,100,by=5), right = FALSE)) %>%
  group_by(Disease, AgeCat) %>%
  add_tally() %>%
  ungroup() %>%
  group_by(Disease, AgeCat, Situation2) %>%
  add_tally(name = "nBySituation") %>%
  ungroup() %>%
  group_by(Disease) %>%
  add_tally(name = "Prevalence") %>%
  ungroup() %>%
  rename(DiseaseNew = Disease) %>%
  left_join(LegendCategories)   %>%
  group_by(Disease, Situation2) %>%
  add_tally(name = "PrevalenceSituation2") %>%
  ungroup()


pdf(paste0(IncludeFigHere,"/DiseaseAgeBySituation.pdf"), 20, 25)
ggplot(DiseasesInTime, aes(AgeCat, nBySituation/PrevalenceSituation2, colour = Situation2,  group = interaction(Disease, Situation2))) +
  geom_point() + 
  geom_line() +  
  theme_bw() +
  #theme(legend.position = "none") + 
  facet_wrap(.~ Disease) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dev.off()

pdf(paste0(IncludeFigHere,"/DiseaseAGeByCat.pdf"))
ggplot(DiseasesInTime, aes(AgeCat, n/Prevalence, colour = Disease, fill = Disease, group = Disease)) +
  geom_point() + 
  geom_line() +  
  theme_bw() +
  theme(legend.position = "none") + 
  facet_wrap(.~ Categories) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

dev.off()


CorrelationDiseases <- function(DFRules,Category, IncludeFigHere) {
  
  # DFRules <- Rules
  #Category <- "Unhealthy_Cross"
  
  PheatmapWide <- DFRules %>% 
    select(-c(Age, Dates)) %>% 
    filter(Situation2 == Category) %>%  #####
  unique() %>%
    pivot_wider(names_from = Disease, values_from = Disease)
  
  FinalDataHMID <- PheatmapWide$eid
  FinalDataHMMatrix <- as.matrix(PheatmapWide %>% 
                                   select(-c(eid:MaxAgeMix)))
  
  FinalDataHMMatrix[which(!is.na(FinalDataHMMatrix) == T)] <- 1
  FinalDataHMMatrix[which(is.na(FinalDataHMMatrix) == T)] <- 0 
  
  class(FinalDataHMMatrix) <- "numeric"
  rownames(FinalDataHMMatrix) = FinalDataHMID
  
  PheatmapCorr <- select(PheatmapWide, -c(eid:MaxAgeMix)) %>%
    mutate_if(is.character, ~replace(., is.na(.), 0)) %>%
    mutate_if(is.character,as.factor) 
  
  #cor <- correlation(PheatmapCorr, method = "biweight") #, redundant  include_factors =   TRUE
  cor <- correlation(PheatmapCorr ,method = "auto", include_factors =   TRUE) #is it working properly?
  
  CorSumRed <- summary(cor, redundant = TRUE)
  
  CorSumRed1 <- CorSumRed %>% 
    select(-contains(".0")) %>% 
    filter(is.na(str_extract(Parameter, ".0")))
  
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  
  CorSumRed1Corr <- as.data.frame(CorSumRed1)
  row.names(CorSumRed1Corr) <- CorSumRed1Corr$Parameter
  CorSumRed1Corr <- CorSumRed1Corr %>%
    select(-Parameter)
  
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd)
    cormat <-cormat[hc$order, hc$order]
  }
  
  
  names(CorSumRed1Corr) <- gsub("\\..*","",names(CorSumRed1Corr))
  rownames(CorSumRed1Corr) <- gsub("\\..*","",rownames(CorSumRed1Corr))
  
  cormat <- reorder_cormat(CorSumRed1Corr) 
  
  
  upper_tri <- get_upper_tri(as.matrix(cormat))
  # Melt the correlation matrix
  melted_cormat <- melt(upper_tri, na.rm = TRUE)
  # Create a ggheatmap
  
  ggheatmap <- ggplot(melted_cormat , aes(Var2, Var1, fill = value))+
    geom_tile(color = "white")+
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                         midpoint = 0, limit = c(-1,1), space = "Lab", 
                         name="Pearson\nCorrelation") +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                     size = 5, hjust = 1)) +
    coord_fixed()
  # Print the heatmap
  
  #pdf(paste0(IncludeFigHere, "/PaiwriseCorrelation2", Category, ".pdf"), 20, 25)
  
  print(ggheatmap)
  
  #dev.off()
  
  return(list(FinalDataHMMatrix,cormat))
  
}




#Relative Risk



#Take out chromosomal abnormalties -- too small 

RelativeRisk <- function(FinalDataHMMatrix, IncludeFigHere,number) {
  
  
  # Hyper <- PheatmapWide %>%
  #  filter(eid %in% filter(PheatmapWide, Disease == "CHRONIC INFECTIOUS DISEASES")$eid)
  #length(unique(Hyper$eid))
  # i <- "Unhealthy_High"
  
  FinalDataHMMatrix <- CorrelationList[[i]][[1]]
  
  
  #StoreFinal
  #disB <- StoreFinal[,disBnm]
  #disA <- StoreFinal[,disAnm]
  
  DataRR <- data.frame(FinalDataHMMatrix) %>%
    # select(-CHROMOSOMAL_ABNORMALITIES) %>%
    add_rownames() %>%
    melt() 
  
  disMat <- DataRR %>% 
    filter(value == 1) %>%
    spread(variable,value,fill=0) %>%
    as.data.frame()
  
  rownames(disMat) <- disMat$rowname
  disMat$rowname <- NULL
  disMat <- as.matrix(disMat)
  
  RRtable <- lapply(colnames(disMat), function(disAnm){
    
    disA <- disMat[,disAnm]
    
    xx <- as.tibble(t(sapply(colnames(disMat), function(disBnm){
      
      disB <- disMat[,disBnm]
      Nab <- sum((disA == 1) & (disB == 1))
      Nnab <- sum((disA == 0) & (disB == 1))
      Ta <- sum(disA == 1)
      Tna <- sum(disA == 0)
      Pexp <- Nab / Ta
      Pnexp <- Nnab / Tna
      RR <- Pexp / Pnexp
      cif <- sqrt((((Ta - Nab) / Nab) / Ta) + (((Tna - Nnab) / Nnab) / Tna))
      ll <- exp(log(RR) - (1.96 * cif))
      ul <- exp(log(RR) + (1.96 * cif))
      c(Nab=Nab,Nnab=Nnab,Ta=Ta,Tna=Tna,Pexp=Pexp,Pnexp=Pnexp,RR=RR,cif=cif,ll=ll,ul=ul)
      
    }))) %>%
      
      mutate(disB = colnames(disMat))
    
    return(xx)
    
  })
  
  names(RRtable) <- colnames(disMat)
  RRtable <- reshape2::melt(RRtable, id.vars = colnames(RRtable[[1]])) %>%
    rename(disA=L1)
  
  RRtableFinal <- RRtable %>%
    select(disA, disB, everything()) %>%
    arrange(-RR) 
  
  RRmat <- RRtableFinal %>%
    #filter(sublevel==F & uplevel == F) %>%
    arrange(-RR) %>%
    select(disA, disB, RR) %>%
    spread(disB,RR,fill = 1) %>%
    as.data.frame()
  
  rownames(RRmat) = RRmat$disA
  RRmat$disA = NULL
  RRmat = as.matrix(RRmat)
  diag(RRmat)=NA
  mymat <- log2(RRmat)
  
  is.na(mymat)<-sapply(mymat, is.infinite)
  
  colnames(mymat) <- gsub( "\\.", " ", colnames(mymat))
  rownames(mymat) <- gsub( "\\.", " ", rownames(mymat))
  
  print(pheatmap(mymat, cellwidth = 15, cellheight = 15, 
                 # filename = paste0(IncludeFigHere, "/RRmat_log2_all", number, ".pdf"), 
                 color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49), 
                 breaks = seq(-5,5,length.out = 50),
                 cutree_rows = 10,
                 #annotation_row = annotdata, 
                 #annotation_col = annotdata,
                 #annotation_colors = annotcols,
                 cutree_cols = 10, cluster_rows = T, cluster_cols = T))
  
  
  RRnet <- RRtable %>%
    filter( ((RR > 1 & ll > 1) | (RR < 1 & ul < 1))) %>% #sublevel==F & uplevel == F &
    arrange(-RR) %>%
    select(disB, disA, RR) %>%
    mutate(RR = log2(RR)) %>%
    filter(RR >= 1  | RR <= -1)%>%
    mutate(sizex = abs(RR) / 5)# %>%
  
  RRnet2 <- RRtable %>%
    filter( ((RR > 1 & ll > 1) | (RR < 1 & ul < 1))) %>% #sublevel==F & uplevel == F &
    arrange(-RR) %>%
    select(disB, disA, RR, cif) %>%
    mutate(RR = log2(RR)) %>%
    filter(RR >= 1  | RR <= -1)%>%
    mutate(sizex = abs(RR) / 5)# %>%
  #filter(!is.infinite(RR)) 
  
  return(list(RRtableFinal, RRnet, RRnet2))
  #filter(!is.infinite(RR)) 
  
  
}


RelativeRisk2 <- function(FinalDataHMMatrix, IncludeFigHere,number) {
  
  
  Store <- list()
  for (i in names(CorrelationList)){
    
    FinalDataHMMatrix <- CorrelationList[[i]][[1]]
    
    DataRR <- data.frame(FinalDataHMMatrix) %>%
      select(-CHROMOSOMAL.ABNORMALITIES) %>%
      add_rownames() %>%
      melt() 
    
    Store[[i]] <- DataRR %>% 
      filter(value == 1) %>%
      spread(variable,value,fill=0) %>%
      as.data.frame()
    
  }
  
  
  StoreFinal <- Store %>%
    bind_rows()
  
  FinalDataHMMatrix <- CorrelationList[[i]][[1]]
  
  
  DataRR <- data.frame(FinalDataHMMatrix) %>%
    select(-CHROMOSOMAL.ABNORMALITIES) %>%
    add_rownames() %>%
    melt() 
  
  disMat <- DataRR %>% 
    filter(value == 1) %>%
    spread(variable,value,fill=0) %>%
    as.data.frame()
  
  rownames(disMat) <- disMat$rowname
  disMat$rowname <- NULL
  disMat <- as.matrix(disMat)
  
  RRtable <- lapply(colnames(disMat), function(disAnm){
    
    #disAnm <- "CHROMOSOMAL.ABNORMALITIES"
    #disBnm <- "OTHER.PSYCHIATRIC.AND.BEHAVIORAL.DISEASES"
    
    
    disA <- disMat[,disAnm]
    disAS <- StoreFinal[,disAnm]
    
    xx <- as.tibble(t(sapply(colnames(disMat), function(disBnm){
      
      disB <- disMat[,disBnm]
      disBS <- StoreFinal[,disBnm]
      Nab <- sum((disA == 1) & (disB == 1))
      Nnab <- sum((disAS == 0) & (disBS == 1))
      Ta <- sum(disA == 1)
      Tna <- sum(disAS == 0)
      Pexp <- Nab / Ta
      Pnexp <- Nnab / Tna
      RR <- Pexp / Pnexp
      cif <- sqrt((((Ta - Nab) / Nab) / Ta) + (((Tna - Nnab) / Nnab) / Tna))
      ll <- exp(log(RR) - (1.96 * cif))
      ul <- exp(log(RR) + (1.96 * cif))
      c(Nab=Nab,Nnab=Nnab,Ta=Ta,Tna=Tna,Pexp=Pexp,Pnexp=Pnexp,RR=RR,cif=cif,ll=ll,ul=ul)
      
    }))) %>%
      
      mutate(disB = colnames(disMat))
    
    return(xx)
    
  })
  
  names(RRtable) <- colnames(disMat)
  RRtable <- reshape2::melt(RRtable, id.vars = colnames(RRtable[[1]])) %>%
    rename(disA=L1)
  
  RRtableFinal <- RRtable %>%
    select(disA, disB, everything()) %>%
    arrange(-RR) 
  
  RRmat <- RRtableFinal %>%
    #filter(sublevel==F & uplevel == F) %>%
    arrange(-RR) %>%
    select(disA, disB, RR) %>%
    spread(disB,RR,fill = 1) %>%
    as.data.frame()
  
  rownames(RRmat) = RRmat$disA
  RRmat$disA = NULL
  RRmat = as.matrix(RRmat)
  diag(RRmat)=NA
  mymat <- log2(RRmat)
  
  is.na(mymat)<-sapply(mymat, is.infinite)
  
  colnames(mymat) <- gsub( "\\.", " ", colnames(mymat))
  rownames(mymat) <- gsub( "\\.", " ", rownames(mymat))
  
  print(pheatmap(mymat, cellwidth = 15, cellheight = 15, 
                 # filename = paste0(IncludeFigHere, "/RRmat_log2_all", number, ".pdf"), 
                 color = colorRampPalette(c(rev(brewer.pal(8,'Blues')),'white',brewer.pal(8,'Reds')))(49), 
                 breaks = seq(-5,5,length.out = 50),
                 cutree_rows = 10,
                 #annotation_row = annotdata, 
                 #annotation_col = annotdata,
                 #annotation_colors = annotcols,
                 cutree_cols = 10, cluster_rows = T, cluster_cols = T))
  
  
  RRnet <- RRtable %>%
    filter( ((RR > 1 & ll > 1) | (RR < 1 & ul < 1))) %>% #sublevel==F & uplevel == F &
    arrange(-RR) %>%
    select(disB, disA, RR) %>%
    mutate(RR = log2(RR)) %>%
    filter(RR >= 1  | RR <= -1)%>%
    mutate(sizex = abs(RR) / 5)# %>%
  #filter(!is.infinite(RR)) 
  
  
  RRnet2 <- RRtable %>%
    filter( ((RR > 1 & ll > 1) | (RR < 1 & ul < 1))) %>% #sublevel==F & uplevel == F &
    arrange(-RR) %>%
    select(disB, disA, RR, cif) %>%
    mutate(RR = log2(RR)) %>%
    filter(RR >= 1  | RR <= -1)%>%
    mutate(sizex = abs(RR) / 5)# %>%
  #filter(!is.infinite(RR)) 
  
  return(list(RRtableFinal, RRnet, RRnet2))
  
}





#Plots 

GraphsRR <- function(RRnet,LegendCategories, number){
  
  RRnet2 <- RRnet %>%
    igraph::graph_from_data_frame(directed = T)
  
  
  AA <- left_join(data.frame(DiseaseL = names((V(RRnet2)))), LegendCategories)
  
  DiseaseB <- inner_join( data.frame(DiseaseL = RRnet$disB), select(AA, DiseaseL, Categories))
  
  
  RRnet1 <- RRnet %>%
    mutate(DiseaseB = inner_join( data.frame(DiseaseL = disB), select(AA, DiseaseL, Categories))$Categories,
           DiseaseA = inner_join( data.frame(DiseaseL = disA), select(AA, DiseaseL, Categories))$Categories, 
           Same = ifelse(DiseaseB == DiseaseA, DiseaseB, "0"), 
           FinalDisease = case_when(Same == 1 ~ as.character(DiseaseB), TRUE ~ "Nothing") 
    ) %>%
    filter(RR >2 | RR < -2)
  
  RRnet2 <- RRnet1 %>%
    igraph::graph_from_data_frame(directed = T)
  
  V(RRnet2)$category = as.character(AA$Categories)
  E(RRnet2)$signx = c('darkred','midnightblue')[1+(E(RRnet2)$RR < 0)]
  
  V(RRnet2)$name <- gsub("\\.", " ", V(RRnet2)$name)
  
  
  print(ggnet2(RRnet2,
               arrow.size = 2,
               arrow.gap = 0.002, 
               edge.size = 'sizex', 
               size = 3, 
               node.color = 'category',
               edge.color = c("color", "grey75"),
               palette = "Paired",
               #palette = c("vowel" = "gold", "consonant" = "grey")
               mode = "circle",
               #edge.color = 'signx',
               edge.alpha = 0.75, 
               color.legend = "ICD10 Categories") + 
          geom_label_repel(label= V(RRnet2)$name, 
                           size=2, 
                           angle = 25, 
                           #alpha = 0.6, 
                           label.padding=.1,
                           fill = "white") + #force_pull   =  -2
          theme(legend.position = 'bottom',
                #legend.title = element_blank(),
                legend.text = element_text(size = 6))) #+
  # scale_color_manual(labels = c("Circulatory"   ,"Digestive"  ,   "Endocrine"    , "Eye/ear"     ,  "Genitourinary", "Infectious"  ,  "Mental"  ,      #"Muscoskeletal" ,"Neoplasm"  , "Nervous" ,      "Respiratory" ,  "Skin"), values = c("red", "red", "red","orange", "pink", "purple", "brown", "black", #"white", "grey", "dark blue", "blue"))
  
  ggsave(filename = paste0(IncludeFigHere, '/CirclePlotColour',number,'.pdf'),useDingbats = F, units = 'cm', width = 19,height = 18)
  
  
  ##edge.color = 'signx', 
  print(ggnet2(RRnet2,
               arrow.size = 2,
               arrow.gap = 0.002, 
               edge.size = 'sizex', 
               size = 3, 
               node.color = 'category',
               edge.color = c("color", "grey75"),
               palette = "Paired",
               mode = "princoord",
               color.legend = "ICD10 Categories",
               #edge.color = 'signx',
               edge.alpha = 0.75) + 
          geom_label_repel(label= V(RRnet2)$name, 
                           size=1.5, 
                           angle = 25, 
                           alpha = 0.6, 
                           label.padding=.1,
                           fill = "white") + #force_pull   =  -2
          theme(legend.position = 'bottom',
                #legend.title = element_blank(),
                legend.text = element_text(size = 6)) 
  )
  ggsave(filename = paste0(IncludeFigHere, '/NoCircle',number,'.pdf'),useDingbats = F, units = 'cm', width = 15,height = 12)
  
  
}


CorrelationList <- list()
RelativeRiskList <- list()
RelativeRiskList2 <- list()
Table <- list()
Rank <- list()



for(i in unique(Rules$Situation2)){ #unique(Rules$Situation2)
  
  # #grid.newpage()
  print(i)
  
  pdf(paste0(IncludeFigHere,"/ALL_2",i,".pdf"), 18,22)
  
  CorrelationList[[i]] <- print(CorrelationDiseases(Rules, i,IncludeFigHere))
  RelativeRiskList[[i]] <- print(RelativeRisk(CorrelationList[[i]][[1]], IncludeFigHere,i))
  RelativeRiskList2[[i]] <- print(RelativeRisk2(CorrelationList[[i]][[1]], IncludeFigHere,i))
  GraphsRR(RelativeRiskList[[i]][[2]],LegendCategories, i)
  
  ######## Clustering
  
  #load("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930/RelativeRiskList2.RData")
  
  #i <- "Unhealthy_High"
  RRnet <- RelativeRiskList2[[i]][[3]] %>% #RelativeRiskList2[[i]][[2]]
    filter(cif < 0.2 & (RR> 2 | RR <  -2)) %>%
    select(-cif)
  
  
  RRnet2 <- RRnet %>%
    igraph::graph_from_data_frame(directed = T)
  
  
  AA <- left_join(data.frame(DiseaseL = names((V(RRnet2)))), LegendCategories)
  
  y = RColorBrewer::brewer.pal(length(unique(LegendCategories$Categories)), "Paired")
  names(y) = (unique(LegendCategories$Categories))
  
  #s <- data.frame(y) %>%
  #  add_rownames(var = "Categories")
  
  #AA<- AA %>% 
  #  inner_join(., s)
  
  DiseaseB <- inner_join( data.frame(DiseaseL = RRnet$disB), select(AA, DiseaseL, Categories))
  
  RRnet1 <- RRnet %>%
    mutate(DiseaseB = inner_join( data.frame(DiseaseL = disB), select(AA, DiseaseL, Categories))$Categories,
           DiseaseA = inner_join( data.frame(DiseaseL = disA), select(AA, DiseaseL, Categories))$Categories, 
           Same = ifelse(DiseaseB == DiseaseA, DiseaseB, "0"), 
           FinalDisease = case_when(Same == 1 ~ as.character(DiseaseB), TRUE ~ "Nothing"))
  
  
  
  RRnet2 <- RRnet1 %>%
    igraph::graph_from_data_frame(directed = T)
  
  
  V(RRnet2)$category = as.character(AA$Categories)
  E(RRnet2)$signx = c('darkred','midnightblue')[1+(E(RRnet2)$RR < 0)]
  
  V(RRnet2)$name <- str_to_sentence(gsub("\\.", " ", V(RRnet2)$name))
  
  saveRDS(RRnet2, paste0(IncludeFigHere,"/RRnet2_Plot_2", i,".RDS"))
  
  print(ggnet2(RRnet2, arrow.size = 2, arrow.gap = 0.02, edge.size = 'sizex', size = 6, edge.color = 'signx', 
               edge.alpha = 0.6, node.color = 'category', palette = y) +  #node.color = 'category', palette = discatcolors
          geom_text_repel(label= V(RRnet2)$name,  box.padding = 0.5,  size = 2.5,
                          max.overlaps = 4) + #max.overlaps = Inf,  fontface = 'bold', fontface = "bold",
          theme(legend.position = 'none')) #bottom
  #       legend.title = element_blank(),
  #       legend.text = element_text(size = 10)
  # ))
  ggsave(filename = paste0(IncludeFigHere,"/RRnet_log2_sig_log2RR1Size2", i,".pdf"),useDingbats = F, units = 'cm', width = 20,height = 16) #15, 12
  
  ggsave(filename = paste0(IncludeFigHere,"/RRnet_log2_sig_log2RR1Size2", i,".pdf"),useDingbats = F, units = 'cm', width = 17,height = 10) #15, 12
  
  
  #ggsave(filename = paste0(IncludeFigHere,"/RRnet_log2_sig_log2RR1.png"), units = 'cm', width = 18,height = 15)
  
  PR <- (page.rank(RRnet2, damping = 0.85)) 
  
  Rank[[i]] <- PR$vector %>%
    as.data.frame() %>%
    add_rownames() %>%
    rename("PR" = ".") 
  
  print(ggplot(Rank[[i]], aes(PR, fct_reorder(rowname, PR))) + geom_bar(stat = "identity") +
          theme_bw()
  ) 
  
  dev.off()
  
  Table[[i]] <- RelativeRiskList2[[i]][[1]] %>%   #RelativeRiskList before!!
    filter(Nnab != 0 ) %>%
    filter(Nab != 0 ) %>%
    #top_n(20, RR) %>% 
    select(disA, disB, RR, cif, ll, ul)
  
  #  write.csv(Table, paste0(IncludeFigHere,"/RRTable", i, ".csv"), row.names = F)
  
  
}


save(RelativeRiskList, file =  paste0(IncludeFigHere,"/RelativeRiskList.RData"))
save(CorrelationList, file =  paste0(IncludeFigHere,"/CorrelationList.RData"))
save(Table, file =  paste0(IncludeFigHere,"/RRTableFINAL.RData"))
save(RelativeRiskList2, file =  paste0(IncludeFigHere,"/RelativeRiskList2.RData"))
save(Rank, file =  paste0(IncludeFigHere,"/Rank.RData"))

#load(paste0(IncludeFigHere,"/CorrelationList.RData"))
#load(paste0(IncludeFigHere,"/RRTableFINAL.RData"))
#load(paste0(IncludeFigHere,"/RelativeRiskList.RData"))
#load(paste0(IncludeFigHere,"/RelativeRiskList2.RData"))
#load(paste0(IncludeFigHere,"/RelativeRiskList2.RData"))

TableFinal <- do.call(rbind, Table) %>%
  add_rownames() %>%
  mutate(rowname = gsub("\\..*","",rowname)) %>%
  #mutate(Diseases = paste0(gsub("\\.", " " , disA), " - ",gsub("\\.", " " , disB)) ) %>%
  #select(-c(disA, disB)) %>%
  pivot_wider( names_from = rowname,  values_from = c(RR, cif, ll, ul), names_glue = "{rowname}_{.value}") %>%
  select(order(colnames(.))) %>% 
  arrange(desc(Unhealthy_High_RR)) %>%
  mutate(across(where(is.numeric), round, 2))


TableFinal$disA <- gsub("\\."," ", TableFinal$disA)
TableFinal$disB <- gsub("\\."," ", TableFinal$disB)

TableFinalGT <- TableFinal %>%
  mutate_all(., as.character) %>%
  mutate(across(everything(), ~ifelse(is.na(.), "-", .))) %>%
  #mutate_all(., list(~na_if(.,"-"))) %>%
  gt() %>%
  tab_header(
    title = "Relative Risk",
    subtitle = "Pairwise associatons between disease categories for all trajectories"
  ) %>%
  tab_spanner(
    label = "Healthy remaining Healthy",
    columns = c(Healthy_Low_RR, Healthy_Low_cif, Healthy_Low_ll, Healthy_Low_ul)
  ) %>% 
  tab_spanner(
    label = "Unhealthy becoming Healthy",
    columns = c(Healthy_Cross_RR, Healthy_Cross_cif, Healthy_Cross_ll, Healthy_Cross_ul)
  ) %>% 
  tab_spanner(
    label = "Unhealthy remaining Unhealthy",
    columns = c(Unhealthy_High_RR, Unhealthy_High_cif, Unhealthy_High_ll, Unhealthy_High_ul)
  ) %>% 
  tab_spanner(
    label = "Healthy becoming Unhealthy",
    columns = c(Unhealthy_Cross_RR, Unhealthy_Cross_cif, Unhealthy_Cross_ll, Unhealthy_Cross_ul)
  ) %>%
  cols_label(
    disA = html("Disease A"),
    disB = html("Disease B"), 
    Unhealthy_Cross_RR  = html("RR"),
    Unhealthy_Cross_cif  = html("cf"),
    Unhealthy_Cross_ll  = html("ll"),
    Unhealthy_Cross_ul  = html("ul"),
    
    Healthy_Cross_RR  = html("RR"),
    Healthy_Cross_cif  = html("cf"),
    Healthy_Cross_ll  = html("ll"),
    Healthy_Cross_ul  = html("ul"),
    
    Unhealthy_High_RR  = html("RR"),
    Unhealthy_High_cif  = html("cf"),
    Unhealthy_High_ll  = html("ll"),
    Unhealthy_High_ul  = html("ul"),
    
    Healthy_Low_RR  = html("RR"),
    Healthy_Low_cif  = html("cf"),
    Healthy_Low_ll  = html("ll"),
    Healthy_Low_ul  = html("ul")
  )


write.csv(TableFinal, paste0(IncludeFigHere,"/RRTableFINAL.csv"), row.names = F)
gtsave(TableFinalGT, paste0(IncludeFigHere,"/RRFINALGT.rtf")) #before TableGT only?

Table <- read.csv(paste0(IncludeFigHere,"/RRTableFINAL.csv"))

#ExtraDiabetes


AssocDiab <- Table %>%
  filter(disA %in% "DIABETES" | disB %in% "DIABETES") %>%
  filter(if_any(c("Healthy_Cross_RR", "Unhealthy_Cross_RR", "Healthy_Low_RR", "Unhealthy_High_RR"), ~ . > 2 )) %>%
  filter(if_any(c("Healthy_Cross_cif", "Unhealthy_Cross_cif", "Healthy_Low_cif", "Unhealthy_High_cif"), ~ . < 0.2 )) %>%
  select(contains("Healthy_Low"),contains("Healthy_High"), contains("dis")) #%>%
#select(contains("dis"), contains("RR"))

write.csv(AssocDiab, paste0(IncludeFigHere,"/AssocDiab.csv"), row.names = F)

#Understanding

Healthy_Cross_Table <-  Table %>%
  filter(Healthy_Cross_RR > 2 & Healthy_Cross_cif < 0.2) %>%
  arrange(-Healthy_Cross_RR) 

Healthy_Cross_TableRepeated <- data.frame(Diseases = c(Healthy_Cross_Table$disA, Healthy_Cross_Table$disB)) %>% 
  count(Diseases) %>%
  arrange(-n)

###

Healthy_Low_Table <-  Table %>%
  filter(Healthy_Low_RR > 2 & Healthy_Low_cif < 0.2) %>%
  arrange(-Healthy_Low_RR)

Healthy_Low_TableRepeated <- data.frame(Diseases = c(Healthy_Low_Table$disA, Healthy_Low_Table$disB)) %>% 
  count(Diseases) %>%
  arrange(-n)

###

Unhealthy_High_Table <-  Table %>%
  filter(Unhealthy_High_RR > 2 & Unhealthy_High_cif < 0.2) %>%
  arrange(-Unhealthy_High_RR)

Unhealthy_High_TableRepeated <- data.frame(Diseases = c(Unhealthy_High_Table$disA, Unhealthy_High_Table$disB)) %>% 
  count(Diseases) %>%
  arrange(-n)

###

Unhealthy_Cross_Table <-  Table %>%
  filter(Unhealthy_Cross_RR > 2 &  Unhealthy_Cross_cif < 0.2) %>%
  arrange(-Unhealthy_Cross_RR)

Unhealthy_Cross_TableRepeated <- data.frame(Diseases = c(Unhealthy_Cross_Table$disA, Unhealthy_Cross_Table$disB)) %>% 
  count(Diseases) %>%
  arrange(-n)

###

PivotTable <- Table %>%
  pivot_longer(
    cols = Healthy_Cross_cif:Unhealthy_High_ul,
    names_to = c("Type", "Type2", "Reading"),
    names_pattern = "(.*)_(.*)_(.*)",
    values_to = "count"
  ) %>%
  mutate(Category = paste0(Type, "_",Type2)) %>%
  pivot_wider(names_from = "Reading", 
              values_from = "count") %>%
  select(-c(Type, Type2)) %>%
  drop_na() %>%
  mutate(TotalRR = paste0(RR, " (", ll, " - ", ul, ")")) %>%
  select(-c(cif, RR, ll, ul)) %>%
  pivot_wider(
    names_from = "Category", 
    values_from = "TotalRR") %>%
  select(c(disA, disB, Healthy_Low, Healthy_Cross, Unhealthy_Cross, Unhealthy_High))


write.csv(PivotTable, paste0(IncludeFigHere,"/PivotTableRRPasted.csv"), row.names = F)

################################################### ################################################### ################################################### 
################################################### Extra Multimorbidity Info  ################################################### ################################################### 
################################################### ################################################### ################################################### 


Rules <-  fread( paste0(IncludeFigHere, "/Rules.csv")) ##################
###########################
###########################
###########################
###########################
#Description of rules dataset - number of diseases per patient etc? 

Stats <- Rules %>%
  group_by(eid) %>%
  count() #%>%
#filter(n != 1) #115,912 with this butnot needed! In the end category does not mean multimorbid, we are interested in types of diseases more than anything else right? 

pdf(paste0(IncludeFigHere, "/MultimorbidityHistogram.pdf"), 4, 3)
ggplot(Stats, aes(n)) + geom_histogram() + theme_bw() + labs(x = "Number of diseases", y = "Frequency")
dev.off()

################################################### ################################################### ################################################### 
################################################### GWAS ################################################### ################################################### 
################################################### ################################################### ################################################### 




Slopes <- fread(paste0(IncludeFigHere,  "/SlopesAndClassifictaion3.csv"))

SelectedAll <-  fread(paste0(IncludeFigHere,  "/Rules.csv")) %>%
  select(eid, MaxAgeMix) %>%
  unique()

AllBiomarkers <- fread("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/AllBiomarkersNames.csv") %>%
  select(eid, genetic_sex_male) %>%
  unique()

InfoCov <- SelectedAll %>%
  inner_join(AllBiomarkers)


set.seed(132) 


GWAS <- Slopes %>% 
  select(eid,  Situation2) %>%
  unique() %>%
  filter(Situation2 != "Unknown") %>%
  filter(eid %in% InfoCov$eid) # why not same?

fwrite(GWAS, paste0(IncludeFigHere, "/GWASData.csv"))

######################### Pheno file 

Pheno <- data.frame(FID = GWAS$eid, IID = GWAS$eid, fastDummies::dummy_cols(GWAS$Situation2 )) 

names(Pheno) <- c("FID", "IID", "Delete", "HealthyCross", "Healthy", "UnhealthyCross", "Unhealthy")

Pheno <- Pheno %>%
  dplyr::select(-Delete) %>%
  drop_na()

write.table(Pheno, file = paste0(IncludeFigHere,"/DummyClust.phen"),  row.names = FALSE, quote = FALSE) #copy and paste to jupyter directly. 


######################## Generate my own covariate file


#include PCAs? - just get full and inner join ( field i ukbiobank)

Cov <- read.table("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/Ethicity_20211026_1222/full_inclusive.cov", sep = "", header = TRUE) %>%
  dplyr::select(FID, pc1:array_used)

Covariate <- InfoCov %>%
  group_by(eid) %>% 
  #slice(which.min(MaxAgeMix)) %>% 
  filter(eid %in% GWAS$eid) %>%
  mutate(eid2 = eid) %>%
  dplyr::select(eid, eid2, genetic_sex_male, MaxAgeMix) 

names(Covariate) <- c("FID", "IID", "sex", "age") 

Covariate2 <- Covariate %>%
  left_join(Cov) %>%
  mutate(age = round(age)) #%>%
# mutate(age2 = age^2)

write.table(Covariate2, file = paste0(IncludeFigHere,"/full_inclusiveMine.cov"),  row.names = FALSE, quote = FALSE) #copy and paste to jupyter directly. 

#Think would be better to get the odlest age diagnosed as covariate??

#Civ <- read.table("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/Ethicity_20211026_1222/full_inclusiveMine.cov", sep = "", header = TRUE) 

############## RUN GWAS ######################

#copy GWAS files? 

system(paste0("cp /castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/WhiteEthnic_20211116_0302/regenie_step1.sh " ,IncludeFigHere,"/regenie_step1.sh"))
system(paste0("cp /castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/Ethicity_20211026_1222/regenie_step2.sh " ,IncludeFigHere,"/regenie_step2.sh"))

#IncludeFigHere <- "/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/AugustManuscript2_20220831_1818"

#setwd(IncludeFigHere)
#system("sbatch regenie_step1.sh") 
#system("sbatch regenie_step2.sh")

################ READ GWAS ###################
#"/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/AugustManuscript2_20220831_1818/regenie_step1.sh"

system(paste0("cp /castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/JanetGWAS/RuthEthnicity2.sh " ,IncludeFigHere,"/RuthEthnicity2.sh"))
system(paste0("cp /castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/JanetGWAS/RuthEthnicityScript3GWAS.R " ,IncludeFigHere,"/RuthEthnicityScript3GWAS.R"))

#system(paste0("sbatch", IncludeFigHere, "/RuthEthnicity2.sh")) #now found here

################################################### ################################################### ################################################### 
################################################### Now To run 2_SummaryStats.R (RuthEthnicityScript3GWAS.R) through sbatch (RuthEthnicity2.sh) ################################################### ################################################### 
################################################### ################################################### ################################################### 


#Now To run RuthEthnicityScript3GWAS.R 
#system("sbatch /castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/ManSlurm_20220801_1138/RuthEthnicity2.sh")
#sbatch /castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/AugustManuscript2_20220831_1818/RuthEthnicity2.sh


#eQTL downloading in proximityAndEQTL.R - GTEx
#SummaryStats analysis Proximity and eQTL in file ProximityAndEQTLFinal3



################################################### ################################################### ################################################### 
################################################### Now To run 2_SummaryStats.R (RuthEthnicityScript3GWAS.R)  ################################################### ################################################### 
################################################### ################################################### ################################################### 

################################################### ################################################### ################################################### 
################################################### With results (ProximityWith8.tsv - follow up with 3_GWASEnrichment.R file)  ################################################### ################################################### 
################################################### ################################################### ################################################### 

################################################### ################################################### ################################################### 
################################################### Ageing Databases  ################################################### ################################################### 
################################################### ################################################### ################################################### 


#ageing databases: https://agingbiotech.info/databases/
#IncludeFigHere <- "/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930"

#https://immunityageing.biomedcentral.com/articles/10.1186/s12979-021-00232-1#Sec14

#ImmAge <-  read_excel("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/ImmAge_FunctinalConserv.xlsx") %>%
#  rename(symbol = external_gene_name) %>%
#  mutate(symbol = toupper(symbol)) %>%
#  add_column(Type = "ImmAge" ) %>%
#  mutate(why = description) %>%
#  select(symbol, why, Type)

# same thing thn orthologs so taking it out. 


ImmAgeOrthologs <- read_excel("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/SharedOrthologs_FunctionalCons.xlsx") %>%
  filter(grepl( "humans" ,Overlap )) %>% 
  rename(symbol = HumanSymbol) %>%
  mutate(symbol = toupper(symbol)) %>%
  add_column(Type = "ImmAge Orthologs" ) %>%
  mutate(why = paste0(Overlap, ": ", HumanDescription)) %>%
  select(symbol, why, Type)


AgingBank <- fread("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/AgingBank_All.txt") %>%
  rename(symbol = name) %>%
  mutate(symbol = toupper(symbol)) %>%
  add_column(Type = "Aging Bank" ) %>%
  rename('pubmedid' = 'pubmed id') %>%
  rename('AntiPro' = 'Anti/Pro') %>%
  mutate(why = paste0(Organism, "_", sample, ": ", AntiPro, " in ", pubmedid)) %>%
  select(symbol, why, Type)

#sig_res <- read.table(paste0(IncludeFigHere,"/AllClust_all_annos.tsv"), header = TRUE) %>%
# rename(symbol = gene) %>%
# select(-c( GENPOS))

#table(sig_res$Type)

#ff <- fread(paste0(IncludeFigHere, "/MergeAll.csv")) %>% 
#  ggplot(., aes(Variable, fill = Type))




sig_res <- fread(paste0("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/FinalAnalysisNov/JanuaryManuscript2_20230106_1930/MergeAll8.csv")) %>% ##result from GWAS!!! 
  rename(symbol = gene) %>%
  rename(ID = SNP) %>%
  rename(Annotation = Type) %>%
  select(-c(ALLELE0, ALLELE1, Common)) %>%
  filter(symbol != "")

sig_res$symbol <- gsub("\\*", "", sig_res$symbol)

#gg <- inner_join(AgingBank, sig_res)

#inflammation related patwhays - blast agaisnt https://www.genome.jp/pathway/hsa04211 (KEGG longevity) - mmm not working? 

cellAge1 <- read_delim("/rds/projects/g/gkoutosg-variant-prediction/Laura2/oldlxb732/BiobankPhenoAge/Datasets/cellAge1.csv",  ";", escape_double = FALSE, trim_ws = TRUE) %>%
  mutate(why = paste0(organism,"_", senescence_effect)) %>% #description out
  select(c(gene_name, why)) %>%
  rename(symbol = gene_name) %>%
  add_column(Type = "Cell Age" )

genage_models <- read_csv("/rds/projects/g/gkoutosg-variant-prediction/Laura2/oldlxb732/BiobankPhenoAge/Datasets/genage_models.csv") %>%
  mutate(why = paste0(organism ,"_", `avg lifespan change (max obsv)` ,"_",  `lifespan effect` ,"_", `longevity influence`)) %>%
  select(c(symbol, why)) %>%
  add_column(Type = "GenAge Models" )

genage_human <- read_csv("/rds/projects/g/gkoutosg-variant-prediction/Laura2/oldlxb732/BiobankPhenoAge/Datasets/genage_human.csv") %>%
  select(c(symbol, why)) %>%
  add_column(Type = "GenAge Human" )


LongevityGenesZenin2019 <- read_excel("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/TreeWAS/LongevityGenesZenin2019.xlsx", skip = 1) %>%
  mutate(Locus = gsub( "CHRNA3/5", "CHRNA3/CHRNA5", Locus)) %>%
  separate_rows(Locus, sep = "/") %>% 
  mutate(why = paste0("Zenin2019: ", trait, "_", DOI)) %>% 
  select(SNP, symbol = Locus, why)

######## https://www.nature.com/articles/s41576-019-0183-6/tables/1 do they come from here? - https://elifesciences.org/articles/39856/figures#content

#snpsNature <- data.frame( SNP = c("rs429358", "rs10455872","rs8042849", "rs142158911", "rs11065979", "rs1556516", "rs34967069", "rs1230666", "rs12924886",
#                "rs1275922", "rs6224", "rs61348208", "rs7844965", "rs4774495", "rs599839", "rs3131621", "rs15285", "rs9872864"),
#                why = "Nature: Logevity", symbol = "None")
## big disease associated SNPs too in same article

#nature <- (intersect(sig_res$ID, snpsNature ))

##########
#https://elifesciences.org/articles/39856/figures#content

eLife <- read_excel("/castles/nr/projects/2017/gkoutosg-variant-prediction/ukbiobank/projects/Multimorbidity/TreeWAS/AlleLIFELifespanGenes.xlsx") %>% 
  rename(symbol = "At or near") %>%
  mutate(symbol = gsub( "CHRNA3/5", "CHRNA3/CHRNA5", symbol)) %>%
  separate_rows(symbol, sep = "/") %>%
  mutate(EffectAge = ifelse(Years < 0, "Decrease", "Increase")) %>% 
  mutate(why = paste0("eLife: ", EffectAge)) %>%
  select(SNP = rsID, symbol, why)

#eLifeSNPs <- (intersect(sig_res$ID, eLife$rsID )) #nothing (strange, not same as above?)
#LongevityGenes <- (intersect(sig_res$ID,LongevityGenesZenin2019$SNP ))

#Check
#https://genomics.senescence.info/longevity/
#http://www.saspatlas.com


WithLongevityStudies <- list( LongevityGenesZenin2019, eLife ) %>%  #snpsNature
  bind_rows() %>% 
  mutate(Type = "Longevity") %>%
  drop_na(symbol)

WithLongevityStudies$SNP <- gsub("\\*", "",  WithLongevityStudies$SNP)

#create a match with SNPs only 

ExtraTable <- WithLongevityStudies %>%
  inner_join(., sig_res %>% rename(SNP = ID) %>% rename(symbol_Annotation = symbol))

write.table(ExtraTable, file = paste0(IncludeFigHere, "/Longevity_Table.tsv"), col.names = TRUE, row.names = FALSE, sep = ",", quote = TRUE) #copy and paste to jupyter directly.

Models <- list(cellAge1  ,genage_models , genage_human  ,AgingBank  , WithLongevityStudies, ImmAgeOrthologs  ) %>% #ImmAge
  bind_rows() %>%
  mutate(symbol = toupper(symbol)) %>%
  select(-SNP)

#require(fuzzyjoin)
#Merge2 <- regex_inner_join(sig_res, Models,by="symbol", ignore_case =TRUE) #interesting
# strange names!!

Merge <- inner_join(Models,sig_res %>% select(-ID) ) %>%
  rename(Databases = Type) %>%
  unique()

write.table(Merge, file = paste0(IncludeFigHere, "/Longevity_Table2.tsv"), col.names = TRUE, row.names = FALSE, sep = ",", quote = TRUE) #copy and paste to jupyter directly.
Merge <- fread(paste0(IncludeFigHere, "/Longevity_Table2.tsv"))


Small <- Merge %>% 
  select(-c(why, Annotation)) %>% 
  unique() %>%
  mutate(Variable = case_when(Variable == "Unhealthy" ~ "NRN", 
                              Variable == "Healthy" ~ "PRP" ))

#options(bitmapType = "cairo")

pdf(paste0(IncludeFigHere, "/AgeingDatabase8.pdf"), 5, 6)

ggplot(Small, aes( Variable,symbol )) + 
  geom_point(aes(colour = Databases, size = Databases)) + 
  geom_point(shape = 1,aes(size = Databases) ,colour = "black") +
  scale_size_manual(values =c(6, 8, 10, 12, 2, 4)) + #4, 6, 8, 10, 2, 12, 14
  scale_colour_manual(values =c("red2", "lightsalmon", "gold2","dodgerblue","chartreuse4", "plum1")) + #4, 6, 8, 10, 2, 12, 14
  # c("Aging Bank", "Cell Age", "GenAge Human", "GenAge Models", "ImmAge Orthologs", "Longevity" )
  theme_bw() + 
  labs(x = "Categories", y = "Genes")

dev.off()



pdf(paste0(IncludeFigHere, "/AgeingDatabase_Legend.pdf"), 5, 6)

ggplot(Small, aes( Variable,symbol )) + 
  geom_point(aes(colour = Databases)) + 
  scale_colour_manual(values =c("red2", "lightsalmon", "gold2","dodgerblue","chartreuse4", "plum1")) + #4, 6, 8, 10, 2, 12, 14
  geom_point(shape = 1,colour = "black") +
  scale_size_manual(values =c(4, 6, 8, 10, 2)) + #c(2, 4, 6, 8, 10)
  theme_bw() + 
  labs(x = "Categories", y = "Genes")

dev.off()

################################################### ################################################### ################################################### 
################################################### Can run now 4_TableExtraction.R to group generated files ################################################### ################################################### 
################################################### ################################################### ################################################### 

