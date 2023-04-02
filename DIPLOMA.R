install.packages('dplyr')
library(dplyr)
FMI = read.csv("FMI variants.csv")
frame = read.csv('FMI variants.csv')
MAVE_VAMP = read.csv('MAVE_VAMP_PTEN_scores.csv')
OncoKB_CKB = read.csv('OncoKB-CKB PTEN annotation.csv')
FMI_samples = read.csv("FMI samples ID's.csv")

MAVE_VAMP = MAVE_VAMP %>%    
  mutate(MAVE_status = ifelse(MAVE_score < -1.1, "LoF", "WT"), VAMP_status = ifelse(VAMP_score < 0.7, "LoF", "WT"))  
 
frame = frame %>% 
  mutate(MAVE_status = ifelse(gene == "PTEN" & amino_acid_change %in% filter(MAVE_VAMP, MAVE_status == "LoF")$amino_acid_change, "LoF", NA), MAVE_status = ifelse(gene == "PTEN" & amino_acid_change %in% filter(MAVE_VAMP, MAVE_status == "WT")$amino_acid_change, "WT", MAVE_status))

frame = frame %>%    
  mutate(VAMP_status = ifelse(gene == "PTEN" & amino_acid_change %in% filter(MAVE_VAMP, VAMP_status == "LoF")$amino_acid_change, "LoF", NA), VAMP_status = ifelse(gene == "PTEN" & amino_acid_change %in% filter(MAVE_VAMP, VAMP_status == "WT")$amino_acid_change, "WT", VAMP_status))


MAVE_VAMPoncoKB <- OncoKB_CKB %>% 
  mutate(oncoKB = ifelse(oncoKB == 'LoF', 'LoF', 'WT'), CKB = ifelse(CKB == 'LoF', 'LoF', 'WT'))

frame = frame %>% 
  mutate(oncoKB = ifelse(gene == 'PTEN' & amino_acid_change %in% filter(OncoKB_CKB, oncoKB == 'LoF')$amino_acid_change, 'LoF', NA), oncoKB = ifelse(gene == 'PTEN' & amino_acid_change %in% filter(OncoKB_CKB, oncoKB == 'WT')$amino_acid_change, 'WT', oncoKB))


frame = frame %>% 
  mutate(CKB = ifelse(gene == 'PTEN' & amino_acid_change %in% filter(OncoKB_CKB, CKB == 'LoF')$amino_acid_change, 'LoF', NA), CKB)


FMI_samples = frame %>%    
  group_by(sample_id) %>%   summarise(PTEN_presence_VAMP_LoF = sum(gene == "PTEN" & VAMP_status == "LoF"), PTEN_presence_VAMP_WT = sum(gene == "PTEN" & VAMP_status == "WT"), PTEN_presence_MAVE_LoF = sum(gene == "PTEN" & VAMP_status == "LoF"), PTEN_presence_MAVE_WT = sum(gene == "PTEN" & VAMP_status == "WT"), PTEN_presence_OncoKB_LoF = sum(gene == "PTEN" & oncoKB == "LoF"), PTEN_presence_CKB_LoF = sum(gene == "PTEN" & CKB == "WT"), PTEN_presence_ANY_LoF = sum(gene == "PTEN" & VAMP_status == "LoF", MAVE_status == "LoF"), PTEN_presence_ANY_WT = sum(gene == "PTEN" & VAMP_status == "WT", MAVE_status == "WT"), MAVE_VAMP_LoF = sum(gene == 'PTEN' & MAVE_status == 'LoF' & VAMP_status == 'LoF'), MAVE_VAMP_WT = sum(gene == 'PTEN' & MAVE_status == 'WT' & VAMP_status == 'WT'), ALL_LoF = sum(MAVE_status == 'LoF' & VAMP_status == 'LoF' & oncoKB == 'LoF' & CKB == 'LoF'), APC_presence = sum(gene == "APC"), SMAD4_presence = sum(gene == "SMAD4"), TP53_presence = sum(gene == "TP53"), KRAS_presence = sum(gene == "KRAS"), NRAS_presence = sum(gene == "NRAS"))
 

#             29 строчка не считает      ANY_presence   если ставишь условие с    oncoKB &   CKB LoF ибо они не перенеслись


table_10 <- data_frame(sample_id = FMI_samples$sample_id)

#доделать
table_10 = FMI_samples %>% 
  mutate(PTEN_presence_VAMP_LoF = ifelse(PTEN_presence_VAMP_LoF > 0, 0, 1), PTEN_presence_VAMP_WT = ifelse(PTEN_presence_VAMP_WT > 0, 0, 1), PTEN_presence_MAVE_LoF = ifelse(PTEN_presence_MAVE_LoF > 0, 0, 1), PTEN_presence_MAVE_WT = ifelse(PTEN_presence_MAVE_WT > 0, 0, 1), PTEN_presence_OncoKB_LoF = ifelse(PTEN_presence_OncoKB_LoF > 0, 0, 1), PTEN_presence_CKB_LoF = ifelse(PTEN_presence_CKB_LoF > 0, 0, 1), PTEN_presence_ANY_LoF = ifelse(PTEN_presence_ANY_LoF > 0, 0, 1), PTEN_presence_ANY_WT = ifelse(PTEN_presence_ANY_WT > 0, 0, 1), APC_presence = ifelse(APC_presence  > 0, 0, 1), SMAD4_presence = ifelse(SMAD4_presence > 0, 0, 1), TP53_presence = ifelse(TP53_presence > 0, 0, 1), KRAS_presence = ifelse(KRAS_presence > 0, 0, 1), NRAS_presence = ifelse(NRAS_presence > 0, 0, 1))


table(table_10$sample_id, table_10$PTEN_presence_VAMP_LoF, table_10$APC_presence)
table(table_10$sample_id, table_10$PTEN_presence_VAMP_WT, table_10$APC_presence)
table(table_10$sample_id, table_10$PTEN_presence_MAVE_LoF, table_10$APC_presence)
table(table_10$sample_id, table_10$PTEN_presence_MAVE_WT, table_10$APC_presence)
table(table_10$sample_id, table_10$PTEN_presence_OncoKB_LoF, table_10$APC_presence)
table(table_10$sample_id, table_10$PTEN_presence_CKB_LoF, table_10$APC_presence)
table(table_10$sample_id, table_10$PTEN_presence_ANY_LoF, table_10$APC_presence)
table(table_10$sample_id, table_10$PTEN_presence_ANY_WT, table_10$APC_presence)

table(table_10$sample_id, table_10$PTEN_presence_VAMP_LoF, table_10$SMAD4_presence)
table(table_10$sample_id, table_10$PTEN_presence_VAMP_LoF, table_10$SMAD4_presence)
table(table_10$sample_id, table_10$PTEN_presence_VAMP_WT, table_10$SMAD4_presence)
table(table_10$sample_id, table_10$PTEN_presence_MAVE_LoF, table_10$SMAD4_presence)
table(table_10$sample_id, table_10$PTEN_presence_MAVE_WT, table_10$SMAD4_presence)
table(table_10$sample_id, table_10$PTEN_presence_OncoKB_LoF, table_10$SMAD4_presence)
table(table_10$sample_id, table_10$PTEN_presence_CKB_LoF, table_10$SMAD4_presence)
table(table_10$sample_id, table_10$PTEN_presence_ANY_LoF, table_10$SMAD4_presence)
table(table_10$sample_id, table_10$PTEN_presence_ANY_WT, table_10$SMAD4_presence)

table(table_10$sample_id, table_10$PTEN_presence_VAMP_LoF, table_10$TP53_presence)
table(table_10$sample_id, table_10$PTEN_presence_VAMP_LoF, table_10$TP53_presence)
table(table_10$sample_id, table_10$PTEN_presence_VAMP_WT, table_10$TP53_presence)
table(table_10$sample_id, table_10$PTEN_presence_MAVE_LoF, table_10$TP53_presence)
table(table_10$sample_id, table_10$PTEN_presence_MAVE_WT, table_10$TP53_presence)
table(table_10$sample_id, table_10$PTEN_presence_OncoKB_LoF, table_10$TP53_presence)
table(table_10$sample_id, table_10$PTEN_presence_CKB_LoF, table_10$TP53_presence)
table(table_10$sample_id, table_10$PTEN_presence_ANY_LoF, table_10$TP53_presence)
table(table_10$sample_id, table_10$PTEN_presence_ANY_WT, table_10$TP53_presence)

table(table_10$sample_id, table_10$PTEN_presence_VAMP_LoF, table_10$KRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_VAMP_LoF, table_10$KRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_VAMP_WT, table_10$KRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_MAVE_LoF, table_10$KRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_MAVE_WT, table_10$KRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_OncoKB_LoF, table_10$KRAS_presence)        #только один столбец??? во всех oncoKB & CKB
table(table_10$sample_id, table_10$PTEN_presence_CKB_LoF, table_10$KRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_ANY_LoF, table_10$KRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_ANY_WT, table_10$KRAS_presence)

table(table_10$sample_id, table_10$PTEN_presence_VAMP_LoF, table_10$NRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_VAMP_LoF, table_10$NRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_VAMP_WT, table_10$NRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_MAVE_LoF, table_10$NRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_MAVE_WT, table_10$NRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_OncoKB_LoF, table_10$NRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_CKB_LoF, table_10$NRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_ANY_LoF, table_10$NRAS_presence)
table(table_10$sample_id, table_10$PTEN_presence_ANY_WT, table_10$NRAS_presence)
