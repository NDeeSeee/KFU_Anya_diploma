# Attach requirement packages and setting WD -----------------------------------
packages_names <-
  c("tidyverse", "data.table", "readxl", "reshape2", "rstudioapi")

lapply(packages_names, require, character.only = TRUE)

setwd(dirname(getActiveDocumentContext()$path))

rename = dplyr::rename
select = dplyr::select
filter = dplyr::filter
group_by = dplyr::group_by
mutate = dplyr::mutate

FMI <- fread("FMI variants.csv") %>% as_tibble()
frame <- fread("FMI variants.csv") %>% as_tibble()
MAVE_VAMP <- fread("MAVE_VAMP_PTEN_scores.csv") %>% as_tibble()
OncoKB_CKB <-
  fread("OncoKB-CKB PTEN annotation.csv") %>% as_tibble()
FMI_samples <- fread("FMI samples ID's.csv") %>% as_tibble()

MAVE_VAMP <- MAVE_VAMP %>%
  mutate(
    MAVE_status = ifelse(MAVE_score < -1.1, "LoF", "WT"),
    VAMP_status = ifelse(VAMP_score < 0.7, "LoF", "WT"),
    oncoKB = ifelse(
      amino_acid_change %in% filter(OncoKB_CKB, oncoKB == "LoF")$amino_acid_change,
      "LoF",
      NA
    ),
    oncoKB = ifelse(
      amino_acid_change %in% filter(OncoKB_CKB, oncoKB == "WT")$amino_acid_change,
      "WT",
      oncoKB
    ),
    CKB = ifelse(
      amino_acid_change %in% filter(OncoKB_CKB, CKB == "LoF")$amino_acid_change,
      "LoF",
      NA
    )
  )


frame <- frame %>%
  mutate(
    MAVE_status = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(MAVE_VAMP, MAVE_status == "LoF")$amino_acid_change,
      "LoF",
      NA
    ),
    MAVE_status = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(MAVE_VAMP, MAVE_status == "WT")$amino_acid_change,
      "WT",
      MAVE_status
    )
  )

frame <- frame %>%
  mutate(
    VAMP_status = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(MAVE_VAMP, VAMP_status == "LoF")$amino_acid_change,
      "LoF",
      NA
    ),
    VAMP_status = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(MAVE_VAMP, VAMP_status == "WT")$amino_acid_change,
      "WT",
      VAMP_status
    )
  )


frame <- frame %>%
  mutate(
    oncoKB = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(OncoKB_CKB, oncoKB == "LoF")$amino_acid_change,
      "LoF",
      NA
    ),
    oncoKB = ifelse(
      gene == "PTEN" &
        amino_acid_change %in% filter(OncoKB_CKB, oncoKB == "WT")$amino_acid_change,
      "WT",
      oncoKB
    )
  )


frame <- frame %>%
  mutate(CKB = ifelse(
    gene == "PTEN" &
      amino_acid_change %in% filter(OncoKB_CKB, CKB == "LoF")$amino_acid_change,
    "LoF",
    NA
  ),
  CKB)


FMI_samples <- frame %>%
  group_by(sample_id) %>%
  summarise(
    PTEN_presence_VAMP_LoF = sum(gene == "PTEN" &
                                   VAMP_status == "LoF"),
    PTEN_presence_VAMP_WT = sum(gene == "PTEN" &
                                  VAMP_status == "WT"),
    PTEN_presence_MAVE_LoF = sum(gene == "PTEN" &
                                   VAMP_status == "LoF"),
    PTEN_presence_MAVE_WT = sum(gene == "PTEN" &
                                  VAMP_status == "WT"),
    PTEN_presence_OncoKB_LoF = sum(gene == "PTEN" &
                                     oncoKB == "LoF"),
    PTEN_presence_CKB_LoF = sum(gene == "PTEN" &
                                  CKB == "LoF"),
    PTEN_presence_ANY_LoF = sum(gene == "PTEN" &
                                  VAMP_status == "LoF", MAVE_status == "LoF"),
    PTEN_presence_ANY_WT = sum(gene == "PTEN" &
                                 VAMP_status == "WT", MAVE_status == "WT"),
    MAVE_VAMP_LoF = sum(gene == "PTEN" &
                          MAVE_status == "LoF" &
                          VAMP_status == "LoF"),
    MAVE_VAMP_WT = sum(gene == "PTEN" &
                         MAVE_status == "WT" &
                         VAMP_status == "WT"),
    ALL_LoF = sum(
      MAVE_status == "LoF" &
        VAMP_status == "LoF" &
        oncoKB == "LoF" &
        CKB == "LoF"
    ),
    APC_presence = sum(gene == "APC"),
    SMAD4_presence = sum(gene == "SMAD4"),
    TP53_presence = sum(gene == "TP53"),
    KRAS_presence = sum(gene == "KRAS"),
    NRAS_presence = sum(gene == "NRAS")
  )



table_10 <- data_frame(sample_id = FMI_samples$sample_id)

table_10 <- FMI_samples %>%
  mutate(
    PTEN_presence_VAMP_LoF = ifelse(PTEN_presence_VAMP_LoF > 0, 1, 0),
    PTEN_presence_VAMP_WT = ifelse(PTEN_presence_VAMP_WT > 0, 1, 0),
    PTEN_presence_MAVE_LoF = ifelse(PTEN_presence_MAVE_LoF > 0, 1, 0),
    PTEN_presence_MAVE_WT = ifelse(PTEN_presence_MAVE_WT > 0, 1, 0),
    PTEN_presence_OncoKB_LoF = ifelse(PTEN_presence_OncoKB_LoF > 0, 1, 0),
    PTEN_presence_CKB_LoF = ifelse(PTEN_presence_CKB_LoF > 0, 1, 0),
    PTEN_presence_ANY_LoF = ifelse(PTEN_presence_ANY_LoF > 0, 1, 0),
    PTEN_presence_ANY_WT = ifelse(PTEN_presence_ANY_WT > 0, 1, 0),
    MAVE_VAMP_LoF = ifelse(MAVE_VAMP_LoF > 0, 1, 0),
    MAVE_VAMP_WT = ifelse(MAVE_VAMP_WT > 0, 1, 0),
    ALL_LoF = ifelse(ALL_LoF > 0, 1, 0),
    APC_presence = ifelse(APC_presence > 0, 1, 0),
    SMAD4_presence = ifelse(SMAD4_presence > 0, 1, 0),
    TP53_presence = ifelse(TP53_presence > 0, 1, 0),
    KRAS_presence = ifelse(KRAS_presence > 0, 1, 0),
    NRAS_presence = ifelse(NRAS_presence > 0, 1, 0)
  )



fisher.test(x = table_10$PTEN_presence_VAMP_LoF,
            y = table_10$APC_presence)
fisher.test(x = table_10$PTEN_presence_VAMP_WT,
            y = table_10$APC_presence)
fisher.test(x = table_10$PTEN_presence_MAVE_LoF,
            y = table_10$APC_presence)
fisher.test(x = table_10$PTEN_presence_MAVE_WT,
            y = table_10$APC_presence)
fisher.test(x = table_10$PTEN_presence_OncoKB_LoF,
            y = table_10$APC_presence)
fisher.test(x = table_10$PTEN_presence_CKB_LoF,
            y = table_10$APC_presence)
fisher.test(x = table_10$PTEN_presence_ANY_LoF,
            y = table_10$APC_presence) # не работает
fisher.test(x = table_10$PTEN_presence_ANY_WT,
            y = table_10$APC_presence) # не работает

fisher.test(x = table_10$PTEN_presence_VAMP_LoF,
            y = table_10$SMAD4_presence)
fisher.test(x = table_10$PTEN_presence_VAMP_WT,
            y = table_10$SMAD4_presence)
fisher.test(x = table_10$PTEN_presence_MAVE_LoF,
            y = table_10$SMAD4_presence)
fisher.test(x = table_10$PTEN_presence_MAVE_WT,
            y = table_10$SMAD4_presence)
fisher.test(x = table_10$PTEN_presence_OncoKB_LoF,
            y = table_10$SMAD4_presence)
fisher.test(x = table_10$PTEN_presence_CKB_LoF,
            y = table_10$SMAD4_presence)
fisher.test(x = table_10$PTEN_presence_ANY_LoF,
            y = table_10$SMAD4_presence)
fisher.test(x = table_10$PTEN_presence_ANY_WT,
            y = table_10$SMAD4_presence)

fisher.test(x = table_10$PTEN_presence_VAMP_LoF,
            y = table_10$TP53_presence)
fisher.test(x = table_10$PTEN_presence_VAMP_WT,
            y = table_10$TP53_presence)
fisher.test(x = table_10$PTEN_presence_MAVE_LoF,
            y = table_10$TP53_presence)
fisher.test(x = table_10$PTEN_presence_MAVE_WT,
            y = table_10$TP53_presence)
fisher.test(x = table_10$PTEN_presence_OncoKB_LoF,
            y = table_10$TP53_presence)
fisher.test(x = table_10$PTEN_presence_CKB_LoF,
            y = table_10$TP53_presence)
fisher.test(x = table_10$PTEN_presence_ANY_LoF,
            y = table_10$TP53_presence)
fisher.test(x = table_10$PTEN_presence_ANY_WT,
            y = table_10$TP53_presence)


fisher.test(x = table_10$PTEN_presence_VAMP_LoF,
            y = table_10$KRAS_presence)
fisher.test(x = table_10$PTEN_presence_VAMP_WT,
            y = table_10$KRAS_presence)
fisher.test(x = table_10$PTEN_presence_MAVE_LoF,
            y = table_10$KRAS_presence)
fisher.test(x = table_10$PTEN_presence_MAVE_WT,
            y = table_10$KRAS_presence)
fisher.test(x = table_10$PTEN_presence_OncoKB_LoF,
            y = table_10$KRAS_presence)
fisher.test(x = table_10$PTEN_presence_CKB_LoF,
            y = table_10$KRAS_presence)
fisher.test(x = table_10$PTEN_presence_ANY_LoF,
            y = table_10$KRAS_presence)
fisher.test(x = table_10$PTEN_presence_ANY_WT,
            y = table_10$KRAS_presence)


fisher.test(x = table_10$PTEN_presence_VAMP_LoF,
            y = table_10$NRAS_presence)
fisher.test(x = table_10$PTEN_presence_VAMP_WT,
            y = table_10$NRAS_presence)
fisher.test(x = table_10$PTEN_presence_MAVE_LoF,
            y = table_10$NRAS_presence)
fisher.test(x = table_10$PTEN_presence_MAVE_WT,
            y = table_10$NRAS_presence)
fisher.test(x = table_10$PTEN_presence_OncoKB_LoF,
            y = table_10$NRAS_presence)
fisher.test(x = table_10$PTEN_presence_CKB_LoF,
            y = table_10$NRAS_presence)
fisher.test(x = table_10$PTEN_presence_ANY_LoF,
            y = table_10$NRAS_presence)
fisher.test(x = table_10$PTEN_presence_ANY_WT,
            y = table_10$NRAS_presence)
