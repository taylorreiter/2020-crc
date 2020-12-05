library(readr)
library(dplyr)
library(purrr)
library(tidyr)
setwd("github/2020-crc")

# gather metadata ---------------------------------------------------------

# Read in accession numbers/ENA metadata from CRC metagenomes
ena <- list.files("inputs/raw_metadata/", pattern = "_tsv.txt$", full.names = T) %>%
  map_dfr(read_tsv, col_types = "ccccccccccccccccccccccccccccccccccccccccccccccccccc")  %>%
  filter(library_strategy == "WGS")

# read in CRC metadata curated by Wirbel et al. 2019 (doi: 10.1038/s41591-019-0406-6)
crc_metadata <- read_tsv("inputs/raw_metadata/41591_2019_406_MOESM3_ESM_TS6.tsv",
                         skip = 19) %>%
  select(id = Sample_ID, diagnosis = Group, age = Age, gender = Gender, bmi = BMI, 
         country = Country, AJCC_stage, TNM_stage, localization = Localization, 
         sampling_rel_to_colonoscopy= Sampling_rel_to_colonoscopy, FOBT, 
         diabetes = Diabetes, vegetarian = Vegetarian, smoking = Smoking) 

# read in CRC metadata for Japan cohort
japan_metadata <- read_tsv("inputs/raw_metadata/41591_2019_458_MOESM3_ESM_TS2-1.tsv", skip = 2) %>%
  mutate(Subject_ID = paste0("human feces_", Subject_ID)) %>%
  filter(Group != "MP") %>% # filters adenomas
  mutate(BMI = as.numeric(BMI)) %>%
  mutate(country = "JAP") %>%
  mutate(diagnosis = Group) %>%
  mutate(diagnosis = gsub("Healthy", "CTR", diagnosis)) %>%
  mutate(diagnosis = gsub("HS", "CTR", diagnosis)) %>%
  mutate(diagnosis = ifelse(diagnosis == "CTR", "CTR", "CRC")) %>%
  select(id = Subject_ID, diagnosis, age = Age,  gender = Gender, bmi = BMI, country,  
         AJCC_stage = Stage, other_stage = Group, localization = `Tumor location`, 
         brinkman_index = `Brinkman Index`, alcohol = Alcohol)

# get CRC metadata for Italian cohorts
# this removes adenoma...decide later if that is the correct thing to do.
library(curatedMetagenomicData)
italy_metadata <- combined_metadata %>% 
  filter(dataset_name %in% c("ThomasAM_2018a", "ThomasAM_2018b")) %>%
  mutate(gender = gsub("female", "F", gender)) %>%
  mutate(gender = gsub("male", "M", gender)) %>%
  filter(study_condition != "adenoma") %>%
  mutate(diagnosis = gsub("control", "CTR", study_condition)) %>%
  select(id = subjectID, age, diagnosis, age, gender, bmi = BMI, country, alcohol)

nrow(italy_metadata) + nrow(crc_metadata) + nrow(japan_metadata)
just_metadata <- bind_rows(italy_metadata, japan_metadata, crc_metadata)

# combine metadata --------------------------------------------------------

# join cohort metadata with ENA metadata

crc1 <- inner_join(ena, crc_metadata, by = c("sample_accession" = "id"))
crc2 <- inner_join(ena, crc_metadata, ena, by = c("run_accession" = "id"))
crc3 <- inner_join(ena, crc_metadata, ena, by = c("sample_alias" = "id")) 
japan <- inner_join(ena, japan_metadata, by = c("sample_title" = "id"))
italy <- inner_join(ena, italy_metadata, by = c("sample_alias" = "id"))


all_metadata <- bind_rows(crc1, crc2, crc3, japan, italy) %>% 
  filter(library_layout == "PAIRED")

all_metadata %>% select(fastq_ftp)

tmp2 <- tmp %>%
  group_by(sample_alias) %>%
  tally()


# format download links ---------------------------------------------------

# some fastq files have 3 download links, where the first one contains what is
# probably a merged or interleaved file. Remove these from fastq download links. 
# fastq 1 download:

all_metadata <- separate(data = all_metadata, 
                         col = fastq_ftp,
                         into = c("fastq_ftp_1", "fastq_ftp_2", "fastq_ftp_3"), 
                         sep = ";", remove = F, extra = "merge", fill = "right")

fastq_1 <- vector()
for(i in 1:length(all_metadata$fastq_ftp)){
  if(grepl(pattern = "_1", x = all_metadata$fastq_ftp_1[i])){
    fastq_1[i] <- all_metadata$fastq_ftp_1[i] 
  } else {
    fastq_1[i] <- all_metadata$fastq_ftp_2[i]
  }
}

fastq_2 <- vector()
for(i in 1:length(all_metadata$fastq_ftp)){
  if(grepl(pattern = "_1", x = all_metadata$fastq_ftp_1[i])){
    fastq_2[i] <- all_metadata$fastq_ftp_2[i] 
  } else {
    fastq_2[i] <- all_metadata$fastq_ftp_3[i]
  }
}

all_metadata$fastq_ftp_1 <- fastq_1
all_metadata$fastq_ftp_2 <- fastq_2

write_tsv(all_metadata, "inputs/metadata.tsv")
