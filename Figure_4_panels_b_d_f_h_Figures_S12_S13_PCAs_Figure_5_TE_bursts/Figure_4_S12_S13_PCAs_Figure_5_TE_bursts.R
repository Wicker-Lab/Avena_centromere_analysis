#### PCAs for FIGURE 4 (b,d,f,h) and FIGURE 5 S12 and S13 ####
### Set path to where you stored the data as string, make sure to keep folder structure intact.
library(ggplot2)
library(tidyr)
library(stringr)
library(hrbrthemes)
library(dplyr)

path <- "Figure_4_panels_b_d_f_h_Figures_S12_S13_PCAs_Figure_5_TE_bursts/"

#### Ava ####
##### Ava all #####

name <- "Ava"
family <- "RLG_Ava"
family_name <- expression(bold("Ava"))

INPUT_DIR <- paste0(path, "dir_Ava_all")
setwd(INPUT_DIR)

pca <- read.table("pca_tab_RLG_Avenas_Ava_consensus_all")
variance <- read.table("pca_variance_RLG_Avenas_Ava_consensus_all")
age <- read.table("../age_distribution_Avenas_Ava_clean_filtered", skip=1)

age$TEnum <- as.numeric(str_extract(age$V1, "(?<=-)[0-9]+(?=_TSD)"))
age$TE_start <- age$TEnum*1000

tab_with_age <- merge(pca, age, by.x = "sample.id", by.y = "V1", all.x = T)
tab_with_age$genome <- str_extract(pca$sample.id, "(?<=chr\\d)[A-Za-z]|(?<=chr)[A-Za-z]+")

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca_tab_RLG_Avenas_Ava_consensus_G"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Ava_", sub(".*_([A-Za-z0-9]+)$", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$group = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}

Ava_all <-rbind(Ava_G1, Ava_G2, Ava_G3)

Ava_all_PCA <- merge(x = Ava_all, y = tab_with_age, by.x = "TE", by.y = "sample.id", all.y = T)

Ava_all_PCA <- Ava_all_PCA %>% replace_na(list(group = "non-selected"))

ggplot(data = Ava_all_PCA, aes(x = EV1, y = EV2, color = group)) +
  geom_point(aes(), size = 0.1) +
  scale_color_manual(values= c('#332288', '#88CCEE', '#44AA99', "gray89"), name = "Subfamily") +
  xlab(paste("PC 1 (",(round(variance[1,1],digits=1)),"%)")) +
  ylab(paste("PC 2 (",(round(variance[2,1],digits=1)),"%)")) +
  ggtitle("Avena", family) +
  theme_ipsum(base_size = 16, plot_title_size = 18, subtitle_size = 16, strip_text_size = 16, caption_margin = 16, axis_title_size = 16)

ggplot(data = Ava_all_PCA, aes(x = EV1, y = EV2)) +
  geom_point(aes(), color = "gray50" , size = 0.1) +
  xlab(paste("PC 1 (",(round(variance[1,1],digits=1)),"%)")) +
  ylab(paste("PC 2 (",(round(variance[2,1],digits=1)),"%)")) +
  ggtitle("Avena", family) +
  facet_grid(.~V9) +
  theme_ipsum(base_size = 16, plot_title_size = 18, subtitle_size = 16, strip_text_size = 16, caption_margin = 16, axis_title_size = 16)

##### G1 #####
name <- "Ava"
family <- "RLG_Ava"
family_name <- expression(bold("Ava"))

INPUT_DIR <- paste0(path, "dir_Ava_G1/")
setwd(INPUT_DIR)

pca <- read.table("pca_tab_RLG_Avenas_Ava_consensus_G1")
variance <- read.table("pca_variance_RLG_Avenas_Ava_consensus_G1")
age <- read.table("../age_distribution_Avenas_Ava_clean_filtered", skip=1)

age$TEnum <- as.numeric(str_extract(age$V1, "(?<=-)[0-9]+(?=_TSD)"))
age$TE_start <- age$TEnum*1000

tab_with_age <- merge(pca, age, by.x = "sample.id", by.y = "V1", all.x = T)
tab_with_age$genome <- str_extract(pca$sample.id, "(?<=chr\\d)[A-Za-z]|(?<=chr)[A-Za-z]+")

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca_tab_RLG_Avenas_Ava_consensus_G"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Ava_G1_",sub(".*_([0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$group = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}

Ava_G1 <-rbind(Ava_G1_1, Ava_G1_2, Ava_G1_3)

Ava_G1_PCA <- merge(x = Ava_G1, y = tab_with_age, by.x = "TE", by.y = "sample.id", all.y = T)
Ava_G1_PCA <- Ava_G1_PCA %>% replace_na(list(group = "non-selected"))
Ava_G1_PCA <- Ava_G1_PCA %>% arrange(group)

ggplot(data = Ava_G1_PCA, aes(x = EV1, y = EV2, color = group)) +
  geom_point(aes(), size = 0.1) +
  scale_color_manual(values= c("#6EAA22", "#008a6d", "#005e70", "gray89"), name = "Subfamily") +
  xlab(paste("PC 1 (",(round(variance[1,1],digits=1)),"%)")) +
  ylab(paste("PC 2 (",(round(variance[2,1],digits=1)),"%)")) +
  ggtitle("G1", family) +
  theme_ipsum(base_size = 16, plot_title_size = 18, subtitle_size = 16, strip_text_size = 16, caption_margin = 16, axis_title_size = 16)

##### G1_1 #####
name <- "Ava"
family <- "RLG_Ava"
family_name <- expression(bold("Ava"))

INPUT_DIR <- paste0(path, "dir_Ava_G1/dir_Ava_G1_1/")
setwd(INPUT_DIR)

pca <- read.table("pca_tab_RLG_Avenas_Ava_consensus_G1_1")
variance <- read.table("pca_variance_RLG_Avenas_Ava_consensus_G1_1")
age <- read.table("../../age_distribution_Avenas_Ava_clean_filtered", skip=1)

age$TEnum <- as.numeric(str_extract(age$V1, "(?<=-)[0-9]+(?=_TSD)"))
age$TE_start <- age$TEnum*1000

tab_with_age <- merge(pca, age, by.x = "sample.id", by.y = "V1", all.x = T)
tab_with_age$genome <- str_extract(pca$sample.id, "(?<=chr\\d)[A-Za-z]|(?<=chr)[A-Za-z]+")

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca_tab_RLG_Avenas_Ava_consensus_G"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Ava_G1_1_",sub(".*_([0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$group = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}
Ava_G1 <-rbind(Ava_G1_1_1, Ava_G1_1_2)

Ava_G1_PCA <- merge(x = Ava_G1, y = tab_with_age, by.x = "TE", by.y = "sample.id", all.y = T)
Ava_G1_PCA <- Ava_G1_PCA %>% replace_na(list(group = "non-selected"))
Ava_G1_PCA <- Ava_G1_PCA %>% arrange(group)

ggplot(data = Ava_G1_PCA, aes(x = EV1, y = EV2, color = group)) +
  geom_point(aes(), size = 0.1) +
  scale_color_manual(values= c("#DDCC77", "#199f47ff", "gray89"), name = "Subfamily") +
  xlab(paste("PC 1 (",(round(variance[1,1],digits=1)),"%)")) +
  ylab(paste("PC 2 (",(round(variance[2,1],digits=1)),"%)")) +
  ggtitle("G1_1", family) +
  theme_ipsum(base_size = 16, plot_title_size = 18, subtitle_size = 16, strip_text_size = 16, caption_margin = 16, axis_title_size = 16)

##### G2 #####
name <- "Ava"
family <- "RLG_Ava"
family_name <- expression(bold("Ava"))

INPUT_DIR <- paste0(path, "dir_Ava_G2/")
setwd(INPUT_DIR)

pca <- read.table("pca_tab_RLG_Avenas_Ava_consensus_G2")
variance <- read.table("pca_variance_RLG_Avenas_Ava_consensus_G2")
age <- read.table("../age_distribution_Avenas_Ava_clean_filtered", skip=1)

age$TEnum <- as.numeric(str_extract(age$V1, "(?<=-)[0-9]+(?=_TSD)"))
age$TE_start <- age$TEnum*1000

tab_with_age <- merge(pca, age, by.x = "sample.id", by.y = "V1", all.x = T)
tab_with_age$genome <- str_extract(pca$sample.id, "(?<=chr\\d)[A-Za-z]|(?<=chr)[A-Za-z]+")

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca_tab_RLG_Avenas_Ava_consensus_G"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Ava_G2_", sub(".*([0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$group = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}

Ava_G2 <-rbind(Ava_G2_1, Ava_G2_2, Ava_G2_3, Ava_G2_4, Ava_G2_5, Ava_G2_6, Ava_G2_7, Ava_G2_8, Ava_G2_9)

Ava_G2_PCA <- merge(x = Ava_G2, y = tab_with_age, by.x = "TE", by.y = "sample.id", all.y = T)
Ava_G2_PCA <- Ava_G2_PCA %>% replace_na(list(group = "non-selected"))
Ava_G2_PCA <- Ava_G2_PCA %>% arrange(group)

ggplot(data = Ava_G2_PCA, aes(x = EV1, y = EV2,color = group)) +
  geom_point(aes(), size = 0.1) +
  scale_color_manual(values=c( "gray20", "gray50", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#80cdc1", "#35978f", "#01665e", "gray89"), name = "Subfamily") +
  xlab(paste("PC 1 (",(round(variance[1,1],digits=1)),"%)")) +
  ylab(paste("PC 2 (",(round(variance[2,1],digits=1)),"%)")) +
  ggtitle("G2", family) +
  theme_ipsum(base_size = 16, plot_title_size = 18, subtitle_size = 16, strip_text_size = 16, caption_margin = 16, axis_title_size = 16)

##### G3 #####
name <- "Ava"
family <- "RLG_Ava"
family_name <- expression(bold("Ava"))

INPUT_DIR <- paste0(path, "dir_Ava_G3/")
setwd(INPUT_DIR)

pca <- read.table("pca_tab_RLG_Avenas_Ava_consensus_G3")
variance <- read.table("pca_variance_RLG_Avenas_Ava_consensus_G3")
age <- read.table("../age_distribution_Avenas_Ava_clean_filtered", skip=1)

age$TEnum <- as.numeric(str_extract(age$V1, "(?<=-)[0-9]+(?=_TSD)"))
age$TE_start <- age$TEnum*1000

tab_with_age <- merge(pca, age, by.x = "sample.id", by.y = "V1", all.x = T)
tab_with_age$genome <- str_extract(pca$sample.id, "(?<=chr\\d)[A-Za-z]|(?<=chr)[A-Za-z]+")

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca_tab_RLG_Avenas_Ava_consensus_G"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Ava_G3_",sub(".*_([0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$group = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}

Ava_G3 <-rbind(Ava_G3_1, Ava_G3_2, Ava_G3_3, Ava_G3_4, Ava_G3_5, Ava_G3_6)

Ava_G3_PCA <- merge(x = Ava_G3, y = tab_with_age, by.x = "TE", by.y = "sample.id", all.y = T)
Ava_G3_PCA <- Ava_G3_PCA %>% replace_na(list(group = "non-selected"))
Ava_G3_PCA <- Ava_G3_PCA %>% arrange(rev(group))

ggplot(data = Ava_G3_PCA, aes(x = EV1, y = EV2, color = group)) +
  geom_point(aes(), size = 0.1) +
  scale_color_manual(values = c("#AA4499", "#CC6677", "gray30", "#882255", "gray60", "black", "gray89"), name = "Subfamily") +
  xlab(paste("PC 1 (",(round(variance[1,1],digits=1)),"%)")) +
  ylab(paste("PC 2 (",(round(variance[2,1],digits=1)),"%)")) +
  ggtitle("G3", family) +
  theme_ipsum(base_size = 16, plot_title_size = 18, subtitle_size = 16, strip_text_size = 16, caption_margin = 16, axis_title_size = 16)

#### Cereba ####
##### Cereba all #####
name <- "Cereba"
family <- "RLG_Cereba"
family_name <- expression(bold("Cereba"))

INPUT_DIR <- paste0(path, "dir_Cereba_all/")
setwd(INPUT_DIR)

pca <- read.table("pca_tab_Avena_Cereba_consensus-1_all")
variance <- read.table("pca_variance_Avena_Cereba_consensus-1_all")
age <- read.table("../age_distribution_Avenas_Cereba_clean_filtered", skip=1)

age$TEnum <- as.numeric(str_extract(age$V1, "(?<=-)[0-9]+(?=_TSD)"))
age$TE_start <- age$TEnum*1000

tab_with_age <- merge(pca, age, by.x = "sample.id", by.y = "V1", all.x = T)
tab_with_age$genome <- str_extract(pca$sample.id, "(?<=chr\\d)[A-Za-z]|(?<=chr)[A-Za-z]+")

Cereba_younger_half <- read.table("PCA_age_pca_tab_Avena_Cereba_consensus-1_G1")
Cereba_younger_half$TE <- rownames(Cereba_younger_half)
Cereba_younger_half$mark <- "Cereba_G1"
Cereba_full_PCA <- merge(tab_with_age, Cereba_younger_half, by.x = "sample.id", by.y = "TE", all.x = T)

ggplot(data = Cereba_full_PCA, aes(x = EV1.x, y = EV2.x, color = mark)) +
  geom_point(aes(), size = 0.5) +
  scale_color_manual(values = c("seagreen", "gray89"), na.value = "gray89", name = "") +
  xlab(paste("PC 1 (",(round(variance[1,1],digits=1)),"%)")) +
  ylab(paste("PC 2 (",(round(variance[2,1],digits=1)),"%)")) +
  ggtitle("Avena", family) +
  theme_ipsum(base_size = 16, plot_title_size = 18, subtitle_size = 16, strip_text_size = 16, caption_margin = 16, axis_title_size = 16)

ggplot(data = Cereba_full_PCA, aes(x = EV1.x, y = EV2.x)) +
  geom_point(aes(), size = 0.5, color = "gray80") +
  xlab(paste("PC 1 (",(round(variance[1,1],digits=1)),"%)")) +
  ylab(paste("PC 2 (",(round(variance[2,1],digits=1)),"%)")) +
  ggtitle("Avena", family) +
  theme_ipsum(base_size = 16, plot_title_size = 18, subtitle_size = 16, strip_text_size = 16, caption_margin = 16, axis_title_size = 16)+
  facet_grid(.~V9)


##### G1 #####
name <- "Cereba"
family <- "RLG_Cereba"
family_name <- expression(bold("Cereba"))

INPUT_DIR <-paste0(path, "dir_Cereba_G1/")
setwd(INPUT_DIR)

pca <- read.table("pca_tab_Avena_Cereba_consensus-1_G1")
variance <- read.table("pca_variance_Avena_Cereba_consensus-1_G1")
age <- read.table("../age_distribution_Avenas_Cereba_clean_filtered", skip=1)

age$TEnum <- as.numeric(str_extract(age$V1, "(?<=-)[0-9]+(?=_TSD)"))
age$TE_start <- age$TEnum*1000

tab_with_age <- merge(pca, age, by.x = "sample.id", by.y = "V1", all.x = T)
tab_with_age$genome <- str_extract(pca$sample.id, "(?<=chr\\d)[A-Za-z]|(?<=chr)[A-Za-z]+")

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca_tab_Avena_Cereba_consensus-1_G"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Cereba_G1_",sub(".*_([0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$group = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}

Cereba <-rbind(Cereba_G1_1, Cereba_G1_2, Cereba_G1_3, Cereba_G1_4, Cereba_G1_5, Cereba_G1_6, Cereba_G1_7)

Cereba_PCA <- merge(x = Cereba, y = tab_with_age, by.x = "TE", by.y = "sample.id", all.y = T)

Cereba_PCA <- Cereba_PCA %>% replace_na(list(group = "non-selected"))

ggplot(data = Cereba_PCA, aes(x = EV1, y = EV2, color = group)) +
  geom_point(aes(), size = 0.1) +
  scale_color_manual(values=c("#332288", "#88CCEE", "#44AA99", "#117733","#999933","steelblue3", "gray25", "gray89"), name = "Subfamily") +
  xlab(paste("PC 1 (",(round(variance[1,1],digits=1)),"%)")) +
  ylab(paste("PC 2 (",(round(variance[2,1],digits=1)),"%)")) +
  ggtitle("Avena", family) +
  theme_ipsum(base_size = 16, plot_title_size = 18, subtitle_size = 16, strip_text_size = 16, caption_margin = 16, axis_title_size = 16)

#### FIGURE 5 ####
ava <- rbind(Ava_G1_1_1, Ava_G2_1, Ava_G2_2, Ava_G3_1, Ava_G3_2, Ava_G3_4)
cereba <- rbind(Cereba_G1_1, Cereba_G1_2, Cereba_G1_3, Cereba_G1_4, Cereba_G1_5, Cereba_G1_6) 

ava_cereba <- rbind(ava, cereba)
ava_cereba$species_Family <- as.factor(str_extract(ava_cereba$TE, "(?<=RLG_)(?:[^_]*_){1}([^_]+)"))
ava_cereba$species <- as.factor(str_extract(ava_cereba$species_Family, "([A-Za-z]{4})"))

##### filter log #####
ava_cereba_filtered <- ava_cereba[!(ava_cereba$group == "Ava_G3_1" & ava_cereba$species == "Asat" & ava_cereba$genome == "A"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G3_1" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "A",]# 1 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Ava_G3_1" & ava_cereba_filtered$species == "Aatl" & ava_cereba_filtered$genome == "A"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G3_1" & ava_cereba_filtered$species == "Aatl" & ava_cereba_filtered$genome == "A",]# 10 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Ava_G3_4" & ava_cereba_filtered$species == "Aatl" & ava_cereba_filtered$genome == "A"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G3_4" & ava_cereba_filtered$species == "Aatl" & ava_cereba_filtered$genome == "A",]#2 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Ava_G3_4" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "C"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G3_4" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "C",]#1 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Ava_G3_4" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "D"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G3_4" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "D",]#1 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Cereba_G1_2" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "A"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Cereba_G1_2" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "A",]#1 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Cereba_G1_2" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "C"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Cereba_G1_2" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "C",]#1 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Cereba_G1_3" & ava_cereba_filtered$species == "Aatl" & ava_cereba_filtered$genome == "A"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Cereba_G1_3" & ava_cereba_filtered$species == "Aatl" & ava_cereba_filtered$genome == "A",]# 6 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Cereba_G1_3" & ava_cereba_filtered$species == "Alon" & ava_cereba_filtered$genome == "A"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Cereba_G1_3" & ava_cereba_filtered$species == "Alon" & ava_cereba_filtered$genome == "A",]#1 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Cereba_G1_3" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "D"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Cereba_G1_3" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "D",]#1 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Cereba_G1_4" & ava_cereba_filtered$species == "Ains" & ava_cereba_filtered$genome == "D"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Cereba_G1_4" & ava_cereba_filtered$species == "Ains" & ava_cereba_filtered$genome == "D",]#4 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Cereba_G1_4" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "D"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Cereba_G1_4" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "D",]#5 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Cereba_G1_5" & ava_cereba_filtered$species == "Aatl" & ava_cereba_filtered$genome == "A"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Cereba_G1_5" & ava_cereba_filtered$species == "Aatl" & ava_cereba_filtered$genome == "A",]#3 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Cereba_G1_5" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "A"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Cereba_G1_5" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "A",]#2 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Ava_G1_1_1" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "C"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G1_1_1" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "C",]#2 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Ava_G2_1" & ava_cereba_filtered$species == "Ains" & ava_cereba_filtered$genome == "D"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G2_1" & ava_cereba_filtered$species == "Ains" & ava_cereba_filtered$genome == "D",]#2 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Ava_G2_1" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "A"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G2_1" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "A",]#1 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Ava_G2_1" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "D"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G2_1" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "D",]#1 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Ava_G2_2" & ava_cereba_filtered$species == "Ains" & ava_cereba_filtered$genome == "D"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G2_2" & ava_cereba_filtered$species == "Ains" & ava_cereba_filtered$genome == "D",]#3 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Ava_G2_2" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "A"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G2_2" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "A",]#1 gone

ava_cereba_filtered <- ava_cereba_filtered[!(ava_cereba_filtered$group == "Ava_G2_2" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "D"),]
ava_cereba_filtered[ava_cereba_filtered$group == "Ava_G2_2" & ava_cereba_filtered$species == "Asat" & ava_cereba_filtered$genome == "D",]#2 gone


ava_cereba_filtered$species_genome <- as.factor(paste0(ava_cereba_filtered$species, ava_cereba_filtered$genome))

subset_valuable_groups <- ava_cereba_filtered %>%
  filter(group %in%c("Cereba_G1_5", "Cereba_G1_1", "Cereba_G1_3", "Ava_G3_1", "Cereba_G1_2", "Ava_G3_4", "Cereba_G1_4", "Ava_G3_2", "Cereba_G1_6", "Ava_G1_1_1", "Ava_G2_1", "Ava_G2_2"))

ggplot(subset_valuable_groups, aes(y=Age, x=group, fill=group)) +
  geom_violin(color="black", stat = "ydensity", position = "dodge", trim = F)  +
  geom_boxplot(width=0.25, color="#e9ecef", alpha=1) +
  ylab("Age [myr]") +
  labs(fill = "Subfamily") +
  xlab("") +
  theme_minimal() +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10), trans = "reverse") +
  facet_grid(species_genome ~ .) +
  coord_flip() +
  theme(text = element_text(size = 30), axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_manual(values=c("#DDCC77", "gray20", "gray50", "#AA4499", "#CC6677", "#882255", '#332288', '#88CCEE', '#44AA99', '#117733','#999933',"steelblue3")) 

