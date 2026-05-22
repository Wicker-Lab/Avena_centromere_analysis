library(tidyverse)
library(stringr)
library(ggplot2)
library(dplyr)

setwd("")
df_CRMs_Centromere_all <- read.csv("df_all_percent_in_centromere_1.csv")
df_CRMs_Centromere_all$cut_off <- "all"
df_CRMs_Centromere_all$surrounding_10Mb <- df_CRMs_Centromere_all$plus10Mb+df_CRMs_Centromere_all$minus10Mb
df_CRMs_Centromere_all <- df_CRMs_Centromere_all[,-c(4:6)]
df_CRMs_Centromere_all$marker1 <- paste0(df_CRMs_Centromere_all$Chromosome, "_", "full")

df_CRMs_Centromere_20 <- read.csv("df_all_percent_in_centromere_02.csv")
df_CRMs_Centromere_20$cut_off <- "0.2"
df_CRMs_Centromere_20$surrounding_10Mb <- df_CRMs_Centromere_20$plus10Mb+df_CRMs_Centromere_20$minus10Mb
df_CRMs_Centromere_20 <- df_CRMs_Centromere_20[,-c(4:6)]
df_CRMs_Centromere_20$marker1 <- paste0(df_CRMs_Centromere_20$Chromosome, "_", "P20")


df <- rbind(df_CRMs_Centromere_all,df_CRMs_Centromere_20)
df$genome <- as.factor(str_extract(df$Chromosome, "(?<=\\d)\\S"))
df <- df %>% arrange(genome)

df_long <- df %>% pivot_longer(cols = c(centromere, surrounding_10Mb), names_to = "Region", values_to = "values")


category_order <- c("1A_full", "1A_P20", "2A_full", "2A_P20", "3A_full", "3A_P20", "4A_full", "4A_P20",
                    "5A_full", "5A_P20", "6A_full", "6A_P20", "7A_full", "7A_P20",
                    "1C_full", "1C_P20", "2C_full", "2C_P20", "3C_full", "3C_P20", "4C_full", "4C_P20",
                    "5C_full", "5C_P20", "6C_full" , "6C_P20", "7C_full", "7C_P20", 
                    "1D_full","1D_P20", "2D_full", "2D_P20", "3D_full","3D_P20", "4D_full", "4D_P20",
                    "5D_full", "5D_P20", "6D_full", "6D_P20", "7D_full", "7D_P20") 


df_long <- df_long %>% 
  arrange(factor(marker1, levels = category_order))

df_long$marker1 <- factor(df_long$marker1 , levels = unique(df_long$marker1))

ggplot(df_long, aes(x = marker1, y = values*100, fill = Region)) +
  geom_bar(stat = "identity", position = position_stack(reverse = T), width = 0.85) +
  scale_y_continuous(limits = c(0, 100)) +
  ylab("% CRM insertions") +
  xlab("") +
  scale_fill_manual(values = c("green4","#7dff7dff")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 14),
        axis.text.y = element_text(size = 14),
        axis.title.y = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(size = 16, face = "bold"),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15, face = "bold"))
