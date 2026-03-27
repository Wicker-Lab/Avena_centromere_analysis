#### FIGURE 3 and FIGURE S11 and S11 ####
### Set path to where you stored the data as string, make sure to keep folder structure intact.
library(stringr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(ggtext)

path <- "Figure_3_Figures_S1_S11_ageXposition"
setwd(path)

#### FIGURE S1 ####
##### Preparation #####
age_tab_Ava <- read.table("age_distribution_Avenas_Ava_clean_filtered", skip=1)
age_tab_Ava$family <- str_extract(age_tab_Ava$V1, "\\S+(?<=_Ava)")
age_tab_Ava$family_short <- "RLG_Ava"
age_tab_Ava$TEnum <- str_extract(age_tab_Ava$V1, "\\d\\S-\\d+")
age_tab_Ava$genome_w_acc <- str_extract(age_tab_Ava$family, "(?<=RLG_)(.{4})")
age_tab_Ava$chr <- str_extract(age_tab_Ava$TEnum, "\\d\\S")
age_tab_Ava$TEstart <- as.numeric(str_extract(age_tab_Ava$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the approximate start position from the ID 
age_tab_Ava$genome <- as.factor(str_extract(age_tab_Ava$chr, "(?<=\\d)\\S"))

age_tab_Ava$TE_ID <- paste(age_tab_Ava$genome_w_acc, age_tab_Ava$TEnum, sep = "_")

df_sorted_Ava <- age_tab_Ava %>% arrange(genome)
df_sorted_Ava$chr <- factor(df_sorted_Ava$chr, levels = rev(unique(df_sorted_Ava$chr)))

age_tab_Cereba <- read.table("age_distribution_Avenas_Cereba_clean_filtered", skip=1)
age_tab_Cereba$family <- str_extract(age_tab_Cereba$V1, "\\S+(?<=_Cereba)")
age_tab_Cereba$family_short <- "RLG_Cereba"
age_tab_Cereba$TEnum <- str_extract(age_tab_Cereba$V1, "\\d\\S-\\d+")
age_tab_Cereba$genome_w_acc <- str_extract(age_tab_Cereba$family, "(?<=RLG_)(.{4})")
age_tab_Cereba$chr <- str_extract(age_tab_Cereba$TEnum, "\\d\\S")
age_tab_Cereba$TEstart <- as.numeric(str_extract(age_tab_Cereba$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the approximate start position from the ID 
age_tab_Cereba$genome <- as.factor(str_extract(age_tab_Cereba$chr, "(?<=\\d)\\S"))

age_tab_Cereba$TE_ID <- paste(age_tab_Cereba$genome_w_acc, age_tab_Cereba$TEnum, sep = "_")

df_sorted_Cereba <- age_tab_Cereba %>% arrange(genome)
df_sorted_Cereba$chr <- factor(df_sorted_Cereba$chr, levels = rev(unique(df_sorted_Cereba$chr)))

age_tab_Beth <- read.table("age_distribution_Avenas_Beth_clean_filtered", skip=1)
age_tab_Beth$family <- str_extract(age_tab_Beth$V1, "\\S+(?<=_Beth)")
age_tab_Beth$family_short <- "RLG_Beth"
age_tab_Beth$TEnum <- str_extract(age_tab_Beth$V1, "\\d\\S-\\d+")
age_tab_Beth$genome_w_acc <- str_extract(age_tab_Beth$family, "(?<=RLX_)(.{4})")
age_tab_Beth$chr <- str_extract(age_tab_Beth$TEnum, "\\d\\S")
age_tab_Beth$TEstart <- as.numeric(str_extract(age_tab_Beth$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the approximate start position from the ID 
age_tab_Beth$genome <- as.factor(str_extract(age_tab_Beth$chr, "(?<=\\d)\\S"))

age_tab_Beth$TE_ID <- paste(age_tab_Beth$genome_w_acc, age_tab_Beth$TEnum, sep = "_")

df_sorted_Beth <- age_tab_Beth %>% arrange(genome)
df_sorted_Beth$chr <- factor(df_sorted_Beth$chr, levels = rev(unique(df_sorted_Beth$chr)))

df_sorted <- rbind(df_sorted_Ava, df_sorted_Beth, df_sorted_Cereba)

##### import sizes of the chromosomes #####
sizeList <- read.table("Asat_OT3098_v2_genome_size_list_chr", header = TRUE)
sizeList <- sizeList %>% arrange(Name)
sizeList$chr <- str_extract(sizeList$Name, "(?<=chr)\\d\\S")
sizeList$chromosome <- str_extract(sizeList$chr, "\\d(?<=\\S)")
sizeList$genome <- as.factor(str_extract(sizeList$chr, "(?<=\\d)\\S"))
sizeList <- sizeList %>% arrange(genome)

sizeList <- sizeList[sizeList$chr %in% unique(df_sorted$chr),]
sizeList <- mutate(sizeList, chr_nr=1:length(sizeList$Name))
sizeList <- sizeList %>% arrange(desc(chr_nr))
sizeList <- mutate(sizeList, chr_nr_sorted=1:length(sizeList$Name))
df_sorted_size <- merge(df_sorted, sizeList, by = "chr")

df_sorted_size <- df_sorted_size %>% arrange(chr_nr_sorted)
df_sorted_size$chr <- factor(df_sorted_size$chr, levels = unique(df_sorted_size$chr))

OT3098_centromere_pos <- read.table("ChiP_inferred_centromere_positions_OT3098", header = T)
OT3098_centromere_pos$chr <- str_extract(OT3098_centromere_pos$OT3098, "(?<=chr)\\d\\S")
df_sorted_size_centromere_pos <-merge(df_sorted_size, OT3098_centromere_pos, by = "chr", all.x = T)


data <- df_sorted_size_centromere_pos[df_sorted_size_centromere_pos$genome_w_acc == "Asat",]
data$V6 <- as.numeric(data$V6)

#### ageXpos plots with ChIP seq data ####
##### gap positions #####
list_of_files <- list.files(path = ".", pattern = paste0("annotation_N_Asat_OT3098_v2"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- sub(".*_(chr[0-9][A-Za-z]).*", "\\1", file)
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  df$group = df_name
  df <- df %>% separate(V1, into = c("Start", "End"), sep = "-")
  df <- df %>% separate(V2, into = c("Name", "Orientation", "Type", "Length", "other"), sep = ";")
  df <- df %>% separate(Length, into = c("delete", "Length"), sep = "=")
  df <- df %>% select(-delete, -other)
  assign(df_name, df, envir = .GlobalEnv)
}

gaps_OT3098 <- rbind(chr1A, chr1C, chr1D, chr2A, chr2C, chr2D, chr3A, chr3C, chr3D, chr4A, chr4C, chr4D, chr5A, chr5C, chr5D, chr6A, chr6C, chr6D, chr7A, chr7C, chr7D)

gaps_OT3098$genome <- str_extract(gaps_OT3098$group, "(?<=chr[0-9])\\S+")
gaps_OT3098$chr <- str_extract(gaps_OT3098$group, "(?<=chr)\\d+\\S+")

gaps_OT3098 <- gaps_OT3098 %>% arrange(group)
gaps_OT3098$group <- factor(gaps_OT3098$group, levels = unique(gaps_OT3098$group))

gaps_OT3098$start <- round(as.numeric(gaps_OT3098$Start)/1000000, 0)
gaps_OT3098$end <- round(as.numeric(gaps_OT3098$End)/1000000, 0)

##### ChIP seq data #####
ChIP_seq_Avena <- read.table("CRR515328_CRR515329_500kb_MQ30.bedgraph")

#### Read in centromere positions ####
OT3098_centromere_pos <- read.table("ChiP_inferred_centromere_positions_OT3098", header = T)
OT3098_centromere_pos$chr <- str_extract(OT3098_centromere_pos$OT3098, "(?<=chr)\\d\\S")
OT3098_centromere_pos$Start <- round(OT3098_centromere_pos$start_full/1000000,2)
OT3098_centromere_pos$End <- round(OT3098_centromere_pos$end_full/1000000,2)

#----CHANGE:
#either
list_of_chr <- c("1A", "2A", "3A", "4A", "5A", "6A", "7A", "1C", "2C", "4C", "5C", "6C", "7C", "2D", "4D", "5D", "6D", "7D")
#or
list_of_chr <- c("1D", "3D", "3C") # run individually with inactivated Line 156 because no gaps in the centromeric region! 

for (chromosome_of_choice_1 in list_of_chr) {
  
  # subset for chromosomes
  start_cent <-  OT3098_centromere_pos[OT3098_centromere_pos$OT3098 == paste0("chr", chromosome_of_choice_1),]$start
  end_cent <-  OT3098_centromere_pos[OT3098_centromere_pos$OT3098 == paste0("chr", chromosome_of_choice_1),]$end
  centromere_mid <- (start_cent+end_cent)/2
  
  chr_choice <- data[data$chr == chromosome_of_choice_1,]
  chr_choice <- chr_choice[(chr_choice$TEstart)>centromere_mid-40000000 & (chr_choice$TEstart)<centromere_mid+40000000,]
  
  ChIP_seq_Avena_choice <- ChIP_seq_Avena[ChIP_seq_Avena$V1 == paste0("Asat_OT3098_v2_chr", chromosome_of_choice_1),]
  ChIP_seq_Avena_choice <- ChIP_seq_Avena_choice[(ChIP_seq_Avena_choice$V2)>centromere_mid-45000000 & (ChIP_seq_Avena_choice$V2)<centromere_mid+45000000,]
  
  # read in centromere positions
  start_cent_Mb <- OT3098_centromere_pos[OT3098_centromere_pos$chr == chromosome_of_choice_1,]$start/1000000
  end_cent_Mb <-  OT3098_centromere_pos[OT3098_centromere_pos$chr == chromosome_of_choice_1,]$end/1000000
  centromere_mid_Mb <- round(((start_cent+end_cent)/2)/1000000,2)
  
  # get gap data
  gaps <- gaps_OT3098[gaps_OT3098$chr == chromosome_of_choice_1 & gaps_OT3098$start > start_cent_Mb-40 & gaps_OT3098$end < end_cent_Mb+35,]
  
  p1 <- ggplot(data = chr_choice, aes(x = TEstart/1000000, y = V6, color = family)) +
    geom_point(shape=20, size=1)  +
    ylab("Insertion age [Myr]") +
    xlab("Position [Mb]") +
    labs(color = "Family") +
    scale_color_manual(values=c("#322480", "#5aae61", "#d058b5ff")) +
    scale_x_continuous(limits = c(centromere_mid_Mb-40, centromere_mid_Mb+40), breaks = scales::pretty_breaks(n = 10)) +
    scale_y_continuous(limits = c(0, 4)) +
    theme_classic(base_size = 16) +
    lapply(1:nrow(gaps), function(i) {annotation_custom(grob = grid::rectGrob(gp = grid::gpar(col = "red", fill = NA, lwd = 1.2)), xmin = gaps$start[i], xmax = gaps$end[i], ymin = -0.28, ymax = -0.2)}) +#do not use if no gaps in the displayed area
    coord_cartesian(clip = "off")
  p1
  
  
  
  p2 <- ggplot(data = ChIP_seq_Avena_choice, aes(x=V2/1000000, y=V4)) +
    geom_line(color = "black",  linewidth = 0.3) +
    theme_classic(base_size = 16) +
    labs(tag = chromosome_of_choice_1) +
    ylab("log2\nCENH3/Input") +
    xlab("") +
    coord_cartesian(xlim = c(start_cent_Mb-40, end_cent_Mb+40),  ylim = c(-2, 6), clip = "off") +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    annotate("rect", xmin = start_cent_Mb, xmax = end_cent_Mb, ymin= -44.5, ymax=7.5, alpha=0.3)
  
  
  p2
  
  layout <- "
B
A
A
A
A
"
  
  plot <- p1+p2+plot_layout(design = layout)
  plot
  
  assign(paste0("p_", chromosome_of_choice_1), plot, envir = .GlobalEnv)
  
  print(plot)

}

(p_1A | p_1C | p_1D)/
  (p_2A | p_2C | p_2D)/
  (p_3A | p_3C | p_3D) + plot_layout(widths = c(1.1,1.1,1.1), heights = c(2.1,2.1,2.1), guides = 'collect') &
  theme(
    panel.spacing.x = unit(5.2, "cm"),
    panel.spacing.y = unit(1.2, "cm"),
    legend.position = "bottom" 
  )

(p_4A | p_4C | p_4D)/  
  (p_5A | p_5C | p_5D)/
  (p_6A | p_6C | p_6D)/
  (p_7A | p_7C | p_7D) + plot_layout(widths = c(1.1,1.1,1.1,1.1), heights = c(2.1,2.1,2.1,2.1), guides = 'collect') &
  theme(
    panel.spacing.x = unit(5.2, "cm"), 
    panel.spacing.y = unit(1.2, "cm"), 
    legend.position = "bottom" 
  )

#### FIGURE 3 ####

custom_theme2 <- theme_classic(base_size = 13) +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(color = "black"),
    axis.line = element_line(color = "black"),
    strip.text = element_text(size = 16),
    panel.grid.major = element_line(color = "lightgray"),
    panel.grid.minor = element_line(color = NA),
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA)
  )

#### Preparation ####

#### Asat OT3098 ####
Asat_Cereba <- read.table("age_distribution_Avenas_Cereba_clean_filtered", skip=1)
Asat_Cereba$family <- str_extract(Asat_Cereba$V1, "\\S+(?<=_Cereba)")

Asat_Ava <-  read.table("age_distribution_Avenas_Ava_clean_filtered", skip=1)
Asat_Ava$family <- str_extract(Asat_Ava$V1, "\\S+(?<=_Ava)")

Asat_Beth <- read.table("age_distribution_Avenas_Beth_clean_filtered", skip=1)
Asat_Beth$family <- str_extract(Asat_Beth$V1, "\\S+(?<=_Beth)")

all_OT3098_AsatOT <- rbind(Asat_Ava, Asat_Beth, Asat_Cereba)

all_OT3098_AsatOT$TEnum <- str_extract(all_OT3098_AsatOT$V1, "\\d\\S-\\d+")
all_OT3098_AsatOT$chr <- str_extract(all_OT3098_AsatOT$TEnum, "\\d\\S")
all_OT3098_AsatOT$TEstart <- as.numeric(str_extract(all_OT3098_AsatOT$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the approximate start position from the ID 
all_OT3098_AsatOT$genome <- as.factor(str_extract(all_OT3098_AsatOT$chr, "(?<=\\d)\\S"))

df_sorted <- all_OT3098_AsatOT %>% arrange(genome)
df_sorted$chr <- factor(df_sorted$chr, levels = rev(unique(df_sorted$chr)))

OT3098_centromere_pos <- read.table("ChiP_inferred_centromere_positions_OT3098", header = T)
OT3098_centromere_pos$chr <- str_extract(OT3098_centromere_pos$OT3098, "(?<=chr)\\d\\S")
OT3098_centromere_pos$Start <- round(OT3098_centromere_pos$start_full/1000000,2)
OT3098_centromere_pos$End <- round(OT3098_centromere_pos$end_full/1000000,2)

##### FIGURE 3 Panel a #####
### gaps ###
list_of_files <- list.files(path = ".", pattern = paste0("annotation_N_Asat_OT3098_v2"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- sub(".*_(chr[0-9][A-Za-z]).*", "\\1", file)
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  df$group = df_name
  df <- df %>% separate(V1, into = c("Start", "End"), sep = "-")
  df <- df %>% separate(V2, into = c("Name", "Orientation", "Type", "Length", "other"), sep = ";")
  df <- df %>% separate(Length, into = c("delete", "Length"), sep = "=")
  df <- df %>% select(-delete, -other)
  assign(df_name, df, envir = .GlobalEnv)
  #df_names <- c(df_names, df_name)
}

gaps_OT3098 <- rbind(chr1A, chr1C, chr1D, chr2A, chr2C, chr2D, chr3A, chr3C, chr3D, chr4A, chr4C, chr4D, chr5A, chr5C, chr5D, chr6A, chr6C, chr6D, chr7A, chr7C, chr7D)

gaps_OT3098$genome <- str_extract(gaps_OT3098$group, "(?<=chr[0-9])\\S+")
gaps_OT3098$chr <- str_extract(gaps_OT3098$group, "(?<=chr)\\d+\\S+")

gaps_OT3098 <- gaps_OT3098 %>% arrange(group)
gaps_OT3098$group <- factor(gaps_OT3098$group, levels = unique(gaps_OT3098$group))

gaps_OT3098$start <- round(as.numeric(gaps_OT3098$Start)/1000000, 2)
gaps_OT3098$end <- round(as.numeric(gaps_OT3098$End)/1000000, 2)

gaps_OT3098 <- gaps_OT3098[gaps_OT3098$start>180 & gaps_OT3098$end<260 & gaps_OT3098$group == "chr5A",]

chr5A_OT3098 <- df_sorted[df_sorted$chr == "5A" & df_sorted$V9 == "AsatOT3098",]
chr5A_OT3098 <- chr5A_OT3098[chr5A_OT3098$TEstart/1000000>180 & chr5A_OT3098$TEstart/1000000<260, ]

Start_cen_5A <- OT3098_centromere_pos[OT3098_centromere_pos$OT3098 == "chr5A",]$Start
End_cen_5A <- OT3098_centromere_pos[OT3098_centromere_pos$OT3098 == "chr5A",]$End

Panel_a <- ggplot(data = chr5A_OT3098, aes(x = TEstart/1000000, y = V6, color = family)) +
  geom_point(shape=20, size=0.8)  +
  ylab("Insertion age [Myr]") +
  xlab("Position on chromosome 5A [Mb]") +
  labs(color = "Family") +
  scale_color_manual(values=c("#322480","#5aae61", "#d058b5ff"), labels = c("*RLG_Ava*","*RLG_Cereba*", "*RLG_Beth*")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  custom_theme2 +
  lapply(1:nrow(gaps_OT3098), function(i) {annotation_custom(grob = grid::rectGrob(gp = grid::gpar(col = "red", fill = NA, lwd = 1.2)), xmin = gaps_OT3098$start[i], xmax = gaps_OT3098$end[i], ymin = -0.20, ymax = -0.15)}) +
  coord_cartesian(xlim = c(180,260), ylim = c(0,3), clip = "off") +
  geom_segment(x = Start_cen_5A, xend = End_cen_5A, y = -0.12, yend = -0.12, color = "black", linewidth = 0.8)+
  theme(legend.text = element_markdown())
Panel_a

##### FIGURE 3 Panel b-d #####
chr5D_OT3098 <- df_sorted[df_sorted$chr == "5D" & df_sorted$V9 == "AsatOT3098",]

### gaps ###
list_of_files <- list.files(path = ".", pattern = paste0("annotation_N_Asat_OT3098_v2"),  full.names = FALSE) 
df_names <- list()
for (file in list_of_files) {
  df_name <- sub(".*_(chr[0-9][A-Za-z]).*", "\\1", file)
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  df$group = df_name
  df <- df %>% separate(V1, into = c("Start", "End"), sep = "-")
  df <- df %>% separate(V2, into = c("Name", "Orientation", "Type", "Length", "other"), sep = ";")
  df <- df %>% separate(Length, into = c("delete", "Length"), sep = "=")
  df <- df %>% select(-delete, -other)
  assign(df_name, df, envir = .GlobalEnv)
}

gaps_OT3098 <- rbind(chr1A, chr1C, chr1D, chr2A, chr2C, chr2D, chr3A, chr3C, chr3D, chr4A, chr4C, chr4D, chr5A, chr5C, chr5D, chr6A, chr6C, chr6D, chr7A, chr7C, chr7D)

gaps_OT3098$genome <- str_extract(gaps_OT3098$group, "(?<=chr[0-9])\\S+")
gaps_OT3098$chr <- str_extract(gaps_OT3098$group, "(?<=chr)\\d+\\S+")

gaps_OT3098 <- gaps_OT3098 %>% arrange(group)
gaps_OT3098$group <- factor(gaps_OT3098$group, levels = unique(gaps_OT3098$group))

gaps_OT3098$start <- round(as.numeric(gaps_OT3098$Start)/1000000, 2)
gaps_OT3098$end <- round(as.numeric(gaps_OT3098$End)/1000000, 2)

gaps_OT3098 <- gaps_OT3098[gaps_OT3098$start>140 & gaps_OT3098$end<240 & gaps_OT3098$group == "chr5D",]

###### FIGURE 3 Panel b ######
OT3098 <- expression(paste(italic("A. sativa"), " OT3098 Position [Mb]"))

chr5D_OT3098 <- chr5D_OT3098[chr5D_OT3098$TEstart/1000000>140 & chr5D_OT3098$TEstart/1000000<240,]

Start_cen_5D <- OT3098_centromere_pos[OT3098_centromere_pos$OT3098 == "chr5D",]$Start
End_cen_5D <- OT3098_centromere_pos[OT3098_centromere_pos$OT3098 == "chr5D",]$End

p_OT3098 <- ggplot(data = chr5D_OT3098, aes(x = TEstart/1000000, y = V6, color = family)) +
  geom_point(shape=20, size=0.8)  +
  ylab("Insertion age [Myr]") +
  xlab(OT3098) +
  labs(color = "Family") +
  scale_color_manual(values=c("#322480","#5aae61", "#d058b5ff"), labels = c("*RLG_Ava*","*RLG_Cereba*", "*RLG_Beth*")) +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10), position = "bottom") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6), position = "left")+
  custom_theme2 +
  lapply(1:nrow(gaps_OT3098), function(i) {annotation_custom(grob = grid::rectGrob(gp = grid::gpar(col = "red", fill = NA, lwd = 1.2)), xmin = gaps_OT3098$start[i], xmax = gaps_OT3098$end[i], ymin = -0.20, ymax = -0.15)}) +
  coord_cartesian(xlim = c(140,240), ylim = c(0,3), clip = "off") +
  annotate("rect", xmin = 192.75, xmax = 203.5, ymin = -4.0, ymax=3, alpha=0.2) +
  geom_segment(x = Start_cen_5D, xend = End_cen_5D, y = -0.12, yend = -0.12, color = "black", linewidth = 0.8)+
  theme(legend.text = element_markdown())
p_OT3098

#### AsatC0648 ####
AsatC0648_Cereba <- read.table("age_tab_Cereba_AsatC0648_clean_filtered", skip=1)
AsatC0648_Cereba$family <- str_extract(AsatC0648_Cereba$V1, "\\S+(?<=_Cereba)")

AsatC0648_Ava <- read.table("age_tab_Ava_AsatC0648_clean_filtered", skip=1)
AsatC0648_Ava$family <- str_extract(AsatC0648_Ava$V1, "\\S+(?<=_Ava)")

AsatC0648_Beth <- read.table("age_tab_Beth_AsatC0648_clean_filtered", skip=1)
AsatC0648_Beth$family <- str_extract(AsatC0648_Beth$V1, "\\S+(?<=_Beth)")

all_AsatC0648 <- rbind(AsatC0648_Beth, AsatC0648_Cereba, AsatC0648_Ava)

all_AsatC0648$TEnum <- str_extract(all_AsatC0648$V1, "\\d\\S-\\d+")
all_AsatC0648$chr <- str_extract(all_AsatC0648$TEnum, "\\d\\S")
all_AsatC0648$TEstart <- as.numeric(str_extract(all_AsatC0648$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the approximate start position from the ID 
all_AsatC0648$genome <- as.factor(str_extract(all_AsatC0648$chr, "(?<=\\d)\\S"))

df_sorted <- all_AsatC0648 %>% arrange(genome)
df_sorted$chr <- factor(df_sorted$chr, levels = rev(unique(df_sorted$chr)))


chr5D_AsatC0648 <- df_sorted[df_sorted$chr == "5D",]

##### gaps #####
list_of_files <- list.files(path = ".", pattern = paste0("annotation_N_Asat_C0648_chr"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- sub(".*_(chr[0-9][A-Za-z]).*", "\\1", file)
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  df$group = df_name
  df <- df %>% separate(V1, into = c("Start", "End"), sep = "-")
  df <- df %>% separate(V2, into = c("Name", "Orientation", "Type", "Length", "other"), sep = ";")
  df <- df %>% separate(Length, into = c("delete", "Length"), sep = "=")
  df <- df %>% select(-delete, -other)
  assign(df_name, df, envir = .GlobalEnv)
  #df_names <- c(df_names, df_name)
}

gaps_AsatC0648 <- rbind(chr1A, chr1C, chr1D, chr2A, chr2C, chr2D, chr3A, chr3C, chr3D, chr4A, chr4C, chr4D, chr5A, chr5C, chr5D, chr6A, chr6C, chr6D, chr7A, chr7C, chr7D)

gaps_AsatC0648$genome <- str_extract(gaps_AsatC0648$group, "(?<=chr[0-9])\\S+")
gaps_AsatC0648$chr <- str_extract(gaps_AsatC0648$group, "(?<=chr)\\d+\\S+")

gaps_AsatC0648 <- gaps_AsatC0648 %>% arrange(group)
gaps_AsatC0648$group <- factor(gaps_AsatC0648$group, levels = unique(gaps_AsatC0648$group))

gaps_AsatC0648$start <- round(as.numeric(gaps_AsatC0648$Start)/1000000, 2)
gaps_AsatC0648$end <- round(as.numeric(gaps_AsatC0648$End)/1000000, 2)

gaps_AsatC0648 <- gaps_AsatC0648[gaps_AsatC0648$start>140 & gaps_AsatC0648$end<240 & gaps_AsatC0648$group == "chr5D",]

###### FIGURE 3 Panel c ###### 
C0648 <- expression(paste(italic("A. sativa"), " C0648 Position [Mb]"))

chr5D_AsatC0648 <- chr5D_AsatC0648[chr5D_AsatC0648$TEstart/1000000>140 & chr5D_AsatC0648$TEstart/1000000<240,]
p_C0648 <- ggplot(data = chr5D_AsatC0648, aes(y = TEstart/1000000, x = V6, color = family)) +
  geom_point(shape=20, size=0.8)  +
  xlab("Insertion age [Myr]") +
  ylab(C0648) +
  labs(color = "Family") +
  scale_color_manual(values=c("#322480","#d058b5ff", "#5aae61"), labels = c("*RLG_Ava*", "*RLG_Beth*" ,"*RLG_Cereba*")) +
  scale_x_reverse(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_reverse(breaks = scales::pretty_breaks(n = 10))+
  custom_theme2 + 
  coord_cartesian(xlim = c(0,3), ylim = c(140,240), clip = "off") + 
  geom_segment(data = gaps_AsatC0648, aes(x = 3.2, xend = 3.15, y = gaps_AsatC0648$start[1], yend = gaps_AsatC0648$end[1]), color = "red", size = 0.3, linetype = "solid", inherit.aes = FALSE) +
  annotate("rect", xmin = -4.05, xmax = 3, ymin = 192.5, ymax=204, alpha=0.2) +
  theme(legend.text = element_markdown())
p_C0648

##### coliniarity plot OT3098 x AsatC0648 #####
comparison <- read.table("log_chr_comp_Asat_OT3098_v2_chr5D__x__Asat_C0648_chr5D.tab", skip=1)
comparison$OT3098 <- (comparison$V2)/1e+06
comparison$AsatC0648 <- (comparison$V3)/1e+06
comparison$Identity <- as.numeric(str_extract(comparison$V5, "\\d+.\\d+|\\d+"))

###### FIGURE 3 Panel d ###### 
p_center <- ggplot(data = comparison, aes(x = OT3098, y = AsatC0648, color = Identity)) +
  geom_point(size=0.1) +
  custom_theme2 +
  xlab(OT3098) +
  ylab(C0648) +
  scale_x_continuous(limits = c(140, 240), position = "bottom", breaks = scales::pretty_breaks(n = 10)) +
  scale_y_reverse(limits = c(240, 140), position = "left", breaks = scales::pretty_breaks(n = 10)) + 
  scale_colour_gradient(low = "gray80", high = "black", space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour")
p_center

(Panel_a | p_OT3098) / 
  (p_C0648 | p_center) +
  plot_layout(widths = c(4, 4), heights = c(4, 4), guides = 'collect')

#### FIGURE S11 ####
#### Preparation Asat OT3098 ####
Asat_Cereba <- read.table("age_distribution_Avenas_Cereba_clean_filtered", skip=1)
Asat_Cereba$family <- str_extract(Asat_Cereba$V1, "\\S+(?<=_Cereba)")

Asat_Ava <- read.table("age_distribution_Avenas_Ava_clean_filtered", skip=1)
Asat_Ava$family <- str_extract(Asat_Ava$V1, "\\S+(?<=_Ava)")

Asat_Beth <- read.table("age_distribution_Avenas_Beth_clean_filtered", skip=1)
Asat_Beth$family <- str_extract(Asat_Beth$V1, "\\S+(?<=_Beth)")

all_OT3098_AsatOT <- rbind(Asat_Ava, Asat_Beth, Asat_Cereba)

all_OT3098_AsatOT$TEnum <- str_extract(all_OT3098_AsatOT$V1, "\\d\\S-\\d+")
all_OT3098_AsatOT$chr <- str_extract(all_OT3098_AsatOT$TEnum, "\\d\\S")
all_OT3098_AsatOT$TEstart <- as.numeric(str_extract(all_OT3098_AsatOT$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the approximate start position from the ID 
all_OT3098_AsatOT$genome <- as.factor(str_extract(all_OT3098_AsatOT$chr, "(?<=\\d)\\S"))

df_sorted <- all_OT3098_AsatOT %>% arrange(genome)
df_sorted$chr <- factor(df_sorted$chr, levels = rev(unique(df_sorted$chr)))


chr4D_OT3098 <- df_sorted[df_sorted$chr == "4D" & df_sorted$V9 == "AsatOT3098",]

### gaps ###
list_of_files <- list.files(path = ".", pattern = paste0("annotation_N_Asat_OT3098_v2"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- sub(".*_(chr[0-9][A-Za-z]).*", "\\1", file)
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  df$group = df_name
  df <- df %>% separate(V1, into = c("Start", "End"), sep = "-")
  df <- df %>% separate(V2, into = c("Name", "Orientation", "Type", "Length", "other"), sep = ";")
  df <- df %>% separate(Length, into = c("delete", "Length"), sep = "=")
  df <- df %>% select(-delete, -other)
  assign(df_name, df, envir = .GlobalEnv)
  #df_names <- c(df_names, df_name)
}

gaps_OT3098 <- rbind(chr1A, chr1C, chr1D, chr2A, chr2C, chr2D, chr3A, chr3C, chr3D, chr4A, chr4C, chr4D, chr5A, chr5C, chr5D, chr6A, chr6C, chr6D, chr7A, chr7C, chr7D)

gaps_OT3098$genome <- str_extract(gaps_OT3098$group, "(?<=chr[0-9])\\S+")
gaps_OT3098$chr <- str_extract(gaps_OT3098$group, "(?<=chr)\\d+\\S+")

gaps_OT3098 <- gaps_OT3098 %>% arrange(group)
gaps_OT3098$group <- factor(gaps_OT3098$group, levels = unique(gaps_OT3098$group))

gaps_OT3098$start <- round(as.numeric(gaps_OT3098$Start)/1000000, 2)
gaps_OT3098$end <- round(as.numeric(gaps_OT3098$End)/1000000, 2)

gaps_OT3098 <- gaps_OT3098[gaps_OT3098$start>110 & gaps_OT3098$end<210 & gaps_OT3098$group == "chr4D",]

chr4D_OT3098 <- chr4D_OT3098[chr4D_OT3098$TEstart/1000000>110 & chr4D_OT3098$TEstart/1000000<210,]

p_right <- ggplot(data = chr4D_OT3098, aes(x = TEstart/1000000, y = V6, color = family)) +
  geom_point(shape=20, size=0.8)  +
  ylab("Insertion age [Myr]") +
  xlab("") +
  labs(color = "Family") +
  scale_color_manual(values=c("#322480", "#5aae61", "#d058b5ff"), labels = c("*RLG_Ava*","*RLG_Cereba*", "*RLG_Beth*")) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 6))+
  custom_theme2 +
  lapply(1:nrow(gaps_OT3098), function(i) {annotation_custom(grob = grid::rectGrob(gp = grid::gpar(col = "red", fill = NA, lwd = 1.2)), xmin = gaps_OT3098$start[i], xmax = gaps_OT3098$end[i], ymin = -0.20, ymax = -0.15)}) +
  coord_cartesian(ylim = c(0,3), xlim = c(110,210), clip = "off") +
  theme(legend.text = element_markdown()) +
  annotate("rect", xmin = 159, xmax = 167.75, ymin = -3.86, ymax=3, alpha=0.2)
p_right

#### Preparation Aste ####
Aste_Cereba <- read.table("age_tab_Cereba_Aste_clean_filtered", skip=1)
Aste_Cereba$family <- str_extract(Aste_Cereba$V1, "\\S+(?<=_Cereba)")

Aste_Ava <- read.table("age_tab_Ava_Aste_clean_filtered", skip=1)
Aste_Ava$family <- str_extract(Aste_Ava$V1, "\\S+(?<=_Ava)")

Aste_Beth <- read.table("age_tab_Beth_Aste_clean_filtered", skip=1)
Aste_Beth$family <- str_extract(Aste_Beth$V1, "\\S+(?<=_Beth)")

all_Aste <- rbind(Aste_Beth, Aste_Cereba, Aste_Ava)

all_Aste$TEnum <- str_extract(all_Aste$V1, "\\d\\S-\\d+")
all_Aste$chr <- str_extract(all_Aste$TEnum, "\\d\\S")
all_Aste$TEstart <- as.numeric(str_extract(all_Aste$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the approximate start position from the ID 
all_Aste$genome <- as.factor(str_extract(all_Aste$chr, "(?<=\\d)\\S"))

df_sorted <- all_Aste %>% arrange(genome)
df_sorted$chr <- factor(df_sorted$chr, levels = rev(unique(df_sorted$chr)))


chr4D_Aste <- df_sorted[df_sorted$chr == "4D",]

### gaps ###
# => no gaps annotated in 4D of A. sterilis

chr4D_Aste <- chr4D_Aste[chr4D_Aste$TEstart/1000000>110 & chr4D_Aste$TEstart/1000000<210,]

p_top <- ggplot(data = chr4D_Aste, aes(y = TEstart/1000000, x = V6, color = family)) +
  geom_point(shape=20, size=0.8)  +
  xlab("Insertion age [Myr]") +
  ylab("") +
  labs(color = "Family") +
  scale_color_manual(values=c("#322480", "#d058b5ff", "#5aae61"), labels = c("*RLG_Ava*", "*RLG_Beth*" ,"*RLG_Cereba*")) + 
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_reverse(breaks = scales::pretty_breaks(n = 10))+
  custom_theme2 +
  coord_cartesian(ylim = c(210, 110), xlim = c(0,3), clip = "off") + 
  theme(legend.text = element_markdown()) +
  annotate("rect", xmin = -3.91, xmax = 3, ymin = 167, ymax = 175.75, alpha=0.2)
p_top


custom_theme_empty <- theme_classic(base_size = 11) +
  theme(
    axis.title = element_text(size = 13),
    axis.text = element_text(color = NA),
    axis.line = element_line(color = NA),
    strip.text = element_text(size = 13),
    panel.grid.major = element_line(color = NA),
    panel.grid.minor = element_line(color = NA),   
    panel.background = element_rect(fill = NA, color = NA),
    plot.background = element_rect(fill = NA, color = NA),
    axis.ticks = element_line(color = NA)
  )

p_top_spacer <- ggplot(data = chr4D_Aste, aes(y = TEstart/1000000, x = V6, color = family)) +
  xlab("") +
  ylab("") +
  labs(color = "Family") +
  scale_x_continuous(breaks = scales::pretty_breaks(n = 6)) +
  scale_y_reverse(breaks = scales::pretty_breaks(n = 10))+
  custom_theme_empty +
  coord_cartesian(ylim = c(210, 110), xlim = c(0,3), clip = "off")
p_top_spacer


#### coliniarity plot OT3098 x Aste ####
comparison <- read.table("log_chr_comp_Asat_OT3098_v2_chr4D__x__Aste_v1_chr4D.tab", skip=1)
comparison$OT3098 <- (comparison$V2)/1e+06
comparison$Aste <- (comparison$V3)/1e+06
comparison$Identity <- as.numeric(str_extract(comparison$V5, "\\d+.\\d+|\\d+"))


p_center <- ggplot(data = comparison, aes(x = OT3098, y = Aste, color = Identity)) +
  geom_point(size=0.1) +
  custom_theme2 +
  xlab("A.sativa OT3098 Position [Mb]") +
  ylab("A.sterilis Position [Mb]") +
  scale_x_continuous(limits = c(110, 210), position = "bottom", breaks = scales::pretty_breaks(n = 10)) +
  scale_y_reverse(limits = c(210, 110), position = "left", breaks = scales::pretty_breaks(n = 10)) +
  scale_colour_gradient(low = "gray80", high = "black", space = "Lab", na.value = "grey50", guide = "colourbar", aesthetics = "colour")
p_center

(p_right | p_top_spacer) /
  (p_center | p_top) +
  plot_layout(widths = c(4, 4), heights = c(4, 4), guides = 'collect')
