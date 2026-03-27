#### Trees for FIGURE 1 (c-f), FIGURE S5, S6, S8 and S10 ####
### Set path to where you stored the data as string, make sure to keep folder structure intact.
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(ggtext)

path <- "Figure_1_Figures_S5_S6_S8_S10/"
setwd(path)

#### FIGURE 1 ####
##### Preparation #####
species <- expression(paste(italic("Avena sativa"), " OT3098"))
family <- "RLG_AsatOT3098"
family_name <- expression(bold("Ava Cereba"))

Ava <- read.table("age_distribution_Avenas_Ava_clean_filtered", skip=1)
AsatOT3098_Ava <- Ava[Ava$V9 == "AsatOT3098",]
AsatOT3098_Ava$family <- paste("Ava")
colnames(AsatOT3098_Ava)[colnames(AsatOT3098_Ava) == "V9"] <- "speciesGenome"
AsatOT3098_Ava$genome <- str_extract(AsatOT3098_Ava$V1, "(?<=chr\\d)[A-Za-z]|(?<=chr)[A-Za-z]+")
AsatOT3098_Ava$species <- str_extract(AsatOT3098_Ava$V1, "(?<=RLG_).*(?=_Ava.*)")

Cereba <- read.table("age_distribution_Avenas_Cereba_clean_filtered", skip=1)
AsatOT3098_Cereba <- Cereba[Cereba$V9 == "AsatOT3098",]
AsatOT3098_Cereba$family <- paste("Cereba")
colnames(AsatOT3098_Cereba)[colnames(AsatOT3098_Cereba) == "V9"] <- "speciesGenome"
AsatOT3098_Cereba$genome <- str_extract(AsatOT3098_Cereba$V1, "(?<=chr\\d)[A-Za-z]|(?<=chr)[A-Za-z]+")
AsatOT3098_Cereba$species <- str_extract(AsatOT3098_Cereba$V1, "(?<=RLG_).*(?=_Cereba.*)")

Beth <- read.table("age_distribution_Avenas_Beth_clean_filtered", skip=1)
AsatOT3098_Beth <- Beth[Beth$V9 == "AsatOT3098",]
AsatOT3098_Beth$family <- paste("Beth")
colnames(AsatOT3098_Beth)[colnames(AsatOT3098_Beth) == "V9"] <- "speciesGenome"
AsatOT3098_Beth$genome <- str_extract(AsatOT3098_Beth$V1, "(?<=chr\\d)[A-Za-z]|(?<=chr)[A-Za-z]+")
AsatOT3098_Beth$species <- str_extract(AsatOT3098_Beth$V1, "(?<=RLX_).*(?=_Beth.*)")

all_OT3098 <- rbind(AsatOT3098_Ava, AsatOT3098_Cereba, AsatOT3098_Beth)

##### import sizes of the chromosomes #####
sizeList <- read.table("Asat_OT3098_v2_genome_size_list_chr", header = TRUE)
sizeList <- sizeList %>% arrange(Name)
sizeList$chr <- str_extract(sizeList$Name, "(?<=chr)\\d\\S")
sizeList$chromosome <- str_extract(sizeList$chr, "\\d(?<=\\S)")
sizeList$genome <- as.factor(str_extract(sizeList$chr, "(?<=\\d)\\S"))
sizeList <- sizeList %>% arrange(genome)

sizeList <- mutate(sizeList, chr_nr=1:length(sizeList$Name))

sizeList <- sizeList %>% arrange(desc(chr_nr))
sizeList <- mutate(sizeList, chr_nr_sorted= rep(1:7, length.out = nrow(sizeList)))

all_OT3098$TEnum <- str_extract(all_OT3098$V1, "\\d\\S-\\d+")
all_OT3098$chr <- str_extract(all_OT3098$TEnum, "\\d\\S")
all_OT3098$TEstart <- as.numeric(str_extract(all_OT3098$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the approximate start position from the ID 
all_OT3098$genome <- as.factor(str_extract(all_OT3098$chr, "(?<=\\d)\\S"))

df_sorted <- all_OT3098 %>% arrange(genome)
df_sorted$chr <- factor(df_sorted$chr, levels = rev(unique(df_sorted$chr)))
df_sorted$chromosome_full <- str_extract(df_sorted$V1, "chr\\d+[A-Za-z]")
df_sorted <- merge(df_sorted, sizeList, by = "chr")

##### import centromere position #####
OT3098_centromere_pos <- read.table("ChiP_inferred_centromere_positions_OT3098", header = T)
OT3098_centromere_pos$chr <- str_extract(OT3098_centromere_pos$OT3098, "(?<=chr)\\d\\S")
df_sorted_centromere_pos <-merge(df_sorted, OT3098_centromere_pos, by = "chr", all.x = T)

##### FIGURE 1 Panel  c #####
df_sorted_centromere_pos_Cereba <- df_sorted_centromere_pos[df_sorted_centromere_pos$family == "Cereba",]
str(df_sorted_centromere_pos_Cereba)

ggplot(data = df_sorted_centromere_pos_Cereba, aes(x = TEstart/1000000, y = chr, color = family)) +
  geom_point(shape=108, size=6.6, alpha= 0.4) +
  labs(color = "Families") +
  xlab("Position [Mb]") +
  ylab("") +
  scale_color_manual(values=c("#5aae61")) +
  scale_x_continuous(limits = c(0, 720), breaks = scales::pretty_breaks(n = 10)) +
  geom_rect(aes(xmin=0, xmax=Length/1000000, ymin=chr_nr_sorted-0.34, ymax=chr_nr_sorted+0.38), color ="black", fill = NA, linewidth=0.1, inherit.aes = F) +
  geom_rect(aes(xmin=start_full/1000000, xmax=end_full/1000000, ymin=chr_nr_sorted+0.44, ymax=chr_nr_sorted+0.48), color = "black", fill = "black", linewidth=0.1, inherit.aes = F) +
  theme_classic(base_size = 12) +
  facet_grid(genome.x~., drop = TRUE, scales = "free_y")+
  theme(strip.text = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18), axis.text.y = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18)) +
  scale_x_continuous(expand = c(0, 0))


##### FIGURE 1 Panel  d #####
df_sorted_centromere_pos_Ava <- df_sorted_centromere_pos[df_sorted_centromere_pos$family == "Ava",]
str(df_sorted_centromere_pos_Ava)

ggplot(data = df_sorted_centromere_pos_Ava, aes(x = TEstart/1000000, y = chr, color = family)) +
  geom_point(shape=108, size=6.6, alpha= 0.4) +
  labs(color = "Families") +
  xlab("Position [Mb]") +
  ylab("") +
  scale_color_manual(values=c("#322480")) +
  scale_x_continuous(limits = c(0, 720), breaks = scales::pretty_breaks(n = 10)) +
  geom_rect(aes(xmin=0, xmax=Length/1000000, ymin=chr_nr_sorted-0.34, ymax=chr_nr_sorted+0.38), color ="black", fill = NA, linewidth=0.1, inherit.aes = F) +
  geom_rect(aes(xmin=start_full/1000000, xmax=end_full/1000000, ymin=chr_nr_sorted+0.44, ymax=chr_nr_sorted+0.48), color = "black", fill = "black", linewidth=0.1, inherit.aes = F) +
  theme_classic(base_size = 12) +
  facet_grid(genome.x~., drop = TRUE, scales = "free_y")+
  theme(strip.text = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18), axis.text.y = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18)) +
  scale_x_continuous(expand = c(0, 0))

##### FIGURE 1 Panel  e #####
genome_labels <- data.frame(genome = c("A", "C", "D"), annotation = c("n = 1138", "n = 53", "n = 1013"), y = -0.8)
ggplot(AsatOT3098_Cereba, aes(y=V6, x=genome, fill=genome)) +
  geom_violin(color="black", stat = "ydensity", position = "dodge", trim = T, alpha = 1, scale = "count")  +
  geom_boxplot(width=0.15, color="white", alpha=1, fill = NA, position = position_dodge(width = 0.9), aes(group = interaction(genome, species))) + # ...age_facet_Groups
  geom_boxplot(width=0.15, color="gray70", alpha=1, fill = NA) +
  ylab("Insertion age [Myr]") +
  labs(fill = "Subgenome") +
  xlab("") +
  theme_minimal() +
  geom_text(data = genome_labels, aes(x = genome, y = y, label = annotation), inherit.aes = FALSE, size = 8, vjust = 1.5) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(text = element_text(size = 30), axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"), axis.text.y = element_text(colour = "black"), legend.position = "none") +
  scale_fill_manual(values=c("#5aae61",'#5aae61',"#5aae61")) +
  coord_cartesian(ylim = c(0, 4.5), clip = "off")

##### FIGURE 1 Panel f #####
genome_labels <- data.frame(genome = c("A", "C", "D"), annotation = c("n = 1103", "n = 613", "n = 827"), y = -0.8)
ggplot(AsatOT3098_Ava, aes(y=V6, x=genome, fill=genome)) +
  geom_violin(color="black", stat = "ydensity", position = "dodge", trim = T, alpha = 1, scale = "count")  +
  geom_boxplot(width=0.15, color="white", alpha=1, fill = NA, position = position_dodge(width = 0.9), aes(group = interaction(speciesGenome, species))) + # ...age_facet_Groups
  geom_boxplot(width=0.15, color="gray70", alpha=1, fill = NA) +
  ylab("Insertion age [Myr]") +
  labs(fill = "Subgenome") +
  xlab("") +
  theme_minimal() +
  geom_text(data = genome_labels, aes(x = genome, y = y, label = annotation), inherit.aes = FALSE, size = 8, vjust = 1.5) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(text = element_text(size = 30), axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"), axis.text.y = element_text(colour = "black"), legend.position = "none") +
  scale_fill_manual(values=c("#322480",'#322480',"#322480")) +
  coord_cartesian(ylim = c(0, 4.5), clip = "off")

#### FIGURE S6 ####
df_sorted_centromere_pos_Beth <- df_sorted_centromere_pos[df_sorted_centromere_pos$family == "Beth",]

##### FIGURE S6 Panel a #####
ggplot(data = df_sorted_centromere_pos_Beth, aes(x = TEstart/1000000, y = chr, color = family)) +
  geom_point(shape=108, size=6.6, alpha= 0.3) +
  labs(color = "Families") +
  xlab("Position [Mb]") +
  ylab("") +
  scale_color_manual(values=c("#d058b5ff")) +
  scale_x_continuous(limits = c(0, 720), breaks = scales::pretty_breaks(n = 10)) +
  geom_rect(aes(xmin=0, xmax=Length/1000000, ymin=chr_nr_sorted-0.34, ymax=chr_nr_sorted+0.38), color ="black", fill = NA, linewidth=0.1, inherit.aes = F) +
  geom_rect(aes(xmin=start_full/1000000, xmax=end_full/1000000, ymin=chr_nr_sorted+0.44, ymax=chr_nr_sorted+0.48), color = "black", fill = "black", linewidth=0.1, inherit.aes = F) +
  theme_classic(base_size = 12) +
  facet_grid(genome.x~., drop = TRUE, scales = "free_y")+
  theme(strip.text = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18), axis.text.y = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18)) +
  scale_x_continuous(expand = c(0, 0))

##### FIGURE S6 Panel b #####
genome_labels <- data.frame(genome = c("A", "C", "D"), annotation = c("n = 1054", "n = 171", "n = 930"), y = -0.8)
ggplot(AsatOT3098_Beth, aes(y=V6, x=genome, fill=genome)) +
  geom_violin(color="black", stat = "ydensity", position = "dodge", trim = T,scale = "count")  +
  geom_boxplot(width=0.15, color="white", alpha=1, fill = NA, position = position_dodge(width = 0.9), aes(group = interaction(genome, species))) + # ...age_facet_Groups
  geom_boxplot(width=0.15, color="#e9ecef", alpha=1, fill = NA) +
  ylab("Insertion age [Myr]") +
  labs(fill = "Subgenome") +
  xlab("") +
  theme_minimal() +
  geom_text(data = genome_labels, aes(x = genome, y = y, label = annotation), inherit.aes = FALSE, size = 8, vjust = 1.5) + 
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10)) +
  theme(text = element_text(size = 30), axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black"), axis.text.y = element_text(colour = "black"), legend.position = "none") +
  scale_fill_manual(values=c("#d058b5ff",'#d058b5ff',"#d058b5ff")) +
  coord_cartesian(ylim = c(min(AsatOT3098_Beth$V6), max(AsatOT3098_Beth$V6)), clip = "off")

#### FIGURE S8 ####
df_sorted_centromere_pos_Cereba$age_group <- ifelse(df_sorted_centromere_pos_Cereba$V6<1, "< 1 Myr", "> 1 Myr")

ggplot(data = df_sorted_centromere_pos_Cereba, aes(x = TEstart/1000000, y = chr, color = age_group)) +
  geom_point(shape=108, size=6.6, alpha= 0.2) +
  labs(color = "") +
  xlab("Position [Mb]") +
  ylab("") +
  scale_color_manual(values=c("#b2182b","#4393c3")) +
  scale_x_continuous(limits = c(0, 720), breaks = scales::pretty_breaks(n = 10)) +
  geom_rect(aes(xmin=0, xmax=Length/1000000, ymin=chr_nr_sorted-0.33, ymax=chr_nr_sorted+0.365), color ="black", fill = NA, linewidth=0.1, inherit.aes = F) +
  geom_rect(aes(xmin=start_full/1000000, xmax=end_full/1000000, ymin=chr_nr_sorted+0.44, ymax=chr_nr_sorted+0.48), color = "black", fill = "black", linewidth=0.1, inherit.aes = F) +
  theme_classic(base_size = 12) +
  facet_grid(genome.x~., drop = TRUE, scales = "free_y")+
  theme(strip.text = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18), axis.text.y = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18)) +
  scale_x_continuous(expand = c(0, 0))

#### FIGURE S10 ####
df_sorted_centromere_pos_Ava$age_group <- ifelse(df_sorted_centromere_pos_Ava$V6<2, "< 2 Myr", "> 2 Myr")

ggplot(data = df_sorted_centromere_pos_Ava, aes(x = TEstart/1000000, y = chr, color = age_group)) +
  geom_point(shape=108, size=6.6, alpha= 0.2) +
  labs(color = "") +
  xlab("Position [Mb]") +
  ylab("") +
  scale_color_manual(values=c("#b2182b","#4393c3")) +
  scale_x_continuous(limits = c(0, 720), breaks = scales::pretty_breaks(n = 10)) +
  geom_rect(aes(xmin=0, xmax=Length/1000000, ymin=chr_nr_sorted-0.33, ymax=chr_nr_sorted+0.365), color ="black", fill = NA, linewidth=0.1, inherit.aes = F) +
  geom_rect(aes(xmin=start_full/1000000, xmax=end_full/1000000, ymin=chr_nr_sorted+0.44, ymax=chr_nr_sorted+0.48), color = "black", fill = "black", linewidth=0.1, inherit.aes = F) +
  theme_classic(base_size = 12) +
  facet_grid(genome.x~., drop = TRUE, scales = "free_y")+
  theme(strip.text = element_blank(), axis.line.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size = 18),
        axis.text.x = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18), axis.text.y = element_text(angle = 0, hjust = 0.5, colour = "black", size = 18)) +
  scale_x_continuous(expand = c(0, 0))

#### FIGURE S5 ####
Aurora <- read.table("age_distribution_Avenas_Aurora_clean_filtered", skip=1)
Ava <- read.table("age_distribution_Avenas_Ava_clean_filtered", skip=1)
Beth <- read.table("age_distribution_Avenas_Beth_clean_filtered", skip=1)
Cereba <- read.table("age_distribution_Avenas_Cereba_clean_filtered", skip=1)

all_TEs <- rbind(Ava, Beth, Cereba, Aurora)
all_TEs$ID_1 <- str_extract(all_TEs$V1, "(?<=RLG_)[A-Za-z]{4}|(?<=RLX_)[A-Za-z]{4}")
all_TEs$ID_2 <- str_extract(all_TEs$V1, "[A-Za-z0-9]+(?=_chr)|[A-Za-z0-9]+(?=_consensus-1_)")
all_TEs$ID_3 <- paste(all_TEs$ID_1, all_TEs$ID_2, all_TEs$V8, sep = "_")

length(all_TEs[all_TEs$ID_3 == "Aatl_Ava_A",]$V1)#485
length(all_TEs[all_TEs$ID_3 == "Alon_Ava_A",]$V1)#885
length(all_TEs[all_TEs$ID_3 == "Asat_Ava_A",]$V1)#1103
length(all_TEs[all_TEs$ID_3 == "Asat_Ava_C",]$V1)#613
length(all_TEs[all_TEs$ID_3 == "Asat_Ava_D",]$V1)#827
length(all_TEs[all_TEs$ID_3 == "Ains_Ava_C",]$V1)#544
length(all_TEs[all_TEs$ID_3 == "Ains_Ava_D",]$V1)#650
length(all_TEs[all_TEs$ID_3 == "Aeri_Ava_C",]$V1)#646

count_df <- all_TEs %>% count(ID_3)
count_df <- count_df %>% separate(ID_3, into = c("species", "family", "genome"), sep = "_")
count_df$family <- paste("RLG", count_df$family, sep = "_")
count_df$ID <- paste(count_df$species, count_df$genome, sep = "_")
str(count_df)

unique(count_df$ID)
count_df$ID <- factor(count_df$ID, levels = rev(c("Aatl_A", "Alon_A", "Asat_A", "Aeri_C", "Ains_C", "Asat_C", "Ains_D", "Asat_D")))

ggplot(data = count_df, aes(x = ID, y = n, fill = family)) +
  geom_bar(stat = "identity", position="dodge") +
  coord_flip() +
  scale_fill_manual(values = c("gold", "#322480", "#d058b5ff", "#5aae61"), name = "TE Family", labels = c("*RLG_Aurora*", "*RLG_Ava*", "*RLG_Beth*", "*RLG_Cereba*")) +
  theme_classic(base_size = 16) +
  theme(legend.text = element_markdown()) +
  xlab("") +
  ylab("Number of full-length copies") +
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))

