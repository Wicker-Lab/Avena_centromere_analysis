#### Trees for FIGURE 4 (a,c,e,g), FIGURE S12 and S13 ####
### Set path to where you stored the data as string, make sure to keep folder structure intact.
library(stringr)
library(dplyr)
library(ape)
library(ggtree)
library(ggnewscale)
library(ggplot2)

#setwd to where you stored the data
path <- "Figure_4_panels_a_c_e_g_Figures_S12_S13_phylogenetic_trees"

#### Ava ####
setwd(path)
#### Import Data for Tree sequences ####
age_tab <- read.table("age_distribution_Avenas_Ava_clean_filtered", skip=1)
age_tab$family <- str_extract(age_tab$V1, "\\S+(?<=_Ava)")
age_tab$TEnum <- str_extract(age_tab$V1, "\\d\\S-\\d+")
age_tab$genome_wo_acc <- str_extract(age_tab$family, "(?<=RLG_)(.{4})")
age_tab$chr <- str_extract(age_tab$TEnum, "\\d\\S")
age_tab$TEstart <- as.numeric(str_extract(age_tab$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the approximate start position from the ID 
age_tab$genome <- as.factor(str_extract(age_tab$chr, "(?<=\\d)\\S"))

age_tab$TE_ID <- paste0(age_tab$genome_wo_acc, age_tab$TEnum, sep = "")
age_tab$TE_ID<- gsub("-", "", age_tab$TE_ID)

df_sorted <- age_tab %>% arrange(genome)
df_sorted$chr <- factor(df_sorted$chr, levels = rev(unique(df_sorted$chr)))
rownames(df_sorted) <- df_sorted[,"TE_ID"]

#### Import Tree file ####
##### Ava all #####
tree <- read.tree(paste0(path, "/", "dir_Ava_RAxML/Ava_all_map_bootstraps.raxml.support"))
plot(tree)
setwd(paste0(path, "/", "dir_Ava_all"))

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Ava_G", sub(".*([0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$subfamily = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}
group_data <- rbind(Ava_G1, Ava_G2, Ava_G3)

data <- merge(df_sorted, group_data, by.x = "V1", by.y = "TE")
data$short_names_1 <- str_extract(data$V1, "(?<=_)[A-Za-z]{4}")
data$short_names_2 <- str_extract(data$V1, "(?<=chr)[0-9][A-Z]-[0-9]+(?=_TSD)")
data$short_names <- paste0(data$short_names_1, "_", data$short_names_2)
data$short_names <-gsub("-", "", data$short_names) 
data$short_names <-gsub("_", "", data$short_names) 

#subset data for tree circle
list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, df_sorted, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,7)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_1 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]


tree_plot <- ggtree(tree, layout = "circular", size = 0.2)

plot_1 <- gheatmap(tree_plot, data_for_tree_sub_1, offset=0.01, width=.1,
                   colnames = F) +
  scale_fill_gradientn(colors = c("yellow", "red", "black"), space= "Lab", na.value = "white", name = "Insertion\nage\n[Myr]") +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_1

list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, data, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,20)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_2 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]

plot_3 <- plot_1 + new_scale_fill()
plot_3 <- gheatmap(plot_3, data_for_tree_sub_2, offset=0.08, width=.1,
                   colnames = F) +
  scale_fill_manual(values = c("#332288", "#88CCEE", "#44AA99", "gray80"), na.value = "gray89", name = "Subfamily") +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_3

##### Ava G1 ####
tree <- read.tree(paste0(path, "/", "dir_Ava_RAxML/Ava_G1_map_bootstraps.raxml.support"))
plot(tree)
setwd(paste0(path, "/", "dir_Ava_G1"))

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Ava_G1_", sub(".*([0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$subfamily = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}
group_data <- rbind(Ava_G1_1,Ava_G1_2,Ava_G1_3)

data <- merge(df_sorted, group_data, by.x = "V1", by.y = "TE")
data$short_names_1 <- str_extract(data$V1, "(?<=_)[A-Za-z]{4}")
data$short_names_2 <- str_extract(data$V1, "(?<=chr)[0-9][A-Z]-[0-9]+(?=_TSD)")
data$short_names <- paste0(data$short_names_1, "_", data$short_names_2)
data$short_names <-gsub("-", "", data$short_names) 
data$short_names <-gsub("_", "", data$short_names) 

#subset data for tree circle
list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, df_sorted, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,7)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_1 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]

tree_plot <- ggtree(tree, layout = "circular", size = 0.2)

plot_1 <- gheatmap(tree_plot, data_for_tree_sub_1, offset=0.01, width=.1,
                   colnames = F) +
  scale_fill_gradientn(colors = c("yellow", "red", "black"), space= "Lab", na.value = "white", name = "Insertion\nage\n[Myr]") +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_1

list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, data, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,20)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_2 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]

plot_3 <- plot_1 + new_scale_fill()
plot_3 <- gheatmap(plot_3, data_for_tree_sub_2, offset=0.05, width=.1,
                   colnames = F) +
  scale_fill_manual(values = c("#6EAA22", "#008a6d", "#005e70", "gray80"), na.value = "gray89", name = "Subfamily") +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_3

##### Ava G1-1 #####
tree <- read.tree(paste0(path, "/", "dir_Ava_RAxML/Ava_G1_1_map_bootstraps.raxml.support"))
plot(tree)
setwd(paste0(path, "/", "dir_Ava_G1/dir_Ava_G1_1/"))

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Ava_G1_", sub(".*([0-9]_[0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$subfamily = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}
group_data <- rbind(Ava_G1_1_1, Ava_G1_1_2)

data <- merge(df_sorted, group_data, by.x = "V1", by.y = "TE")
data$short_names_1 <- str_extract(data$V1, "(?<=_)[A-Za-z]{4}")
data$short_names_2 <- str_extract(data$V1, "(?<=chr)[0-9][A-Z]-[0-9]+(?=_TSD)")
data$short_names <- paste0(data$short_names_1, "_", data$short_names_2)
data$short_names <-gsub("-", "", data$short_names) 
data$short_names <-gsub("_", "", data$short_names) 

#subset data for tree circle
list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, df_sorted, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,7)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_1 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]


tree_plot <- ggtree(tree, layout = "circular", size = 0.2)
plot_1 <- gheatmap(tree_plot, data_for_tree_sub_1, offset=0.01, width=.1,
                   colnames = F) +
  scale_fill_gradientn(colors = c("yellow", "red", "black"), space= "Lab", na.value = "white", name = "Insertion\nage\n[Myr]") +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_1

list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, data, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,20)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_2 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]

plot_3 <- plot_1 + new_scale_fill()
plot_3 <- gheatmap(plot_3, data_for_tree_sub_2, offset=0.022, width=.1,
                   colnames = F) +
  scale_fill_manual(values = c("#DDCC77", "#199f47ff"), na.value = "gray89", name = "Subfamily") +
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_3

##### Ava G2 #####
tree <- read.tree(paste0(path, "/", "dir_Ava_RAxML/Ava_G2_map_bootstraps.raxml.support"))
plot(tree)
setwd(paste0(path, "/", "dir_Ava_G2/"))

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca"),  full.names = FALSE) 
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Ava_G2_", sub(".*([0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$subfamily = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}
group_data <- rbind(Ava_G2_1, Ava_G2_2, Ava_G2_3, Ava_G2_4, Ava_G2_5, Ava_G2_6, Ava_G2_7, Ava_G2_8, Ava_G2_9)

data <- merge(df_sorted, group_data, by.x = "V1", by.y = "TE")
data$short_names_1 <- str_extract(data$V1, "(?<=_)[A-Za-z]{4}")
data$short_names_2 <- str_extract(data$V1, "(?<=chr)[0-9][A-Z]-[0-9]+(?=_TSD)")
data$short_names <- paste0(data$short_names_1, "_", data$short_names_2)
data$short_names <-gsub("-", "", data$short_names) 
data$short_names <-gsub("_", "", data$short_names) 

#subset data for tree circle
list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, df_sorted, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,7)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_1 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]


tree_plot <- ggtree(tree, layout = "circular", size = 0.2)

plot_1 <- gheatmap(tree_plot, data_for_tree_sub_1, offset=0.01, width=.1,
                   colnames = F) +
  scale_fill_gradientn(colors = c("yellow", "red", "black"), space= "Lab", na.value = "white", name = "Insertion\nage\n[Myr]")+
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_1


list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, data, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,20)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_2 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]

plot_3 <- plot_1 + new_scale_fill()
plot_3 <- gheatmap(plot_3, data_for_tree_sub_2, offset=0.06, width=.1,
                   colnames = F) +
  scale_fill_manual(values = c("gray20", "gray50", "#8c510a", "#bf812d", "#dfc27d", "#f6e8c3", "#80cdc1", "#35978f", "#01665e"), na.value = "gray89", name = "Subfamily")+
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_3 

##### Ava G3 #####
tree <- read.tree(paste0(path, "/", "dir_Ava_RAxML/Ava_G3_map_bootstraps.raxml.support"))
plot(tree)
setwd(paste0(path, "/", "dir_Ava_G3/"))

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Ava_G3_", sub(".*([0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$subfamily = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}
group_data <- rbind(Ava_G3_1, Ava_G3_2, Ava_G3_3, Ava_G3_4, Ava_G3_5, Ava_G3_6)

data <- merge(df_sorted, group_data, by.x = "V1", by.y = "TE")
data$short_names_1 <- str_extract(data$V1, "(?<=_)[A-Za-z]{4}")
data$short_names_2 <- str_extract(data$V1, "(?<=chr)[0-9][A-Z]-[0-9]+(?=_TSD)")
data$short_names <- paste0(data$short_names_1, "_", data$short_names_2)
data$short_names <-gsub("-", "", data$short_names) 
data$short_names <-gsub("_", "", data$short_names) 

#subset data for tree circle
list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, df_sorted, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,7)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_1 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]


tree_plot <- ggtree(tree, layout = "circular", size = 0.2)

plot_1 <- gheatmap(tree_plot, data_for_tree_sub_1, offset=0.01, width=.1,
                   colnames = F) +
  scale_fill_gradientn(colors = c("yellow", "red", "black"), space= "Lab", na.value = "white", name = "Insertion\nage\n[Myr]")+
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_1

list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, data, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,20)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_2 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]

plot_3 <- plot_1 + new_scale_fill()
plot_3 <- gheatmap(plot_3, data_for_tree_sub_2, offset=0.07, width=.1,
                   colnames = F) +
  scale_fill_manual(values = c("#AA4499", "#CC6677", "gray30", "#882255", "gray60", "black"), na.value = "gray89", name = "Subfamily")+
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_3 

#### Cereba ####
setwd(path)
#### Import Data for Tree sequences ####
age_tab <- read.table("age_distribution_Avenas_Cereba_clean_filtered", skip=1)
age_tab$family <- str_extract(age_tab$V1, "\\S+(?<=_Cereba)")
age_tab$TEnum <- str_extract(age_tab$V1, "\\d\\S-\\d+")
age_tab$genome_wo_acc <- str_extract(age_tab$family, "(?<=RLG_)(.{4})")
age_tab$chr <- str_extract(age_tab$TEnum, "\\d\\S")
age_tab$TEstart <- as.numeric(str_extract(age_tab$TEnum, "(?<=-)\\d+"))*1000 #time 1000 to get the start position from the ID 
age_tab$genome <- as.factor(str_extract(age_tab$chr, "(?<=\\d)\\S"))

age_tab$TE_ID <- paste0(age_tab$genome_wo_acc, age_tab$TEnum, sep = "")
age_tab$TE_ID<- gsub("-", "", age_tab$TE_ID)

df_sorted <- age_tab %>% arrange(genome)
df_sorted$chr <- factor(df_sorted$chr, levels = rev(unique(df_sorted$chr)))
rownames(df_sorted) <- df_sorted[,"TE_ID"]

##### Cereba all #####
tree <- read.tree(paste0(path, "/", "dir_Cereba_RAxML/Cereba_all_250_initial_run.raxml.bestTree"))
plot(tree)
setwd(paste0(path, "/", "dir_Cereba_all/"))

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Cereba_G", sub(".*([0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$subfamily = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}

Cereba_G1 <- rbind(Cereba_G1)

data <- merge(df_sorted, Cereba_G1, by.x = "V1", by.y = "TE")
data$short_names_1 <- str_extract(data$V1, "(?<=_)[A-Za-z]{4}")
data$short_names_2 <- str_extract(data$V1, "(?<=chr)[0-9][A-Z]-[0-9]+(?=_TSD)")
data$short_names <- paste0(data$short_names_1, "_chr", data$short_names_2)
data$short_names <-gsub("-", "", data$short_names) 
data$short_names <-gsub("_", "", data$short_names) 

#subset data for tree circle
list_for_tree <- as.data.frame(tree$tip.label)

data_for_tree <- merge(list_for_tree, df_sorted, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,7)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_1 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]


tree_plot <- ggtree(tree, layout = "circular", size = 0.2)

plot_1 <- gheatmap(tree_plot, data_for_tree_sub_1, offset=0.01, width=.1,
                   colnames = F) +
  scale_fill_gradientn(colors = c("yellow", "red", "black"), space= "Lab", name = "Insertion\nage\n[Myr]")+
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_1

list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, data, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,20)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_2 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]

plot_3 <- plot_1 + new_scale_fill()
plot_3 <- gheatmap(plot_3, data_for_tree_sub_2, offset=0.08, width=.1,
                   colnames = F) +
  scale_fill_manual(values = c("seagreen"), na.value = "gray89", name = "Subfamily")+
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_3 


##### Cereba G1 #####
tree <- read.tree(paste0(path, "/", "dir_Cereba_RAxML/G1_1_250_initial_run.raxml.bestTree"))
plot(tree)
setwd(paste0(path, "/", "dir_Cereba_G1/"))

list_of_files <- list.files(path = ".", pattern = paste0("PCA_age_pca"),  full.names = FALSE)
df_names <- list()
for (file in list_of_files) {
  df_name <- paste0("Cereba_G1_", sub(".*([0-9]).*", "\\1", file))
  df <- read.table(file, header = FALSE, sep = "\t") 
  cat("Processing file:", file, "and assigning to variable:", df_name, "\n")
  names(df)[names(df) == "V1"] <- "TE"
  names(df)[names(df) == "V2"] <- "Age"
  df$subfamily = df_name
  df$genome = as.factor(str_extract(df$TE, "(?<=chr[1-7])\\S"))
  assign(df_name, df, envir = .GlobalEnv)
}

Cereba_G1 <- rbind(Cereba_G1_1, Cereba_G1_2, Cereba_G1_3, Cereba_G1_4, Cereba_G1_5, Cereba_G1_6, Cereba_G1_7)

data <- merge(df_sorted, Cereba_G1, by.x = "V1", by.y = "TE")
data$short_names_1 <- str_extract(data$V1, "(?<=_)[A-Za-z]{4}")
data$short_names_2 <- str_extract(data$V1, "(?<=chr)[0-9][A-Z]-[0-9]+(?=_TSD)")
data$short_names <- paste0(data$short_names_1, "_chr", data$short_names_2)
data$short_names <-gsub("-", "", data$short_names) 
data$short_names <-gsub("_", "", data$short_names) 

#subset data for tree circle
list_for_tree <- as.data.frame(tree$tip.label)

data_for_tree <- merge(list_for_tree, df_sorted, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,7)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_1 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]


tree_plot <- ggtree(tree, layout = "circular", size = 0.2)

plot_1 <- gheatmap(tree_plot, data_for_tree_sub_1, offset=0.01, width=.1,
                   colnames = F) +
  scale_fill_gradientn(colors = c("yellow", "red", "black"), space= "Lab", name = "Insertion\nage\n[Myr]")+
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_1

list_for_tree <- as.data.frame(tree$tip.label)
data_for_tree <- merge(list_for_tree, data, by.x = "tree$tip.label", by.y = "TE_ID")
rownames(list_for_tree) <- tree$tip.label

data_for_tree_sub <- data_for_tree[, c(1,20)]
rownames(data_for_tree_sub) <- data_for_tree_sub$`tree$tip.label`
data_for_tree_sub_2 <- data_for_tree_sub[,as.numeric(c(2)), drop = FALSE]

plot_3 <- plot_1 + new_scale_fill()
plot_3 <- gheatmap(plot_3, data_for_tree_sub_2, offset=0.029, width=.1,
                   colnames = F) +
  scale_fill_manual(values = c('#332288', '#88CCEE', '#44AA99', '#117733','#999933',"steelblue3", "gray25"), na.value = "gray89", name = "Subfamily")+
  theme(legend.title = element_text(size = 14), legend.text = element_text(size = 12))
plot_3 