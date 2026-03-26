
library("ggplot2")
library("ggrepel")
library("dplyr")
library("tidyverse")
library("ggtext")


setwd("~/data/dir_Avena_genome/dir_R_scripts")

getwd()

# this script was used for Fig S7, S9  and S16




# TE insertion age distribution box plots Fig S7 and S14 --------------------------

infile <- "age_distribution_clean_filtered_all"

df0 <-read.table(infile, sep='\t', header=T)
df0 <-read.table(infile, sep='\t', stringsAsFactors = FALSE, header=T)
head(df0)

# select subset for family to plot
fam <- "Aurora"
fam <- "Cereba"
fam <- "Ava"
fam <- "Beth"

fam2 <- paste("RLG",fam,sep='_')

df <- df0[grep(fam,df0$TE),]
head(df)

# make counts for groups in genome
n_labels <- df %>%
  group_by(Genome) %>%
  summarise(n = n(), .groups = "drop")

n_labels

# allow re-ordering
df$Genome <- as.factor(df$Genome)

title <- paste("Insertion age distribution",fam2,sep=' ')
colors <- c("#0072B2",'#009E73',"#E69F00")

# boxplot insertion ages
p <- ggplot(df, aes(x = age, y = Genome)) +
  geom_violin(alpha=0.5,width=0.8,scale = "count",aes(fill=Subgenome,color=Subgenome))+
  scale_fill_manual(values=colors)+
  scale_color_manual(values=colors)+
  geom_boxplot(alpha=0.5) +
  geom_jitter(width=0, height=0.2,size=0.2,alpha=0.3)+
  xlab("Insertion age [MYA]")+
  ylab("SpeciesGenome")+
  ggtitle(title)+
  geom_text(data = n_labels, aes(x = 2.5,y = Genome,label = paste0("Copies = ",n)),hjust = 0, vjust=-1,size = 3) +
  xlim(0,3) +
  theme_classic()
p


outfile <- paste("age_distribution_bx_plot",fam,sep='_')
path <- getwd()
out_png <- paste(outfile,".png",sep='')
wide <- 7
high <- 5
ggsave(p, filename=out_png, path = path, width=wide,height=high,dpi=600)





# compare insertion ages of TEs inside centromere vs. outside Fig. S9------------------

# the input file was produced with the perl script TEpop_find_centromere_copies
# which uses  age_distribution_clean_filtered_all as input

infile <- "TEpop_centrormere_dist_all_cent_v3"

df <-read.table(infile, sep='\t', header=T)
head(df)

# make counts for groups (needed for figure later)
n_labels <- df %>%
  group_by(flag) %>%
  summarise(n = n(), .groups = "drop")

n_labels
out_table <- "Cent_TEs_copy_numbers"
write.table(n_labels,file=out_table,sep="\t",quote = FALSE)


xmin <- 0
xmax <- 5

# boxplot 
p <- ggplot(df, aes(x = age, y = group)) +
  geom_boxplot(data=df,aes(color=rem)) +
  geom_violin(alpha=0.5,scale=("count"),aes(fill=rem))+
  geom_jitter(width=0, height=0.2,size=0.1,alpha=0.2)+
  xlab("Insertion age [MYA]")+
  ylab("TE family/location")+
  xlim(xmin, xmax) +
  theme_classic() +
  facet_grid(genome~.)

p


# write plot to output file 
outfile <- "box_blot_cent_vs_outside"
path <- getwd()
out_png <- paste(outfile,".png",sep='')
wide <- 10
high <- 6
ggsave(p, filename=out_png, path = path, width=wide,height=high,dpi=600)





# do shapiro test for nowmal distribution 

qqnorm(df[df$group == "RLG_Beth_IN" & df$genome == "A",]$age)
qqline(df[df$group == "RLG_Beth_IN" & df$genome == "A",]$age)
hist(df[df$group == "RLG_Beth_IN" & df$genome == "A",]$age)
shapiro.test(df[df$group == "RLG_Beth_IN" & df$genome == "A",]$age)

# ==> gives relatively reasonable values, but in any case the wilcox test is used


# wilcox test for differences in age of TEs inside vs. out side of centromeres ----------------------
geno <- c("A","C","D")
fami <- c("Ava","Cereba","Beth")

# loop over genomes and families 
for (g in geno) {
  
  for (i in fami) {
    
    df_comp <- df[df$family == i & df$genome == g,]  
    wilc <- wilcox.test(age~rem, data=df_comp, alternative = "less")
    pval <- wilc$p.value
    # round p-value
    pv2 <-  format(pval, digits = 3, scientific = TRUE)
    print(paste(g,i,pv2, sep=' '))
    
  }
  
}


# p-values and n-counts were added manually to the figure using Gimp 




# plot distance to centromere vs insertion age, used for plots in Fig S9b----------------

infile <- "TEpop_centromere_dist_vs_age"
df <-read.table(infile, sep='\t', header=T)
head(df)

fam <-"RLG_Ava"
fam <-"RLG_Beth"
fam <-"RLG_Cereba"

title <- paste("Distance to centrome\nvs. insertion age ",fam,sep='')

p <- ggplot(df[df$fam==fam,],aes(x=Cent_dist/1000,y=age)) +
  geom_bin2d() +
  scale_fill_gradient(low = "blue", high = "red") +
  geom_point(size=0.1,aplha=0.1,color="gray60")+
  labs(fill = "Copies") +
  ggtitle(title)+
  xlab("Distance to Centromere[kb]") +
  ylab("Insertion age [Myr]") +
  xlim(0,15)+
  theme_classic()
p


# write plot to output file
outfile <- paste(infile,"_",fam,sep='')
path <- getwd()
out_png <- paste(outfile,".png",sep='')
wide <- 5
high <- 5
ggsave(p, filename=out_png, path = path, width=wide,height=high,dpi=600)












