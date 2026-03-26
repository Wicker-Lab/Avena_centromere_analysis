
library("ggplot2")
library("ggrepel")
library("dplyr")
library("tidyverse")

setwd("~/data/dir_Avena_genome/dir_R_scripts")

getwd()




# plot for fraction of TE consensus aligned with individual TE copies, used for Fig S15 ------------------------

# run for all families and consensus lengths individually for the panels in Fig S15
fam <- "Ava"
len <- 8710

fam <- "Aurora"
len <- 14747

fam <- "Beth"
len <- 6861 # Beth

fam <- "Cereba"
len <- 7751 # Cereba


infile <- paste("evaluate_pair_RLG_",fam,"_all_species",sep='')
df <-read.table(infile, sep='\t', header=T)
head(df)


# make counts for groups in genome
n_labels <- df %>%
  group_by(query) %>%
  summarise(n = n(), .groups = "drop")

xmin <-70
xmax <- 100

title <- substitute("Full-length alignment TE copies vs. consensus " * italic(RLG_ * fam),list(fam = fam))

df$aln_rel <- df$len_aln/len*100
head(df)

xmin <-99
xmax <- 100

# boxplot for aligned bases 
p <- ggplot(df, aes(x = aln_rel, y = query)) +
  geom_boxplot() +
  geom_jitter(width=0, height=0.2,size=0.1,alpha=0.3)+
  geom_violin(alpha=0.5, scale=("count"))+
  ggtitle(title)+
  scale_y_discrete(limits = rev)+
  #theme_classic()
  
  geom_text(data = n_labels, aes(x = xmin,y = query,label = paste0("Copies = ",n)),hjust = 0, vjust=0,size = 3) +
  xlim(xmin, xmax) +
  xlab("Fraction laigned with consensus")+
  ylab("Species")+
  theme_classic()

p

# write output png 
outfile <- infile
path <- getwd()
out_png <- paste(outfile,"_len_aln.png",sep='')
wide <- 7
high <- 5
ggsave(p, filename=out_png, path = path, width=wide,height=high,dpi=600)


# Relationship between numbers of SNPs between LTRs of a full-length TE copy and insertion age, Fig S16 -------------------------------------
# the plot were produced with the provided Perl script png_TEpop_eval_insertion_age
# which uses age_distribution_clean_filtered_all as input 
# the script needs the following packages from ubuntu repositories:
#libgd-graph-perl
#libgd-perl
#libgd-text-perl



# Within-species pairwise comparisons of TE copies, used for Fig S17a --------------------------
# Full-length retrotransposon copies with estimated insertion ages younger than 50,000 years 
# were aligned all vs. all across their entire length

infile <- "evaluate_pair_all_x_all_v2_copies"

df <-read.table(infile, sep='\t', header=T)
head(df)

# get table wiuth sample sizes
infile <- "evaluate_pair_all_x_all_v2_sample_size"
labels <-read.table(infile, sep='\t', header=T)
head(labels)


# make counts for groups in seq1
n_labels <- df %>%
  group_by(seq1) %>%
  summarise(n = n(), .groups = "drop")


# boxplot with n-values
p <- ggplot(df, aes(x = sim, y = seq1)) +
  geom_boxplot() +
  geom_violin(alpha=0.5)+
  geom_jitter(data=df[df$sim<=99.2,],width=0, height=0.2,size=0.2,alpha=0.3)+
  geom_jitter(data=df[df$sim>99.2,], width=0, height=0.2,size=0.7, alpha=0.4,color='red')+
  xlab("Pairwise sequence indetity of TE copies [%]")+
  ylab("Species/TE family")+
  
  geom_text(
    data = labels,
    aes(x = 95.2,           # fixed x position for all labels
        y = seq1, 
        label = paste0("Copies = ",copies,"\nPairs = ", pairs)),
    hjust = 0,              # align text left of x=95.5
    size = 3
  ) +
  xlim(95, 100) +           # keep your desired x range
  theme_classic()

p


# write plot to output file, half page (8 cm) 
outfile <- "box_blot_copies_all_vs_all_sim_v2"
path <- getwd()
out_png <- paste(outfile,".png",sep='')
wide <- 6
high <- 3.5
ggsave(p, filename=out_png, path = path, width=wide,height=high,dpi=600)


# Proportions of C-to-T substitutions derived from variant calling of copies <100,000 years, used for Fig17b -----------
# the provided perl script NGS_extract_met_SNPs_from_vcf extracts C->T and G->A substitutions from the vcf files:
# RLG_Ava_age_less_than_100ky.vcf
# RLG_Cereba_age_less_than_100ky.vcf



