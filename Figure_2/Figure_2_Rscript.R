### Contact: matthias.heuberger@botinst.uzh.ch or matthias.heuberger@hotmail.com
### Script to make figure 2

library(ggplot2)
library(dplyr)



ava <- read.table("~/data2/Asat_Ava_bases_covered_2Mb")

ava <- mutate(ava, TE_family="Ava")

cereba <- read.table("~/data2/Asat_Cereba_bases_covered_2Mb")

cereba <- mutate(cereba, TE_family="Cereba")

aurora <- read.table("~/data2/Asat_Aurora_bases_covered_2Mb")

aurora <- mutate(aurora, TE_family="Aurora")

beth <- read.table("~/data2/Asat_Beth_bases_covered_2Mb")

beth <- mutate(beth, TE_family="Beth")

Cen48 <- read.table("~/data2/Asat_Cen48_bases_covered_2Mb")

Cen48 <- mutate(Cen48, TE_family="Cen48")

Cen87 <- read.table("~/data2/Asat_Cen87_bases_covered")

Cen87 <- mutate(Cen87, TE_family="Cen87")


df <- rbind(ava, cereba, aurora, beth, Cen48, Cen87)

df$TE_family <- factor(df$TE_family, levels = c("Aurora", "Ava", "Beth", "Cereba", "Cen48", "Cen87"))

colorz <- c("#322480","#d058b5ff", "#5aae61", "orange", "red", "black")

chromosome_order <- c("Asat_OT3098_v2_chr1A","Asat_OT3098_v2_chr2A","Asat_OT3098_v2_chr3A","Asat_OT3098_v2_chr4A","Asat_OT3098_v2_chr5A","Asat_OT3098_v2_chr6A","Asat_OT3098_v2_chr7A","Asat_OT3098_v2_chr1C","Asat_OT3098_v2_chr2C","Asat_OT3098_v2_chr3C",
                      "Asat_OT3098_v2_chr4C",
                      "Asat_OT3098_v2_chr5C",
                      "Asat_OT3098_v2_chr6C",
                      "Asat_OT3098_v2_chr7C",
                      "Asat_OT3098_v2_chr1D",
                      "Asat_OT3098_v2_chr2D",
                      "Asat_OT3098_v2_chr3D",
                      "Asat_OT3098_v2_chr4D",
                      "Asat_OT3098_v2_chr5D",
                      "Asat_OT3098_v2_chr6D",
                      "Asat_OT3098_v2_chr7D")


df$V1 <- factor(df$V1, levels = rev(chromosome_order))

centromere_positions <- read.table("~/Downloads/centromere_positions_OT3098")

centromere_positions$V1 <- factor(centromere_positions$V1, levels = rev(chromosome_order))

a <- ggplot(df[df$V6>0&df$TE_family!="Aurora",], aes(x=V2/1000000, y=V6, color=TE_family)) +
  geom_point(size=0.5, alpha=0.5) +
  facet_wrap(~V1, ncol=1, as.table=FALSE) +
  #geom_line(data=df[df$TE_family=="Cen48"&df$TE_family=="Cen87",],size=0.75, linetype="dashed") +
  geom_point(data = centromere_positions, aes(x=V4/1000000, y=0.5), inherit.aes = F, size=0.25) +
  geom_point(data = centromere_positions, aes(x=V5/1000000, y=0.5), inherit.aes = F, size=0.25) +
  scale_color_manual(values = colorz) +
  scale_y_continuous(name = "Fraction of genome (%)") +
  scale_x_continuous(name = "Position on Chromosome (Mb)") +
  theme_bw() +
  theme(strip.text = element_blank(),
        legend.position = "top",
        legend.justification = "left",
        legend.title = element_blank(),
        text = element_text(size=15),
        axis.text.y = element_text(size=8)) 


taberu <- read.table("~/Downloads/Wunschtabelle.csv", header = T)

cen_order <- c("Asat_OT3098_v2_chr1A_cen","Asat_OT3098_v2_chr2A_cen","Asat_OT3098_v2_chr3A_cen","Asat_OT3098_v2_chr4A_cen","Asat_OT3098_v2_chr5A_cen","Asat_OT3098_v2_chr6A_cen","Asat_OT3098_v2_chr7A_cen","Asat_OT3098_v2_chr1C_cen","Asat_OT3098_v2_chr2C_cen","Asat_OT3098_v2_chr3C_cen",
               "Asat_OT3098_v2_chr4C_cen",
               "Asat_OT3098_v2_chr5C_cen",
               "Asat_OT3098_v2_chr6C_cen",
               "Asat_OT3098_v2_chr7C_cen",
               "Asat_OT3098_v2_chr1D_cen",
               "Asat_OT3098_v2_chr2D_cen",
               "Asat_OT3098_v2_chr3D_cen",
               "Asat_OT3098_v2_chr4D_cen",
               "Asat_OT3098_v2_chr5D_cen",
               "Asat_OT3098_v2_chr6D_cen",
               "Asat_OT3098_v2_chr7D_cen")


taberu$ID <- factor(taberu$ID, levels = rev(cen_order))

colorz_2 <- c("#322480","#d058b5ff", "#5aae61","lightblue" ,"darkorange", "grey70")

b <- ggplot(taberu, aes(fill=cat, y=value/total_len*100, x=ID)) +
  geom_bar(position = "stack", stat = "identity", width =  0.75, color="grey5", linewidth=0.25) +
  coord_flip() +
  theme_classic() +
  scale_fill_manual(values = colorz_2) +
  scale_y_continuous(name="Fraction of centromere (%)") +
  theme(legend.position = "top",
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank(),
        legend.title = element_blank(),
        text = element_text(size=15)) 

library(patchwork)

a+b