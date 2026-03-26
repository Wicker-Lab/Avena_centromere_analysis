### Contact: matthias.heuberger@botinst.uzh.ch or matthias.heuberger@hotmail.com
### Script to make figure for supplementary figure S4

library(ggplot2)
library(patchwork)

df <- read.table("fl_copies_table.csv", header = T)

df$subgenome <- factor(df$subgenome, levels = rev(c("A","C","D")))

colorz <- c("#322480", "#d058b5ff", "#5aae61")

a <- ggplot(df[df$fam=="Ava",], aes(x=subgenome, y=percentage, fill=class)) +
  geom_bar(stat = "identity", color="grey20", width = 0.75) + 
  coord_flip() +
  scale_x_discrete(name="Ava") +
  scale_fill_manual(values = rev(c("#322480","#A99bf6"))) +
  scale_y_continuous(name="Proportion of centromere (%)",limits = c(0,50)) +
  #labs(title = "Ava") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.justification = "top",
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size=15))
b <- ggplot(df[df$fam=="Beth",], aes(x=subgenome, y=percentage, fill=class)) +
  geom_bar(stat = "identity", color="grey20", width = 0.75) + 
  coord_flip() +
  scale_x_discrete(name="Beth") +
  scale_fill_manual(values = rev(c("#d058b5ff","#e8a8ddff"))) +
  scale_y_continuous(name="Proportion of centromere (%)",limits = c(0,50)) +
  #labs(title = "Ava") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.justification = "top",
        axis.title.x = element_blank(),
        axis.line.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        text = element_text(size=15))

c <- ggplot(df[df$fam=="Cereba",], aes(x=subgenome, y=percentage, fill=class)) +
  geom_bar(stat = "identity", color="grey20", width = 0.75) + 
  coord_flip() +
  scale_x_discrete(name="Cereba") +
  scale_fill_manual(values = rev(c("#5aae61","#9dd5a5"))) +
  scale_y_continuous(name="Proportion of centromere (%)",limits = c(0,50)) +
  #labs(title = "Ava") +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.justification = "top",
        text = element_text(size=15))

a/b/c