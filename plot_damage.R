#!/usr/bin/env Rscript   
library(ggplot2)


args <- commandArgs(TRUE)




if (length(args)<2) {
  stop("At least two argument must be supplied [1] imput_file (output of estimate_damage.pl) and [2] output file (for example figure1.png)", call.=FALSE)
}


mutation <- read.table(args[1], header=FALSE, sep="")
colnames(mutation) <- c("abs","type","experiment","count","family", "damage");

typ <-c("G_T", "C_A", "C_T", "G_A", "T_A","A_T","A_G", "T_C","C_G","G_C","T_G","A_C")
new_mutation = subset(mutation,  type %in% typ)
new_mutation$type <- factor(new_mutation$type, level=typ)



#coloring scheme (feel free to change)
local_color <- c("cornflowerblue", "royalblue4","grey1", "grey10","grey20", "grey30", "grey40", "grey50", "grey60", "grey70", "grey80", "grey90", "grey100")

d<-ggplot(new_mutation, aes(x = reorder(type, damage), y = log2(damage), color=experiment)) + 
 #geom_point(position=position_jitter(height = 0.1), alpha = 0.6, size=1.5) +
  geom_point( alpha = 0.6, size=1.5) +    
    scale_colour_manual(values = local_color) +
  geom_hline(yintercept= log2(1.5), color = "#990000", linetype="dashed") +
annotate("text", x= 3, y = log2(1.6), label = "Above this line one third of variants is due to damage",  color = "#990000") +
 geom_hline(yintercept= 0, color = "grey") +	
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
        legend.position="bottom", axis.text.x = element_text(angle = 90, hjust = 1, size=11)) +
  #ylim(0, 1) +
  ggtitle("GIV scores for variant types")
ggsave(args[2], d, width=9, height=6) 


