#!/usr/bin/env Rscript   
library(ggplot2)


args <- commandArgs(TRUE)




if (length(args)<2) {
  stop("At least two argument must be supplied [1] imput_file (output of estimate_damage.pl) and [2] output file (for example figure1.png)", call.=FALSE)
}

argument2 = args[2]
mutation <- read.table(args[1], header=FALSE, sep="")
colnames(mutation) <- c("abs","type","experiment","count","family", "damage");

#typ <-c("G_T", "C_A", "C_T", "G_A", "T_A","A_T","A_G", "T_C","C_G","G_C","T_G","A_C")
typ <- unique(mutation$type)


new_mutation = subset(mutation,  type %in% typ)
new_mutation$type <- factor(new_mutation$type, level=typ)

keep_only <- unique(new_mutation$experiment)

#coloring scheme (feel free to change)
l = length(keep_only)
local_color <- rainbow(l)


for (selected_type in typ)
{
  
  out = paste(selected_type, argument2, sep = "_")

  d<-ggplot(subset(new_mutation,type %in% c(selected_type)), aes(count, reorder(experiment, count, mean), colour=experiment)) + 
             
  geom_point(position="jitter", alpha = 0.6, size=1.5) +
   
    scale_colour_manual(values = local_color) +
 
  theme(panel.background = element_rect(fill = 'white', colour = 'white'), 
        legend.position="none", axis.text.x = element_text(angle = 90, hjust = 1, size=11)) +
  ggtitle("error rate for variant types")
  
  
  ggsave(out, d, width=9, height=6) 

}
