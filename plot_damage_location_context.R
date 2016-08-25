#!/usr/bin/env Rscript
library(ggplot2)
library(reshape2)

args <- commandArgs(TRUE)
if (length(args)<2) {
  stop("At least two argument must be supplied [1] imput_file (output of estimate_damage.pl) and [2] output file (for example figure1.png)", call.=FALSE)
}
argument2 = args[2]



mutation <- read.table(args[1], header=FALSE, sep="")
colnames(mutation) <- c("experiment","type","read","count","loc","context", "abs");

typ = unique(mutation$type)

keep_only <-unique(mutation$experiment)

new_mutation = subset(mutation,  experiment %in% keep_only)
new_mutation$experiment <- factor(new_mutation$experiment, level=keep_only)

exp = unique(new_mutation$experiment)
l = length(keep_only)
local_color <- c("red", "brown","royalblue4","green","orange","slateblue2" ,"purple","brown4","brown1", "orange","royalblue4", "orange", "royalblue1", "yellow", "slateblue4","red","brown4")


for (selected_type in typ)
{

    out = paste(selected_type, argument2, sep = "_")
       
    d<-ggplot(subset(new_mutation,type %in% c(selected_type))) + 
      geom_point(aes(x=loc, y=count, group=experiment, color=experiment, shape=read)) +
      facet_grid(context~read, scales = "free") +
#      scale_colour_manual(values = local_color) + 
      theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
      #theme(legend.position = "none") +      
      #scale_x_continuous(limits = c(0, 100)) +
     # scale_y_continuous(limits = c(0, 0.013)) +
      ggtitle(selected_type)
    ggsave(out, d, width=15, height=20) 
  
   
}


