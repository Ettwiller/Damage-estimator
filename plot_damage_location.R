library(ggplot2)




mutation <- read.table("Damage_estimation_for_Run_NoPreCRBQ5BC5_S2_L001_R1_001_val_1.fq", header=FALSE, sep="")
colnames(mutation) <- c("experiment","type","read","count","abs","loc");

typ = unique(mutation$type)

keep_only <-unique(mutation$experiment)

new_mutation = subset(mutation,  experiment %in% keep_only)
new_mutation$experiment <- factor(new_mutation$experiment, level=keep_only)

exp = unique(new_mutation$experiment)
l = length(keep_only)
local_color <- c("red", "brown","royalblue4","green","orange","slateblue2" ,"purple","brown4","brown1", "orange","royalblue4", "orange", "royalblue1", "yellow", "slateblue4","red","brown4")


for (selected_type in typ)
{
    file_name = paste(selected_type, "_UDG_effect.png" ,sep="")
    
    d<-ggplot(subset(new_mutation,type %in% c(selected_type))) + 
      geom_point(aes(x=loc, y=count, group=experiment, color=experiment, shape=read)) +
      facet_grid(~read, scales = "fixed") +
      scale_colour_manual(values = local_color) + 
      theme(panel.background = element_rect(fill = 'white', colour = 'white')) +
      #theme(legend.position = "none") +      
      #scale_x_continuous(limits = c(0, 100)) +
     # scale_y_continuous(limits = c(0, 0.013)) +
      ggtitle(selected_type)
    ggsave(file_name, d, width=15, height=10) 
  
   
}


