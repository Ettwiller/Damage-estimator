library(ggplot2)

mutation <- read.table("all_for_R", header=FALSE, sep="")
colnames(mutation) <- c("abs","type","experiment","count","family", "damage");
typ = unique(mutation$type)
#typ <- c("A_C","A_T","T_G","T_A","G_T","G_C","C_A","C_G","A_G","T_C","G_A","C_T","G_G","C_C","A_A","T_T")




keep_only <- unique(mutation$experiment)
#keep_only <- c()

new_mutation = subset(mutation,  experiment %in% keep_only)
new_mutation$experiment <- factor(new_mutation$experiment, level=keep_only)


exp = unique(new_mutation$experiment)
l = length(keep_only)

#local_color <- c("blue", "orange", "red", "brown4", "blue", "orange", "red", "brown4","royalblue4","slateblue4","slateblue2" ,"black")
local_color <- c("darkblue", "blue", "slateblue4", "orange","orange", "orange","orange", "red","blue", "orange","orange", "orange","orange", "orange", "orange","orange", "blue", "red","orange","blue", "red","orange","blue", "orange","blue", "orange")

vector = c()
for (selected_type in typ)
{
  file_name = paste(selected_type, "_all.png")
  
  sub <- mutation[ which(mutation$type==selected_type), ]
  
  tiff(file_name,  width = 1000, height = 700, units = "px") 
  par(mar=c(2,2,2,2), oma=c(20,1,1,1))
  stripchart(count~experiment, subset = type ==selected_type, data = new_mutation,pch = 21, bg = "#00000050",vertical = TRUE, las =3, at = 1:l, method = "jitter",jitter =0.4, cex=1, col = local_color,  ylab="enrichment")
  
  
  #    boxplot(count~experiment, subset = type ==selected_type, at = 1:l +0.25, boxwex = 0.25, data = new_mutation, cex = 1, add = TRUE, outline = FALSE, col="gold", bty="n", xaxt="n")
  mtext(selected_type, cex=1.5)
  dev.off()

  
}
