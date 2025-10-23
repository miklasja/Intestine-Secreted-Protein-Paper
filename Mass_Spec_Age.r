#Install the following package

is_available <- require("EnhancedVolcano")
if(!is_available) {
    install.packages("EnhancedVolcano")
}

#Load in your data
gene_data <- read.csv("OvY_Cyto_Volcano.csv")

#Load in custom labels
Up_genes <- readLines("OvY_Cyto_Up.txt")
Down_genes <- readLines("OvY_Cyto_Down.txt")

#Remove zero-value hits
gene_data <- gene_data[!is.na(gene_data[,3]),] ##empty FC values.

#Reformat
sig_hits <- gene_data[, 1:3]
colnames(sig_hits) <- c("protein", "padj", "log2FoldChange")

#Color list below
#Young colour #00AEEF
#Middle-age colour #CCC9FF
#Old colour #D6AAB6


colCustom = rep("grey", dim(sig_hits)[1])
colCustom[which(sig_hits$padj < 0.1)] <- "#0000C0"
colCustom[which(sig_hits$padj > 0.1)] <- "grey"

names(colCustom)[colCustom == "grey"] <- "Not significant"
names(colCustom)[colCustom == "#0000C0"] <- "Significant"

#For any set of genes, add a new color
add_new_color <- function(custom_colour, geneset, desired_col, label) {
	custom_colour[which(sig_hits$protein %in% geneset)] <- desired_col
	names(custom_colour)[custom_colour == desired_col] <- label
	return(custom_colour)
}

colCustom <- add_new_color(colCustom, Up_genes, "#D6AAB6", "Identified old significant")
colCustom <- add_new_color(colCustom, Down_genes, "#00AEEF", "Identified young significant")

sig_hits_ordered <- sig_hits[order(sig_hits$padj),]

top_genes <- sig_hits_ordered$protein[1:10]

top_genes <- ""

out_title <- "example volcanoplot"

library(EnhancedVolcano)
library(EnhancedVolcano)

axis_tick_size <- 35 ##size of axis ticks
axis_label_size <- 40 ##size of axis title

other_volc <- EnhancedVolcano(sig_hits, lab = sig_hits$protein,
                              x = "log2FoldChange", y = "padj", 
                              
                              pCutoff = 0.1,
                              
                              selectLab = top_genes, 
                              
                              FCcutoff = 0, 
                              
                              ylab = "P-value",
                              title = NULL, subtitle = out_title, caption = "",  
                              
                              colCustom = colCustom, 
                              
                              gridlines.major = FALSE, gridlines.minor = FALSE,
                              labSize = 6.0, axisLabSize = axis_tick_size, pointSize = 4) #legendPosition = "topright",

DESIRED_YMAX <- NA #Change this if you have a particular -log10 threshold to use, if you want to hardcode replace NA with number.

if(!is.na(DESIRED_YMAX)) {
  y_max <- c(0, DESIRED_YMAX)
} else {
  y_max <- c(0, max(-log10(sig_hits$padj), na.rm = TRUE))
}

other_volc <- other_volc + ylim(y_max) + theme(
  axis.title.x = element_text(size = axis_label_size),
  axis.title.y = element_text(size = axis_label_size)
)

#Display in your RStudio

x_min = min(sig_hits$log2FoldChange)
x_max = max(sig_hits$log2FoldChange)
MAX_Y = max(-log10(sig_hits$padj))

print(other_volc)

#save to PDF

pdf("volcano_out.pdf", width = 12, height = 12)
print(other_volc)

dev.off()
