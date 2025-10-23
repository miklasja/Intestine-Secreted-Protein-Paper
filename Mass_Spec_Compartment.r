#Install the following package

is_available <- require("EnhancedVolcano")
if(!is_available) {
    install.packages("EnhancedVolcano")
}

#Load in your data
gene_data <- read.csv("young_Compartment_Volcano.csv")

#Load in custom labels
target_genes <- readLines("Y_Cyto_genes.txt")
ER_genes <- readLines("Y_ER_genes.txt")

#Remove zero-value hits
gene_data <- gene_data[!is.na(gene_data[,3]),] ##empty FC values.

#Reformat
sig_hits <- gene_data[, 1:3]
colnames(sig_hits) <- c("protein", "padj", "log2FoldChange")

colCustom = rep("grey", dim(sig_hits)[1])
colCustom[which(abs(sig_hits$log2FoldChange) < 2)] <- "grey"
colCustom[which(sig_hits$padj > 0.1)] <- "grey"

names(colCustom)[colCustom == "grey"] <- "Not significant"

add_new_color <- function(custom_colour, geneset, desired_col, label) {
	##quite simply
	custom_colour[which(sig_hits$protein %in% geneset)] <- desired_col
	names(custom_colour)[custom_colour == desired_col] <- label
	return(custom_colour)
}

colCustom <- add_new_color(colCustom, target_genes, "#F39958", "Cytoplasm & Significant")
colCustom <- add_new_color(colCustom, ER_genes, "#0000C0", "ER & Significant")

sig_hits_ordered <- sig_hits[order(sig_hits$padj),]

top_genes <- sig_hits_ordered$protein[1:10]

top_genes <- ""

out_title <- "example volcanoplot"

library(EnhancedVolcano)
library(EnhancedVolcano)

axis_tick_size <- 30 ##size of axis ticks
axis_label_size <- 40 ##size of axis title

other_volc <- EnhancedVolcano(sig_hits, lab = sig_hits$protein,
                              x = "log2FoldChange", y = "padj", pCutoff = 0, selectLab = top_genes, FCcutoff = 0, ylab = "P-value",
                              title = NULL, subtitle = out_title, caption = "",  colCustom = colCustom, gridlines.major = FALSE, gridlines.minor = FALSE,
                              labSize = 6.0, axisLabSize = axis_tick_size, pointSize = 3) #legendPosition = "topright",

DESIRED_YMAX <- NA

if(!is.na(DESIRED_YMAX)) {
  y_max <- c(0, DESIRED_YMAX)
} else {
  y_max <- c(0, max(-log10(sig_hits$padj), na.rm = TRUE))
}

other_volc <- other_volc + ylim(y_max) + theme(
  axis.title.x = element_text(size = axis_label_size),
  axis.title.y = element_text(size = axis_label_size)
)

#Display in your rstudio

x_min = min(sig_hits$log2FoldChange)
x_max = max(sig_hits$log2FoldChange)
MAX_Y = max(-log10(sig_hits$padj))

other_volc <- other_volc + geom_segment(aes(x = -1, xend = -1, y = 1, yend = MAX_Y, colour = "segment"), color = "#F39958", lty = 2)
other_volc <- other_volc + geom_segment(aes(x = -16, xend = -1, y = 1, yend = 1, colour = "segment"), color = "#F39958", lty = 2)
other_volc <- other_volc + geom_segment(aes(x = 2, xend = 2, y = .3, yend = MAX_Y, colour = "segment"), color = "#0000C0", lty = 2)
other_volc <- other_volc + geom_segment(aes(x = 2, xend = 16, y = .3, yend = .3, colour = "segment"), color = "#0000C0", lty = 2)

print(other_volc)

#Save to a pdf

pdf("volcano_out.pdf", width = 12, height = 12)
print(other_volc)

dev.off()
