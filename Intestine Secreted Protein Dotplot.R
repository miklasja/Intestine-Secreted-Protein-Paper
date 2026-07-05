library(ggplot2)
library(RColorBrewer)

# Set the working directory to the directory of the currently open script
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
d <- read.csv("Lifespan_Screen_Summary.csv")

ggplot(data = d, aes(x = Adjusted_P_value, y = reorder(Gene, Adjusted_P_value), 
                        color = Lifespan_percent_change, size = Adjusted_P_value)) + 
  geom_point() +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu")))+
  theme_bw() + 
  ylab("Gene") + 
  xlab("-Log(P-value)") + 
  ggtitle("Lifespan screen summary")
ggsave("Adjusted_P_value_ordered3.pdf", width = 5.5, height = 6.5)
