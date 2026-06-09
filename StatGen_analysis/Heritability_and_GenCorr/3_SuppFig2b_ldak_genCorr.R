library(tidyverse)
library(ggpubr)
library(data.table)
library(RColorBrewer)

# Read in genetic correlation results
gencor <- read.table("all_bcx2_genetic_correlation.tsv", sep="\t", header=T)

# Plot genetic correlation results
p_new <- ggplot(gencor, aes(x = Trait.2, y = rg, color = p.value)) +
  geom_pointrange(aes(ymin = rg - 1.96*se, ymax = rg + 1.96*se)) + 
  scale_color_gradient2(low="red", mid="grey", high = "black", 
                        midpoint=0.05, transform="log",
                        breaks=c(0.01, 0.02, 0.05, 0.1, 0.5)) +
  labs(x="Blood Trait", y="Genetic Correlation (95% C.I.)") +
  coord_flip() +
  theme_minimal() 
  #scale_color_gradient2(low = "blue", mid = "gray", high = "red") + #, midpoint = 0.05
p_new <- p_new + geom_hline(aes(yintercept=0), color="black", linetype=4, linewidth=0.5)

# Save the plot
ggsave("SuppFig2b_ldak_genCorr.pdf", plot=p_new, height = 4, width = 6, units="in")