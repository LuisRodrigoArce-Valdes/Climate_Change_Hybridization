rm(list = ls())
# With this script we will do a bar plot of the hybridization cases per year
library(ggplot2)

# Reading and formating
Years <- read.delim("../data/Years.txt", header = T)
Years[,-c(1,2)] -> Years

# Tableing
Years <- as.data.frame(table(Years$Year, Years$Type))
Years$Var1 <- as.numeric(as.character(Years$Var1))
Years$Var2 <- factor(Years$Var2, levels = unique(as.character(Years$Var2)),
                     labels = c("Climate-induced hybridization between\nallopatrically diverged species",
                                "Predicted range overlap and hybridization\nbetween allopatrically diverged species\ndue to climate change",
                                "Altered natural hybrid zones\ndue to climate change"))

# Plotting
png("../results/06_Bars.png", width = 32, height = 16, units = "cm", res = 300)
ggplot(Years) +
  geom_col(aes(x=Var1, y=Freq, fill=Var2), color="black") +
  scale_x_continuous(breaks = 2002:2021) +
  scale_fill_manual(values = c('#7fc97f','#beaed4','#fdc086')) +
  labs(y="Frequency") +
  theme_classic() +
  theme(axis.title.x = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank(),
        text = element_text(family = "serif", size = 18),
        legend.text = element_text(size = 12),
        panel.grid.major.y = element_line(color = "grey"))
dev.off()
