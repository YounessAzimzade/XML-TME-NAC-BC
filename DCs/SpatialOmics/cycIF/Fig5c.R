# Fig5c
# Author: Youness Azimzade
# Email: younessazimzade@gmail.com
# Date: 10/08/2023


# Load necessary libraries for data visualization and manipulation
library(ggplot2)
library(reshape2)

# Read the CSV file containing distance and frequency data
Distance <- read.csv("~/CumFreq.csv")

# Convert the 'Radius' column to character type
Distance$Radius <- as.character(Distance$Radius)

# Melt the data to long format for better visualization
Distance.m <- melt(Distance)

# Convert 'Radius' and 'variable' columns to appropriate numeric and character types
Distance.m$Radius <- as.numeric(Distance.m$Radius)
Distance.m$variable <- as.character(Distance.m$variable)

# Rename specific categories in the 'variable' column for clarity
Distance.m$variable[Distance.m$variable == "ERn58Epi"] <- "ER-"
Distance.m$variable[Distance.m$variable == "ERn86Epi"] <- "ER-"
Distance.m$variable[Distance.m$variable == "ERp3019Epi"] <- "ER+"
Distance.m$variable[Distance.m$variable == "ERp73Epi"] <- "ER+"

# Rename column names for improved clarity in the plot
colnames(Distance.m) <- c("Distance from DC (micrometer)", "Sample", "Averaged number of ECs")

# Generate a PDF plot for visualization
pdf("~/Fig5c.pdf", width = 10, height = 7)

# Create a scatter plot using ggplot2
ggplot(Distance.m, aes(x = `Distance from DC (micrometer)`, y = `Averaged number of ECs`, color = Sample)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(values = c("blue", "red")) +
  theme(text = element_text(size = 30)) +
  xlim(0, 200) +
  ylim(0, 200) +
  labs(x = "Distance from DC (micrometer)", y = "Averaged number of ECs")

# Save the PDF plot and close the PDF device
dev.off()
