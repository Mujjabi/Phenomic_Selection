
PCA <- data <- read.table("PCA_TASSEL.txt", header = TRUE, sep = ",")  
write.csv(PCA, "PCA.csv", row.names = FALSE)

PCA2  <- read.csv("PCA.csv", header = TRUE, stringsAsFactors = FALSE)

# Load necessary library
library(ggplot2)
library(viridis)
library(RColorBrewer)



custom_colors <- c("red", "blue",  "green",  "yellow",  "magenta",  "cyan",  "orange",  
                   "purple",  "brown",  "gray",  "black",  "pink",  "saddlebrown",  "chartreuse",  
                   "olivedrab",  "steelblue",  "chocolate",  "seagreen",  "gold", "antiquewhite", "orchid4")
  



# If your TASSEL data has 'Sample', 'PC1', 'PC2', and 'Group' columns:
highlight_sample <- "MOB318"  # Replace with the name of the sample you want to highlight

p <- ggplot(PCA2, aes(x = PC1, y = PC2, color = Group, label = Taxa)) +
  geom_point(size = 4, alpha = 0.8) +  # Plot all points
  geom_text(vjust = -0.5, size = 3, show.legend = FALSE) +  # Add labels
  geom_point(data = PCA2[PCA2$Taxa == highlight_sample, ], 
             aes(x = PC1, y = PC2), 
             color = "red",  # Highlight color
             size = 10, 
             shape = 1,   # Circle shape
             stroke = 1.5) +  # Line thickness for the circle
  labs(
    title = "PCA Plot",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = custom_colors)

print(p)




group_centroids <- PCA2 %>%
  group_by(Group) %>%
  summarise(
    PC1_mean = mean(PC1),
    PC2_mean = mean(PC2)
  )

p <- ggplot(PCA2, aes(x = PC1, y = PC2, color = Group, label = Taxa)) +
  geom_point(size = 4, alpha = 0.8) +  # Plot all points
  geom_text(vjust = -0.5, size = 3, show.legend = FALSE) +  # Add labels for samples
  # Add group centroids and labels
  geom_text(data = group_centroids, 
            aes(x = PC1_mean, y = PC2_mean, label = Group), 
            color = "black", size = 5, fontface = "bold") +
  labs(
    title = "PCA Plot",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = custom_colors)

print(p)


p <- ggplot(PCA2, aes(x = PC1, y = PC2, color = Group, label = Taxa)) +
  geom_point(size = 4, alpha = 0.8) +  # Plot all points
  geom_text(vjust = -0.5, size = 2, show.legend = FALSE) +  # Add labels for samples
  # Add group centroids and labels
  #geom_text(data = group_centroids, 
           #aes(x = PC1_mean, y = PC2_mean, label = Group), 
           #color = "black", size = 4, fontface = "bold") +
  labs(
    title = "PCA Plot",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",   # Move the legend to the left
    legend.title = element_blank(),  # Optional: Remove legend title
    legend.text = element_text(size = 12),  # Adjust the size of the legend text
    legend.key.size = unit(0.5, "cm"),  # Adjust the size of the legend keys
    legend.spacing.y = unit(0.5, "cm"),  # Adjust the space between legend items vertically
    legend.box.margin = margin(0, 10, 0, 0)  # Add margin to the right of the legend for spacing
  ) +
  scale_color_manual(values = custom_colors) +
  guides(color = guide_legend(ncol = 1))  # Ensure the legend is in a single column

print(p)





p <- ggplot(PCA2, aes(x = PC1, y = PC2, color = Group, label = Taxa)) +
  geom_point(size = 4, alpha = 0.8) +  # Plot all points
  geom_text(vjust = -0.5, size = 2, aes(color = Group), show.legend = FALSE) +  # Add labels with group color
  # Add group centroids and labels
  geom_text(data = group_centroids, 
            aes(x = PC1_mean, y = PC2_mean, label = Group, color = Group), 
            size = 5, fontface = "bold") +
  labs(
    title = "PCA Plot",
    x = "PC1",
    y = "PC2"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",   # Remove the legend
    legend.title = element_blank(),  # Optional: Remove legend title if you don't need it
    legend.text = element_text(size = 12),  # Adjust the size of the legend text (though it won't show now)
    legend.key.size = unit(0.5, "cm"),  # Adjust the size of the legend keys
    legend.spacing.y = unit(0.5, "cm"),  # Adjust the space between legend items vertically
    legend.box.margin = margin(0, 10, 0, 0)  # Add margin to the right of the legend for spacing
  ) +
  scale_color_manual(values = custom_colors)  # Ensure the colors match the custom colors for the groups

print(p)







