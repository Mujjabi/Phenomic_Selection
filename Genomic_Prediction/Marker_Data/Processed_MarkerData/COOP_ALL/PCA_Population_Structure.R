
PCA <- data <- read.table("PCA_COOP.txt", header = TRUE, sep = ",")  
write.csv(PCA, "PCA_COOP.csv", row.names = FALSE)

PCA2  <- read.csv("PCA_COOP_Groups.csv", header = TRUE, stringsAsFactors = FALSE)

# Load necessary library
library(ggplot2)
library(viridis)
library(RColorBrewer)



## Structure across breeding programs

custom_colors <- c("red", "blue",  "green",  "yellow")

# Define the samples to highlight and their custom names
highlight_sample <- c("MOB318", "PS030", "MOB752", "MOB323")  # Samples to highlight
custom_labels    <- c("B73", "A427", "PHG35", "Mo17")  # Custom names for the highlighted samples

# Create a subset of the data for the highlighted samples
highlight_data <- PCA2[PCA2$Taxa %in% highlight_sample, ]
highlight_data$CustomLabel <- custom_labels[match(highlight_data$Taxa, highlight_sample)]

# Define colors for the highlighted points
highlight_colors <- c("PHG35" = "steelblue", "B73" = "chocolate", "PS030" = "seagreen", "Mo17" = "blue")

p <- ggplot(PCA2, aes(x = PC1, y = PC2, color = Program)) +
  geom_point(size = 4, alpha = 0.8) +  # Plot all points
  geom_point(data = highlight_data, 
             aes(x = PC1, y = PC2),  # Remove color mapping for Taxa to exclude from legend
             size = 10, 
             shape = 1,   # Circle shape
             stroke = 1.5,
             color = "black",  # Use a fixed color to avoid legend entry
             show.legend = FALSE) +  # Remove from legend
  geom_text(data = highlight_data, 
            aes(x = PC1, y = PC2, label = CustomLabel), 
            vjust = -1.5, size = 5, color = "black") +  # Add custom labels
  labs(
    title = "PCA Plot",
    x = "PC1 (8.36%)",
    y = "PC2 (5.35%)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = custom_colors)  # Maintain custom colors for Program

print(p)


### 3 dimensional wtih 3 PCs

library(plotly)

# Create a 3D scatter plot
p <- plot_ly(PCA2, 
             x = ~PC1, y = ~PC2, z = ~PC3, 
             color = ~Program, colors = custom_colors,
             type = "scatter3d", mode = "markers",
             marker = list(size = 4, opacity = 0.8)) %>%
  
  # Add highlighted points with black circles
  add_trace(data = highlight_data, 
            x = ~PC1, y = ~PC2, z = ~PC3,
            type = "scatter3d", mode = "markers",
            marker = list(size = 10, color = "rgba(0,0,0,0)", 
                          line = list(color = "black", width = 2)),
            showlegend = FALSE) %>%
  
  # Add labels for highlighted points
  add_trace(data = highlight_data, 
            x = ~PC1, y = ~PC2, z = ~PC3, 
            type = "scatter3d", mode = "text",
            text = ~CustomLabel, textposition = "top center",
            textfont = list(size = 14, color = "black"),
            showlegend = FALSE) %>%
  
  layout(
    title = "Population Structure within COOP Materials",
    scene = list(
      xaxis = list(title = "PC1 (8.36%)"),
      yaxis = list(title = "PC2 (5.35%)"),
      zaxis = list(title = "PC3 (3.13%)")
    )
  )

p



## We can now look at subgroups within each breeding program




custom_colors2 <- c("#ff7f00",  "yellow", "green1", "magenta", 
                   "cyan" ,"#6A3D9A", "red4", "#1f78b4","black","#33a02c",
                   "#fb9a99","red","blue","deeppink1",
                   "purple1","yellow4","#b15928")


c25 <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black", "gold1",
  "skyblue2", "#FB9A99", # lt pink
  "palegreen2",
  "#CAB2D6", # lt purple
  "#FDBF6F", # lt orange
  "gray70", "khaki2",
  "maroon", "orchid1", "deeppink1", "blue1", "steelblue4",
  "darkturquoise", "green1", "yellow4", "yellow3",
  "darkorange4", "brown"
)



# Define the samples to highlight and their custom names
highlight_sample <- c("MOB318", "PS030", "MOB752", "MOB323")  # Samples to highlight
custom_labels    <- c("B73", "A427", "PHG35", "Mo17")  # Custom names for the highlighted samples

# Create a subset of the data for the highlighted samples
highlight_data <- PCA2[PCA2$Taxa %in% highlight_sample, ]
highlight_data$CustomLabel <- custom_labels[match(highlight_data$Taxa, highlight_sample)]

# Define colors for the highlighted points
highlight_colors <- c("PHG35" = "steelblue", "B73" = "chocolate", "PS030" = "seagreen", "Mo17" = "blue")

p <- ggplot(PCA2, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 4, alpha = 0.8) +  # Plot all points
  geom_point(data = highlight_data, 
             aes(x = PC1, y = PC2),  # Remove color mapping for Taxa to exclude from legend
             size = 10, 
             shape = 1,   # Circle shape
             stroke = 1.5,
             color = "black",  # Use a fixed color to avoid legend entry
             show.legend = FALSE) +  # Remove from legend
  geom_text(data = highlight_data, 
            aes(x = PC1, y = PC2, label = CustomLabel), 
            vjust = -1.5, size = 5, color = "black") +  # Add custom labels
  labs(
    title = "PCA Plot",
    x = "PC1 (8.36%)",
    y = "PC2 (5.35%)"
  ) +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  ) +
  scale_color_manual(values = custom_colors2)  # Maintain custom colors for Program

print(p)







