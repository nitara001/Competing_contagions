# Load the igraph library
library(igraph)

# Define the top-level folder containing the subfolders
top_level_dir <- "C:/Users/s2607536/OneDrive - University of Edinburgh/Data/asnr-master/asnr-master/Networks"

# List all subfolders and sub-subfolders
all_subfolders <- list.dirs(top_level_dir, full.names = TRUE, recursive = TRUE)

# Create a list to store subfolder names and their respective graphs
all_graphs <- list()

# Loop through each subdirectory
for (sub in all_subfolders) {
  subfolder_name <- basename(sub)
  graphml_files <- list.files(sub, pattern = ".graphml", full.names = TRUE)
  subfolder_graphs <- list()
  
  for (graphml_file in graphml_files) {
    graph <- read_graph(graphml_file, format = "graphml")
    subfolder_graphs <- c(subfolder_graphs, list(graph))
  }
  
  all_graphs[[subfolder_name]] <- subfolder_graphs
}

