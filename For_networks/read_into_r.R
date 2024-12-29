# Load the igraph library
library(igraph)

# Define the top-level folder containing the subfolders
top_level_dir <- "C:\\Users\\s2607536\\Documents\\GitHub\\Competing_contagions\\Networks"


# List all subfolders and sub-subfolders
all_subfolders <- list.dirs(top_level_dir, full.names = TRUE, recursive = TRUE)
all_subfolders

# Create a list to store subfolder names and their respective graphs
all_graphs <- list()

# Loop through each subdirectory
for (sub in all_subfolders) {
  subfolder_name <- basename(sub)
  graphml_files <- list.files(sub, pattern = ".graphml", full.names = TRUE)
  subfolder_graphs <- list()
  
  for (graphml_file in graphml_files) {
    # Try to load the graph, and handle any errors
    graph <- tryCatch({
      read_graph(graphml_file, format = "graphml")
    }, error = function(e) {
      message(paste("Error reading:", graphml_file))
      return(NULL)  # Return NULL if there’s an error
    })
    
    # Only add non-null graphs to the list
    if (!is.null(graph)) {
      subfolder_graphs <- c(subfolder_graphs, list(graph))
    }
  }
  
  all_graphs[[subfolder_name]] <- subfolder_graphs
}

# Extract networks from all_graphs and apply filtering
extracted_networks <- list()
for (i in seq_along(all_graphs)) {
  for (n in seq_along(all_graphs[[i]])) {
    graph <- all_graphs[[i]][[n]]
    if (igraph::vcount(graph) > 10 && igraph::is_connected(graph)) {
      # Add it to the extracted_networks list only if the graph is connected
      extracted_networks[[paste0(names(all_graphs)[i], "_", n)]] <- graph
    }
  }
}

# Check the extracted networks
extracted_networks
