 library(igraph)
#  Detect communities
g= hens
communities <- cluster_edge_betweenness(g)
membership <- membership(communities)

# Preserve Modular Configuration
edges_within_communities <- vector(mode = "list", length = length(communities))

for (i in 1:length(communities)) {
  vertices_in_community <- which(membership == i)
  edges_within_community <- E(g)[from(vertices_in_community) %--% vertices_in_community]
  edges_within_communities[[i]] <- edges_within_community
}

# Randomise Between-Subgroup Connections
modular_null_network <- g

for (i in 1:length(communities)) {
  vertices_in_community <- which(membership == i)
  vertices_outside_community <- setdiff(V(g), vertices_in_community)
  
  # Get the edges between the current community and outside vertices
  edges_between_community_outside <- E(modular_null_network)[from(vertices_in_community) %--% vertices_outside_community]
num_edges_to_sample <- floor(length(edges_between_community_outside) / 2) * 2

# Randomie the order of edges between the current community and outside vertices
edges_between_community_outside_randomized <- sample(edges_between_community_outside)

  # Delete existing edges between the current community and outside vertices
  delete_edges(modular_null_network, edges_between_community_outside)
  
  # Add randomised edges between the current community and outside vertices
  add_edges(modular_null_network, edges_between_community_outside_randomized)
}

calculate_modularity_and_avg_module_size(modular_null_network)
calculate_modularity_and_avg_module_size(hens)

