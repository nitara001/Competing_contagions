 library(igraph)
#  Detect communities
g= baboon
communities <- igraph::cluster_louvain(g)
membership <- membership(communities)

# Preserve Modular Configuration
for (i in 1:length(communities)) {
  vertices_in_community <- which(membership == i)
  # Use combn to generate all possible combinations of vertices within the community
  all_vertex_pairs <- combn(vertices_in_community, 2)
  edges_within_community <- vector("list", length = ncol(all_vertex_pairs))
  for (j in 1:ncol(all_vertex_pairs)) {
    edges_within_community[[j]] <- E(g)[from(all_vertex_pairs[1, j]) %--% to(all_vertex_pairs[2, j])]
  }
  edges_within_community_combined <- do.call(c, edges_within_community)
  # Store the edges within the community 
  edges_within_communities[[i]] <- edges_within_community_combined
}

# Randomise Between-Subgroup Connections
modular_null_network <- g
for (i in 1:length(communities)) {
  vertices_in_community <- which(membership == i)
  vertices_outside_community <- setdiff(V(g), vertices_in_community)
  
  # Use combn to generate all possible combinations of nodes between the current community and outside nodes
  all_vertex_pairs <- combn(vertices_in_community, vertices_outside_community, simplify = FALSE)
  edges_between_community_outside <- vector("list", length = length(all_vertex_pairs))
  
  for (j in 1:length(all_vertex_pairs)) {
    edges_between_community_outside[[j]] <- E(modular_null_network)[from(all_vertex_pairs[[j]][1]) %--% all_vertex_pairs[[j]][2]]
  }
  edges_between_community_outside_combined <- do.call(c, edges_between_community_outside)
  edges_between_community_outside_randomized <- sample(edges_between_community_outside_combined)
  
  # Delete existing edges between the current community and outside nodes
  delete_edges(modular_null_network, edges_between_community_outside_combined)
  
  # Add randomized edges between the current community and outside nodes
  add_edges(modular_null_network, edges_between_community_outside_randomized)
}


calculate_modularity_and_avg_module_size(modular_null_network)
calculate_modularity_and_avg_module_size(baboon)

