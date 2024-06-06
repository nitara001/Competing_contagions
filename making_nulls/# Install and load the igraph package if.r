 library(igraph)
# Step 1: Detect communities using edge-betweenness method
g= hens
communities <- cluster_edge_betweenness(g)

# Get membership vector
membership <- membership(communities)

# Step 2: Preserve Modular Configuration
# Create a list to store the edges within each community
edges_within_communities <- vector(mode = "list", length = length(communities))

# Iterate through each community
for (i in 1:length(communities)) {
  # Get the vertices in the community
  vertices_in_community <- which(membership == i)
  # Extract the edges within the community
  edges_within_community <- E(g)[from(vertices_in_community) %--% vertices_in_community]
  # Store the edges within the community
  edges_within_communities[[i]] <- edges_within_community
}

# Step 3: Randomize Between-Subgroup Connections
# Create a copy of the original network
modular_null_network <- g

for (i in 1:length(communities)) {
  # Get the vertices in the current community
  vertices_in_community <- which(membership == i)
  
  # Get the vertices outside the current community
  vertices_outside_community <- setdiff(V(g), vertices_in_community)
  
  # Get the edges between the current community and outside vertices
  edges_between_community_outside <- E(modular_null_network)[from(vertices_in_community) %--% vertices_outside_community]
  
# Round down the length of edges_between_community_outside to the nearest even number
num_edges_to_sample <- floor(length(edges_between_community_outside) / 2) * 2

# Randomize the order of edges between the current community and outside vertices
edges_between_community_outside_randomized <- sample(edges_between_community_outside)

  # Delete existing edges between the current community and outside vertices
  delete_edges(modular_null_network, edges_between_community_outside)
  
  # Add randomized edges between the current community and outside vertices
  add_edges(modular_null_network, edges_between_community_outside_randomized)
}

calculate_modularity_and_avg_module_size(modular_null_network)
calculate_modularity_and_avg_module_size(dolphins)



library(igraph)

# Step 1: Detect communities using Louvain method
modular_null_networks <- list()

# List of network names in the sampled_networks list
network_names <- names(sampled_networks)

for (network_name in network_names) {
  # Load the current network
  g <- sampled_networks[[network_name]]
  communities <- cluster_louvain(g)
  # Get membership vector
  membership <- membership(communities)

  # Step 2: Preserve Modular Configuration
  # Create a list to store the edges within each community
  edges_within_communities <- vector(mode = "list", length = length(communities))

  # Iterate through each community
  for (i in 1:length(communities)) {
    # Get the vertices in the community
    vertices_in_community <- which(membership == i)
    # Extract the edges within the community
    edges_within_community <- E(g)[from(vertices_in_community) %--% vertices_in_community]
    # Store the edges within the community
    edges_within_communities[[i]] <- edges_within_community
  }

  # Step 3: Randomize Between-Subgroup Connections
  # Create a copy of the original network
  modular_null_network <- g

  for (i in 1:length(communities)) {
    # Get the vertices in the current community
    vertices_in_community <- which(membership == i)

    # Get the vertices outside the current community
    vertices_outside_community <- setdiff(V(g), vertices_in_community)

    # Get the edges between the current community and outside vertices
    edges_between_community_outside <- E(modular_null_network)[from(vertices_in_community) %--% vertices_outside_community]

    # Round down the length of edges_between_community_outside to the nearest even number
    num_edges_to_sample <- floor(length(edges_between_community_outside) / 2) * 2

    # Randomize the order of edges between the current community and outside vertices
    edges_between_community_outside_randomized <- sample(edges_between_community_outside, size = num_edges_to_sample, replace = TRUE)

    # Delete existing edges between the current community and outside vertices
    delete_edges(modular_null_network, edges_between_community_outside)

    # Add randomized edges between the current community and outside vertices
    add_edges(modular_null_network, edges_between_community_outside_randomized)
  }

  # Store the modular null network with its name
  modular_null_networks[[network_name]] <- modular_null_network
}
