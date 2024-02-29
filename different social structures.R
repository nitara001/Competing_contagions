#sealions
sealions <- all_graphs[["sealion_proximity_weighted"]][[1]]

edge_weights <- get.edge.attribute(sealions, "weight")
weighted_sealions <- sealions
E(weighted_sealions)$weight <- edge_weights
plot(
  weighted_sealions,
  vertex.label=NA,
  main = "Sea lions"
)

#dolphins
dolphins<- all_graphs[["dolphin_association_weighted"]][[1]]
dolphin_edge_weights <- get.edge.attribute(dolphins, "weight")
weighted_dolphins <- dolphins
E(weighted_dolphins)$weight <- dolphin_edge_weights
plot(
  weighted_dolphins,
  layout = layout_nicely(weighted_dolphins),
  vertex.label= NA,
  main = "Dolphins")

#ants
ants<- all_graphs[["ants_proximity_weighted"]][[1]]
ants_edge_weights<- get.edge.attribute(ants, "weight")
weighted_ants <- ants
E(weighted_ants)$weight <- ants_edge_weights
plot(
  weighted_ants,
  layout = layout_nicely(weighted_ants),
  vertex.label= NA,
  main = "Ants")

#songbirds
songbirds<- all_graphs[["songbird_association_weighted"]][[1]]
songbirds_edge_weights<- get.edge.attribute(songbirds, "weight")
weighted_songbirds <- songbirds
E(weighted_songbirds)$weight <- songbirds_edge_weights
plot(
  weighted_songbirds,
  layout = layout_nicely(weighted_songbirds),
  vertex.label= NA,
  main = "songbirds")


#voles
voles<- all_graphs[["voles_social_projection_bipartite_weighted"]][[1]]
voles_edge_weights<- get.edge.attribute(voles, "weight")
weighted_voles <- voles
E(weighted_voles)$weight <- voles_edge_weights
plot(
  weighted_voles,
  layout = layout_nicely(weighted_voles),
  vertex.label= NA,
  main = "voles")

#rhesus mac
rhesus<- all_graphs[["rhesusmacaque_association_weighted"]][[1]]
rhesus_edge_weights<- get.edge.attribute(rhesus, "weight")
weighted_rhesus <- rhesus
E(weighted_rhesus)$weight <- rhesus_edge_weights
plot(
  weighted_rhesus,
  layout = layout_nicely(weighted_rhesus),
  vertex.label= NA,
  main = "rhesus macaque")

#sticklebacks
stickleback<- all_graphs[["fishstickleback_proximity_weighted"]][[1]]
stickleback_edge_weights<- get.edge.attribute(stickleback, "weight")
weighted_stickleback <- stickleback
E(weighted_stickleback)$weight <- stickleback_edge_weights
plot(
  weighted_stickleback,
  layout = layout_nicely(weighted_stickleback),
  vertex.label= NA,
  main = "Fish sticklebacks")

#cowbird
cowbird<- all_graphs[["cowbird_yokel_sexual_weighted"]][[1]]
cowbird_edge_weights<- get.edge.attribute(cowbird, "weight")
weighted_cowbird <- cowbird
E(weighted_cowbird)$weight <- cowbird_edge_weights
plot(
  weighted_cowbird,
  layout = layout_nicely(weighted_cowbird),
  vertex.label= NA,
  main = "cowbirds")

#baboon
baboon<- all_graphs[["baboon_association_weighted"]][[9]]
baboon_edge_weights<- get.edge.attribute(baboon, "weight")
weighted_baboon <- baboon
E(weighted_baboon)$weight <- baboon_edge_weights
plot(
  weighted_baboon,
  layout = layout_nicely(weighted_baboon),
  vertex.label= NA,
  main = "baboon
")

#ants
ants <- all_graphs[["ants_trophallaxis_weighted"]][[1]]

edge_weights <- get.edge.attribute(ants, "weight")
weighted_ants <- ants
E(weighted_ants)$weight <- edge_weights
plot(
  weighted_ants,
  vertex.label=NA,
  main = " weighted_ants"
)

#hens
hens <- all_graphs[["hens_dominance_weighted"]][[1]]




zebra <- all_graphs[["zebra_groupmembership_weighted"]][[1]]

edge_weights <- get.edge.attribute(zebra, "weight")
weighted_zebra <- zebra
E(weighted_zebra)$weight <- edge_weights
plot(
  weighted_zebra,
  vertex.label=NA,
  main = " weighted_zebra"
)
