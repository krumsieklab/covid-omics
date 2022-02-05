# figure 3
library(igraph)
library(magrittr)
load(glue::glue("{root}/2c_networks/results/igraph_networks.Rdata"))

plot_networks <- function(){
  par(mfrow = c(1,2))
  # plot overall network 
  plot(
    net1, 
    edge.color = adjustcolor( "black", alpha.f = 0.23),
    mark.col = adjustcolor( "black", alpha.f = 0.05),
    mark.border = NA,
    vertex.border = NA,
    vertex.shape = (shapes()[c(1,9)])[as.numeric(V(net1)$shape)],
    layout = LO_net1,
    vertex.color = V(net1)$color,
    vertex.label =NA, 
    vertex.label.cex = 0.7,
    vertex.label.color = "black",
    vertex.size = 2
  )
  
  # plot annotated subnetwork 
  plot(net2, 
       edge.color = adjustcolor( "black", alpha.f = 0.7),
       vertex.border = NA,
       vertex.shape = (shapes()[c(1,9)])[V(net2)$shape],
       vertex.color = V(net2)$color,
       mark.groups = cluster_louvain(net2)$membership %>% 
         {lapply(seq(max(.)), function(i) which(.==i))},  
       mark.col = adjustcolor( "black", alpha.f = 0.05),
       mark.border = NA,
       layout = LO_net2,
       vertex.label.cex = 0.1,
       vertex.label.color = "black",
       vertex.size=2 )
}

# save the plottable networks
save(file = "figure3.fig", plot_networks, net1, net2, LO_net1, LO_net2)

