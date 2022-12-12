co_expression_network <- function(DEG, tresholdCor, method = "pearson"){
    
    #compute the correlation matrix Genes x Genes
    cor.mat <- cor(t(DEG), method = method)
    # remove the self-loop
    diag(cor.mat) <- 0
    # |ρ| < threshold
    adj.mat1 <- ifelse(cor.mat >= tresholdCor, 1,
                         ifelse(cor.mat <= -tresholdCor,-1, 0 ))
    
    # correlation test 
    cor.padj <- corr.p(cor.mat, nrow(cor.mat), adjust="fdr",ci=FALSE)$p
    cor.padj[lower.tri(cor.padj)] <- t(cor.padj)[lower.tri(cor.padj)]
    
    # significantly adjustment  
    # if p.value < 0.05 we confirm correlation =/= 0 
    adj.mat2 <- ifelse(abs(cor.padj) > 0.05, 0, 1) 
    adj.mat <- adj.mat1 * adj.mat2
    
    return(list(cor.mat = as.matrix(cor.mat), adj.mat = as.matrix(adj.mat)))
}


getHubs <- function(net, measure = "degree"){
  if (measure == "degree") {
    degs <- igraph::degree(net, mode="all")
  }
  else if (measure == "betweenness") {
    degs <- igraph::betweenness(net, normalized = T)
  }
  else {
    stop("NO measure selected")
  }
  q_treshold = quantile(degs[degs>0], 0.95)
  hubs <- degs[degs>=q_treshold]
  hubs <- sort(hubs, decreasing = TRUE)
  
  return(hubs)
}


getName <- function(genes.raw.name){
  genes.clean.name <- (str_split(genes.raw.name, "[.]",simplify = TRUE))[,1]
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes.names <- getBM(filters= "ensembl_gene_id", 
                       attributes= c("ensembl_gene_id", "hgnc_symbol"),
                       values=genes.clean.name, mart= mart)
  return(genes.names)
}

getHubSubnetwork <- function(hubs.c, hubs.n, adj.mat.c, net.c){
  # find the common hubs and remove take only cancer related node
  common.hubs = intersect(names(hubs.c), names(hubs.n))
  condition = which(names(hubs.c) %in% setdiff(names(hubs.c), common.hubs))
  hubs.c.only = hubs.c[condition]
  
  # map name -> ids
  hubs.c.ids <- vector("integer", length(hubs.c.only))
  for (i in 1:length(hubs.c.ids)){hubs.c.ids[i] <- match(names(hubs.c.only)[i],
                                                         rownames(adj.mat.c))}
  #identifying the neighborhood
  hubs.c.neigh <- c()
  for (f in hubs.c.ids){
    hubs.c.neigh <- append(hubs.c.neigh, neighbors(net.c, f, mode = "all"))
  }
  
  hubs.c.neigh <- unique(hubs.c.neigh)
  hubs.c.neigh.names <- rownames(adj.mat.c[hubs.c.neigh,])
  subnet.name  <- unique(c(names(hubs.c), hubs.c.neigh.names))
  
  #creating the subnetwork
  hub.c.adj <- adj.mat.c[subnet.name, subnet.name]
  genes.name <- getName(rownames(hub.c.adj))$hgnc_symbol
  rownames(hub.c.adj) <- genes.name
  colnames(hub.c.adj) <- genes.name
  
  
  return(hub.c.adj)
}

getHubSubnetworkZ <- function(hubs.z, adj.mat.z, net.z){

  # map name -> ids
  hubs.z.ids <- vector("integer", length(hubs.z))
  for (i in 1:length(hubs.z.ids)){hubs.z.ids[i] <- match(names(hubs.z)[i],
                                                         rownames(adj.mat.z))}
  #identifying the neighborhood
  hubs.z.neigh <- c()
  for (f in hubs.z.ids){
    hubs.z.neigh <- append(hubs.z.neigh, neighbors(net.z, f, mode = "all"))
  }
  
  hubs.z.neigh <- unique(hubs.z.neigh)
  hubs.z.neigh.names <- rownames(adj.mat.z[hubs.z.neigh,])
  subnet.name  <- unique(c(names(hubs.z), hubs.z.neigh.names))
  
  #creating the subnetwork
  hub.z.adj <- adj.mat.z[subnet.name, subnet.name]
  genes.name <- getName(rownames(hub.z.adj))$hgnc_symbol
  rownames(hub.z.adj) <- genes.name
  colnames(hub.z.adj) <- genes.name
  
  return(hub.z.adj)
}




powerlaw_test <- function(data){
  # testing power law
  data_pl <- displ$new(data)
  est <- estimate_xmin(data_pl)
  data_pl$xmin <- est
  data_pl  
  
  fitting = fit_power_law(data, xmin = data_pl$xmin)
  return(list(KS.p = fitting$KS.p,
              KS.stats = fitting$KS.stat, 
              alpha = fitting$alpha))
  
}

getFishernet <- function(cor.mat.c, cor.mat.n, tresholdZ, n.patient){
  # fisher scores
  z.c = fisherz(cor.mat.c)
  z.n = fisherz(cor.mat.n)
  
  n.c = n.n = n.patient
  Z =  (z.c - z.n)/sqrt( (1/(n.c-3)+(1/(n.n-3))) )
  # build the adjacency matrix 
  adj.mat.Z = ifelse(abs(Z) < tresholdZ, 0, 1)
  diag(adj.mat.Z) <- 0
  
  # to graph
  net.Z = graph_from_adjacency_matrix(adj.mat.Z, mode = "undirected",
                                      weighted = TRUE)
  
  return(list(adj.mat = adj.mat.Z, net = net.Z))
}

testingVS <- function(data){
  data_pl <- displ$new(frequencies)
  est <- estimate_xmin(data_pl)
  data_pl$xmin <- est
  data_pl
  
  # testing power law vs Log-Normal()
  data_alt <- dislnorm$new(frequencies)
  data_alt$xmin <- est$xmin
  data_alt$pars <- estimate_pars(data_alt)
  LogNcomp <- compare_distributions(data_pl, data_alt)
  
  # testing power law vs Exp()
  data_alt <- disexp$new(frequencies)
  data_alt$xmin <- est$xmin
  data_alt$pars <- estimate_pars(data_alt)
  Expcomp <- compare_distributions(data_pl, data_alt)
  
  # testing power law vs Pois()
  data_alt <- dispois$new(frequencies)
  data_alt$xmin <- est$xmin
  data_alt$pars <- estimate_pars(data_alt)
  Poiscomp <- compare_distributions(data_pl, data_alt)
  
  list(LogNormal = list(test_stats = LogNcomp$test_statistic, 
                        p_value = LogNcomp$p_two_sided),
       Exp = list(test_stats = Expcomp$test_statistic,
                  p_value = Expcomp$p_two_sided),
       Pois = list(test_stats = Poiscomp$test_statistic,
                   p_value = Poiscomp$p_two_sided)
  )
}


applyHardThreshold <- function(cor.mat, adj.mat, tau=0.5) {
  idx = as.vector(which(cor.mat >= tau))
  hard.adj.mat = adj.mat
  hard.adj.mat[idx] = 1
  hard.adj.mat[-idx] = 0
  return(as.matrix(hard.adj.mat))
}

sigmoid <- function(alpha, tau0, x) {
  exponent = -alpha*(x-tau0)
  res = 1/(1-exp(exponent))
  return(res)
}

applySoftThreshold <- function(cor.mat, adj.mat, alpha=0.1, tau0=0.1) { #
  soft.adj.mat = adj.mat
  apply(soft.adj.mat, 1:2, function(x) sigmoid(alpha, tau0,x))
  return(as.matrix(soft.adj.mat))
}

normalize <- function(data){
  normalized = (data-min(data))/(max(data)-min(data))
  return(normalized)
}

getInfoNetworkTopology <- function(g) {
  
  infos = matrix(NA, nrow = 1, ncol = 10)
  
  w = abs(E(g)$weight)
  
  infos[1] = igraph::global_efficiency(g, weights = w)
  infos[2] = igraph::average_local_efficiency(g, weights = w, mode = "in")
  infos[3] = igraph::average_local_efficiency(g, weights = w, mode = "out")
  infos[4] = igraph::diameter(g, weights = w)
  infos[5] = max(igraph::degree(g, normalized = T))
  #n of connected components
  infos[6] = igraph::count_components (g, mode = c("weak","strong"))
  #mean Subgraph centrality of a vertex measures the number of subgraphs a vertex participates in, 
  #weighting them according to their size.
  clos_array = igraph::closeness(g, weights = w)[!is.na(igraph::closeness(g, weights = w))]
  infos[7] = mean(normalize(clos_array))
  infos[8] = igraph::centr_betw(g, directed=F, normalized=T)$centralization
  infos[9] = mean(igraph::hub_score(g, weights = w, scale = TRUE)$vector)
  infos[10] = mean(igraph::authority_score(g, weights = w, scale = TRUE)$vector)

  return(infos)
  
}

getCO.EXPR <- function(method, DEG.c, DEG.n){
  tresholdCor = 0.65
  # CANCER ADJ MATRIX
  # get adj matrix
  data = co_expression_network(DEG.c, tresholdCor, method = method)
  cor.mat.c = data$cor.mat
  adj.mat.c = data$adj.mat
  
  # NORMAL ADJ MATRIX
  # get adj matrix
  data = co_expression_network(DEG.n, tresholdCor, method = method)
  cor.mat.n = data$cor.mat
  adj.mat.n = data$adj.mat
  
  # Build the networks
  net.c <- graph_from_adjacency_matrix(adj.mat.c, mode = "undirected", 
                                       weighted = TRUE)
  net.n <- graph_from_adjacency_matrix(adj.mat.n, mode = "undirected", 
                                       weighted = TRUE)
  
  # Find the hubs (5% of the nodes with highest degree values)
  # sorted by degree
  hubs.c = getHubs(net.c) 
  hubs.n = getHubs(net.n) 
  
  # graph nodes w/ hubs 
  hubs.c.genes.name <- getName(names(hubs.c))
  hubs.n.genes.name <- getName(names(hubs.n))
  
  l.c.hubs = length(hubs.c.genes.name$hgnc_symbol)
  l.n.hubs = length(hubs.n.genes.name$hgnc_symbol)
  
  # compare hubs sets related to the two condition (cancer, normal)
  common.hubs = getName(intersect(names(hubs.c), names(hubs.n)))$hgnc_symbol
  
  return(list(common.hubs = common.hubs, l.c.hubs = l.c.hubs, 
              l.n.hubs = l.n.hubs))
}
