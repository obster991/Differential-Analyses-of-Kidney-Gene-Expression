# DIGITAL EPIDEMIOLOGY AND PRECISION MEDICINE Homework

library(ggplot2);  
library(stats); library(psych); library(stringr)
library(igraph); library(GGally); library(intergraph); 
library(sna); library(scales); library(lsa)
library('biomaRt')
library(poweRlaw)
library(kableExtra);

getwd()

source("code_clean/graph_functions.R")
source("code_clean/networks_functions.R")

# create figure directory
dir.create(file.path(getwd(), "figure"), showWarnings = FALSE)

# 1. Download data ---------------------------------------------------
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DT)

proj <- "TCGA-KICH" 
dir.create(file.path(proj))

rna.query.C <- GDCquery(project = proj, data.category = "Transcriptome Profiling",
                        data.type = "Gene Expression Quantification",
                        workflow.type = "STAR - Counts",
                        sample.type = "Primary Tumor")

GDCdownload(query = rna.query.C, directory = "GDCdata", method = "api")
rna.data.C <- GDCprepare(rna.query.C)
rna.expr.data.C <- assay(rna.data.C)
norm.data.C <- rna.data.C@assays@data@listData[["fpkm_unstrand"]]
rownames(norm.data.C) <- rownames(rna.expr.data.C)
colnames(norm.data.C) <- colnames(rna.expr.data.C)
write.csv(norm.data.C, "rna_expr_data_C.csv",
          sep = ",", row.names=TRUE, col.names=TRUE, quote = FALSE)



rna.query.N <- GDCquery(project = proj, data.category = "Transcriptome Profiling", 
                        data.type = "Gene Expression Quantification", 
                        workflow.type = "STAR - Counts", 
                        sample.type = "Solid Tissue Normal")

GDCdownload(query = rna.query.N, directory = "GDCdata", method = "api")
rna.data.N <- GDCprepare(rna.query.N)
rna.expr.data.N <- assay(rna.data.N)
norm.data.N <- rna.data.N@assays@data@listData[["fpkm_unstrand"]]
rownames(norm.data.N) <- rownames(rna.expr.data.N)
colnames(norm.data.N) <- colnames(rna.expr.data.N)
write.csv(norm.data.N, "rna_expr_data_N.csv",sep = ",", row.names=TRUE, col.names=TRUE, quote = FALSE)

#clinical info
clinical.query <- GDCquery_clinic(project = proj, type = "clinical", save.csv = FALSE)
write.csv(clinical.query, file = file.path(proj,paste(proj, "_clinical_data.txt",sep="")), row.names = FALSE, quote = FALSE)


# 2. DEGs -------------------------------------------------------------

# load checkpoint data
load("data_R/data.RData")

# Get log2FC and adjusted-pvals
data = getFC_pvals(filtr.expr.c, filtr.expr.n, method = "fdr")
log2FC = data$log2FC
pvals_adj = data$pvals_adj

tresholdFC = 3.0

# Get the subset of the data frame C and N
tresholdPval = 0.001
condition = abs(log2FC) >= tresholdFC & pvals_adj <= tresholdPval
mean(condition)  # keep 2% of genes

mask = which(condition)
length(mask)

# get the subset of hundreds of genes
DEG.c = filtr.expr.c[mask,]
DEG.n = filtr.expr.n[mask,]

# get volcano plot
getFCHistogram(log2FC, tresholdFC)
getVolcanoPlot(log2FC, pvals_adj, tresholdFC, tresholdPval)

# 3: Co-expression network ---------------------------------------

# **Computation**. Using only DEGs, compute the gene co-expression networks 
# related to the 2 conditions (cancer, normal) considering:
#   1. Pearson’s correlation (or another measure of similarity);
#   2. Binary adjacency matrix where aij=0 if |ρ|< threshold.


tresholdCor = 0.65
# CANCER ADJ MATRIX
# get adj matrix
data = co_expression_network(DEG.c, tresholdCor, method = "pearson")
cor.mat.c = data$cor.mat
adj.mat.c = data$adj.mat

# NORMAL ADJ MATRIX
# get adj matrix
data = co_expression_network(DEG.n, tresholdCor, method = "pearson")
cor.mat.n = data$cor.mat
adj.mat.n = data$adj.mat

# get corr histogram
getCorrHistogram(cor.mat.n, tresholdCor, output_file = "figure/Corr_histogram_Normal.pdf")
getCorrHistogram(cor.mat.c, tresholdCor, output_file = "figure/Corr_histogram_Cancer.pdf")



# Build the networks
net.c <- graph_from_adjacency_matrix(adj.mat.c, mode = "undirected", 
                                    weighted = TRUE)
net.n <- graph_from_adjacency_matrix(adj.mat.n, mode = "undirected", 
                                     weighted = TRUE)

ggnet2(net.c, size = 1)
edge_density(simplify(net.c), loops = FALSE)

# analysis on graph: CANCER
deg.c <- igraph::degree(net.c)
deg.c <- deg.c[deg.c > 0]
C.frequencies <- as.data.frame(table(deg.c))[, 2]

getDegreedistribution(net.c, title = "Cancer", output_file = "figure/C_deg_distr.pdf")
getHistoDistro(C.frequencies, title = "Cancer", output_file = "figure/C_histo.pdf")

# Testing Fit Power Law
powerlaw_test(C.frequencies)

# analysis on graph: NORMAL
deg.n <- igraph::degree(net.n)
N.frequencies <- as.data.frame(table(deg.n))[, 2]

getDegreedistribution(net.n, title = "Normal", output_file = "figure/N_deg_distr.pdf")
getHistoDistro(N.frequencies, title = "Cancer", bw = 2, output_file = "figure/N_histo.pdf")

# Fit Power Law
powerlaw_test(N.frequencies)

# Find the hubs (5% of the nodes with highest degree values)
# sorted by degree
hubs.c = getHubs(net.c) 
hubs.n = getHubs(net.n) 

# graph nodes w/ hubs 
hubs.c.genes.name <- getName(names(hubs.c))
hubs.n.genes.name <- getName(names(hubs.n))


length(hubs.c.genes.name$hgnc_symbol)
length(hubs.n.genes.name$hgnc_symbol)

hubs.test = sort(hub_score(net.c)$vector, decreasing = T)
hubs.test = hubs.test[hubs.test > 0]
names.hubs.test = names(hubs.test[1:15])
getName(names.hubs.test)

# visualize hubs name as list
cat(paste(hubs.c.genes.name[,2], collapse='\n' ) )
cat(paste(hubs.n.genes.name[,2], collapse='\n' ) )

# 3.3: compare hubs sets related to the two condition (cancer, normal) and 
# identify the hubs characterizing only cancer tissue

# compare hubs sets related to the two condition (cancer, normal)
common.hubs = getName(intersect(names(hubs.c), names(hubs.n)))
common.hubs$hgnc_symbol

hub.c.adj <- getHubSubnetwork(hubs.c, hubs.n, adj.mat.c, net.c)
net.hubs.c <- graph_from_adjacency_matrix(hub.c.adj, mode = "undirected", 
                                            weighted = TRUE)

suppressWarnings(getSubnetworkGraph(net.hubs.c, hubs.c.genes.name))

## Info on hubs 
common.hubs = intersect(names(hubs.c), names(hubs.n))
condition = which(names(hubs.c) %in% setdiff(names(hubs.c), common.hubs))
hubs.c.only = hubs.c[condition]
genes_hubs.c.only = getName(names(hubs.c.only))

# high betweenness genes
betweenness_c_hubs = igraph::betweenness(net.hubs.c, directed=F, weights=NA)
sort(betweenness_c_hubs, decreasing = T)
names_bet = names(sort(betweenness_c_hubs, decreasing = T))
head(names_bet[names_bet %in% genes_hubs.c.only$hgnc_symbol])
# high centrality 
closeness_c_hubs = igraph::closeness(net.hubs.c, mode="all", weights=NA)
sort(closeness_c_hubs, decreasing = T)
names_clos = names(sort(closeness_c_hubs, decreasing = T))
head(names_clos[names_clos %in% genes_hubs.c.only$hgnc_symbol])
# eigen centrality
eigen_cent = igraph::eigen_centrality(net.hubs.c, directed=F, weights=NA)$vector
sort(eigen_cent, decreasing = T)
names_eig = names(sort(closeness_c_hubs, decreasing = T))
head(names_eig[names_eig %in% genes_hubs.c.only$hgnc_symbol])

# 4: Differential Co-expression network --------------------------------------

#Computation. Using only DEGs, compute the differential co-expression 
#network (Cancer vs. Normal) following the procedure described in the slides 
#(DEPM_5). Binary adjacency matrix with aij=0 if |Z|<3

tresholdZ = 4
data <- getFishernet(cor.mat.c, cor.mat.n, tresholdZ, n.patient = dim(filtr.expr.c)[2])
adj.mat.Z <- data$adj.mat
net.Z <- data$net

ggnet2(net.Z, size = 1)
edge_density(simplify(net.Z), loops = FALSE)

deg.Z <- igraph::degree(net.Z)
Z.frequencies <- as.data.frame(table(deg.Z))[, 2]

getDegreedistribution(net.Z, title = "DiffNetwork", output_file = "figure/Z_degree_distr.pdf")
getHistoDistro(Z.frequencies, title = "z", bw = 10, output_file = "figure/Z_histo.pdf")


# Fit Power Law
powerlaw_test(Z.frequencies)

# hubs 
hubs.Z = getHubs(net.Z) 
hubs.Z.genes.name <- getName(names(hubs.Z))

length(hubs.Z.genes.name$ensembl_gene_id)
# find the common hubs 
common.hubs.zn = getName(intersect(names(hubs.Z), names(hubs.n)))
common.hubs.zc = getName(intersect(names(hubs.Z), names(hubs.c)))

intersect(common.hubs.zn$hgnc_symbol, common.hubs.zc$hgnc_symbol)


hubs.Z.adj = getHubSubnetworkZ(hubs.Z, adj.mat.Z, net.Z)
net.hubs.z <- graph_from_adjacency_matrix(hubs.Z.adj, mode = "undirected", 
                                          weighted = TRUE)


suppressWarnings(getSubnetworkGraph(net.hubs.z, hubs.Z.genes.name))


# 5: PSN and community detection ------------------------------

#Compute the Patient Similarity Network using cancer gene expression 
#profile

# compute the correlation matrix Patient x Patient
cor.mat.PSN <- cor(filtr.expr.c, method = "pearson")
# remove the self-loop
diag(cor.mat.PSN) <- 0
#Perform the community detection (e.g. apply Louvain algorithm to the PSN)
net.PSN = graph_from_adjacency_matrix(cor.mat.PSN, 
                                      mode = "undirected",
                                      weighted = TRUE)

ggnet2(net.PSN, size = 8, mode = "circle", edge.size = E(net.PSN)$weight, 
       edge.alpha = 0.5)
edge_density(net.PSN)

clusters = cluster_louvain(net.PSN, weights = E(net.PSN)$weight, resolution = 1)
clusters$memberships

coords = layout_with_fr(net.PSN)
plot(clusters, net.PSN, layout=coords)


# PSN with clinical data ----------------------------------------------

# extract only patient in common C vs N 
common.patients = colnames(filtr.expr.c)
clinical_data = clinical.query[clinical.query$submitter_id %in% common.patients, ]
# remove columns with NA
clinical_data = clinical_data[ , colSums(is.na(clinical_data))==0]
length(colnames(clinical_data)) # 41
colnames(clinical_data)
# remove non informative columns where all the patients have the same information
clinical_data = clinical_data[vapply(clinical_data, function(x) length(unique(x)) > 1, logical(1L))]
length(colnames(clinical_data)) # 23
colnames(clinical_data)
# remove non informative columns like datetime and bar code
keeping_col = c("submitter_id", 
                "ajcc_pathologic_stage", 
                "prior_malignancy", 
                "ajcc_staging_system_edition", 
                "ajcc_pathologic_t", 
                "ajcc_pathologic_n", 
                "race", 
                "gender",
                "treatments_pharmaceutical_treatment_or_therapy")
clinical_data = clinical_data[, keeping_col]
patients.name <- clinical_data$submitter_id
clinical_str_data = clinical_data
# to factor
clinical_data = sapply(clinical_data[, 2:length(colnames(clinical_data))], 
                        function(x) as.numeric(factor(x)))
rownames(clinical_data) <- patients.name


# build the PSN w clinical data
# compute the correlation matrix Patient x Patient
cor.mat.PSN.clinical <- cor(t(clinical_data), method = "pearson")
# remove the self-loop
diag(cor.mat.PSN.clinical) <- 0
#Perform the community detection (e.g. apply Louvain algorithm to the PSN)
net.PSN = igraph::graph_from_adjacency_matrix(cor.mat.PSN.clinical, 
                                      mode = "undirected",
                                      weighted = TRUE)

clusters = cluster_louvain(net.PSN,
                           weights = abs(E(net.PSN)$weight),
                           resolution = 1)

net.PSN = set_vertex_attr(net.PSN, "x", index = V(net.PSN), coords[, 1])
net.PSN = set_vertex_attr(net.PSN, "y", index = V(net.PSN), coords[, 2])

coords = layout_with_fr(net.PSN, weights = abs(E(net.PSN)$weight))
pdf(file="PSN_clinical.pdf")
plot(clusters, net.PSN, 
     layout=coords,
     vertex.label = NA,
     edge.color = "darkgray")
dev.off()


see <- c("TCGA-KN-8432", "TCGA-KN-8431", "TCGA-KN-8425", "TCGA-KN-8428", 
         "TCGA-KL-8336","TCGA-KL-8329" )
rownames(clinical_str_data) <- patients.name

c1 = clinical_str_data[clusters[[1]],]
c2 = clinical_str_data[clusters[[2]],]
c3 = clinical_str_data[clusters[[3]],]

c1
c2
c3


# BONUS 1 ------------

# Hard and soft thresholding

hard.adj.mat.c = applyHardThreshold(cor.mat.c, adj.mat.c)
hard.adj.mat.n = applyHardThreshold(cor.mat.n, adj.mat.n)

soft.adj.mat.c = applySoftThreshold(cor.mat.c, adj.mat.c)
soft.adj.mat.n = applySoftThreshold(cor.mat.n, adj.mat.n)

# ADJUST PARAMETERS TO HAVE SCALE-FREE
net.c.hard <- igraph::graph_from_adjacency_matrix(hard.adj.mat.c, mode = "undirected", 
                                                  weighted = TRUE)
net.n.hard <- igraph::graph_from_adjacency_matrix(hard.adj.mat.n, mode = "undirected", 
                                                  weighted = TRUE)

net.c.soft <- igraph::graph_from_adjacency_matrix(soft.adj.mat.c, mode = "undirected", 
                                                  weighted = TRUE)
net.n.soft <- igraph::graph_from_adjacency_matrix(soft.adj.mat.n, mode = "undirected", 
                                                  weighted = TRUE)

#Compare network topology


info1 = getInfoNetworkTopology(net.c.hard)
info2 = getInfoNetworkTopology(net.n.hard)
info3 = getInfoNetworkTopology(net.c.soft)
info4 = getInfoNetworkTopology(net.n.soft)

topology.infos = rbind(info1, info2, info3, info4)
dim(topology.infos)

colnames(topology.infos) = c("Global efficiency", "Avg local efficiency (IN)",
                             "Avg local efficiency (OUT)", "Diameter", "Max degree",
                             "Number of connected components",
                             "Avg Closeness Centrality", "Avg Betweeness Centrality",
                             "Hub score", "Authority score")
rownames(topology.infos) = c("Hard Tresholding Cancer Network", "Hard Tresholding Normal Network", 
                             "Soft Tresholding Cancer Network", "Soft Tresholding Normal Network")

kbl(topology.infos) %>% 
  kable_styling(bootstrap_options = c("striped", "hover", "condensed", "responsive"))

# BONUS 2 ----------

# We used betweeness centrality
net.c.abs = net.c

# if weight is <0 get abs(w)
E(net.c.abs)$weight[which(E(net.c.abs)$weight < 0)] = abs(E(net.c.abs)$weight[which(E(net.c.abs)$weight < 0)])

# hubs per betw
hubs.c.betw = getHubs(net.c.abs, measure = "betweenness")
hubs.n.betw = getHubs(net.c.abs, measure = "betweenness")

hubs.c.betw = getName(names(hubs.c.betw))
hubs.n.betw = getName(names(hubs.n.betw))

#sort(igraph::betweenness(net.c.abs, normalized = T), decreasing = T)[1]

# Names and compare sets
cat(sort(paste(hubs.c.betw[,2])))
cat(sort(paste(hubs.c.genes.name[,2])))

common.c.hubs = intersect(hubs.c.betw$hgnc_symbol, hubs.c.genes.name$hgnc_symbol)
common.c.hubs
# "AFMID"
common.n.hubs = intersect(hubs.n.betw$hgnc_symbol, hubs.n.genes.name$hgnc_symbol)
common.n.hubs
# " AK4 "

# BONUS 3 ---------------

# Perform the study using different similarity measure
#1. Building co-expression networks
#with ::: biweight midcorrelation ::
# BiocManager::install("WGCNA", force = T)
library(WGCNA)

bi.corr.c.obj = bicorAndPvalue(t(DEG.c))
bi.corr.n.obj = bicorAndPvalue(t(DEG.n))

cor.mat <- bi.corr.c.obj$bicor
diag(cor.mat) <- 0
adj.mat1 <- ifelse(cor.mat >= 0.5, 1,ifelse(cor.mat <= -0.5,-1, 0 ))
cor.padj <- corr.p(cor.mat, nrow(cor.mat), adjust="fdr",ci=FALSE)$p
cor.padj[lower.tri(cor.padj)] <- t(cor.padj)[lower.tri(cor.padj)]
adj.mat2 <- ifelse(abs(cor.padj) > 0.05, 0, 1) 
bi.corr.c <- adj.mat1 * adj.mat2

sum(bi.corr.c < 0) / 2
sum(bi.corr.c > 0) / 2
sum(bi.corr.c == 0) / 2

cor.mat <- bi.corr.n.obj$bicor
diag(cor.mat) <- 0
adj.mat1 <- ifelse(cor.mat >= 0.6, 1,ifelse(cor.mat <= -0.6,-1, 0 ))
cor.padj <- corr.p(cor.mat, nrow(cor.mat), adjust="fdr",ci=FALSE)$p
cor.padj[lower.tri(cor.padj)] <- t(cor.padj)[lower.tri(cor.padj)]
adj.mat2 <- ifelse(abs(cor.padj) > 0.05, 0, 1) 
bi.corr.n <- adj.mat1 * adj.mat2

# Build the networks
net.c.bicor <- graph_from_adjacency_matrix(bi.corr.c, mode = "undirected", 
                                     weighted = TRUE)
net.n.bicor <- graph_from_adjacency_matrix(bi.corr.n, mode = "undirected", 
                                     weighted = TRUE)

ggnet2(net.c.bicor, node.size = 0.2)

# are scale free?
deg.c <- igraph::degree(net.c.bicor)
deg.c <- deg.c[deg.c > 0]
C.frequencies <- as.data.frame(table(deg.c))[, 2]

getDegreedistribution(net.c.bicor, title = "CancerBicORR", output_file = "figure/BiCorrC_distr.pdf")
# Testing Fit Power Law
powerlaw_test(C.frequencies)


# analysis on graph: NORMAL
deg.n <- igraph::degree(net.n.bicor)
N.frequencies <- as.data.frame(table(deg.n))[, 2]
getDegreedistribution(net.n.bicor, title = "NormalBicORR", output_file = "figure/BiCorrN_distr.pdf")

# Fit Power Law
powerlaw_test(N.frequencies)

# Find the hubs (5% of the nodes with highest degree values)
# sorted by degree
hubs.c.bicor = getHubs(net.c.bicor) 
hubs.n.bicor = getHubs(net.n.bicor) 

# graph nodes w/ hubs 
hubs.c.genes.name.bicor <- getName(names(hubs.c.bicor))
hubs.n.genes.name.bicor <- getName(names(hubs.n.bicor))


length(hubs.c.genes.name.bicor$hgnc_symbol)
length(hubs.n.genes.name.bicor$hgnc_symbol)

# compare hubs sets related to the two condition (cancer, normal)
common.hubs.bicor = getName(intersect(names(hubs.c.bicor), names(hubs.n.bicor)))
common.hubs.bicor$hgnc_symbol

hub.c.adj.bicor <- getHubSubnetwork(hubs.c.bicor, hubs.n.bicor, bi.corr.c, net.c.bicor)
net.hubs.c.bicor <- graph_from_adjacency_matrix(hub.c.adj.bicor, mode = "undirected", 
                                          weighted = TRUE)

suppressWarnings(getSubnetworkGraph(net.hubs.c.bicor, hubs.c.genes.name.bicor,
                                    output_file = "bicorr_hubs.pdf"))

## Info on hubs 
common.hubs.bicor = intersect(names(hubs.c.bicor), names(hubs.n.bicor))
condition = which(names(hubs.c.bicor) %in% setdiff(names(hubs.c.bicor), common.hubs.bicor))
hubs.c.only.bicor = hubs.c.bicor[condition]
genes_hubs.c.only.bicor = getName(names(hubs.c.only.bicor))

# using other methods
method = c("pearson", "kendall", "spearman")

# run and collect
data1 <- getCO.EXPR(method[1],DEG.c, DEG.n)
data2 <- getCO.EXPR(method[1],DEG.c, DEG.n)
data3 <- getCO.EXPR(method[1],DEG.c, DEG.n)

other.corr.analysis = rbind(data1$common.hubs, data2$common.hubs, data3$common.hubs)

rownames(other.corr.analysis) = c("Pearson", "Kendall", "Spearman")
other.corr.analysis

# BONUS 4 ---------------
cor.mat.PSN.n <- cor(filtr.expr.n, method = "pearson")
diag(cor.mat.PSN.n) <- 0
net.PSN.n = graph_from_adjacency_matrix(cor.mat.PSN.n, 
                                        mode = "undirected",
                                        weighted = TRUE)

ggnet2(net.PSN.n, size = 8, mode = "circle", 
       edge.size = E(net.PSN.n)$weight, 
       edge.alpha = 0.5)

edge_density(net.PSN.n)

clusters.n = cluster_louvain(net.PSN.n, weights = E(net.PSN.n)$weight, resolution = 1)
clusters.n$memberships

coords.n = layout_with_fr(net.PSN.n)
pdf(file="PSN_genexpr_n.pdf")
plot(clusters.n, net.PSN.n, layout=coords.n, vertex.label = NA, edge.color="darkgrey")
dev.off()





