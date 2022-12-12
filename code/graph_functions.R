getFC_pvals <- function(cancer.data, normal.data, method = "fdr"){
  #Fold change:
  mean_C = rowMeans(cancer.data)
  mean_N = rowMeans(normal.data)
  log2FC = log2(mean_C)-log2(mean_N)
  
  #paired t.test
  data <- cbind(filtr.expr.c, filtr.expr.n)
  pvals <- apply(data, 1, function(data){
    t.test(data[1:24], data[25:48], paired = TRUE)$p.value
  })
  
  # add FDR corrections for multiple test
  pvals_adj = p.adjust(pvals, method = method)
  
  return(data.frame(log2FC = log2FC,
                    pvals_adj = pvals_adj))
}

getFCHistogram <- function(data, treshold, output_file = "figure/FC_histogram.pdf"){
  
  df <- data.frame(log2FC = data)
  df$legend <- ifelse(abs(data) <= treshold, "ignore", "accept")
  p  <- ggplot(df, aes(x=log2FC, fill = legend)) + 
    geom_histogram(color="black", bins = 50) +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0), breaks = scales::pretty_breaks(n = 10)) +
    scale_fill_manual(values = c("ignore" = "darkgray", "accept" = "red")) +
    labs(x = "Log2 Fold-change (FC)", 
         y = "Frequency") +
    guides(colour = guide_legend(override.aes = list(size=3))) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
          legend.title = element_blank(), 
          legend.position="top",
          legend.text=element_text(size=9, face = "bold"),
          axis.text = element_text(size = 15, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          aspect.ratio=1,
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
  
  print(p)
  ggsave(output_file, p, width = 5, height = 5, dpi = 300)
  return(p)
}


# Volcano plot
getVolcanoPlot <- function(log2FC, pvals_adj, tresholdFC, tresholdPval,
                           output_file = "figure/volcanoPlot.pdf"){
  
  data <- data.frame(FC = log2FC, pval = -log10(pvals_adj))
  up_regulated_genes   <- (pvals_adj <= tresholdPval) & (log2FC > tresholdFC)
  down_regulated_genes <- (pvals_adj <= tresholdPval) & (log2FC < -tresholdFC)
  
  data$legend <- ifelse(up_regulated_genes, "up-regulated genes", 
                        ifelse(down_regulated_genes, "down-regulated genes", "ignore")) 
  
  p <- ggplot(data, aes(x = FC, y = pval, color = legend)) +
    geom_point() +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0)) +
    scale_color_manual(values = c("ignore" = "darkgrey", 
                                  "up-regulated genes" = "red", 
                                  "down-regulated genes" = "green")) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20), 
          legend.title = element_text(colour = "black", size = 10, face = "bold"), 
          legend.position="top", 
          text = element_text(size = 15, face = "bold")) +
    guides(colour = guide_legend(override.aes = list(size=4))) +
    labs(x = "Log2 Fold-change (FC)", 
         y = "-log10(adjusted p-value)") +
    geom_hline(yintercept = -log10(tresholdPval), linetype = "dashed", 
               color = "black") +
    labs(colour = NULL) +
    geom_vline(xintercept = tresholdFC, linetype = "dashed", color = "black") +
    geom_vline(xintercept = -tresholdFC, linetype = "dashed", color = "black")
  
  print(p)
  ggsave(output_file, p, width = 5, height = 5, dpi = 300)
  return(p)
  
}


getCorrHistogram <- function(corr.matrix, treshold, 
                             output_file = "figure/Corr_histogram.pdf"){
  
  array_of_corr = c(corr.matrix[lower.tri(corr.matrix)])
  df <- data.frame(corr = array_of_corr)
  df$legend <- ifelse(abs(array_of_corr) <= treshold, "ignore", "accept")
  p  <- ggplot(df, aes(x=corr, fill = legend)) + 
    geom_histogram(color="black", bins = 50) +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0), breaks = scales::pretty_breaks(n = 10)) +
    scale_fill_manual(values = c("ignore" = "darkgray", "accept" = "red")) +
    labs(title = "Histogram of correlation", x = "correlation", 
         y = "Frequency") +
    #guides(colour = guide_legend(override.aes = list(size=3))) +
    guides(colour = "none") + 
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
          legend.title = element_blank(), 
          legend.position = c(0.87, 0.9), 
          legend.text=element_text(size=9, face = "bold"),
          axis.text = element_text(size = 15, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          aspect.ratio=1,
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
  
  print(p)
  ggsave(output_file, p, width = 5, height = 5, dpi = 300)
}

getDegreedistribution <- function(net, title, output_file = "figure/Cancer_deg_distr.pdf"){
  
  G.degrees <- igraph::degree(net)
  G.degree.histogram <- as.data.frame(table(G.degrees))
  G.degree.histogram[,1] <- as.numeric(G.degree.histogram[,1])
  
  log10(G.degree.histogram[,2])
  G.degree.histogram[,2]
  
  p <- ggplot(G.degree.histogram, aes(x = G.degrees, y = Freq)) +
    geom_point(size = 3) +
    scale_x_continuous("Degree",
                       trans = 'log2',
                       breaks = trans_breaks('log2', function(x) 2^x),
                       labels = trans_format('log2', math_format(2^.x))) +
    scale_y_continuous("Frequency",
                       trans = 'log2',
                       breaks = trans_breaks('log2', function(x) 2^x),
                       labels = trans_format('log2', math_format(2^.x))) +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 15), 
          legend.title = element_blank(), 
          legend.position="bottom",
          legend.text=element_text(size=9, face = "bold"),
          axis.text = element_text(size = 15, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          aspect.ratio=1,
          axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
  print(p)
  ggsave(output_file, p, width = 5, height = 5, dpi = 300)
}

getSubnetworkGraph <- function(net.hubs.c, hubs.c.genes.name,
                               output_file = "figure/Hubs_c.pdf"){
  
  V(net.hubs.c)$label <- names(V(net.hubs.c))
  V(net.hubs.c)$label[!(names(V(net.hubs.c)) %in% hubs.c.genes.name[,2])] <- NA
  
  V(net.hubs.c)$size  <- 1
  V(net.hubs.c)$size[!(names(V(net.hubs.c)) %in% hubs.c.genes.name[,2])] <- 0.5
  
  V(net.hubs.c)$color <- ifelse(names(V(net.hubs.c)) %in% hubs.c.genes.name[,2],
                                "hubs", "neigh")
  
  p <- ggnet2(net.hubs.c, palette = "Set1", 
              color = V(net.hubs.c)$color,
              label = V(net.hubs.c)$label, 
              size =  V(net.hubs.c)$size,
              label.size = 2, 
              edge.alpha = 0.2,
              legend.size = 12, legend.position = "bottom",
              mode = "fruchtermanreingold", layout.par = list(cell.jitter = 0.75)) +
    theme(legend.title = element_text(colour = "black", size = 10, face = "bold")) +
    guides(size = "none") +
    theme(aspect.ratio=1)
  
  suppressWarnings(print(p))
  ggsave(output_file, suppressWarnings(print(p)), width = 5, height = 5, dpi = 300)
  
}

getHistoDistro <- function(frequencies, bw = 10, title, output_file){
  
  df <- data.frame(freq = frequencies)
  p  <- ggplot(df, aes(x=freq)) + 
    geom_histogram(aes(y=after_stat(count)/sum(after_stat(count))), fill = "gray", 
                   color = "black", binwidth = bw) +
    scale_x_continuous(expand = c(0,0)) + 
    scale_y_continuous(expand = c(0,0), breaks = scales::pretty_breaks(n = 5)) +
    labs(x = "Degree", 
         y = "Frequency") +
    theme_classic() +
    theme(axis.text = element_text(size = 15, face = "bold"),
          axis.title = element_text(size = 15, face = "bold"),
          aspect.ratio=1,
          legend.position = "none",
          axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
  print(p)
  ggsave(output_file, p, width = 5, height = 5, dpi = 300)
}
