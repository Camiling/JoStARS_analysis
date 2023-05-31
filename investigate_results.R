rm(list=ls())
library(tailoredGlasso) # From Github/Camiling/tailoredGlasso
load('data/PanCancer_data.RData')
load('data/res_stabJGL_PanCancer.RData')
source('useful_functions.R')

save.results=T

# FIND PRECISION MATRICES AND CHECK EDGE AGREEMENT 

p = ncol(res.joint$opt.fit[[1]]) # 131

res.joint$opt.lambda1
# 0.3226316
res.joint$opt.lambda2
# 0.007894737

# Get precision matrices from joint estimates
theta.est.brca  = cov2cor(res.joint$opt.fit[[1]])
theta.est.brca[which(abs(theta.est.brca) < 1e-5, arr.ind = T)] = 0
theta.est.OVCA  = cov2cor(res.joint$opt.fit[[2]])
theta.est.OVCA[which(abs(theta.est.OVCA) < 1e-5, arr.ind = T)] = 0
theta.est.ucec  = cov2cor(res.joint$opt.fit[[3]])
theta.est.ucec[which(abs(theta.est.ucec) < 1e-5, arr.ind = T)] = 0

# Sparsity
sparsity(theta.est.brca!=0)
# 0.04885496
sparsity(theta.est.OVCA!=0)
# 0.03581914
sparsity(theta.est.ucec!=0)
# 0.0386377

# How many edges do they agree on?
confusion.matrix(theta.est.brca!=0, theta.est.OVCA!=0)
#      [,1] [,2]
#[1,]  199  154
#[2,]  400 7762

confusion.matrix(theta.est.brca!=0, theta.est.ucec!=0)
#    [,1] [,2]
#[1,]  173  302
#[2,]  426 7614

confusion.matrix(theta.est.ucec!=0, theta.est.OVCA!=0)
#     [,1] [,2]
#[1,]  147  206
#[2,]  328 7834

(sum(theta.est.brca!=0 & theta.est.ucec!=0 & theta.est.OVCA!=0)-p)/2
# 111

# Also MCC
MCC(theta.est.brca!=0, theta.est.OVCA!=0)
# 0.4012364

MCC(theta.est.brca!=0, theta.est.ucec!=0)
# 0.2793082

MCC(theta.est.ucec!=0, theta.est.OVCA!=0)
# 0.3267991




# CHECK TOP HUBS -----------------------------

g.brca = igraph::graph.adjacency(theta.est.brca!=0, mode='undirected', diag=F)
g.ucec = igraph::graph.adjacency(theta.est.ucec!=0, mode='undirected', diag=F)
g.OVCA = igraph::graph.adjacency(theta.est.OVCA!=0, mode='undirected', diag=F)

df.degree = data.frame(degree=c(igraph::degree(g.brca), igraph::degree(g.ucec),igraph::degree(g.OVCA)), 
                       Group=factor(c(rep('BRCA', p), rep('UCEC', p),rep('OVCA', p))))

df.degree.brca = data.frame(protein=colnames(brca_dat),gene=mapping.frame$gene[match(colnames(brca_dat),mapping.frame$protein)], 
                            degree=df.degree[which(df.degree$Group=='BRCA'),1])
df.degree.OVCA = data.frame(protein=colnames(ovac_dat),gene=mapping.frame$gene[match(colnames(ovac_dat),mapping.frame$protein)], 
                            degree=df.degree[which(df.degree$Group=='OVCA'),1])
df.degree.ucec = data.frame(protein=colnames(ucec_dat),
                            gene=mapping.frame$gene[match(colnames(ucec_dat),mapping.frame$protein)], 
                            degree=df.degree[which(df.degree$Group=='UCEC'),1])

df.degree.brca.ordered = df.degree.brca[rev(order(df.degree.brca$degree)),]
df.degree.OVCA.ordered = df.degree.OVCA[rev(order(df.degree.OVCA$degree)),]
df.degree.ucec.ordered = df.degree.ucec[rev(order(df.degree.ucec$degree)),]



# Print top table ----------------------
top.frac=0.9

df.degree.brca.ordered.top = df.degree.brca.ordered[which(df.degree.brca.ordered$degree>quantile(df.degree.brca.ordered$degree,top.frac)),]
df.degree.OVCA.ordered.top = df.degree.OVCA.ordered[which(df.degree.OVCA.ordered$degree>quantile(df.degree.OVCA.ordered$degree,top.frac)),]
df.degree.ucec.ordered.top = df.degree.ucec.ordered[which(df.degree.ucec.ordered$degree>quantile(df.degree.ucec.ordered$degree,top.frac)),]

which.unique.brca = which(! df.degree.brca.ordered.top$protein %in% c(df.degree.OVCA.ordered.top$protein,df.degree.ucec.ordered.top$protein))
which.unique.ucec = which(! df.degree.ucec.ordered.top$protein %in% c(df.degree.OVCA.ordered.top$protein,df.degree.brca.ordered.top$protein))
which.unique.OVCA = which(! df.degree.OVCA.ordered.top$protein %in% c(df.degree.ucec.ordered.top$protein,df.degree.brca.ordered.top$protein))
which.common.brca = which(df.degree.brca.ordered.top$protein %in% df.degree.OVCA.ordered.top$protein & df.degree.brca.ordered.top$protein %in% df.degree.ucec.ordered.top$protein)
which.common.ucec = which(df.degree.ucec.ordered.top$protein %in% df.degree.brca.ordered.top$protein & df.degree.ucec.ordered.top$protein %in% df.degree.OVCA.ordered.top$protein)
which.common.OVCA= which(df.degree.OVCA.ordered.top$protein %in% df.degree.brca.ordered.top$protein & df.degree.OVCA.ordered.top$protein %in% df.degree.ucec.ordered.top$protein)


dim.neighs.top = max(c(length(df.degree.OVCA.ordered.top$protein),length(df.degree.ucec.ordered.top$protein),length(df.degree.OVCA.ordered.top$protein)))

# Also print gene name
for(i in 1:dim.neighs.top){
  if(is.na(df.degree.brca.ordered.top[i,1])) { cat(' & & && ') }
  else {
    if(i %in% which.common.brca){
      cat(paste0('\\textbf{',df.degree.brca.ordered.top$protein[i], '} & \\emph{ ',df.degree.brca.ordered.top$gene[i], '} & ', df.degree.brca.ordered.top$degree[i], ' && '))
    }
    else if(i %in% which.unique.brca){
      cat(paste0('\\color{red}{',df.degree.brca.ordered.top$protein[i], '} & \\emph{ ',df.degree.brca.ordered.top$gene[i], '} & ', df.degree.brca.ordered.top$degree[i], ' && '))
    }
    else { cat(paste0(df.degree.brca.ordered.top$protein[i], ' & \\emph{',df.degree.brca.ordered.top$gene[i], '} & ', df.degree.brca.ordered.top$degree[i], ' && ')) }
  }
  if(is.na(df.degree.ucec.ordered.top[i,1])) { cat(' & & && ') }
  else {
    if(i %in% which.common.ucec){
      cat(paste0('\\textbf{',df.degree.ucec.ordered.top$protein[i], '} & \\emph{ ',df.degree.ucec.ordered.top$gene[i], '} & ', df.degree.ucec.ordered.top$degree[i], ' && '))
    }
    else if(i %in% which.unique.ucec){
      cat(paste0('\\color{red}{',df.degree.ucec.ordered.top$protein[i], '} & \\emph{ ',df.degree.ucec.ordered.top$gene[i], '} & ', df.degree.ucec.ordered.top$degree[i], ' && '))
    }
    else { cat(paste0(df.degree.ucec.ordered.top$protein[i], ' & \\emph{',df.degree.ucec.ordered.top$gene[i], '} & ', df.degree.ucec.ordered.top$degree[i], ' && ')) }
  }
  if(is.na(df.degree.OVCA.ordered.top[i,1])) { cat(' & & \\\\ \n') }
  else {
    if(i %in% which.common.OVCA){
      cat(paste0('\\textbf{',df.degree.OVCA.ordered.top$protein[i], '} & \\emph{ ',df.degree.OVCA.ordered.top$gene[i], '} & ', df.degree.OVCA.ordered.top$degree[i], ' \\\\ \n  '))
    }
    else if(i %in% which.unique.OVCA){
      cat(paste0('\\color{red}{',df.degree.OVCA.ordered.top$protein[i], '} & \\emph{ ',df.degree.OVCA.ordered.top$gene[i], '} & ', df.degree.OVCA.ordered.top$degree[i], ' \\\\ \n  '))
    }
    else { cat(paste0(df.degree.OVCA.ordered.top$protein[i], ' & \\emph{',df.degree.OVCA.ordered.top$gene[i], '} & ', df.degree.OVCA.ordered.top$degree[i], ' \\\\ \n ')) }
  }
}

# MAKE EDGE LIST IN FILE --------------------------------

edges.brca = get_and_print_edges(theta.est.brca!=0,colnames(brca_dat),theta.est.brca)
colnames(edges.brca) =  c("Protein1", "Protein2", "PartialCor")
edges.ucec = get_and_print_edges(theta.est.ucec!=0,colnames(ucec_dat),theta.est.ucec)
colnames(edges.ucec) =  c("Protein1", "Protein2", "PartialCor")
edges.OVCA = get_and_print_edges(theta.est.OVCA!=0,colnames(ovac_dat),theta.est.OVCA)
colnames(edges.OVCA) =  c("Protein1", "Protein2", "PartialCor")
write.csv(edges.brca, file = "edge_lists/edgesBRCA.csv", row.names = F,quote=F)
write.csv(edges.ucec, file = "edge_lists/edgesUCEC.csv", row.names = F,quote=F)
write.csv(edges.OVCA, file = "edge_lists/edgesOVCA.csv", row.names = F,quote=F)


# PLOT RESULTING NETWORKS ------------------------------

# Plot with common edges marked
unique.list = list()
unique.list[[1]] = (theta.est.brca!=0) + 0  # If edge in brca
unique.list[[1]][which(theta.est.brca!=0 & theta.est.OVCA!=0 & theta.est.ucec!=0)] = 3 # Present in all
unique.list[[2]] = (theta.est.OVCA!=0) + 0  # If edge in OVCA
unique.list[[2]][which(theta.est.brca!=0 & theta.est.OVCA!=0 & theta.est.ucec!=0)] = 3 # Present in all
unique.list[[3]] = (theta.est.ucec!=0) + 0  # If edge in ucec
unique.list[[3]][which(theta.est.brca!=0 & theta.est.OVCA!=0 & theta.est.ucec!=0)] = 3 # Present in all

# Get layout for Responders
set.seed(123)
net.layout = network::network(unique.list[[1]],directed=F, ignore.eval=F,names.eval='weights')
x.layout = sna::gplot.layout.fruchtermanreingold(net.layout, NULL)

nets.layout=list()
net.layout.brca = network::network(unique.list[[1]],directed=F, ignore.eval=F,names.eval='weights')
network::set.edge.attribute(net.layout.brca, "color", c("black", "grey75","red", "aquamarine4")[(net.layout.brca %e% "weights")+1])
net.layout.brca %v% "x" = x.layout[, 1]
net.layout.brca %v% "y" = x.layout[, 2]
nets.layout[[1]] = GGally::ggnet2(net.layout.brca,node.size = 1.5, edge.size = 0.25,alpha=0.9,mode = c('x','y'),color = 'royalblue1', edge.color = 'color')+
  ggplot2::ggtitle('BRCA') + ggplot2::theme(plot.title=ggplot2::element_text(size=15,hjust=0.5))
net.layout.OVCA = network::network(unique.list[[2]],directed=F, ignore.eval=F,names.eval='weights')
network::set.edge.attribute(net.layout.OVCA, "color", c("black", "grey75","red", "aquamarine4")[(net.layout.OVCA %e% "weights")+1])
net.layout.OVCA %v% "x" = x.layout[, 1]
net.layout.OVCA %v% "y" = x.layout[, 2]
nets.layout[[2]] = GGally::ggnet2(net.layout.OVCA,node.size = 1.5, edge.size = 0.25,alpha=0.9,mode = c('x','y'),color = 'royalblue1', edge.color = 'color')+
  ggplot2::ggtitle('OVCA') + ggplot2::theme(plot.title=ggplot2::element_text(size=15,hjust=0.5))
net.layout.ucec = network::network(unique.list[[3]],directed=F, ignore.eval=F,names.eval='weights')
network::set.edge.attribute(net.layout.ucec, "color", c("black", "grey75","red", "aquamarine4")[(net.layout.ucec %e% "weights")+1])
net.layout.ucec %v% "x" = x.layout[, 1]
net.layout.ucec %v% "y" = x.layout[, 2]
nets.layout[[3]] = GGally::ggnet2(net.layout.ucec,node.size = 1.5, edge.size = 0.25,alpha=0.9,mode = c('x','y'),color = 'royalblue1', edge.color = 'color')+
  ggplot2::ggtitle('UCEC') + ggplot2::theme(plot.title=ggplot2::element_text(size=15,hjust=0.5))

if(save.results){
  pdf("plots/nets_edgesmarked_samelayout.pdf", 10,4)
  print(gridExtra::grid.arrange(grobs=nets.layout,ncol=3))
  dev.off()
}


# CHECK EVIDENCE IN STRING  -----------------------------

# Get edges with gene names
edges.brca.gene.all = get_and_print_edges(theta.est.brca!=0,mapping.frame$gene[match(colnames(brca_dat), mapping.frame$protein)])
edges.OVCA.gene.all = get_and_print_edges(theta.est.OVCA!=0,mapping.frame$gene[match(colnames(ovac_dat), mapping.frame$protein)])
edges.ucec.gene.all = get_and_print_edges(theta.est.ucec!=0,mapping.frame$gene[match(colnames(ucec_dat), mapping.frame$protein)])
edges.brca.gene = edges.brca.gene.all[-which(duplicated(edges.brca.gene.all)),]
edges.OVCA.gene = edges.OVCA.gene.all[-which(duplicated(edges.OVCA.gene.all)),]
edges.ucec.gene = edges.ucec.gene.all[-which(duplicated(edges.ucec.gene.all)),]
colnames(edges.brca.gene) =  c("Gene1", "Gene2")
colnames(edges.ucec.gene) =  c("Gene1", "Gene2")
colnames(edges.OVCA.gene) =  c("Gene1", "Gene2")


# Print unique names for STRING
cat(unique(mapping.frame$gene), sep='\n')
# We get the network https://string-db.org/cgi/network?taskId=bb7kLoAahcIx&sessionId=bZmv5BTkUpEx
string.output = 'https://string-db.org/cgi/generatetaskspecificdownloadfile?taskId=bb7kLoAahcIx&downloadDataFormat=tsv_short&cpnonce=b91LPDBj3ACO&downloadFileName=string_interactions_short.tsv'
string.rppa = utils::read.csv(string.output, sep = "\t")
# Sort so that gene pairs are given in alphabetical order.
string.rppa[, 1:2] = t(apply(string.rppa[, 1:2], 1, sort)) # First gene in alphabet is always in col 1
colnames(string.rppa)[1:2] =  c("Gene1", "Gene2")

# Combine edges in STRING and stabJGL into one table, so that we can check how many edges occur twice in the table 
# BRCA first
df.all.edges.brca =  rbind(edges.brca.gene, string.rppa[, 1:2])
ids.in.string = which(duplicated(df.all.edges.brca)) - nrow(edges.brca.gene) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with these gene pairs XX HERE
edges.brca.gene.evidence = string.rppa[ids.in.string,1:2]
n.extra.pairs.evidence.brca = 0
for(i in 1:nrow(edges.brca.gene.evidence)){
  n.rep = 0
  for(j in 1:nrow(edges.brca.gene.all)){
    if(all(edges.brca.gene.all[j,] == edges.brca.gene.evidence[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.brca = n.extra.pairs.evidence.brca + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string)+n.extra.pairs.evidence.brca # 68
(length(ids.in.string)+n.extra.pairs.evidence.brca )/nrow(edges.brca.gene.all) # 0.1243144
# Also the amount of unique edges
length(ids.in.string)/nrow(edges.brca.gene) # 0.1233141
# Then UCEC
df.all.edges.ucec =  rbind(edges.ucec.gene, string.rppa[, 1:2])
ids.in.string = which(duplicated(df.all.edges.ucec)) - nrow(edges.ucec.gene) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with these gene pairs 
edges.ucec.gene.evidence = string.rppa[ids.in.string,1:2]
n.extra.pairs.evidence.ucec = 0
for(i in 1:nrow(edges.ucec.gene.evidence)){
  n.rep = 0
  for(j in 1:nrow(edges.ucec.gene.all)){
    if(all(edges.ucec.gene.all[j,] == edges.ucec.gene.evidence[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.ucec = n.extra.pairs.evidence.ucec + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string)+n.extra.pairs.evidence.ucec # 44
(length(ids.in.string)+n.extra.pairs.evidence.ucec)/nrow(edges.ucec.gene.all) # 0.108642
# Also the amount of unique edges
length(ids.in.string)/nrow(edges.ucec.gene) # 0.09974425
# Then OVCA
df.all.edges.OVCA =  rbind(edges.OVCA.gene, string.rppa[, 1:2])
ids.in.string = which(duplicated(df.all.edges.OVCA)) - nrow(edges.OVCA.gene) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with these gene pairs XX HERE
edges.OVCA.gene.evidence = string.rppa[ids.in.string,1:2]
n.extra.pairs.evidence.OVCA = 0
for(i in 1:nrow(edges.OVCA.gene.evidence)){
  n.rep = 0
  for(j in 1:nrow(edges.OVCA.gene.all)){
    if(all(edges.OVCA.gene.all[j,] == edges.OVCA.gene.evidence[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.OVCA = n.extra.pairs.evidence.OVCA + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string) + n.extra.pairs.evidence.OVCA # 42
(length(ids.in.string)+n.extra.pairs.evidence.OVCA)/nrow(edges.OVCA.gene.all) # 0.1257485
# Also the amount of unique edges
length(ids.in.string)/nrow(edges.OVCA.gene) # 0.1238095



# COMPARE TO RESULTS FROM FGL -----------------------------

load('data/res_FGL_PanCancer.RData')
res.fgl$opt.lambda1 # 0.01
res.fgl$opt.lambda2 # 0

# Get precision matrices from joint estimates
theta.est.brca.fgl  = cov2cor(res.fgl$opt.fit[[1]])
theta.est.brca.fgl[which(abs(theta.est.brca.fgl) < 1e-5, arr.ind = T)] = 0
theta.est.OVCA.fgl  = cov2cor(res.fgl$opt.fit[[2]])
theta.est.OVCA.fgl[which(abs(theta.est.OVCA.fgl) < 1e-5, arr.ind = T)] = 0
theta.est.ucec.fgl  = cov2cor(res.fgl$opt.fit[[3]])
theta.est.ucec.fgl[which(abs(theta.est.ucec.fgl) < 1e-5, arr.ind = T)] = 0

# Sparsity
sparsity(theta.est.brca.fgl!=0)
# 0.6892543
sparsity(theta.est.OVCA.fgl!=0)
# 0.7086318
sparsity(theta.est.ucec.fgl!=0)
# 0.6793893

# Also MCC
MCC(theta.est.brca.fgl!=0, theta.est.OVCA.fgl!=0)
# 0.03297134
MCC(theta.est.brca.fgl!=0, theta.est.ucec.fgl!=0)
# 0.02863449
MCC(theta.est.ucec.fgl!=0, theta.est.OVCA.fgl!=0)
# 0.04184616

# Density
g.brca.fgl = igraph::graph.adjacency(theta.est.brca.fgl!=0, mode='undirected', diag=F)
g.ucec.fgl = igraph::graph.adjacency(theta.est.ucec.fgl!=0, mode='undirected', diag=F)
g.OVCA.fgl = igraph::graph.adjacency(theta.est.OVCA.fgl!=0, mode='undirected', diag=F)
df.degree.fgl = data.frame(degree=c(igraph::degree(g.brca.fgl), igraph::degree(g.ucec.fgl),igraph::degree(g.OVCA.fgl)), 
                           Group=factor(c(rep('BRCA', p), rep('UCEC', p),rep('OVCA', p))))


# Printing top list for FGL
top.frac=0.9
df.degree.brca.fgl = data.frame(protein=colnames(brca_dat),gene=mapping.frame$gene[match(colnames(brca_dat),mapping.frame$protein)], 
                            degree=df.degree.fgl[which(df.degree.fgl$Group=='BRCA'),1])
df.degree.OVCA.fgl = data.frame(protein=colnames(ovac_dat),gene=mapping.frame$gene[match(colnames(ovac_dat),mapping.frame$protein)], 
                            degree=df.degree.fgl[which(df.degree.fgl$Group=='OVCA'),1])
df.degree.ucec.fgl = data.frame(protein=colnames(ucec_dat),
                            gene=mapping.frame$gene[match(colnames(ucec_dat),mapping.frame$protein)], 
                            degree=df.degree.fgl[which(df.degree.fgl$Group=='UCEC'),1])

df.degree.brca.ordered.fgl = df.degree.brca.fgl[rev(order(df.degree.brca.fgl$degree)),]
df.degree.OVCA.ordered.fgl = df.degree.OVCA.fgl[rev(order(df.degree.OVCA.fgl$degree)),]
df.degree.ucec.ordered.fgl = df.degree.ucec.fgl[rev(order(df.degree.ucec.fgl$degree)),]
df.degree.brca.ordered.top.fgl = df.degree.brca.ordered.fgl[which(df.degree.brca.ordered.fgl$degree>quantile(df.degree.brca.ordered.fgl$degree,top.frac)),]
df.degree.OVCA.ordered.top.fgl = df.degree.OVCA.ordered.fgl[which(df.degree.OVCA.ordered.fgl$degree>quantile(df.degree.OVCA.ordered.fgl$degree,top.frac)),]
df.degree.ucec.ordered.top.fgl = df.degree.ucec.ordered.fgl[which(df.degree.ucec.ordered.fgl$degree>quantile(df.degree.ucec.ordered.fgl$degree,top.frac)),]

which.unique.brca.fgl = which(! df.degree.brca.ordered.top.fgl$protein %in% c(df.degree.OVCA.ordered.top.fgl$protein,df.degree.ucec.ordered.top.fgl$protein))
which.unique.ucec.fgl = which(! df.degree.ucec.ordered.top.fgl$protein %in% c(df.degree.OVCA.ordered.top.fgl$protein,df.degree.brca.ordered.top.fgl$protein))
which.unique.OVCA.fgl = which(! df.degree.OVCA.ordered.top.fgl$protein %in% c(df.degree.ucec.ordered.top.fgl$protein,df.degree.brca.ordered.top.fgl$protein))
which.common.brca.fgl = which(df.degree.brca.ordered.top.fgl$protein %in% df.degree.OVCA.ordered.top.fgl$protein & df.degree.brca.ordered.top.fgl$protein %in% df.degree.ucec.ordered.top.fgl$protein)
which.common.ucec.fgl = which(df.degree.ucec.ordered.top.fgl$protein %in% df.degree.brca.ordered.top.fgl$protein & df.degree.ucec.ordered.top.fgl$protein %in% df.degree.OVCA.ordered.top.fgl$protein)
which.common.OVCA.fgl = which(df.degree.OVCA.ordered.top.fgl$protein %in% df.degree.brca.ordered.top.fgl$protein & df.degree.OVCA.ordered.top.fgl$protein %in% df.degree.ucec.ordered.top.fgl$protein)


dim.neighs.top.fgl = max(c(length(df.degree.OVCA.ordered.top.fgl$protein),length(df.degree.ucec.ordered.top.fgl$protein),length(df.degree.OVCA.ordered.top.fgl$protein)))

# Also print gene name
for(i in 1:dim.neighs.top.fgl){
  if(is.na(df.degree.brca.ordered.top.fgl[i,1])) { cat(' & & && ') }
  else {
    if(i %in% which.common.brca.fgl){
      cat(paste0('\\textbf{',df.degree.brca.ordered.top.fgl$protein[i], '} & \\emph{',df.degree.brca.ordered.top.fgl$gene[i], '} & ', df.degree.brca.ordered.top.fgl$degree[i], ' && '))
    }
    else if(i %in% which.unique.brca.fgl){
      cat(paste0('\\color{red}{',df.degree.brca.ordered.top.fgl$protein[i], '} & \\emph{',df.degree.brca.ordered.top.fgl$gene[i], '} & ', df.degree.brca.ordered.top.fgl$degree[i], ' && '))
    }
    else { cat(paste0(df.degree.brca.ordered.top.fgl$protein[i], ' & \\emph{',df.degree.brca.ordered.top.fgl$gene[i], '} & ', df.degree.brca.ordered.top.fgl$degree[i], ' && ')) }
  }
  if(is.na(df.degree.ucec.ordered.top.fgl[i,1])) { cat(' & & && ') }
  else {
    if(i %in% which.common.ucec.fgl){
      cat(paste0('\\textbf{',df.degree.ucec.ordered.top.fgl$protein[i], '} & \\emph{',df.degree.ucec.ordered.top.fgl$gene[i], '} & ', df.degree.ucec.ordered.top.fgl$degree[i], ' && '))
    }
    else if(i %in% which.unique.ucec.fgl){
      cat(paste0('\\color{red}{',df.degree.ucec.ordered.top.fgl$protein[i], '} & \\emph{',df.degree.ucec.ordered.top.fgl$gene[i], '} & ', df.degree.ucec.ordered.top.fgl$degree[i], ' && '))
    }
    else { cat(paste0(df.degree.ucec.ordered.top.fgl$protein[i], ' & \\emph{',df.degree.ucec.ordered.top.fgl$gene[i], '} & ', df.degree.ucec.ordered.top.fgl$degree[i], ' && ')) }
  }
  if(is.na(df.degree.OVCA.ordered.top.fgl[i,1])) { cat(' & & \\\\ \n') }
  else {
    if(i %in% which.common.OVCA.fgl){
      cat(paste0('\\textbf{',df.degree.OVCA.ordered.top.fgl$protein[i], '} & \\emph{',df.degree.OVCA.ordered.top.fgl$gene[i], '} & ', df.degree.OVCA.ordered.top.fgl$degree[i], ' \\\\ \n  '))
    }
    else if(i %in% which.unique.OVCA.fgl){
      cat(paste0('\\color{red}{',df.degree.OVCA.ordered.top.fgl$protein[i], '} & \\emph{',df.degree.OVCA.ordered.top.fgl$gene[i], '} & ', df.degree.OVCA.ordered.top.fgl$degree[i], ' \\\\ \n  '))
    }
    else { cat(paste0(df.degree.OVCA.ordered.top.fgl$protein[i], ' & \\emph{',df.degree.OVCA.ordered.top.fgl$gene[i], '} & ', df.degree.OVCA.ordered.top.fgl$degree[i], ' \\\\ \n ')) }
  }
}

# 0.04184616


# Get edges with gene names
edges.brca.gene.fgl.all = get_and_print_edges(theta.est.brca.fgl!=0,mapping.frame$gene[match(colnames(brca_dat), mapping.frame$protein)])
edges.OVCA.gene.fgl.all = get_and_print_edges(theta.est.OVCA.fgl!=0,mapping.frame$gene[match(colnames(ovac_dat), mapping.frame$protein)])
edges.ucec.gene.fgl.all = get_and_print_edges(theta.est.ucec.fgl!=0,mapping.frame$gene[match(colnames(ucec_dat), mapping.frame$protein)])
edges.brca.gene.fgl = edges.brca.gene.fgl.all[-which(duplicated(edges.brca.gene.fgl.all)),]
edges.OVCA.gene.fgl = edges.OVCA.gene.fgl.all[-which(duplicated(edges.OVCA.gene.fgl.all)),]
edges.ucec.gene.fgl = edges.ucec.gene.fgl.all[-which(duplicated(edges.ucec.gene.fgl.all)),]
colnames(edges.brca.gene.fgl) =  c("Gene1", "Gene2")
colnames(edges.ucec.gene.fgl) =  c("Gene1", "Gene2")
colnames(edges.OVCA.gene.fgl) =  c("Gene1", "Gene2")
# Check with STRING
df.all.edges.brca.fgl =  rbind(edges.brca.gene.fgl, string.rppa[, 1:2])
ids.in.string.fgl = which(duplicated(df.all.edges.brca.fgl)) - nrow(edges.brca.gene.fgl) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with evidence
edges.brca.gene.evidence.fgl = string.rppa[ids.in.string.fgl,1:2]
n.extra.pairs.evidence.brca.fgl = 0
for(i in 1:nrow(edges.brca.gene.evidence.fgl)){
  n.rep = 0
  for(j in 1:nrow(edges.brca.gene.fgl.all)){
    if(all(edges.brca.gene.fgl.all[j,] == edges.brca.gene.evidence.fgl[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.brca.fgl = n.extra.pairs.evidence.brca.fgl + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string.fgl) + n.extra.pairs.evidence.brca.fgl # 408
(length(ids.in.string.fgl)+ n.extra.pairs.evidence.brca.fgl)/nrow(edges.brca.gene.fgl.all) # 0.06951781 for stabJGL: 0.1243144
# Also the amount of unique edges
length(ids.in.string.fgl)/nrow(edges.brca.gene.fgl) # 0.05414471 for stabJGL: 0.1233141
# Then UCEC
df.all.edges.ucec.fgl =  rbind(edges.ucec.gene.fgl, string.rppa[, 1:2])
ids.in.string.fgl = which(duplicated(df.all.edges.ucec.fgl)) - nrow(edges.ucec.gene.fgl) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with evidence
edges.ucec.gene.evidence.fgl = string.rppa[ids.in.string.fgl,1:2]
n.extra.pairs.evidence.ucec.fgl = 0
for(i in 1:nrow(edges.ucec.gene.evidence.fgl)){
  n.rep = 0
  for(j in 1:nrow(edges.ucec.gene.fgl.all)){
    if(all(edges.ucec.gene.fgl.all[j,] == edges.ucec.gene.evidence.fgl[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.ucec.fgl = n.extra.pairs.evidence.ucec.fgl + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string.fgl)+n.extra.pairs.evidence.ucec.fgl # 398
(length(ids.in.string.fgl)+n.extra.pairs.evidence.ucec.fgl)/nrow(edges.ucec.gene.fgl.all) # 0.06879862 for stabJGL: 0.108642
# Also the amount of unique edges
length(ids.in.string.fgl)/nrow(edges.ucec.gene.fgl) # 0.05635492 for stabJGL: 0.09974425
# Then OVCA
df.all.edges.OVCA.fgl =  rbind(edges.OVCA.gene.fgl, string.rppa[, 1:2])
ids.in.string.fgl = which(duplicated(df.all.edges.OVCA.fgl)) - nrow(edges.OVCA.gene.fgl) # True the second time an edge occurs (which is the duplicate). We subtract nrow to get the index of the edge in STRING df.
# Count all pairs with evidence
edges.OVCA.gene.evidence.fgl = string.rppa[ids.in.string.fgl,1:2]
n.extra.pairs.evidence.OVCA.fgl = 0
for(i in 1:nrow(edges.OVCA.gene.evidence.fgl)){
  n.rep = 0
  for(j in 1:nrow(edges.OVCA.gene.fgl.all)){
    if(all(edges.OVCA.gene.fgl.all[j,] == edges.OVCA.gene.evidence.fgl[i,])){
      n.rep = n.rep + 1
    }
  }
  n.extra.pairs.evidence.OVCA.fgl = n.extra.pairs.evidence.OVCA.fgl + (n.rep-1)
}
# Number of edges with evidence in STRING database
length(ids.in.string.fgl)+n.extra.pairs.evidence.OVCA.fgl # 424
(length(ids.in.string.fgl)+n.extra.pairs.evidence.OVCA.fgl)/nrow(edges.OVCA.gene.fgl.all) # 0.07026848 for stabJGL: 0.1257485
# Also the amount of unique edges
length(ids.in.string.fgl)/nrow(edges.OVCA.gene.fgl) # 0.05656849 for stabJGL: 0.1238095



# INVESTIGATE COMMON SUBNETWORK --------------------------------

subnet.common = (theta.est.brca!=0 & theta.est.OVCA!=0 & theta.est.ucec!=0) + 0
#colnames(subnet.common) = rownames(subnet.common) = mapping.frame$gene[match(colnames(subnet.common),mapping.frame$protein)]
sparsity(subnet.common) # 0.01855549

set.seed(111)
net.common = network::network(subnet.common ,directed=F)
net.plot = GGally::ggnet2(net.common, edge.size = 0.3,alpha=0.8,mode = "fruchtermanreingold",color = 'royalblue1',label=T,
                          size='degree',size.min = 1)+
            ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5, size=0.6))


pdf("plots/net_common_PanCan.pdf", 12, 12)
net.plot
dev.off()

# Can discuss top hubs a bit

# Focus on some relevant interactions: MIG6 (ERRFI1) and EGFR: 



# What about FGL? Can check for interactions
subnet.common.fgl = (theta.est.brca.fgl!=0 & theta.est.OVCA.fgl!=0 & theta.est.ucec.fgl!=0) + 0



# Plot similarity matrices  -------------------------------

# MCC of stabJGL estimates
tumor.names=c('BRCA','UCEC','OVCA')
thetas.est.all = list(theta.est.brca, theta.est.ucec,theta.est.OVCA)
thetas.est.all.fgl = list(theta.est.brca.fgl, theta.est.ucec.fgl,theta.est.OVCA.fgl)
MCC_vals_all = sapply(thetas.est.all, function(x) sapply(thetas.est.all, function(y) MCC(x!=0,y!=0))) # A matrix of pairwise MCC scores
MCC_vals_all[!lower.tri(MCC_vals_all)]=NA
#diag(MCC_vals_all)=NA
rownames(MCC_vals_all)=colnames(MCC_vals_all) =tumor.names

# MCC of FGL
MCC_vals_all_fgl= sapply(thetas.est.all.fgl, function(x) sapply(thetas.est.all.fgl, function(y) MCC(x!=0,y!=0))) # A matrix of pairwise MCC scores
MCC_vals_all_fgl[!lower.tri(MCC_vals_all_fgl)]=NA
#diag(MCC_vals_all_fgl)=NA
rownames(MCC_vals_all_fgl)=colnames(MCC_vals_all_fgl) =tumor.names

library(reshape2)
data_melt <- melt(MCC_vals_all,na.rm=T)  
ggp <- ggplot(data_melt, aes(Var1, Var2)) + geom_tile(aes(fill = value))+
  scale_fill_gradient2(low = "darkblue", high = "red", midpoint = 0.25, limit = c(0.01,0.53), space = "Lab",name="MCC") +
  theme_minimal()+ #theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  coord_fixed()+xlab('')+ylab('')+theme(text = element_text(size = 20),plot.title=ggplot2::element_text(size=25,hjust=0.5))+ggtitle('stabJGL')
ggp 

data_melt.fgl <- melt(MCC_vals_all_fgl ,na.rm=T)  
ggp.fgl  <- ggplot(data_melt.fgl , aes(Var1, Var2)) + geom_tile(aes(fill = value))+
  scale_fill_gradient2(low = "darkblue", high = "red", midpoint = 0.25, limit = c(0.01,0.53), space = "Lab",name="MCC") +
  theme_minimal()+ #theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1))+
  coord_fixed()+xlab('')+ylab('')+theme(text = element_text(size = 20),plot.title=ggplot2::element_text(size=25,hjust=0.5))+ggtitle('FGL')
ggp.fgl  

library(patchwork)
pdf("plots/MCC_plot.pdf", 8,4)
(ggp.fgl +theme(legend.position = 'none')) +ggp
dev.off()


# Plot degree distribution ----------------------------------------------



# Histogram
gg.hist = ggplot2::ggplot(df.degree, aes(degree, group=Group, colour=Group))+ geom_histogram(aes(fill=Group),stat='count')+ggtitle('(a)')+
  scale_color_manual(values=c("cornflowerblue","hotpink" ,"darkorchid"))+
  facet_wrap(~ Group)+theme_minimal()+scale_fill_manual(values=c("cornflowerblue","hotpink"  ,"darkorchid"))+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5),text = element_text(size = 11)) +xlim(-1,110)+ylim(-1,17)
gg.hist.fgl = ggplot2::ggplot(df.degree.fgl, aes(degree, group=Group, colour=Group))+ geom_histogram(aes(fill=Group),stat='count')+ggtitle('(b)')+
   scale_color_manual(values=c("cornflowerblue","hotpink" ,"darkorchid"))+
  facet_wrap(~ Group)+theme_minimal()+scale_fill_manual(values=c("cornflowerblue","hotpink"  ,"darkorchid"))+
  ggplot2::theme(plot.title=ggplot2::element_text(hjust=0.5),text = element_text(size = 11))+xlim(0,110)+ylim(-1,17)

pdf(paste0("plots/degreehist_PanCan_stabJGL.pdf"), 8, 3)
print((gg.hist+ theme(legend.position = 'none')) +gg.hist.fgl)
dev.off()



# Find edges FGL could not identify -----------------------------------


mm=which(((theta.est.ucec.fgl!=0) + (theta.est.brca.fgl!=0) +  (theta.est.OVCA.fgl!=0))<3 & ((theta.est.ucec!=0) +  (theta.est.brca!=0) +  (theta.est.OVCA!=0)) == 3,arr.ind=T)
kk=cbind(colnames(brca_dat)[mm[,1]],colnames(brca_dat)[mm[,2]])[,]
kk
cbind(mapping.frame$gene[match(kk[,1],mapping.frame$protein)],mapping.frame$gene[match(kk[,2],mapping.frame$protein)])

mm=which(theta.est.OVCA!=0 & theta.est.OVCA.fgl==0,arr.ind=T)
kk=cbind(colnames(brca_dat)[mm[,1]],colnames(brca_dat)[mm[,2]])[,]
kk
cbind(mapping.frame$gene[match(kk[,1],mapping.frame$protein)],mapping.frame$gene[match(kk[,2],mapping.frame$protein)])











