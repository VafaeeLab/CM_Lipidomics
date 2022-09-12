


# Read data ---------------------------------------------------------------

lipids <- read.table("./data/lipids", header = TRUE, sep = "\t")
raw    <- read.table("./data/raw_data", header = TRUE, sep = "\t", check.names = FALSE, fill = TRUE)
pheno  <- read.table("./data/metadata", header = TRUE, sep = "\t")
c.chain<- read.table("./data/carbon.chain", header = FALSE)


# Aggregate identical lipid ions ------------------------------------------

# PI(18:0_20:3)-H has duplicate row, with correlation 0.77, aggregated by mean
raw.m  <- aggregate(raw, by = list(raw$LipidIon), FUN = "mean")
rownames(raw.m) <- raw.m[,1]
raw.m <- raw.m[,c(-1,-2)]

# raw.m <- raw
# Normalisation -----------------------------------------------------------

library(limma)

# Meadian or quantile normalisation after log-2 transformation, adding an epsilon to avoid log(0)
norm <- normalizeMedianAbsValues(log2(raw.m + 1e-5))
# norm <- normalizeQuantiles(log2(as.matrix(raw.m + 1e-5)))

par(mfrow=c(2,2))
boxplot(log2(raw.m + 1e-5), main = "Not Normalised", ylab = "Log2 (ratio)")
boxplot(norm, main = "Normalised", ylab = "Log2 (ratio)")


# Dimentionality Reduction and Visualisation ------------------------------

# library(umap)
# norm.umap <- as.data.frame(umap(t(norm))$layout)
# norm.umap$class <- pheno$Class
# ggplot(norm.umap, aes(x=V1, y=V2, color = class)) + geom_point() + geom_text(label=rownames(norm.umap))

library(ggfortify)

# norm.pca <- as.data.frame(prcomp(t(norm))$x)
# norm.pca$class <- pheno$Class
# ggplot(norm.pca, aes(x=PC1, y=PC2, color = class)) + geom_point() + geom_text(label=rownames(norm.umap))

autoplot(prcomp(t(norm)), data = pheno, colour = "Class", frame = TRUE)+
  scale_color_manual(values=c("#00BA83", "#F8766D", "#619CFF"))+
  theme_light()


# library(tsne)
# norm.tsne <- as.data.frame(tsne(t(norm)))
# norm.tsne$class <- pheno$Class
# row.names(norm.tsne) <- pheno$Sample
# ggplot(norm.tsne, aes(x=V1, y=V2, color = class)) + geom_point() + geom_text(label=rownames(norm.umap))

# Heatmap of most variable lipid ions -------------------------------------

# subset profile of the top 50 most variable lipidions
norm.var <- norm[names(sort(apply(norm, 1, var), decreasing=TRUE))[1:300],]

require(RColorBrewer)
require(gplots)

# Get some nicer colours
morecols <- colorRampPalette(brewer.pal(11,"RdYlBu"))

# Set up colour vector for celltype variable
col.class <- c("skyblue","tomato", "orange")[pheno$Class]

# Plot the heatmap
heatmap.2(norm.var,col=rev(morecols(50)),trace="none", 
          main="Top 50 most variable lipid ions across samples",
          ColSideColors=col.class, scale="row", labRow = rownames(norm.var), Colv = NA, dendrogram = "row" )


# DE analysis -------------------------------------------------------------
require(gridExtra)
require(limma)
require(reshape2)

source("./VolcanoPlot.R")


design <- model.matrix(~ 0 + pheno$Class)
colnames(design) <- levels(pheno$Class)

cont.matrix <- makeContrasts(PbA - Control, Py - Control, PbA - Py,levels=design)

pval_cutoff = 0.01
logFC_cutoff = 1

fit <- lmFit(norm, design = design)
fit.cont <- contrasts.fit(fit, cont.matrix)
fit.cont <- eBayes(fit.cont)
summa.fit <- decideTests(fit.cont, lfc = logFC_cutoff, p.value = pval_cutoff)
summary(summa.fit)
vennDiagram(summa.fit, include=c("up", "down"), 
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3"))


par(mfrow=c(2,3))

res <- topTable(fit.cont, coef = 1, n = Inf)
VolcanoPlot(res, title = "PbA vs Control", abs_logFC_cutoff = logFC_cutoff, fdr_cutoff = pval_cutoff) #volcanoplot(fit.cont) #limma
A <- row.names(subset(res, adj.P.Val < pval_cutoff & abs(logFC)>logFC_cutoff))

res <- topTable(fit.cont, coef = 2, n = Inf)
VolcanoPlot(res, title = "Py vs Control", abs_logFC_cutoff = logFC_cutoff, fdr_cutoff = pval_cutoff)
B <- row.names(subset(res, adj.P.Val < pval_cutoff & abs(logFC)>logFC_cutoff))

res <- topTable(fit.cont, coef = 3, n = Inf)
VolcanoPlot(res, title = "PbA vs Py", abs_logFC_cutoff = logFC_cutoff, fdr_cutoff = pval_cutoff)
C <- row.names(subset(res, adj.P.Val < pval_cutoff & abs(logFC)>logFC_cutoff))


vennDiagram(summa.fit, include=c("up", "down"), 
            counts.col=c("red", "blue"),
            circle.col = c("red", "blue", "green3"))


venn(list(PbA_Control=A, Py_Control=B, PbA_Py=C))


bp <- function(ions, title){
  if(length(ions)>1){
    df = t(norm[ions,])
    row.names(df) = pheno$Class
    df = melt(df[,sort(colnames(df))])
    
  }else{
    df = as.data.frame(norm[ions,])
    colnames(df) = "value"
    df$Var2 <- rep(ions,nrow(df))
    df$Var1 <- pheno$Class
  }
     
   p<- ggplot(data = df, aes(x=Var2, y=value)) + geom_boxplot(aes(fill=Var1), outlier.shape=NA)+# + geom_jitter(width = 0.1) +  
    ggtitle(title) +  xlab("DE lipid Ions") + ylab("Value")+
  theme_light() + theme(axis.text.x = element_text(angle = 90), panel.grid.major = element_blank(), panel.grid.minor = element_blank())
   
  return(p)
}


for(t in sort(unique(lipids$LipidType))){
  print(t)
  p <- list()
    i = 1
  z = A[grepl(t, A)]
  if(length (z)>0){
    p[[i]] <- bp(z, paste(t, ", PbA vs Control"))
    i = i+1
  }
  z = B[grepl(t, B)]
  if(length (z)>0){
    p[[i]] <- bp(z, paste(t, ", Py vs Control"))
    i = i+1
  }
  z = C[grepl(t, C)]
  if(length (z)>0){
    p[[i]] <- bp(z, paste(t, ", PbA vs Py"))
    i = i+1
  }
  if(length(p)>0){do.call(grid.arrange,p)}
  
}


p = list(bp(A, "PbA vs Control"), bp(B, "Py vs Control"), bp(C, "PbA vs Py"))
do.call(grid.arrange,p)


source("./carbonEnrich.R")
out <- cbind(carbonEnrich(A, c.chain,lipids), carbonEnrich(B, c.chain,lipids), carbonEnrich(C, c.chain,lipids))
par(mfrow=c(1,1))
write.csv(out, file = "./results/carbon_chain_pvalues.csv")

# Lipid Classes -----------------------------------------------------------


df = t(norm)
row.names(df) = pheno$Class
df = melt(df[,sort(colnames(df))])
type <- sub("\\(.*", "", df[,2])
df$Type <- type
colnames(df) <- c("Class", "Ions", "value", "Type")

pval <- c()
p <- list()
i = 1
for(t in sort(unique(df$Type))){
  z <- subset(df, Type %in% t)
  pval <- rbind(pval,c(t.test(z[which(z$Class == "PbA"), "value"],z[which(z$Class == "Control"), "value"] )$p.value,
            t.test(z[which(z$Class == "Py"), "value"],z[which(z$Class == "Control"), "value"] )$p.value,
            t.test(z[which(z$Class == "PbA"), "value"],z[which(z$Class == "Py"    ), "value"] )$p.value))
  p[[i]] <- ggplot(data = z, aes(x=Class, y=value)) + geom_boxplot(width=0.3, aes(fill=Class), outlier.colour = NA) +
    scale_fill_manual(values=c("#00BA83", "#F8766D", "#619CFF"))+
    ggtitle(t) + xlab("")+ geom_jitter(width = 0.05, size = 0.3)+ 
    theme_light() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none") 
  i = i+1
}
do.call(grid.arrange,p)

rownames(pval) <- sort(unique(df$Type))
colnames(pval) <- c("PbA vs Control", "Py vs Control", "PbA vs Py")
z = melt(pval[,sort(colnames(pval))])
ggplot(data = z, aes(x=Var1, y=-log10(value))) + geom_bar(stat = "identity", aes(fill=Var2), position="dodge")+xlab("")+ylab("-l0g10(p-value)")+
  geom_hline(yintercept = -log10(0.05))
# ggplot(data = df, aes(fill=Type, y=value, x=Class)) + geom_bar(stat="identity", position="fill")

# Correlation analysis ----------------------------------------------------
iqr <- apply(norm,1, IQR)
hist(iqr, breaks = 50)
abline(v=1, col="blue")
abline(v=1.5, col="red")

keep <- iqr>1 #roughly equivalent to top quartile #quantile(iqr,0.75)
varied <- norm[keep,]
c <- cor(t(varied))#, method = "spearman")
diag(c) <- 0
require(corrplot)
pval <- cor.mtest(t(varied), conf.level = .95)

library(corrplot)
# corrplot(c, type = "full", p.mat = pval$p, order = "hclust", addrect = 8, tl.pos='n', method = "color", sig.level = .05, insig = "blank")
z<-corrplot(c, type = "full",  order = "hclust", method = "color", tl.col = "black", tl.cex = 0.45)# tl.srt = 45)
tmp = z
res <- topTable(fit.cont, coef = 1, n = Inf)
fc1 <- res[rownames(tmp),1]

res <- topTable(fit.cont, coef = 2, n = Inf)
fc2 <- res[rownames(tmp),1]

res <- topTable(fit.cont, coef = 3, n = Inf)
fc3 <- res[rownames(tmp),1]

fc <- cbind(fc1,fc3,fc2)
rownames(fc) <- rownames(z)
colnames(fc) <- c("PbA vs Control", "PbA vs Py", "Py vs Control")
fc[fc>1] = 1
fc[fc<1*-1] = -1

morecols <- colorRampPalette(colors = c("chartreuse3", "white", "firebrick2"))

# morecols <- colorRampPalette(colors = c("red", "white", "blue"))

heatmap.2(fc,trace="none", col=rev(morecols(100)), Rowv = FALSE, 
          main="FC",
          scale="none", labRow = rownames(fc), Colv = NA, dendrogram = "none" )

library(pvclust)
result <- pvclust(t(varied), method.dist="correlation", method.hclust="average", nboot=1000)
plot(result)
pvrect(result, alpha = 0.95)
seplot(result, identify = TRUE)


#Cytoscape:
net = c()
nodes = c()
for(i in 1:nrow(c)){
  for(j in i:ncol(c)){
    if(abs(c[i,j])>0.7){
      ion1 = gsub("\\+.*","",rownames(c)[i])
      ion2 = gsub("\\+.*","",rownames(c)[j])
      net = rbind(net,c(ion1,ion2,c[i,j]))
      nodes = c(nodes, ion1)
      nodes = c(nodes,ion2)
    }
  }
}
colnames(net) = c("Ion1", "Ion2", "Corr")
nodes = unique(nodes)
type = gsub("\\(.*","",nodes)
nodes = cbind(nodes,type)
colnames(nodes) = c("Ion","type")

write.table(net, file = "./results/net-0.7.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(nodes, file = "./results/nodes-0.7.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# past --------------------------------------------------------------------


rm.norm <- norm[!grepl("TG", rownames(norm)),]
norm.var <- rm.norm[names(sort(apply(rm.norm, 1, var), decreasing=TRUE))[1:50],]

norm.var <- norm[names(sort(apply(norm, 1, var), decreasing=TRUE))[1:50],]

c <- cor(t(norm.var))
diag(c) <- 0

thr = 0.7


library(corrplot)
corrplot(c, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

edges <- c()
for (i in 1:nrow(c)) {
  for(j in which(abs(c[i,])>thr)){
    edges <- rbind(edges, c(names(c[i,])[i], names(c[i,])[j], c[i,j]))
  }
}
colnames(edges) <- c("Lipid1", "Lipid2", "Correlation")
write.table(edges, file = paste("./results/net.corr_",thr,"_.txt"), row.names = FALSE, quote = FALSE, sep = "\t")


require(igraph)
library(network)
library(ggnetwork)

# graph representation of highly correlated lipids

delete.isolates <- function(graph, mode = 'all') {
  isolates <- which(degree(graph, mode = mode) == 0) 
  delete.vertices(graph, isolates)
}


c <-abs(c)
c[which(c<thr)] <- 0
g <- graph.adjacency(c, mode="undirected", weighted=TRUE) %>% 
  set_vertex_attr("name", value = rownames(norm))#%>%set_vertex_attr("label", value = d)
plot(g)
g2 <- delete.isolates(g)
plot(g2)

# Functional analysis ----------------------------------------------------


# Investigation of 22:6 ---------------------------------------------------

chain22.6 <- norm[grep(".22:6.",rownames(norm)),]


