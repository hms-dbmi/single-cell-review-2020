library(rgl)
library(pagoda2)
require(ggplot2)

cam <- readobj::read.obj('elephant.obj')
colnames(cam$shapes$polySurface1$positions) <- paste0('p',1:ncol(cam$shapes$polySurface1$positions))

# surface view
caml <- readobj::tinyobj2shapelist3d(cam)
shade3d(caml)


dotplot <- function(mat,color='lightgray',size=2,clear=T) {
  rgl.clear()
  rgl.points(mat[1,],mat[2,],mat[3,],color=color,size=size)
}

graphplot <- function(graph,mat,edge.color='gray50', ...) {
  # intersect with the vertices
  mat <- mat[,colnames(mat) %in% V(graph)$name]
  dotplot(mat, ...)
  x <- as.vector(t(as_edgelist(graph)))
  segments3d(x=mat[1,x],y=mat[2,x],z=mat[3,x],axes=FALSE, color=edge.color, plot=TRUE)
}




# downsampling
set.seed(0)
n.vertices <- 5000;
cs <- list(shapes=list(elephant=list(positions=cam$shapes$polySurface1$positions[,sample(1:ncol(cam$shapes$polySurface1$positions),n.vertices)])))

camp <- cs$shapes[[1]]$positions



# whole perspective
rgl.bg(color = "white")
rgl.viewpoint(theta=55,phi=05,zoom=0.8)
dotplot(camp)

# pagoda2

r <- Pagoda2$new(camp*1e3)
r$reductions <- list(PCA=t(camp))
r$makeKnnGraph(k=5,type='PCA',distance='L2');


r$getEmbedding(type='PCA',embeddingType='UMAP',perplexity=3,distance='L2')

M<- 10; emb <- r$getEmbedding(type='PCA',embeddingType='largeVis',distance='L2',gamma=1/10,perplexity=100,M=M)
set.seed(1)

set.seed(0)
emb <- r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=100,distance='L2')

emb <- r$embeddings$PCA$UMAP
emb <- r$embeddings$PCA$largeVis
emb <- r$embeddings$PCA$tSNE


conos::embeddingPlot(emb,colors=col,size=0.1,alpha=0.2,raster=T,raster.height = 3,raster.width=3)+ theme_void() + theme(panel.border = element_rect(color = 1, size=0.2,fill=NA,linetype=1),axis.line=element_blank())

set.seed(0)
r$getKnnClusters(method=function(x) conos:::leiden.community(x,r=0.04),type='PCA',min.cluster.size=20)
r$getKnnClusters(method=function(x) conos:::leiden.community(x,r=1),type='PCA',min.cluster.size=20)
f1 <- r$clusters$PCA$community



plot(emb[,1],emb[,2],col=sccore:::fac2col(f1)[rownames(emb)],cex=0.2,pch=19)
r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='clusters (lV)')

dotplot(camp,sccore:::fac2col(f1)[colnames(camp)])



# trunk view

trunk.points <- names(f1)[f1 %in% c('1')]
tcamp <- camp[,trunk.points]

# figure out a perspective for the trunk
nat:::pan3d(2)


t.zoom<-par3d()$zoom
t.userMatrix<-par3d()$userMatrix
t.windowRect<-par3d()$windowRect

# draw trunk wireframe

# subsample again
set.seed(0)
tscamp <- tcamp[,sample(1:ncol(tcamp),300)]

# sparsely connectedn graph
tr <- Pagoda2$new(tscamp*1e3,min.cells.per.gene=-1,min.transcripts.per.cell=-1,trim=0)
tr$reductions <- list(PCA=t(tscamp))
tr$makeKnnGraph(k=5,type='PCA',distance='L2');

graphplot(tr$graphs$PCA,camp,color=1,edge.color='gray50',size=3)
rgl.viewpoint(theta=55,phi=05,zoom=0.8)
rgl.postscript('graph1.pdf',fmt='pdf')


# whole elephant
trunk.points <- names(f1)[f1 %in% c('1','3','12','5','8','2','8')]
tcamp <- camp[,trunk.points]
tcamp <- camp

set.seed(4)
tscamp <- tcamp[,sample(1:ncol(tcamp),1000)]

# sparsely connectedn graph
tr <- Pagoda2$new(tscamp*1e3,min.cells.per.gene=-1,min.transcripts.per.cell=-1,trim=0)
tr$reductions <- list(PCA=t(tscamp))
tr$makeKnnGraph(k=5,type='PCA',distance='L2');

graphplot(tr$graphs$PCA,camp,color=1,edge.color='gray50',size=3)
rgl.viewpoint(theta=55,phi=05,zoom=0.8)

rgl.postscript('graph2.pdf',fmt='pdf')


# cluster views
## res1
set.seed(0)
f1 <- r$getKnnClusters(method=function(x) conos:::leiden.community(x,r=0.04),type='PCA',min.cluster.size=20)
f1 <- r$clusters$PCA$community

set.seed(6)
col <- sccore:::fac2col(f1,s=0.8,v=0.7,shuffle=T)[colnames(camp)]
dotplot(camp,col)
rgl.postscript('clust1.pdf',fmt='pdf')

#r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,groups=f1)
a1 <- conos::embeddingPlot(emb,colors=col,size=0.5,alpha=0.2,raster=T,raster.height = 3,raster.width=3)+ theme_void() + theme(panel.border = element_rect(color = 1, size=0.2,fill=NA,linetype=1),axis.line=element_blank())
a1
pdf(file='clust1.tsne.pdf',width=2,height=2); print(a1); dev.off()

## res2
set.seed(0)
f1 <- r$getKnnClusters(method=function(x) conos:::leiden.community(x,r=0.5),type='PCA',min.cluster.size=20)
f1 <- r$clusters$PCA$community

set.seed(6)
col <- sccore:::fac2col(f1,s=0.8,v=0.7,shuffle=T)[colnames(camp)]
dotplot(camp,col)
rgl.postscript('clust2.pdf',fmt='pdf')

#r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,groups=f1)
a1 <- conos::embeddingPlot(emb,colors=col,size=0.5,alpha=0.2,raster=T,raster.height = 3,raster.width=3)+ theme_void() + theme(panel.border = element_rect(color = 1, size=0.2,fill=NA,linetype=1),axis.line=element_blank())
a1
pdf(file='clust2.tsne.pdf',width=2,height=2); print(a1); dev.off()


## res3
set.seed(0)
f1 <- r$getKnnClusters(method=function(x) conos:::leiden.community(x,r=3),type='PCA',min.cluster.size=20)
f1 <- r$clusters$PCA$community

set.seed(6)
col <- sccore:::fac2col(f1,s=0.8,v=0.7,shuffle=T)[colnames(camp)]
dotplot(camp,col)
rgl.postscript('clust3.pdf',fmt='pdf')

#r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,groups=f1)
a1 <- conos::embeddingPlot(emb,colors=col,size=0.5,alpha=0.2,raster=T,raster.height = 3,raster.width=3)+ theme_void() + theme(panel.border = element_rect(color = 1, size=0.2,fill=NA,linetype=1),axis.line=element_blank())
a1
pdf(file='clust3.tsne.pdf',width=2,height=2); print(a1); dev.off()

