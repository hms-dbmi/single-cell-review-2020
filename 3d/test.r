library(rgl)

xd <- read.table("camel_b.obj.txt",colClasses = c('NULL','numeric','numeric','numeric'),header=F,sep=' ')

cam <- readobj::read.obj('camel_b.obj.txt')

cam <- readobj::read.obj('elephant.obj')


caml <- readobj::tinyobj2shapelist3d(cam)
shade3d(caml)

camp <- cam$shapes[[1]]$positions
rgl.open()
rgl.bg(color = "white")
rgl.points(camp[1,], camp[2,], camp[3,], color ="lightgray",size=2)

rgl.bg(color = "white") # Setup the background color
rgl.points(x, y, z, color = "blue", size = 5) # Scatter plot


library(pagoda2)


cam <- readobj::read.obj('elephant.obj')

colnames(cam$shapes$polySurface1$positions) <- paste0('p',1:ncol(cam$shapes$polySurface1$positions))


# create a downsampled set

set.seed(0)
n.vertices <- 5000;
cs <- list(shapes=list(elephant=list(positions=cam$shapes$polySurface1$positions[,sample(1:ncol(cam$shapes$polySurface1$positions),n.vertices)])))

camp <- cs$shapes[[1]]$positions

# view
rgl.open()
rgl.bg(color = "white")
rgl.points(camp[1,], camp[2,], camp[3,], color ="lightgray",size=2)

rgl.viewpoint(theta=55,phi=05,zoom=0.8)
rgl.viewpoint(theta=55,phi=05)
rgl.clear()
rgl.points(camp[1,], camp[2,], camp[3,], color =sccore:::fac2col(f1)[colnames(camp)],size=2)

rgl.postscript('proj1.pdf',fmt='pdf')


# zoom in view

# clustering
## make p2 object
require(pagoda2)
require(conos)
r <- Pagoda2$new(camp*1e3)
r$reductions <- list(PCA=t(camp))



r$makeKnnGraph(k=50,type='PCA',distance='L2');
r$getKnnClusters(method=function(x) leiden.community(x,r=1),type='PCA')

f1 <- r$clusters$PCA$community
table(f1)
tf1 <- table(f1)
f1[f1 %in% names(tf1)[tf1<5]] <- NA;
f1 <- droplevels(f1)
table(f1)

r$getEmbedding(type='PCA',embeddingType='tSNE',perplexity=200,verbose=F,distance='L2')
r$getEmbedding(type='PCA',embeddingType='largeVis',distance='L2')

emb <- r$embeddings$PCA$tSNE
emb <- r$embeddings$PCA$largeVis
plot(emb[,1],emb[,2],col=sccore:::fac2col(f1)[rownames(emb)],cex=0.2,pch=19)

r$plotEmbedding(type='PCA',show.legend=F,mark.clusters=T,min.group.size=1,shuffle.colors=F,mark.cluster.cex=1,alpha=0.1,main='clusters (lV)')
## different resolutions

# tSNE



# take trunk (clusters 5,2)
trunk.points <- names(f1)[f1 %in% c('5','2')]
tcamp <- camp[,trunk.points]
tpos <- rowMeans(tcamp)
tcamp <- tcamp-rowMeans(tcamp)

tcam <- cam;
tcam$shapes <- lapply(tcam$shapes, function(x) {
  x$positions <- x$positions-tpos;
  x
})
str(tcam)


tcaml <- readobj::tinyobj2shapelist3d(tcam)
rgl.clear()
shade3d(tcaml)
rgl.viewpoint(theta=55,phi=05,zoom=0.2)


nat:::pan3d(2)

t.zoom<-par3d()$zoom
t.userMatrix<-par3d()$userMatrix
t.windowRect<-par3d()$windowRect

open3d(zoom = t.zoom, userMatrix = t.userMatrix, windowRect=t.windowRect)

perps3d(x=x,y=y,z=z) # plus whatever arguments you need, but ignoring all perspective arguments
rgl.snapshot( filename, fmt="png", top=TRUE)

rgl.viewpoint(zoom=t.zoom,userMatrix=t.userMatrix)

rgl.points(camp[1,], camp[2,], camp[3,], color =sccore:::fac2col(f1)[colnames(camp)],size=2)
rgl.points(tcamp[1,], tcamp[2,], tcamp[3,], color =sccore:::fac2col(f1)[colnames(tcamp)],size=2)
