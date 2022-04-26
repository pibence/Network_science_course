library(quantmod)
library(dplyr)
library(igraph)
library(ggplot2)
library(matrixcalc)


##############
##Exercise 1##
##############

#reading data
msg <- read.table(gzfile("CollegeMsg.txt.gz"))
colnames(msg) <- c("source","target","sec_since_epoch")
msg$day <- (msg$sec_since_epoch - min(msg$sec_since_epoch)) %/% (60*60*24) + 1

#making sure data is in increasing order in time
msg <- msg[order(msg$day,decreasing = F),]

#creating graph from all edges and counting all edges
total <- graph.data.frame(msg,directed = T)
Ne <- gsize(total)
Ve <- gorder(total)

#creating graph from first Ne/2 edges in edgelist
static <- graph.data.frame(msg[1:(Ne %/% 2 + 1),],directed = T)
#performing link prediction methods

#common neighbors
adj <- as.matrix(as_adjacency_matrix(static))
colnames(adj)[1:5]
adj2 <- matrix.power(adj,2)

CN <- get.data.frame(graph.adjacency(adj2,weighted = T))
CN <- CN[order(CN$weight,decreasing = T),]
head(CN)
nrow(CN)

#Adamic-Adar
AAinmatrix <- similarity(static,method = "invlogweighted",mode = "in",loops = F)
colnames(AAinmatrix) <- colnames(adj)
rownames(AAinmatrix) <- rownames(adj)

AAoutmatrix <- similarity(static,method = "invlogweighted",mode = "out",loops = F)
colnames(AAoutmatrix) <- colnames(adj)
rownames(AAoutmatrix) <- rownames(adj)

AAmatrix <- AAinmatrix + AAoutmatrix

AA <- get.data.frame(graph.adjacency(AAmatrix,weighted = T,diag = F))
AA <- AA[order(AA$weight,decreasing = T),]

#Local path index with parameter ß=0.1 ir order to give higher weight for shorter paths. source: https://www.nature.com/articles/s41598-020-76860-2
beta <- 0.1
LPmatrix <- adj2 + beta * matrix.power(adj,3)
LP <- get.data.frame(graph.adjacency(LPmatrix,weighted = T))
LP <- LP[order(LP$weight,decreasing = T),]
head(LP)

#evaluating link prediction methods
#since i created the static network from Ne%/%2+1 edges to make sure Ne/2 edges appear, now I am going to check for the residuum, Ne-(Ne%/%2+1)
Nehalf <- Ne - (Ne %/% 2 + 1)
#creating vector for all measures
fraction <- c()
#creating benchmark graph from original graph
BM <- graph.data.frame(msg[(Nehalf + 2):Ne,],directed = T)

#common neighbors
fraction[1] <- gsize(intersection(graph.data.frame(CN[1:Nehalf,1:2],directed = T),BM))/Nehalf

#Adamic-Adar
fraction[2] <- gsize(intersection(graph.data.frame(AA[1:Nehalf,1:2],directed = T),BM))/Nehalf

#Local path index
fraction[3] <- gsize(intersection(graph.data.frame(LP[1:Nehalf,1:2],directed = T),BM))/Nehalf

result <- data.frame(method = c("Common neighbors","Adamic-Adar", "Local path index"),fraction = fraction)

##############
##Exercise 2##
##############

date_start <- "2015-01-01"
date_end <- "2020-11-01"

components <- as.data.frame(read.csv("constituents.csv",sep = ";"))
samp_size <- 50
samp <- sample_n(components[,1:3],samp_size)
top <- components[1:samp_size,1:3]
data <- c()

#downloading necessary data
data <- as.data.frame(do.call(cbind,sapply(seq(1,samp_size), function(x) getSymbols(as.character(top$Symbol[x]), src = "yahoo", from = date_start, to = date_end, auto.assign = FALSE )[,6])))

#500 stock calculation block
#data_full <- data #for total sp500 calculation. once downloaded and saved
#top_full <- top #for total sp500 calculation. once downloaded and saved
#data <- data_full
#top <- top_full
#firstdates <- c("2000-01-01","2001-01-01","2002-01-01", "2003-01-01", "2004-01-01","2005-01-01","2006-01-01","2007-01-01","2008-01-01","2009-01-01","2010-01-01","2011-01-01","2012-01-01", "2013-01-01", "2014-01-01","2015-01-01","2016-01-01","2017-01-01","2018-01-01","2019-01-01", "2020-01-01", "2021-01-01")
####

ncol(data)
#return calculation
returns <- as.data.frame(sapply(data, function(x) diff(log(x))))

colnames(returns) <- top$Name
rownames(returns) <- as.Date(rownames(data[2:nrow(data),]))

firstdates <- c("2015-01-01","2016-01-01","2017-01-01","2018-01-01","2019-01-01", "2020-01-01", "2021-01-01")
#creating palette for visualization
cc <- c("orange","red","green3","blue","cyan","magenta","yellow","gray","brown","purple","beige")
color <- data.frame(factor=levels(top$Sector),color=as.character(cc))
#graph creation
mg <- list()
tree <-list()
palette(cc)

for (i in (1:(length(firstdates) - 1))){
  #calculating correlation and distance
  correlation <- cor(returns[rownames(returns) > firstdates[i] & rownames(returns) < firstdates[i+1],])
  distance <- 1 - correlation
  #creating graph from distance matrix ir order to obtain edgelist
  g  <- graph.adjacency(distance,weighted = TRUE,diag = F,mode = "upper")
  #getting edgelist from graph
  el <- get.data.frame(g)
  #creating increasing order in distance. 
  el <- el[order(el$weight,decreasing = F),]
  #creating new graph
  mg[[i]] <- graph.data.frame(el,directed = F)
  #creating attribute to nodes in order to draw them in differenct colors
  industry <- sapply(seq(1,50),function(x) top[top$Name == V(mg[[i]])$name[x],3])
  mg[[i]] <- set_vertex_attr(mg[[i]],name = "sector",value = industry)
  V(mg[[i]])$color <- sapply(1:gorder(mg[[i]]), function(x) color[color$factor == industry[x],2])
  #calculating minimum spanning tree  
  tree[[i]] <- mst(mg[[i]],algorithm = "prim")
  
}

#plotting minimum spanning trees
par(bg = "black")
for (i in (1:(length(firstdates) - 1))){
  plot(tree[[i]], 
       #general
       palette = cc,
       layout = layout_with_fr,
       #vertex attribuites
       vertex.size = 5,
       vertex.label.font = 1,
       vertex.label.color = V(mg[[i]])$color,
       vertex.label.dist = 1,
       vertex.label.cex = 0.75,
       
       #edge attributes
       edge.color = "orange",
       edge.curved = 0.3)
  title(paste("Minimum spanning tree of top",samp_size, "SP500 companies in",
              unlist(strsplit(firstdates[i],"-"))[1]), col.main = "azure", font.main = 2)
  legend("topleft", legend = color$factor, pch = 21, pt.bg = color$color,
         text.col = "white", cex = 0.75, bg = "transparent", pt.cex = 2,bty = "o",
         title = 'Sector')
}

#comparing minimum spanning trees by calculating Jaccard measure
jaccard <- c()
jaccard$date <- sapply(seq(1,(length(firstdates) - 2)), 
                       function(x) paste(unlist(strsplit(firstdates[x],"-"))[1],"-", unlist(strsplit(firstdates[x+1],"-"))[1]))
jaccard$jaccard <- sapply(seq(1,(length(firstdates)-2)), 
                          function(x) gsize(intersection(tree[[x]],tree[[x+1]]))/gsize(union(tree[[x]],tree[[x+1]])))
jaccard <- as.data.frame(jaccard)
#the higher the jaccard measure the higher the similarity

#plotting jaccard measure over the years
dev.off()
gg <- ggplot(data = jaccard, aes(x = date)) + geom_point(aes(y = jaccard), col = "red",size = 3) + 
  ggtitle("Jaccard similarity measure over the years") + theme(plot.title = element_text(hjust = 0.5))
#gg_full<-gg #once calculated and saved
gg_full 
gg2 <- ggplot(data = jaccard, aes(x = seq(1:20))) + geom_line(aes(y = jaccard), col = "red", size = 2) + 
  ggtitle("Jaccard similarity measure over the years") + theme(plot.title = element_text(hjust = 0.5))
gg2
#gg2_full <- gg2 #once calculated and saved

