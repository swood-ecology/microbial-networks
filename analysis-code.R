# ----------------------------------------------------------------------
# Contingency in soil microbial response to tropical forest conversion
# Analaysis by Steve Wood
# Date Started: 12 Feb 2016
# ----------------------------------------------------------------------

# MAJOR ANALYSIS COMPONENTS
# 1. How does bacterial diversity vary by land cover type?
# 2. Which phyla dominate in different land use types?
# 3. How does cross-site similarity differ between bacterial and fungal communities? 
# 4. Are there tight couplings between particular bacterial types?
# 5. Does land use disrupt co-occurence patterns of bacteria and fungi?


# NEEDED PACKAGES
source("~/Documents/Work/Statistics/Code/setRows.R")  # for setting rows
source("~/Documents/Work/Statistics/Code/CAst.R")     # for abundance weighted C score
library(vegan)      # for diversity indices
library(ggplot2)    # for plotting
library(MASS)       # for box-cox transformation
library(reshape)    # for shaping data sets
library(netassoc)   # network analyses
library(igraph)     # graphing networks
library(network)    # for network analysis
library(rgexf)      # converts igraph objects to gexf for export to Gephi
library(picante)    # for phylogeny-based analyses
# source("https://bioconductor.org/biocLite.R")
# biocLite("ggtree")


# READ IN DATA
setwd("~/Dropbox/Pasoh_Bacteria/Bacterial_analyses/Steve_analyses")

OTU <- read.csv("OTU.csv",header=T)
OTU <- set.rownames(OTU)
OTU <- t(OTU)
Phyla <- read.csv("PhylaRA.csv",header=T)
Phyla <- set.rownames(Phyla)
Phyla <- t(Phyla)

mapping <- read.csv("Mapping.csv",header=T)
mapping <- set.rownames(mapping)


# SECTION 1
# How does bacterial diversity vary by land cover type?

# 1.1: NMDS (unconstrained ordination)

# ANOSIM by soil depth to assess assumption of pooling depths
plot(anosim(OTU,grouping=mapping$Horizon))  # no sig. diff. b/t depths

# ANOSIM and NMDS by forest and forest type
plot(anosim(OTU,grouping=mapping$Forest_Oilpalm))  
plot(anosim(OTU,grouping=mapping$Site.Name))

nmds <- metaMDS(OTU,distance="bray",k=2,trymax=100)

data.scores <- as.data.frame(scores(nmds))  
data.scores$site <- rownames(data.scores) 
data.scores$grp <- mapping$Forest_Oilpalm  
species.scores <- as.data.frame(scores(nmds, "species"))  
species.scores$species <- rownames(species.scores)

ggplot() + 
  geom_point(data=data.scores,aes(x=NMDS1,y=NMDS2,shape=grp,colour=grp),size=4) + # add the point markers
  coord_equal() +
  theme_bw() +
  theme(axis.title.x=element_text(size=18),
        axis.title.y=element_text(size=18),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.text=element_text(size=16),
        legend.title=element_blank(),
        legend.position="top",
        panel.border=element_rect(color="black",size=0.75))


# 1.2: Diversity Analysis

# Species accumulation curves
SAC <- specaccum(OTU,"random")
summary(SAC)
plot(SAC, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue")
boxplot(SAC, col="yellow",add=T,pch="+")

# Diversity calculations
S <- specnumber(OTU)
H <- diversity(OTU)
simp <- diversity(OTU, "simpson")
invsimp <- diversity(OTU, "inv")
alpha <- fisher.alpha(OTU)
J <- H/log(S)
pairs(cbind(H,simp,invsimp,alpha,S,J),pch="+",col="blue")

div.data <- data.frame(S,H,simp,invsimp,alpha,J,mapping)

shapiro.test(S)        # non-normal
shapiro.test(H)        # normal
shapiro.test(simp)     # non-normal
shapiro.test(invsimp)  # non-normal
shapiro.test(alpha)    # non-normal
shapiro.test(J)        # normal

# code for boxcox transformation
b.c <- function(var,lambda){
  return(((var^lambda)-1)/lambda)
}

H.model <- lm(H~Forest_Oilpalm,data=div.data)
summary(H.model)
J.model <- lm(J~Forest_Oilpalm,data=div.data)
summary(J.model)

boxplot(H~Forest_Oilpalm,data=div.data)


boxcox(lm(S~Forest_Oilpalm,data=div.data),plotit=FALSE)
S.model <- lm(b.c(S,lambda=-0.9)~Forest_Oilpalm,data=div.data)
summary(S.model)

boxcox(lm(simp~Forest_Oilpalm,data=div.data),plotit=FALSE)
simp.model <- lm(b.c(simp,lambda=-2.0)~Forest_Oilpalm,data=div.data)
summary(simp.model)

boxcox(lm(invsimp~Forest_Oilpalm,data=div.data),plotit=FALSE)
invsimp.model <- lm(b.c(invsimp,lambda=-0.2)~Forest_Oilpalm,data=div.data)
summary(invsimp.model)

boxcox(lm(alpha~Forest_Oilpalm,data=div.data),plotit=FALSE)
alpha.model <- lm(b.c(alpha,lambda=-0.7)~Forest_Oilpalm,data=div.data)
summary(alpha.model)




# SECTION 2
# Which phyla dominate in different land use types?

Phyla.LM <- read.csv("PhylaRA.csv",header=T)
Phyla.LM <- melt(Phyla.LM,id.vars='Phylum')
Phyla.LM <- cast(Phyla.LM,variable~Phylum)
names(Phyla.LM)[1] <- "Sample.ID"
mapping <- read.csv("Mapping.csv",header=T)
Phyla.LM <- merge(Phyla.LM,mapping,by='Sample.ID')

summary(lm(Phyla.LM$Chlorobi ~ Phyla.LM$Forest_Oilpalm))

# Significantly impacted phyla: Acidobacteria(-),Actinobacteria(+),BRC1(+),
# Chlamydiae(+),Chlorobi(+),Chloroflexi(+),Crenarchaeota(+),Cyanobacteria(+),
# Elusimicrobia(+),FBP(+),FCPU426(+),Firmicutes(+),NC10(+),Nitrospirae(+),OD1(+),
# OP11(+),OP3(+),Spirochaetes(+),Tenericutes(+),TM6(+),TM7(+),Verrucomicrobia(+),
# ZB3(+)


Phyla.Sub.1 <- Phyla[,c(4,5,13,21,39)]
Phyla.Sub.2 <- Phyla[,c(11,14:16,19,27,35)]
                     
Phyla.Plot.1 <- melt(Phyla.Sub.1)
names(Phyla.Plot.1)[c(1,2)] <- c("Sample.ID","Taxon")
mapping <- read.csv("Mapping.csv",header=T)
Phyla.New.1 <- merge(Phyla.Plot.1,mapping,by='Sample.ID')
Phyla.New.1 <- Phyla.New.1[,c(2:4)]

Phyla.Ag.1 <- aggregate(. ~ Taxon + Site.Name, Phyla.New.1, function(x) c(mean = mean(x), se = sqrt(var(x)/length(x))))
Phyla.Ag.1 <- cbind(Phyla.Ag.1,Phyla.Ag.1[,3][,2])
Phyla.Ag.1[3] <- Phyla.Ag.1[,3][,1]
names(Phyla.Ag.1)[c(3,4)] <- c("value.mean","value.se")

Phyla.Plot.2 <- melt(Phyla.Sub.2)
names(Phyla.Plot.2)[c(1,2)] <- c("Sample.ID","Taxon")
mapping <- read.csv("Mapping.csv",header=T)
Phyla.New.2 <- merge(Phyla.Plot.2,mapping,by='Sample.ID')
Phyla.New.2 <- Phyla.New.2[,c(2:4)]

Phyla.Ag.2 <- aggregate(. ~ Taxon + Site.Name, Phyla.New.2, function(x) c(mean = mean(x), se = sqrt(var(x)/length(x))))
Phyla.Ag.2 <- cbind(Phyla.Ag.2,Phyla.Ag.2[,3][,2])
Phyla.Ag.2[3] <- Phyla.Ag.2[,3][,1]
names(Phyla.Ag.2)[c(3,4)] <- c("value.mean","value.se")


#plot the stacked bar plot with polar coordinates
# ggplot(Phyla.New.2, aes(x = Site.Name)) + 
#   geom_bar(aes(weight=value, fill = Taxon), position = 'fill') + 
#   scale_y_continuous("", breaks=NULL) + coord_polar() + scale_fill_hue(l=30, c=120) +
#   theme_bw() +
#   theme(axis.title.x=element_blank(),
#         axis.title.y=element_blank(),
#         panel.grid.major=element_blank(),
#         panel.grid.minor=element_blank(),
#         legend.text=element_text(size=13),
#         legend.title=element_blank(),
#         panel.border=element_rect(color="black",size=0.75))

# CLUSTERED BAR PLOT
cbbPalette <- c("#000000", "#56B4E9", "#009E73", "#0072B2", "#D55E00")

ggplot(Phyla.Ag.1, aes(x=Taxon,y=value.mean,fill=Site.Name)) + 
  geom_bar(position=position_dodge(),stat="identity") +
  geom_errorbar(aes(ymin=value.mean-value.se, ymax=value.mean+value.se),
                position=position_dodge(.9),width=0.2,color="grey") +
  theme_bw() + ylab("Relative abundance\n") + 
  scale_fill_manual(values=cbbPalette) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=14),
        axis.text.y=element_text(size=12),
        axis.ticks = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        legend.position="top",
        panel.border=element_rect(color="black",size=0.75))

#cbbPalette2 <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
ggplot(Phyla.Ag.2, aes(x=Taxon,y=value.mean,fill=Site.Name)) + 
  geom_bar(position=position_dodge(),stat="identity") +
  geom_errorbar(aes(ymin=value.mean-value.se, ymax=value.mean+value.se),
                position=position_dodge(.9),width=0.2,color="grey") +
  theme_bw() + ylab("Relative abundance\n") + 
  scale_fill_manual(values=cbbPalette) +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size=14),
        axis.text.x=element_text(size=14,angle=45, hjust=1),
        axis.text.y=element_text(size=12),
        axis.ticks = element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        legend.text=element_text(size=12),
        legend.title=element_blank(),
        legend.position="top",
        panel.border=element_rect(color="black",size=0.75))





# SECTION 4
# 4. Are there tight couplings between particular bacterial types?

# another approach
network <- read.csv("OTU.Phyla.csv",header=T)
network <- aggregate(. ~ Phylum,data=network[,-1],sum) 
network <- set.rownames(network)

network.all <- make_netassoc_network(network,numnulls=1000)
network.all <- network.all$network_all

network.primary <- network[,c(16:24)]
network.logged <- network[,c(1:8)]
network.forest <- network[,c(1:8,16:24)]
network.palm <- network[,9:15]

net.forest <- make_netassoc_network(network.forest[c(-8,-33,-35),],numnulls=1000)
net.forest <- net.forest$network_all

net.logged <- make_netassoc_network(network.logged[c(-8,-24,-33,-35,-42),],numnulls=1000)
logged.heat <- net.logged$matrix_spsp_ses_thresholded
net.logged <- net.logged$network_all

net.primary <- make_netassoc_network(network.primary[c(-6,-8,-15,-16,-17,-33,-35),],numnulls=1000)
primary.heat <- net.primary$matrix_spsp_ses_thresholded
net.primary <- net.primary$network_all

E(net.forest)$width <- E(net.forest)$weight
E(net.forest)$edge.color <- "gray80"
V(net.forest)$label.cex=.75
V(net.forest)$label.font=1.5

net.palm <- make_netassoc_network(network.palm[-26,],numnulls=1000)
palm.heat <- net.palm$matrix_spsp_ses_thresholded
net.palm <- net.palm$network_all

E(net.palm)$width <- E(net.palm)$weight
E(net.palm)$edge.color <- "gray80"
V(net.palm)$label.cex=.65
V(net.palm)$label.font=2
deg.palm <- degree(net.palm, mode="all")
V(net.palm)$size <- deg.palm

E(net.logged)$width <- E(net.logged)$weight
E(net.logged)$edge.color <- "gray80"
V(net.logged)$label.cex=.65
V(net.logged)$label.font=2
deg.logged <- degree(net.logged, mode="all")
V(net.logged)$size <- deg.logged

E(net.primary)$width <- E(net.primary)$weight
E(net.primary)$edge.color <- "gray80"
V(net.primary)$label.cex=.65
V(net.primary)$label.font=2
V(net.primary)$size <- V(net.primary)
deg.primary <- degree(net.primary, mode="all")
V(net.primary)$size <- deg.primary


wtc <- cluster_walktrap(net.palm)
plot(wtc,net.palm,edge.arrow.size=0,edge.curved=.1,vertex.frame.color="red",
     vertex.label.color="black",vertex.color="white",main="Palm Oil")

plot(simplify(net.palm),edge.arrow.size=0,vertex.frame.color="red",
     vertex.label.color="black",vertex.color="white",main="Palm Oil",layout=layout_nicely)

plot(simplify(net.logged),edge.arrow.size=0,vertex.frame.color="red",
     vertex.label.color="black",vertex.color="white",main="Logged Forest",layout=layout_nicely)

plot(simplify(net.primary),edge.arrow.size=0,vertex.frame.color="red",
     vertex.label.color="black",vertex.color="white",main="Primary Forest",layout=layout_nicely)



# print(igraph.to.gexf(net.palm),"PalmNetwork.gexf",replace=T)
# print(igraph.to.gexf(net.forest),"ForestNetwork.gexf",replace=T)


# network metrics
# explained in: http://download.springer.com/static/pdf/321/chp%253A10.1007%252F978-1-4939-0983-4_4.pdf?originUrl=http%3A%2F%2Flink.springer.com%2Fchapter%2F10.1007%2F978-1-4939-0983-4_4&token2=exp=1457466968~acl=%2Fstatic%2Fpdf%2F321%2Fchp%25253A10.1007%25252F978-1-4939-0983-4_4.pdf%3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com%252Fchapter%252F10.1007%252F978-1-4939-0983-4_4*~hmac=56c17ec2c1c5f5840bb018841b48f261e4085ae6cffe619f42179733dd2f024b


max_cliques(net.palm)
cliques(net.logged, min=4)
cliques(net.primary, min=4)
 
largest_cliques(net.palm)
largest_cliques(net.logged)
largest_cliques(net.primary)


edge_density(net.palm,loops=T)
edge_density(net.logged,loops=T)
edge_density(net.primary,loops=T)

transitivity(net.palm)
transitivity(net.logged)
transitivity(net.primary)

vertex_connectivity(net.palm)
vertex_connectivity(net.logged)
vertex_connectivity(net.primary)

assortativity_degree(net.palm)
assortativity_degree(net.logged)
assortativity_degree(net.primary)

wtc <- cluster_walktrap(net.palm)
modularity(net.palm,membership(wtc))
plot(wtc,net.palm,edge.arrow.size=0)

wtc <- cluster_walktrap(net.logged)
modularity(net.logged,membership(wtc))
plot(wtc,net.logged,edge.arrow.size=0)

wtc <- cluster_walktrap(net.primary)
modularity(net.primary,membership(wtc))
plot(wtc,net.primary,edge.arrow.size=0)

plot_dendrogram(fit_hrg(net.forest),mode="phylo")
plot_dendrogram(fit_hrg(net.palm),mode="phylo")

betweenness(net.palm,weights=abs(E(net.palm)$weight))
betweenness(net.logged,weights=abs(E(net.logged)$weight))
betweenness(net.primary,weights=abs(E(net.primary)$weight))



# Compare network statistics to randomized networks

mod <- vector()
ass <- vector()
trans <- vector()
edge <- vector()
vert <- vector()

# generate 10,000 random networks based on number of edges and vertices
# in the regional network
for(i in 1:10000){
  x <- erdos.renyi.game(gorder(network.all),ecount(network.all)/2,type="gnm")
  mod[i] <- modularity(x,membership(cluster_walktrap(x)))
  ass[i] <- assortativity_degree(x)
  trans[i] <- transitivity(x)
  vert[i] <- vertex_connectivity(x)
  edge[i] <- edge_density(x)
}

# generate 10,000 random networks based on number of edges and vertices
# in the palm network
for(i in 1:10000){
  x <- erdos.renyi.game(gorder(net.palm),ecount(net.palm)/2,type="gnm")
  mod[i] <- modularity(x,membership(cluster_walktrap(x)))
  ass[i] <- assortativity_degree(x)
  trans[i] <- transitivity(x)
  vert[i] <- vertex_connectivity(x)
  edge[i] <- edge_density(x)
}

# generate 10,000 random networks based on number of edges and vertices
# in the logged network
for(i in 1:10000){
  x <- erdos.renyi.game(gorder(net.logged),ecount(net.logged)/2,type="gnm")
  mod[i] <- modularity(x,membership(cluster_walktrap(x)))
  ass[i] <- assortativity_degree(x)
  trans[i] <- transitivity(x)
  vert[i] <- vertex_connectivity(x)
  edge[i] <- edge_density(x)
}

# generate 10,000 random networks based on number of edges and vertices
# in the primary network
for(i in 1:10000){
  x <- erdos.renyi.game(gorder(net.primary),ecount(net.primary)/2,type="gnm")
  mod[i] <- modularity(x,membership(cluster_walktrap(x)))
  ass[i] <- assortativity_degree(x)
  trans[i] <- transitivity(x)
  vert[i] <- vertex_connectivity(x)
  edge[i] <- edge_density(x,loops=T)
}

mean(mod) #mean modularity index of random networks
sd(mod) #standard deviation
mean(ass)
sd(ass)
mean(trans)
sd(trans)
mean(vert)
sd(vert)
mean(edge)
sd(edge)

z.mod.palm <- (modularity(net.palm,membership(cluster_walktrap(net.palm)))-mean(mod))/sd(mod)
z.mod.logged <- (modularity(net.logged,membership(cluster_walktrap(net.logged)))-mean(mod))/sd(mod)


# #determine how many of the random networks have modularity scores #lower than observed network
# ebc <- cluster_edge_betweenness(network.all,weights=abs(E(network.all)$weight),merges=T)
# sum(mod<=max(ebc$modularity))


# Look at pairings
network.all <- make_netassoc_network(network,numnulls=1000)
net.all.assoc <- network.all$matrix_spsp_ses_thresholded
net.all.assoc <- net.all.assoc[,]
net.all.assoc


# z<-lm(Phyla.LM$ZB3~Phyla.LM$Actinobacteria)
# plot(Phyla.LM$Actinobacteria,Phyla.LM$ZB3,xlab="Actinobacteria",ylab="ZB3")
# abline(z)
# cor(Phyla.LM$Actinobacteria,Phyla.LM$ZB3)
# text(0.0005,0.00025,"r = 0.19")
# 
# z<-lm(Phyla.LM$ZB3~Phyla.LM$Planctomycetes)
# plot(Phyla.LM$Planctomycetes,Phyla.LM$ZB3,xlab="Planctomycetes",ylab="ZB3")
# abline(z)
# cor(Phyla.LM$Planctomycetes,Phyla.LM$ZB3)
# text(0.045,0.00025,"r = 0.42")
# 
# z<-lm(Phyla.LM$NKB19~Phyla.LM$Bacteroidetes)
# plot(Phyla.LM$Bacteroidetes,Phyla.LM$NKB19,xlab="Bacteroidetes",ylab="NKB19")
# abline(z)
# cor(Phyla.LM$Bacteroidetes,Phyla.LM$NKB19)
# text(0.07,0.0006,"r = 0.35")
# 
# z<-lm(Phyla.LM$NKB19~Phyla.LM$Fibrobacteres)
# plot(Phyla.LM$Fibrobacteres,Phyla.LM$NKB19,xlab="Fibrobacteres",ylab="NKB19")
# abline(z)
# cor(Phyla.LM$Fibrobacteres,Phyla.LM$NKB19)
# text(0.0009,0.0006,"r = 0.91")
# 
# z<-lm(Phyla.LM$NKB19~Phyla.LM$`[Parvarchaeota]`)
# plot(Phyla.LM$`[Parvarchaeota]`,Phyla.LM$NKB19,xlab="Parvarchaeota",ylab="NKB19")
# abline(z)
# cor(Phyla.LM$`[Parvarchaeota]`,Phyla.LM$NKB19)
# text(0.0009,0.0006,"r = 0.77")
# 
# z<-lm(Phyla.LM$BRC1~Phyla.LM$`[Caldithrix]`)
# plot(Phyla.LM$`[Caldithrix]`,Phyla.LM$BRC1,xlab="Tenericutes",ylab="NKB19")
# abline(z)
# cor(Phyla.LM$BRC1,Phyla.LM$`[Caldithrix]`)
# text(0.0009,0.0006,"r = 0.18")


# look at correlations
library(reshape2)
library(Hmisc)

corr.matrix <- rcorr(as.matrix(t(network)),type="pearson")[[1]]
network.all <- make_netassoc_network(network,numnulls=1000)
assoc.matrix <- network.all$matrix_spsp_ses_all[,]

z <- lm(melt(assoc.matrix)$value~melt(corr.matrix)$value)
plot(melt(corr.matrix)$value,melt(assoc.matrix)$value,col="gray",xlab="Correlation",ylab="Co-occurence")
abline(z)
text(0.5,150,"y = -0.46 + 9.76x")
text(0.5,135,"Adj. R2 = 0.19")
text(0.5,120,"P < 0.00")


# heat maps
ggplot(data=melt(logged.heat[,]),aes(x=X1,y=X2)) + geom_tile(aes(fill=value)) +
  scale_fill_gradient2(na.value="white") + 
    theme_bw() + 
    theme(legend.title=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1),
        panel.border=element_rect(color="black",size=0.75))

ggplot(data=melt(primary.heat[,]),aes(x=X1,y=X2)) + geom_tile(aes(fill=value)) +
  scale_fill_gradient2(na.value="white") + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1),
        panel.border=element_rect(color="black",size=0.75))

ggplot(data=melt(palm.heat[,]),aes(x=X1,y=X2)) + geom_tile(aes(fill=value)) +
  scale_fill_gradient2(na.value="white") + 
  theme_bw() + 
  theme(legend.title=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1),
        panel.border=element_rect(color="black",size=0.75))
