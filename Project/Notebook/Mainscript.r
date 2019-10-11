#-----------------------------------------------------------------------------------------------------#
#		Block 00		GENERAL INFORMATION
#-----------------------------------------------------------------------------------------------------#

# Copyright statement comment:
#   All rights reserved.
#
# Author comment:
#   Rick Reijnders
#   Script version: 11-10-2019

# File description:
#	Name
#	  Mainscript.R
#
  
rm(list=ls()) # clean workspace


#-----------------------------------------------------------------------------------------------------#
#		Block 01		(Install &) Load packages
#-----------------------------------------------------------------------------------------------------#
# Install packages if needed and load, or just load packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", ask = F)

# Insert all packages in requiredpackages
requiredpackages <-
  c("GRENITS",
    "igraph",
	"corrr",
	"GENIE3",
	"parmigene",
    "RCy3")
	
for (i in requiredpackages) {
	if (!requireNamespace(i, quietly = TRUE))
		BiocManager::install(i, ask = F, dependencies = TRUE)
	require(as.character(i), character.only = TRUE)
	print(i)
}


#-----------------------------------------------------------------------------------------------------#
#		Block 02a		GENERAL Settings
#-----------------------------------------------------------------------------------------------------#
options(stringsAsFactors 	= F)
motifsize 					= 4

#-----------------------------------------------------------------------------------------------------#
#		Block 02b		INPUT Settings
#-----------------------------------------------------------------------------------------------------#
DataFolderLocation 		= "/Github/MSB1014-Network_Biology-Project/Project/Notebook/" 	# Must be in format "/location...DatalocationFolder/"


#-----------------------------------------------------------------------------------------------------#
#		Block 03		Prepare folders and load data
#-----------------------------------------------------------------------------------------------------#
# Data_dir should contain the data file!
if(Sys.info()[6]=="RR"){
	data_dir <- paste("D:",DataFolderLocation,sep="")  # Rick's link (laptop)
}else if(Sys.info()[6]=="Admin"){
	data_dir <- paste("E:",DataFolderLocation,sep="")  # Rick's other link (desktop)
}

# Get folder and move to one level above, then make new folder containing timestamp
setwd(data_dir)


#-----------------------------------------------------------------------------------------------------#
#							Load data
#-----------------------------------------------------------------------------------------------------#
data <- read.delim("insilico-data.txt", row.names=1)
set.seed(1) # for reproducibility of results

#-----------------------------------------------------------------------------------------------------#
#							True network definition and analalysis
#-----------------------------------------------------------------------------------------------------#
trueGraph = graph_from_literal(ErbB3.P--+Ras.GTP,ErbB3.P--+Gab1.P,ErbB3.P--+ErbB1.P,ErbB3.P--+PIP2,ErbB3.P--+SHC.P,ErbB3.P--+ErbB2.P,ErbB3.P--+PIP3,ERK.P.P--+ErbB3.P,ERK.P.P--+Gab1.P,ERK.P.P--+SOS.P,ERK.P.P--+ErbB1.P,ERK.P.P--+ErbB2.P,ERK.P.P--+ERK.P,Ras.GTP--+ErbB3.P,Ras.GTP--+ErbB1.P,Ras.GTP--+ErbB4.P,Ras.GTP--+SHC.P,Ras.GTP--+RAF.P,Ras.GTP--+ErbB2.P,MEK.P.P--+ERK.P.P,MEK.P.P--+MEK.P,MEK.P.P--+ERK.P,Gab1.P--+ErbB3.P,Gab1.P--+Ras.GTP,Gab1.P--+ErbB1.P,Gab1.P--+ErbB4.P,Gab1.P--+PIP2,Gab1.P--+ErbB2.P,Gab1.P--+PIP3,ErbB1.P--+ErbB3.P,ErbB1.P--+Ras.GTP,ErbB1.P--+Gab1.P,ErbB1.P--+SOS.P,ErbB1.P--+ErbB4.P,ErbB1.P--+PIP2,ErbB1.P--+SHC.P,ErbB1.P--+ErbB2.P,ErbB1.P--+PIP3,AKT.P--+AKT.P.P,ErbB4.P--+Ras.GTP,ErbB4.P--+Gab1.P,ErbB4.P--+ErbB1.P,ErbB4.P--+PIP2,ErbB4.P--+SHC.P,ErbB4.P--+ErbB2.P,ErbB4.P--+PIP3,PIP2--+Gab1.P,PIP2--+ErbB4.P,PIP2--+ErbB2.P,PIP2--+PIP3,MEK.P--+MEK.P.P,SHC.P--+ErbB3.P,SHC.P--+Ras.GTP,SHC.P--+SOS.P,SHC.P--+ErbB1.P,SHC.P--+ErbB4.P,SHC.P--+ErbB2.P,RAF.P--+Ras.GTP,RAF.P--+MEK.P.P,RAF.P--+MEK.P,ErbB2.P--+ErbB3.P,ErbB2.P--+Ras.GTP,ErbB2.P--+Gab1.P,ErbB2.P--+ErbB1.P,ErbB2.P--+ErbB4.P,ErbB2.P--+PIP2,ErbB2.P--+SHC.P,ErbB2.P--+PIP3,ERK.P--+ERK.P.P,PIP3--+AKT.P,PIP3--+PIP2,PIP3--+AKT.P.P,AKT.P.P--+AKT.P,AKT.P.P--+RAF.P)

#nodes:
V(trueGraph)
length(V(trueGraph))

#links:
E(trueGraph)
length(E(trueGraph))


#-----------------------------------------------------------------------------------------------------#
#		Script block 3: FUNCTIONS
#-----------------------------------------------------------------------------------------------------#

# Get motifs fromgraph and plot the motifs with frequency
GetMotifsFromGraph = function(graph.name, motifsize = 4, directed=TRUE, amountRandomNets = 1000){

	graph = get(graph.name)
	if(exists("randomNetMotifsTotal")){rm(randomNetMotifsTotal)}
	
	# Find network motifs in the graph "graph":
	mygraphmotifs <- graph.motifs(graph, motifsize)
	# Find which motifs occur:

	for (i in 1:length(mygraphmotifs))
	{
			if(is.na(mygraphmotifs[i])){next}
			

			motif <- mygraphmotifs[i]
			if (motif > 0) # There are some occurrences of this motif
			{
			   # print(i)
			   # Find out what the motif looks like:
			   motifgraph <- graph.isocreate(size=motifsize, number=i-1, directed=directed)
			   edges <- E(motifgraph)
			   plot(motifgraph)
			   text(-0.5,-1.2,paste("This motif(",i-1,") occurs ",motif," times:",sep=""))
			   print(paste("Motif",i-1,"occurs",motif,"times."))
			   #cat(edges)
			   
			}
	}
	
	print(mygraphmotifs)
    mygraphmotifs[is.na(mygraphmotifs)]=0

	
	# start statistical testing
	graph = delete.vertices(graph, degree(graph)==0)
	
	meanInDegree = mean(degree(graph,mode = "in"))
	meanOutDegree = mean(degree(graph,mode = "out"))
	
	# start random net generation; create distribution
	for(p in 1:amountRandomNets){
		g <- barabasi.game(n=length(V(graph)),directed = directed)
		randomNetMotifs = graph.motifs(g, motifsize)
		randomNetMotifs[is.na(randomNetMotifs)]=0
		
		if(!exists("randomNetMotifsTotal")){
			randomNetMotifsTotal = randomNetMotifs
			}else{
			randomNetMotifsTotal = randomNetMotifsTotal + randomNetMotifs
		}
	}
	randomNetMotifsTotal = randomNetMotifsTotal/amountRandomNets

	# Which occured more than random:
	impMotifIndex = which(mygraphmotifs>(randomNetMotifsTotal))-1
	print(impMotifIndex)
}
#-----------------------------------------------------------------------------------------------------#
#							Bayesian networks
#-----------------------------------------------------------------------------------------------------#
# directed
# 20N 16L
# 
output.folder <- paste(data_dir,"Bayes",sep="") 
LinearNet(output.folder, data)
analyse.output(output.folder)


prob.file <- paste(output.folder, "/NetworkProbability_Matrix.txt", sep = "")
prob.mat <- read.table(prob.file)
prob.mat <- as.matrix(prob.mat)
threshold = 0.08
prob.mat[prob.mat < threshold] <- 0

graph.Bayes <- igraph::graph_from_adjacency_matrix(prob.mat,weighted=TRUE)
plot(graph.Bayes)
#nodes: 20
V(graph.Bayes)
length(V(graph.Bayes))

#links: 16
E(graph.Bayes)
length(E(graph.Bayes))
GetMotifsFromGraph("graph.Bayes")
#createNetworkFromIgraph(graph.Bayes, title="GRENITS",collection="Bayesian")
#copyVisualStyle("default","GRENITS")
#setVisualStyle("GRENITS")
#setEdgeTargetArrowShapeDefault("Arrow", style.name = "GRENITS")


#-----------------------------------------------------------------------------------------------------#
#							Correlation
#-----------------------------------------------------------------------------------------------------#
#undirected ->directed..?
# 20N 72L
res.cor <- correlate(t(data), diagonal = 0)
row.names(res.cor) <- res.cor$rowname
res.cor[1] <- NULL

res.cor.filtered <- res.cor
res.cor.filtered[res.cor.filtered < 0.3] <- 0

graph.Cor <- igraph::graph_from_adjacency_matrix(as.matrix(res.cor.filtered), weighted=TRUE)
plot(graph.Cor)
#nodes: 20
V(graph.Cor)
length(V(graph.Cor))

#links: 72
E(graph.Cor)
length(E(graph.Cor))
GetMotifsFromGraph("graph.Cor")
#createNetworkFromIgraph(graph.Cor,title="Correlation_Network",collection="Correlation")


#-----------------------------------------------------------------------------------------------------#
#							Regression
#-----------------------------------------------------------------------------------------------------#
# undirected ->directed..?
# 20N 40L
weightMat <- GENIE3(as.matrix(data))

linkList <- getLinkList(weightMat)

linkList <- getLinkList(weightMat, reportMax=5)
linkList <- getLinkList(weightMat, threshold=0.1)

graph.Reg <- igraph::graph_from_data_frame(linkList, directed=TRUE, vertices=rownames(data))
plot(graph.Reg)
#nodes: 20
V(graph.Reg)
length(V(graph.Reg))

#links: 40
E(graph.Reg)
length(E(graph.Reg))
GetMotifsFromGraph("graph.Reg")

#createNetworkFromIgraph(graph.Reg, title="GENIE3_Network",collection="GENIE3")
#copyVisualStyle("default","GENIE3")
#setVisualStyle("GENIE3")
#setEdgeTargetArrowShapeDefault("Arrow", style.name = "GENIE3")


#-----------------------------------------------------------------------------------------------------#
#							Mutual information
#-----------------------------------------------------------------------------------------------------#
mi <- knnmi.all(data)
grn1 <- aracne.a(mi)
grn2 <- aracne.m(mi)
grn3 <- clr(mi)
grn4 <- mrnet(mi)

if(F){ #LOCKED
	graph.MU1 <- igraph::graph_from_adjacency_matrix(grn1, weighted=TRUE)
	plot(graph.MU1)
	#nodes: 20
	V(graph.MU1)
	length(V(graph.MU1))

	#links: 152
	E(graph.MU1)
	length(E(graph.MU1))

	#createNetworkFromIgraph(graph.MU1, title="ARACNE.A_Network",collection="MI")
	#setVisualStyle("default")

	graph.MU3 <- igraph::graph_from_adjacency_matrix(grn3, weighted=TRUE)
	plot(graph.MU3)
	#nodes: 20
	V(graph.MU3)
	length(V(graph.MU3))

	#links: 232
	E(graph.MU3)
	length(E(graph.MU3))

	#createNetworkFromIgraph(graph.MU3, title="CLR_Network",collection="MI")
	#setVisualStyle("default")

	graph.MU4 <- igraph::graph_from_adjacency_matrix(grn4, weighted=TRUE)
	plot(graph.MU4)
	#nodes: 20
	V(graph.MU4)
	length(V(graph.MU4))

	#links: 234
	E(graph.MU4)
	length(E(graph.MU4))

	#createNetworkFromIgraph(graph.MU4, title="MRNET_Network",collection="MI")
	#setVisualStyle("default")
}

graph.MU2 <- igraph::graph_from_adjacency_matrix(grn2, weighted=TRUE)
plot(graph.MU2)
#nodes: 20
V(graph.MU2)
length(V(graph.MU2))

#links: 82
E(graph.MU2)
length(E(graph.MU2))
GetMotifsFromGraph("graph.MU2")
#createNetworkFromIgraph(graph.MU2, title="ARACNE.M_Network",collection="MI")
#setVisualStyle("default")


#-----------------------------------------------------------------------------------------------------#
#							motif
#-----------------------------------------------------------------------------------------------------#
 motifsize = 4
 g <- barabasi.game(n=15,directed = T)
 #erdos.renyi.game # randomn
plot(g)
motifs(g, motifsize)
which(motifs(g, motifsize)>=1)


# create the motif MINUS ONE BECAUSE THESE FACKERS START AT 0; THIS IS R NOT FUCKING C++
#plot(graph.isocreate(size=3, number=3-1,directed = T))

graph = graph.Cor

# Find network motifs in the graph "graph":
	mygraphmotifs <- graph.motifs(graph, motifsize)
	# Find which motifs occur:

	for (i in 1:length(mygraphmotifs))
	{
			if(is.na(mygraphmotifs[i])){next}
			

			motif <- mygraphmotifs[i]
			if (motif > 0) # There are some occurrences of this motif
			{
			   print(i)
			   # Find out what the motif looks like:
			   motifgraph <- graph.isocreate(size=motifsize, number=i-1, directed=TRUE)
			   edges <- E(motifgraph)
			   plot(motifgraph)
			   text(-1,-1,paste("This motif occurs",motif,"times."))
			   print(paste("This motif occurs",motif,"times:"))
			   print(edges)
			}
	}
	
	











































#-----------------------------------------------------------------------------------------------------#
#							Ensembl
#-----------------------------------------------------------------------------------------------------#
ensembl= data.frame(matrix(0,20,20))

#selectedlist = c("graph.Bayes","graph.Cor","graph.Reg")
selectedlist = ls()[grep(pattern = "graph",x = ls())]

SelectedlistLength = length(selectedlist)
#for( i  in 1:length(ls()[grep(pattern = "graph",x = ls())])){

for( i  in 1:SelectedlistLength){


a=get(selectedlist[i])[]
ensembl = ensembl + as.matrix(a)
}
ensembl= as.data.frame(ensembl)

rownames(ensembl) = rownames(a)
colnames(ensembl) = colnames(a)



ensembl_edit = ensembl

ensembl_edit = (ensembl_edit-min(ensembl_edit))/(max(ensembl_edit)-min(ensembl_edit))

#ensembl_edit2 = sign(ensembl_edit) * abs(ensembl_edit)^(1/3)

ensembl_edit[is.na(ensembl_edit)] = 0
hist(quantile(unlist(as.matrix(ensembl_edit[])),probs=(1:100)/100),breaks=100)


cutoff = 0.2

ensembl_edit[ensembl_edit<cutoff] = 0
ensembl_edit[ensembl_edit>cutoff] = 1
graphR.Ensembl = graph_from_adjacency_matrix(adjmatrix = as.matrix(ensembl_edit), mode ="directed",diag = T)
sum(graphR.Ensembl[])



plot(graphR.Ensembl)

createNetworkFromIgraph(graphR.Ensembl, title="graphR.Ensembl",collection="Ensembl")
copyVisualStyle("default","graphR.Ensembl")
setVisualStyle("graphR.Ensembl")
setEdgeTargetArrowShapeDefault("Arrow", style.name = "graphR.Ensembl")

# predicted 12 correct from 30