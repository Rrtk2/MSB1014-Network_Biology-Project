# skills 2 - week2 - networkinference


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
# Insert settings
if(!exists("Sourcemode")){
options(stringsAsFactors 	= F)
Verbose 					= 3		#0-3:	0= no feedback  1= prints results  2=prints results + feedback  3= prints all
ErrorReporter 				= 0 # if error, report where and stop script
}
if(Verbose>=3) cat(Warning.style("Done running block 02a\n"))


#-----------------------------------------------------------------------------------------------------#
#		Block 02b		INPUT Settings
#-----------------------------------------------------------------------------------------------------#
if(!exists("Sourcemode")){
DataFolderLocation 		= "/Dropbox/MU Systems Biology/Courses/MSB1014 Network biology/Skills/Week 2 - Network inference/" 	# Must be in format "/location...DatalocationFolder/"
}
if(Verbose>=3) cat(Warning.style("Done running block 02b\n"))


#-----------------------------------------------------------------------------------------------------#
#		Block 03		Prepare folders and load data
#-----------------------------------------------------------------------------------------------------#
if(ErrorReporter==0){ # start reporter
tryCatch({ #Start functional reporter block 
			
	# Data_dir should contain the data file!
	if(Sys.info()[6]=="RR"){
		data_dir <- paste("D:",DataFolderLocation,sep="")  # Rick's link (laptop)
	}else if(Sys.info()[6]=="Admin"){
		data_dir <- paste("E:",DataFolderLocation,sep="")  # Rick's other link (desktop)
	}

	# Get folder and move to one level above, then make new folder containing timestamp
	setwd(data_dir)

	  
	if(Verbose>=3) cat(Warning.style("Done running block 03\n"))

}, #End functional reporter block
error=function(cond) {
assign("ErrorReporter","Error in block 03",envir = .GlobalEnv)
assign("cond",cond,envir = .GlobalEnv)
})} # end reporter


#-----------------------------------------------------------------------------------------------------#
#							Load data
#-----------------------------------------------------------------------------------------------------#
data <- read.delim("insilico-data.txt", row.names=1)
set.seed(123) # for reproducibility of results


#-----------------------------------------------------------------------------------------------------#
#							Bayesian networks
#-----------------------------------------------------------------------------------------------------#


output.folder <- data_dir 
LinearNet(output.folder, data)
analyse.output(output.folder)


prob.file <- paste(output.folder, "/NetworkProbability_Matrix.txt", sep = "")
prob.mat <- read.table(prob.file)
prob.mat <- as.matrix(prob.mat)
threshold = 0.08
prob.mat[prob.mat < threshold] <- 0

graph.Bayes <- igraph::graph_from_adjacency_matrix(prob.mat,weighted=TRUE)
plot(graph.Bayes)
createNetworkFromIgraph(graph.Bayes, title="GRENITS",collection="Bayesian")
copyVisualStyle("default","GRENITS")
setVisualStyle("GRENITS")
setEdgeTargetArrowShapeDefault("Arrow", style.name = "GRENITS")


#-----------------------------------------------------------------------------------------------------#
#							Correlation
#-----------------------------------------------------------------------------------------------------#


res.cor <- correlate(t(data), diagonal = 0)
row.names(res.cor) <- res.cor$rowname
res.cor[1] <- NULL

res.cor.filtered <- res.cor
res.cor.filtered[res.cor.filtered < 0.3] <- 0

graph.Cor <- igraph::graph_from_adjacency_matrix(as.matrix(res.cor.filtered), weighted=TRUE)
createNetworkFromIgraph(graph.Cor,title="Correlation_Network",collection="Correlation")


#-----------------------------------------------------------------------------------------------------#
#							Regression
#-----------------------------------------------------------------------------------------------------#
 

weightMat <- GENIE3(as.matrix(data))

linkList <- getLinkList(weightMat)

linkList <- getLinkList(weightMat, reportMax=5)
linkList <- getLinkList(weightMat, threshold=0.1)

graph.Reg <- igraph::graph_from_data_frame(linkList, directed=TRUE, vertices=rownames(data))
createNetworkFromIgraph(graph.Reg, title="GENIE3_Network",collection="GENIE3")
copyVisualStyle("default","GENIE3")
setVisualStyle("GENIE3")
setEdgeTargetArrowShapeDefault("Arrow", style.name = "GENIE3")


#-----------------------------------------------------------------------------------------------------#
#							Mutual information
#-----------------------------------------------------------------------------------------------------#

mi <- knnmi.all(data)
grn1 <- aracne.a(mi)
grn2 <- aracne.m(mi)
grn3 <- clr(mi)
grn4 <- mrnet(mi)

graph.MU1 <- igraph::graph_from_adjacency_matrix(grn1, weighted=TRUE)
createNetworkFromIgraph(graph.MU1, title="ARACNE.A_Network",collection="MI")
setVisualStyle("default")

graph.MU2 <- igraph::graph_from_adjacency_matrix(grn2, weighted=TRUE)
createNetworkFromIgraph(graph.MU2, title="ARACNE.M_Network",collection="MI")
setVisualStyle("default")

graph.MU3 <- igraph::graph_from_adjacency_matrix(grn3, weighted=TRUE)
createNetworkFromIgraph(graph.MU3, title="CLR_Network",collection="MI")
setVisualStyle("default")

graph.MU4 <- igraph::graph_from_adjacency_matrix(grn4, weighted=TRUE)
createNetworkFromIgraph(graph.MU4, title="MRNET_Network",collection="MI")
setVisualStyle("default")




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

graph = g

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