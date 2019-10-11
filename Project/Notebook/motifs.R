

 
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
	
	


count_motifs(g, 3) # tot number of motifs
sample_motifs(g, 3) #??????


########### test below discard

g = barabasi.game(10, 
    m = 5,
    power = 0.6, 
    out.pref = TRUE,
    zero.appeal = 0.5,
    directed = TRUE)

# Label nodes to more easily keep track during subsets/deletions
V(g)$name = c('one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine', 'ten')

subGraph = graph.neighborhood(g, order = 1, V(g)[1], mode = 'all')[[1]]
allMotifs = triad.census(subGraph)
removeNode = delete.vertices(subGraph, 'one')
node1Motifs = allMotifs - triad.census(removeNode)

# see graph.isoclass