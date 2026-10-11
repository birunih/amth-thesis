library(redist)
library(igraph)

data(iowa)
iowa_map <- redist_map(iowa, existing_plan = cd_2010, total_pop = pop)
adj <- get_adj(iowa_map)

#Redist adjacency is a list for which each entry is an integer list of nodes which are adjacent to a given node It is also 0-indexed. The igraph package requires an adjacency matrix as input (there are other options, but list->matrix is most convenient).
list_to_matrix <- function(adj_list) {
  mat <- matrix(0, nrow = length(adj_list), ncol = length(adj_list))
  
  for (i in 1:length(adj_list)) {
    for(j in 1:length(adj_list[[i]])) {
      mat[i,adj_list[[i]][j]+1] <- 1
    }
  }
  return(mat)
}

adj_graph <- graph_from_adjacency_matrix(list_to_matrix(adj), mode = 'undirected', weighted = NULL) #creates a graph object from the adjacency matrix

random_edges <- sample_spanning_tree(adj_graph) #samples a spanning tree at random over the adjacency graph and returns its edges

(random_tree <- subgraph_from_edges(adj_graph, random_edges)) #creates subgraph of spanning tree with these edges

V(random_tree)$name <- as.character(V(random_tree)) #was running into indexing errors, consultation with ChatGPT revealed I should do this to lock in node names

#Demonstrating what it would look like to cut 2 edges of the spanning tree.
cut_edges <- E(random_tree)[c(1,61)]

cut_ends <- ends(random_tree, cut_edges)

(subgraphs <- decompose(delete_edges(random_tree, cut_edges)))

#Identifying which of the 3 resulting trees is the middle one (containing 1 node from each edge cut).
middle <- subgraphs[[which(sapply(subgraphs, function(g) {
  sum(as.character(cut_ends[1,]) %in% V(g)$name) == 1 &&
  sum(as.character(cut_ends[2,]) %in% V(g)$name) == 1
}))]]

sum(iowa$pop[-as.numeric(V(middle)$name)]) #combined population of the two outer trees, i.e. the district created

#I want to calculate what this newly defined district's population would be for every possible pair of edges. This function takes in a spanning tree and list of node populations and outputs a matrix which shows resulting new-district populations.
new_district_pop <- function(spanning_tree, pop) {
  e <- length(E(spanning_tree))
  mat <- matrix(0, nrow = e, ncol = e)
  
  for(i in 1:e) {
    for(j in 1:e) {
      if(j > i) {
        cut_edges <- E(spanning_tree)[c(i,j)]
        
        cut_ends <- ends(spanning_tree, cut_edges)
        
        subgraphs <- decompose(delete_edges(spanning_tree, cut_edges))
        
        middle <- subgraphs[[which(sapply(subgraphs, function(g) {
          sum(as.character(cut_ends[1,]) %in% V(g)$name) == 1 &&
            sum(as.character(cut_ends[2,]) %in% V(g)$name) == 1
        }))]]
        
        mat[i,j] <- sum(pop[-as.numeric(V(middle)$name)])
      }
    }
  }
  
  return(mat)
}

#Here I apply the function using the spanning tree we generated earlier.
new_district_pop_mat <- new_district_pop(random_tree,iowa$pop)

#Say we want 30 districts total
n <- 30
ideal_pop <- sum(iowa$pop)/n #ideal district population, assuming population balance is a goal

dev_from_ideal_mat <- abs(new_district_pop_mat - ideal_pop)

k <- 3
sorted_indices <- order(dev_from_ideal_mat)[1:k]

k_best_edges <- arrayInd(sorted_indices, dim(dev_from_ideal_mat))

(edges_to_cut <- k_best_edges[sample(1:k,1),]) #picking an edge pair at random from the k best pairs ('best' meaning resulting in a new district closest to the intended population) 

#To make a redist_plans object, we need a vector of precinct assignments. The following code will create a map showing what happens when we cut edges 1 & 61 of whichever random spanning tree was generated earlier.
precinct_assignments <- rep(1, 99) 

for (i in 1:length(V(middle)$name)) {
  precinct_assignments[as.integer(V(middle)$name)] <- 2
}

redist_plans(precinct_assignments, iowa_map, algorithm = 'custom')

redist.plot.plans(redist_plans(precinct_assignments, iowa_map, algorithm = 'custom'), draws = 1, shp = iowa_map)

#This function allows you to visualize multiple plans side-by-side.
noncontig_plans <- function(map, n, k, runs) {
  adj <- get_adj(map)
  
  adj_graph <- graph_from_adjacency_matrix(list_to_matrix(adj), mode = 'undirected', weighted = NULL) 
  
  precinct_assignments <- matrix(1, nrow = nrow(map), ncol = runs)
  
  ideal_pop <- sum(map$pop)/n 
  
  for (i in 1:runs) {
  
    random_edges <- sample_spanning_tree(adj_graph) 
    
    random_tree <- subgraph_from_edges(adj_graph, random_edges)
    
    V(random_tree)$name <- as.character(V(random_tree))
    
    dev_from_ideal_mat <- abs(new_district_pop(random_tree, map$pop) - ideal_pop)
  
    sorted_indices <- order(dev_from_ideal_mat)[1:k]
    
    k_best_edges <- arrayInd(sorted_indices, dim(dev_from_ideal_mat))
    
    edges_to_cut <- k_best_edges[sample(1:k,1),]
    
    cut_edges <- E(random_tree)[edges_to_cut]
    
    cut_ends <- ends(random_tree, cut_edges)
    
    (subgraphs <- decompose(delete_edges(random_tree, cut_edges)))
    
    middle <- subgraphs[[which(sapply(subgraphs, function(g) {
      sum(as.character(cut_ends[1,]) %in% V(g)$name) == 1 &&
      sum(as.character(cut_ends[2,]) %in% V(g)$name) == 1
    }))]]
    
    
    for (j in 1:length(V(middle)$name)) {
        precinct_assignments[as.integer(V(middle)$name),i] <- 2
    }
  }
  
  redist_plans(precinct_assignments, map, algorithm = 'custom')
  
  redist.plot.plans(redist_plans(precinct_assignments, map, algorithm = 'custom'), draws = 1:runs, shp = map)
}

noncontig_plans(iowa_map, 10, 5, 6)
noncontig_plans(iowa_map, 3, 5, 100)
