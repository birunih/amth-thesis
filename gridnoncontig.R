library(foreign)
library(igraph)
library(redist)
library(sf)

#Importing 3 functions from iowanoncontig.R
list_to_matrix <- function(adj_list) {
  mat <- matrix(0, nrow = length(adj_list), ncol = length(adj_list))
  
  for (i in 1:length(adj_list)) {
    for(j in 1:length(adj_list[[i]])) {
      mat[i,adj_list[[i]][j]+1] <- 1
    }
  }
  return(mat)
}

new_district_pop <- function(spanning_tree, pop) {
  e <- length(E(spanning_tree))
  mat <- matrix(0, nrow = e, ncol = e)
  
  for(i in 1:e) {
    for(j in 1:e) {
      if(j > i) {
        cut_edges <- E(spanning_tree)[c(i,j)]
        
        cut_ends <- ends(spanning_tree, cut_edges)
        
        subgraphs <- decompose(delete_edges(spanning_tree, cut_edges), mode = 'weak')
        
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
    
    subgraphs <- decompose(delete_edges(random_tree, cut_edges), mode = 'weak')
    
    middle <- subgraphs[[which(sapply(subgraphs, function(g) {
      sum(as.character(cut_ends[1,]) %in% V(g)$name) == 1 &&
        sum(as.character(cut_ends[2,]) %in% V(g)$name) == 1
    }))]]
    
    
    for (j in 1:length(V(middle)$name)) {
      precinct_assignments[as.integer(V(middle)$name),i] <- 2
    }
  }
  
  redist_plans(precinct_assignments, map, algorithm = 'custom')
  
  print(redist.plot.plans(redist_plans(precinct_assignments, map, algorithm = 'custom'), draws = 1:runs, shp = map))
  
  return(list(
    adj_graph = adj_graph,
    precinct_assignments = precinct_assignments
  )) #new, explained later
}

#Inquired with https://r-spatial.github.io/sf/#cheatsheet & ChatGPT as to how to set up a 10x10 grid geometry
grid <- st_make_grid(st_bbox(c(xmin = 0, ymin = 0, xmax = 10, ymax = 10), crs = NA), cellsize = c(1,1))

#Making dataframe of a 100-precinct region (10x10 grid), each precinct contains the same population (of 1)
biruniville <- data.frame(
  precinctid = 1:100,
  pop = rep(1, 100),
  geometry = grid
)

bville_map <- redist_map(biruniville, ndists = 10)

noncontig_plans(bville_map, 10, 5, 3)
fiftyruns <- noncontig_plans(bville_map, 10, 5, 50)
hundredruns <- noncontig_plans(bville_map, 10, 5, 100)

#Wanted to make a function to check if the newly formed district is contiguous, so I added a line to noncontig_plans that returns adj_graph & precinct_assignments.
check_contiguous <- function(assignment_matrix, adj){
  
  check_contiguous_run <- function(assignment_vector, adj) {
    subg <- subgraph(adj, vids = which(assignment_vector == 1))
    is_connected(subg)
  }
    
  contig <- rep(0, ncol(assignment_matrix))
  
  for(i in 1:ncol(assignment_matrix)) {
    contig[i] <- check_contiguous_run(assignment_matrix[,i], adj)
  }
  
  return(contig)
}

#Want to get some sense of how often contiguity occurs using this adapted algorithm. I assume it will differ with map, amount of precincts, & intended district size. When we want 1/10 of the 10x10 grid, it's about 29% contiguity.
sum(check_contiguous(fiftyruns$precinct_assignment, fiftyruns$adj_graph))
sum(check_contiguous(hundredruns$precinct_assignment, fiftyruns$adj_graph))

#I expect that creating a district with half the total population will result in more contiguity.
fiftyrunshalf <- noncontig_plans(bville_map, 2, 5, 50)
sum(check_contiguous(fiftyrunshalf$precinct_assignment, fiftyruns$adj_graph)) #I am getting 64%

#Want to count contiguity in the Iowa examples too.
data(iowa)
iowa_map <- redist_map(iowa, existing_plan = cd_2010, total_pop = pop)

fiftyrunshalfiowa <- noncontig_plans(iowa_map, 2, 5, 50)
sum(check_contiguous(fiftyrunshalfiowa$precinct_assignment, fiftyrunshalfiowa$adj_graph)) #I am getting 36%

montgomery <- read.dbf('MD_Montgomery_points.dbf')

#pseudocoding this out - a function which applies the modified algorithm to generate and plot a single complete redistricting plan (i.e. 1 run)
#I want a top-level function which will repeat a step function (sorry, I know I shouldn't say step function since that's already a thing) that makes one more district out of the remaining unclaimed area n-1 times, where n is the desired number of districts.
#
#step function
#inputs: adjacency graph of remaining area, desired population, k (the larger k is, the more variable the sizes of the districts)
#
#top-level function
#inputs: map, desired # districts, k

noncontig_district <- function(adj_graph, pop, ideal_pop, k) {
  
  current_nodes <- as.integer(V(adj_graph)$name)
  
  current_pop <- pop[current_nodes]
  
  random_edges <- sample_spanning_tree(adj_graph) 
  
  random_tree <- subgraph_from_edges(adj_graph, random_edges)
  
  dev_from_ideal_mat <- abs(new_district_pop(random_tree, pop) - ideal_pop)
  
  sorted_indices <- order(dev_from_ideal_mat)[1:k]
  
  k_best_edges <- arrayInd(sorted_indices, dim(dev_from_ideal_mat))
  
  edges_to_cut <- k_best_edges[sample(1:k,1),]
  
  cut_edges <- E(random_tree)[edges_to_cut]
  
  cut_ends <- ends(random_tree, cut_edges)
  
  subgraphs <- decompose(delete_edges(random_tree, cut_edges), mode = 'weak')
  
  middle <- subgraphs[[which(sapply(subgraphs, function(g) {
    sum(as.character(cut_ends[1,]) %in% V(g)$name) == 1 &&
      sum(as.character(cut_ends[2,]) %in% V(g)$name) == 1
  }))]]
  
  remaindernodes <- V(middle)$name
  
  remainder <- induced_subgraph(
    adj_graph,
    vids = which(V(adj_graph)$name %in% remaindernodes))
  
  assigned <- setdiff(as.integer(V(adj_graph)$name), remaindernodes)
  
  return(list(
    remainder = remainder,
    assigned = assigned))
}

noncontig_plans1 <- function(map, n, k) {
  
  adj <- get_adj(map)
  
  adj_graph <- graph_from_adjacency_matrix(list_to_matrix(adj), mode = 'undirected', weighted = NULL)
  
  V(adj_graph)$name <- as.character(V(adj_graph))
  
  precinct_assignments <- rep(n, nrow(map))
  
  ideal_pop <- sum(map$pop)/n
  
  for (i in 1:(n-1)) {
    
    step <- noncontig_district(adj_graph, map$pop, ideal_pop, k)
    
    assigned <- step$assigned
    
    for (j in 1:length(assigned)) {
      precinct_assignments[assigned] <- i
    }
    
    adj_graph <- step$remainder
  }
  
  redist_plans(precinct_assignments, map, algorithm = 'custom')
  
  print(redist.plot.plans(redist_plans(precinct_assignments, map, algorithm = 'custom'), draws = 1, shp = map))
  
  return(list(
    adj_graph = adj_graph,
    precinct_assignments = precinct_assignments))
}