## Script to define function to calculate some graph measures

############################################
########              ########  
#######   SINGLELAYER  #########
########              ##########
###########################################



#### GLOBAL MEASURES #########################

# GLOBAL EFFICIENCY
calculate_global_efficiency <- function(graph) {
  global_efficiency(graph, directed = FALSE, weights = E(graph)$weight)
}



# MEAN PATH LENGTH
calculate_mean_path_length <- function(graph) {
  mean_distance(graph, directed = FALSE, weights = E(graph)$weight)
}


# TRANSITIVITY, GLOBAL CLUSTERING

calculate_global_transitivity <- function(graph) {
  transitivity(graph, type = "undirected", weights = E(graph)$weight)
}

# DIAMETER

calculate_diameter <- function(graph) {
  diameter(graph, directed = FALSE, weights = E(graph)$weight)
}

# MODULARITY

calculate_modularity <- function(graph) {
  communities <- cluster_louvain(graph) #fisrt we need to find the communities
  modularity(communities)
}


# ASSORTATIVITY for the moment being we do not use it as it is not defined for weighted networks

calculate_assortativity <- function(graph){
  assortativity(graph, directed = FALSE)
}

# GLOBAL STRENGTH is calculated as the mean of individual strentgh
# we do not really need a function as we will apply mean to local strength

######## LOCAL MEASURES ########################################

# LOCAL EFFICIENCY
calculate_local_efficiency <- function(graph){
  local_efficiency(graph, weights = E(graph)$weight, directed = FALSE)
}

# LOCAL TRANSITIVITY, LOCAL CLUSTERING (we will not use it as it does not show statistically significant results)
calculate_local_transitivity <- function(graph) {
  transitivity(graph, type = "localundirected", weights = E(graph)$weight)
}

# Closeness centrality of vertices

calculate_closeness_centrality <- function(graph) {
  closeness(graph, weights = E(graph)$weight)
}

# STRENGTH
#We do no need function as we only need to lappy
#Example
# Graph_strength <- lapply(graph_list, strength)



##### APPLY FUNCTION TO GRAPH LIST ###################3

apply_function_to_graph_list <- function(graph_list, func) {
  result_list <- lapply(graph_list, func)
  return(result_list)
}



#### add global measurements to clinical dataframe
## If we add or remove measures we have to modify the function
## df --> dataframe with clinical data, where we will add results of measures
## graphs --> list of graphs
add_global_graph_measures <- function(df, graphs, matrix_type = c("FA", "GM", "fMRI")){
  
  result <- apply_function_to_graph_list(graphs, calculate_diameter)
  df[paste0(matrix_type,"_diameter")] <- unlist(result)
  
  result <- apply_function_to_graph_list(graphs, calculate_global_efficiency)
  df[paste0(matrix_type,"_gl-efficiency")] <- unlist(result)
  
  result <- apply_function_to_graph_list(graphs, calculate_mean_path_length)
  df[paste0(matrix_type,"_mean_path-length")] <- unlist(result)
  
  result <- apply_function_to_graph_list(graphs, calculate_global_transitivity)
  df[paste0(matrix_type,"_gl-transitivity")] <- unlist(result)
  
  result <- apply_function_to_graph_list(graphs, calculate_modularity)
  df[paste0(matrix_type,"_gl-modularity")] <- unlist(result)
  
  result <- lapply(graphs, strength)
  df[paste0(matrix_type,"_gl-strength")] <- unlist(lapply(result, mean))
  
  df[paste0(matrix_type,"_density")] <- unlist(lapply(graphs, edge_density))
  
  
  return(df)
  
}


### add local measurements to dataframe

add_local_graph_measures <- function(df, graphs, matrix_type = c("FA", "GM", "fMRI")) {
  result <- apply_function_to_graph_list(graphs, calculate_local_efficiency)
  df_temp <- do.call(rbind, result)
  colnames(df_temp) <- paste0(matrix_type,"_loc-eff_", colnames(df_temp))
  df <- cbind(df, df_temp)
  
  
  result <- apply_function_to_graph_list(graphs, calculate_closeness_centrality)
  df_temp <- do.call(rbind, result)
  colnames(df_temp) <- paste0(matrix_type,"_loc-clos-cent_", colnames(df_temp))
  df <- cbind(df, df_temp)  
  
  result <- lapply(graphs, strength)
  df_temp <- do.call(rbind, result)
  colnames(df_temp) <- paste0(matrix_type,"_loc-strength_", colnames(df_temp))
  df <- cbind(df, df_temp)
  
  return(df)
}


############################################
########              ########  
#######   MULTILAYER #########
########             ##########
###########################################
  
#### add global measurements to clinical dataframe
## If we add or remove measures we have to modify the function

############ DISCARDED ###############
### Global measures does not produce significant results #####
add_global_multi_graph_measures <- function(df, multilayer_graphs){
  
  result <- lapply(multilayer_graphs, GetAverageGlobalClustering, Layers, number_nodes)
  df["MultiGlobalClustering"] <- unlist(result)
  
  result <- lapply(multilayer_graphs, GetMultiStrength, Layers, number_nodes, FALSE)
  df["MultiGlobalStrength"] <- unlist(lapply(result, mean))
  
  result <- lapply(multilayer_graphs, GetMultiClosenessCentrality, Layers, number_nodes)
  result <- lapply(result, function(x) x$closeness)
  df["MultiClosenessCentrality"] <- unlist(lapply(result, mean))
  
  
  return(df)
  
}
#######################################



add_multi_graph_measures <- function(df, multilayer_graphs) {
  result <- lapply(multilayer_graphs, GetMultiPageRankCentrality, Layers, Nodes)
  df_temp <- do.call(rbind, lapply(result, function(matrix){
    t(matrix)
  }))
  colnames(df_temp) <- paste0("MultiPRCent_", vertex_names)
  df <- cbind(df, df_temp)
  
  result <- lapply(multilayer_graphs, GetMultiStrengthSum, Layers, Nodes, FALSE)
  df_temp <- do.call(rbind, lapply(result, function(matrix){
    t(matrix)
  }))
  colnames(df_temp) <- paste0("MultiSTSum_", vertex_names)
  df <- cbind(df, df_temp)  
  
  result <-lapply(multilayer_graphs, GetMultiPathStatistics, Layers, Nodes)
  
  result <- sapply(result, function(data){data$closeness})
  result <- t(result)
  df_temp <- as.data.frame(result)
  colnames(df_temp) <- paste0("MultiPath_", vertex_names)
  df <- cbind(df, df_temp)
  
  return(df)
}


