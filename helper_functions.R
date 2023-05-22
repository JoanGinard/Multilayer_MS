## Helper functions

## Function to filter matrices without corresponding patient or conversely patients without corresponding patient
## INPUT  orig_matrices_dir -> matrices folder name 
## INPUT csv_append -> append to the patient id in matrices file name EX: "_r_matrix.csv" --> "patientID_r_matrix.csv"
## INPUT clinical_data -> dataframe with, at least, an id variable with patients id

filter_patients <- function(orig_matrices_dir, csv_append, clinical_data, verbose = T){
  list_files <- list_of_files(dataDIR, orig_matrices_dir) #file list from folder --> function defined in this script
  # cReate dataframe with patients IDs and filenames associated with id
  df <- data.frame(clinical_data$id)
  df$name <- paste0(df$clinical_data.id, csv_append)
 
  if(verbose){
    cat(paste("Cheking data...\n"))
    cat(paste("Matrices without patients...\n"))  
    discarded_matrices <- setdiff(list_files, df$name) #Files that have no patients
    print(length(discarded_matrices))
    cat(paste("Number of patients without matrices...\n"))
    discarded_patients <- setdiff(df$name, list_files) #patients that have no file
    print(length(discarded_patients))
    discarded_data <- c(discarded_matrices, discarded_patients)
  }else{
    discarded_matrices <- setdiff(list_files, df$name)
    discarded_patients <- setdiff(df$name, list_files)
    discarded_data <- c(discarded_matrices, discarded_patients)
  }
  
  list_files <- list_files[!(list_files %in% discarded_data)] # final file list
  return(list_files)
}

####################################################################################################################



###### function to obtain a list of csv files from DIR
#INPUT --> dataDIR and matrices folder name
#Return file list of csv files
list_of_files <- function(dataDIR, matricesDIR){
  file_list <-list.files(file.path(dataDIR, matricesDIR), pattern = "\\.csv$")
  
   return(file_list)
}

########################################################################################################################

###### function to load data from DIR of matrices
#INPUT --> dataDIR and matrices folder name
# RETURN a list of matrices

load_data <- function(dataDIR, matricesDIR){
  # Obtain list of files if from original folders we filters, else directly from folder
  if(matricesDIR == "subjects_FA/"){
    file_list <- filter_patients(FADIR, FA_append, clinical, verbose = F)
  }else if(matricesDIR == "subjects_fMRI/"){
    file_list <- filter_patients(fMRIDIR, fMRI_append, clinical, verbose = F)
  }else if(matricesDIR == "subjects_GM/"){
    file_list <- filter_patients(GMDIR, GM_append, clinical, verbose = F)
  }else {
    file_list <- list_of_files(dataDIR, matricesDIR) #obtain file list
    
  }

  tables <- vector("list", length = length(file_list)) #allocate space in a vector list
  
  #Load matrices
  for(i in seq_along(tables)){
    tables[[i]] <- read.csv(file.path(dataDIR, matricesDIR, file_list[i]), header = FALSE)
  }
  
  return(tables)
  
}


##########################################################################################################################

#### function to load data, but specifically as matrix
#Same function as previous but returns a matrix 

load_as_matrix <- function(dataDIR, matricesDIR){
  dfs <- load_data(dataDIR, matricesDIR)
  matrices_list <- lapply(dfs, as.matrix)

  return(matrices_list)
}

###########################################################################################################################

###### function to save matrices after step of preprocessing
#old_append --> final part of the original name of data. EX: 001FA_matrices_orig.csv --> _orig.csv
#new_append --> new final part of the name. EX: _corrected.csv --> 001FA_matrices_corrected.csv
#newMatricesDIR --> folder to save new matrices
#file_list --> list of files inside original folder of matrices
#tables --> list of matrices to store

save_data <- function(dataDIR, newMatricesDIR, file_list, tables, new_append, old_append = ".csv"){
  #create new folder if not exists
  if(!dir.exists(file.path(dataDIR, newMatricesDIR))){
    dir.create(file.path(dataDIR, newMatricesDIR)) 
  }
  
  #Iterate over matrix list
  for(k in seq_along(file_list)){
    new_name <- strsplit(file_list[k], old_append)[[1]] #strip file names old_append
    new_name <- paste0(new_name, new_append) #create new name appending new_append
    
    write.table(tables[[k]], file.path(dataDIR, newMatricesDIR, new_name),
                sep = ",", quote = FALSE, col.names = FALSE, row.names = FALSE)
    
    
  }
  
}

#############################################################################################################

###### function to create a dataframe with HS ids and the name of the loaded files
## csv_append --> -> final part of the original name of data. EX: 001FA_matrices_orig.csv --> _orig.csv
HS_dataframe <- function(dataDIR, csv_append){

  #Load clinical data
  clinical_data <- read.csv(file.path(dataDIR, "clinic.csv"))
  
  #Create dataframe with column Id and file name
  HS_id <- clinical_data$id[clinical_data$controls_ms == 0]
  HS_df <- data.frame(HS_id)
  HS_df$file_name <- paste0(HS_df$HS_id, csv_append)
  
  return(HS_df)
}

#############################################################################################################

###### function to obtain a dataframe with the number of subjects where a connection is present
## Returns a matrix like dataframe where if (i,j) = 13, 13 HS subject have this connection present
#INPUT --> dataDIR and matrices folder name
## csv_append --> -> final part of the original name of data. EX: 001FA_matrices_orig.csv --> _orig.csv
HS_60connections <- function(dataDIR, matricesDIR, csv_append){
  HS_df <- HS_dataframe(dataDIR, csv_append) #function defined in this script
  
  HS_tables <- vector("list", length = length(HS_df$HS_id)) # we create vector-list to store matrices of Healthy subjects
  
  #Load matrices,
  for(i in seq_along(HS_df$file_name)){
    HS_tables[[i]] <- read.csv(file.path(dataDIR, matricesDIR, HS_df$file_name[i]), header = FALSE)
  }
  
  
  # We obtain a list of dataframes, each dataframe with 1 and 0; 1 where we have values greater than 0, 0 otherwise.
  HS_tables_filtered <- lapply(HS_tables, function(x) {data.frame(x) %>%
      mutate(across(1:ncol(data.frame(x)), ~ifelse(.x > 0, 1, 0)))})
  
  rows <- dim(HS_tables[[1]])[1]
  columns <- rows #matrices are simmetric
  
  #We sum all dataframes
  HS_all <- data.frame(matrix(nrow = rows, ncol = columns))  #Create empty dataframe
  HS_all[,][is.na(HS_all[,])] <- 0 #fill NAs in empty dataframe
  
  
  for(i in seq_along(HS_tables_filtered)){
    HS_all <- HS_all +HS_tables_filtered[[i]] #sum to all 0,1 matrix like dataframes
  }
  #Now HS_all is a dataframe indicating in each position how many healthy subjects have a weighted connection of more than 0
  return(HS_all)
    
}

########################################################################################################################

###### function to check if files in folder and id in clinical data are the same. 
## OUTPUTS a message
## INPUT --> list of files and 
## csv_append --> -> final part of the original name of data. EX: 001FA_matrices_orig.csv --> _orig.csv

same_order <- function(list_files, csv_append, clinical_data){
  df <- data.frame(clinical_data$id, list_files)
  df$name <- paste0(df$clinical_data.id, csv_append)  #Create variable to compare
  
  value <- sum(df$list_files == df$name) #How many there are in the same order
  
  if(value == length(list_files)){
    print("CORRECT. SAME ORDER")
  }else{
    print("PROBLEM. REVIEW PROCESS")
  }
  

}

########################################################################################################

###### function to correct matrices by sex and age (linear correction)
#### INPUT: tables, a list where matrices are stored, 
###         new_tables, another list where resulting matrices will be stored
###         HS_indices: indices that correspond to HS
###         age: vector of subjects ages
###         sex: vector with subject gender
### OUTPUT: list of tables (matrices) with corrected matrices
age_sex_correction <- function(tables, new_tables, HS_indices, sex, age){
  
  rows <- dim(tables[[1]])[1]
  columns <- rows #squared matrix
  
  #We iterate only over upper diagonal and then apply same numbers to lower diagonal
  # skip main diagonal as it should be zero
  for(i in 1:(rows-1)){
    for(j in (i+1):columns){
      #We create a vector with values in the same position in all matrices, 
      #i.e: A vector containing all [i,j] values from all matrices
  
      y<- sapply(tables, function(x)  x[i,j]) 
      

      fit <- lm(y ~ sex + age) #mind age and sex are in correct order
      res <- summary(fit)$residuals #get residuals of our fit
      
      #Mean of the healthy subjects
      control_means <- mean(y[HS_indices])
      res <- res + control_means #add healthy subjects mean to our residuals
      for(k in seq_along(tables)){
        new_tables[[k]][i,j] <- ifelse(res[k] > 0, res[k], 0) #after correction it is possible some case will have a negative value
        new_tables[[k]][j,i] <- ifelse(res[k] > 0, res[k], 0) #It must be symmetric
      }
    }
    
  }
  
  return(new_tables)
  
  
}

#######################################################################################################################################

###### function to construct graphs from matrices
## It constructs UNDIRECTED and WEIGHTED graphs
## Returns a list of graphs, same order as patients id in clinical data_frame
## IMPORTANT --> "vertex_names": names of the vertex must be loaded in the script before executing the function
## REQUIRES ---> igraph

graph_constructor <- function(matrices) {
  matrices <- lapply(matrices, as.matrix) # Defined function
  #From matrix to graph
  networks <- lapply(matrices, graph_from_adjacency_matrix, 
                      mode = "undirected", 
                      weighted = TRUE) 
  
  
  # Set the names
  
  networks <- lapply(networks, function(graph){
    V(graph)$name <- vertex_names
    graph
  })
  return(networks)
}


#####################################################################################################################################

###### function to load graphs from matrices
#INPUT matrices dir, clinical data
# Normalization flag, default TRUE, to indicate if svd normalization is required
# loads data from folder, checks if order is the same as clinical_data
#converts matrices to graphs
# RETURNS a list of graphs
## IMPORTANT --> vertex names must be loaded, also dataDIR (folder of data)

load_graphs <- function(matricesDIR, clinical_data,  normalization = TRUE) {
  list_files <- list_of_files(dataDIR, matricesDIR) #Files in data
  
  if(matricesDIR == FADIR){
    csv_append <- "_FA_factor_corrected.csv"
  }else if(matricesDIR == GMDIR){
    csv_append <- "_GM_matrix_corrected.csv"
  }else if(matricesDIR == fMRIDIR){
    csv_append <- "_r_matrix_corrected.csv"
  }
  
  same_order(list_files, csv_append, clinical_data) #Check files are in the same order as in clinical data
  
  
  
  ### Load matrices
  
  matrices <- vector(mode = "list", length = length(list_files)) #Vector to store data
  
  matrices <- load_data(dataDIR, matricesDIR) #We use a helper function
  
  if(normalization){
    matrices <- lapply(matrices, svd_normalization) #Apply normalization if required
  }
  
  
  graph_list <- graph_constructor(matrices)
  
  return(graph_list)
}

##########################################################################################################################

######## function to obtain range and histogram of values
#INPUT list of matrices dataframes
# returns a histogram graph while prints out range of values of matrices

range_values <- function(table_list){
  all_values <- unlist(table_list) 
  range_val <- range(all_values) #Obtain range
  
  print(paste("Minimun values is", range_val[1], "and maximum is", range_val[2]))
  
  #Histogram
  network_name <- deparse(substitute(table_list))
  network_name <- strsplit(network_name, "_")[[1]][1]
  
  df_all_values <- data.frame(values = all_values)
  
  p <- ggplot(df_all_values, aes(x = values)) + 
    geom_histogram(bins = 30, color = "black", fill = "navy") +
    labs(x = "Values", y = "Count", title = paste("Histogram of", network_name,"connections"))
  
  return(p)
}

##########################################################################################################################

### Function to apply SVD normalization
## INPUT matrix to normlaize
## OUTPUTS normalized matrices

svd_normalization <- function(matrix){
  
  data.svd <- svd(matrix)
  result<- data.svd$u %*% (10*diag(data.svd$d)/data.svd$d[1]) %*% t(data.svd$v) #Normalization according to bibliography
  
  #Set to zero negative values and extremely small values, 
  #most of them are in fact zero (~1e-16)
  threshold <- 1e-6
  
  result <- ifelse(result > threshold, result, 0)
  return(result)
  
}

#########################################################################################################################

### Function to make statistical tests and compare two variables
### If normal variables performs t.test if not Shapiro-Wilk
### INPUT df -> dataframe with clinical data and variables, including the one to test
### INPUT name -> variable name to check statistical differences
### Returns p value
###


statistical_test <- function(df, name){

  PwMS <- df[df$controls_ms == "PwMS", ]
  HS <- df[df$controls_ms == "HS", ]
  val <- df[, c(name)]
  
  #Check for normality
  
  #In some case all values are the same, like all 1s, in those case shapiro.test will not work
  # we set p value to zero, it does not have a normal distribution
  unique_values_HS <- length(unique(HS[, c(name)]))
  unique_values_PwMS <- length(unique(PwMS[, c(name)]))
  
  if(unique_values_HS == 1 | unique_values_PwMS == 1){
    HS_p <- 0
    PwMS_p <- 0
  }else{
    HS_test <- shapiro.test(HS[, c(name)])
    PwMS_test <- shapiro.test(PwMS[, c(name)])
    HS_p <- HS_test$p
    PwMS_p <- PwMS_test$p
  }
  

  if(PwMS_p < 0.05 | HS_p < 0.05){
    #we cannot assume normal distribution, for name
    #Perform Whitney U test
    test <- wilcox.test(PwMS[, c(name)], HS[, c(name)], PAIRED = FALSE)
      return(test$p.value)

  }else{
    #Can assume normality for name
    #Check for homoscedasticity (we could use a different test)
    
    bart <-bartlett.test( val ~ controls_ms, data = df)
    var <- FALSE
    if(bart$p.value > 0.05){
    var <- TRUE  
    }
    #Perform t test
    test <- t.test(PwMS[, c(name)], HS[, c(name)], PAIRED = FALSE, var.equal = var)
      return(test$p.value)

    
  } 
  
}




