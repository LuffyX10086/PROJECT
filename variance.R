install.packages("ape")
library(ape)
library(dplyr)
# Load libraries ----------------------------------------------------------
library(getopt)
library(cluster)
suppressMessages(require(getopt))
suppressMessages(require(ape))

##extract the tip names of backbone tree------------------------------------------
# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

# Load libraries ----------------------------------------------------------


# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

spec_backbone <- matrix(c(
  'help'   , 'h', 0, "logical",
  'phy'    , 'p', 1, "character",
  'nodes'  , 'n', 2, "boolean"
), byrow = T, ncol = 4)

# Read options and do help -----------------------------------------------

opt_backbone <- getopt(spec_backbone)

if ( !is.null(opt_backbone$help) ){
  cat(getopt(spec_backbone, usage = T))
  q(status = 1)
}

# Set defaults -----------------------------------------------------------

if ( is.null( opt_backbone$nodes   ) ) { opt_backbone$nodes  = F   }

#testing the file
phy_backbone <- read.tree('9_backbone_v3_multifurcatingconsensus_2022-02-16_TJCauto_GBMID_RN (2).tre')

# Load in backbone data -----------------------------------------------------------
#phy_backbone <- read.tree(opt_backbone$phy)

if( opt_backbone$nodes){
  out_backbone <- phy_backbone$node.label
} else {
  out_backbone <- phy_backbone$tip.label
}

#write.csv(out_backbone, file = "output_backbone.csv", row.names = FALSE)


##Load the sequence data
spec <- matrix(c(
  'help'     , 'h', 0, "logical",
  'alignment', 'a', 1, "character",
  'model'    , 'm', 2, "character",
  'distmax'  , 'd', 2, "integer",
  'cores'    , 'c', 2, "integer",
  'ngroups','n',1,'integer'
), byrow = TRUE, ncol = 4)

# Read options
opt <- getopt(spec)

# Set defaults ------------------------------------------------------------

if (is.null(opt$alignment)) {
  stop("Error: path to alignment is required")
}

if (is.null(opt$model)) opt$model <- "F84"
if (is.null(opt$distmax)) opt$distmax <- 65536
if (is.null(opt$cores)) opt$cores <- 1

if (opt$distmax < 2 | opt$distmax > 65536) {
  stop("Error: -d/--distmax must be greater than 1 and less than 65,537")
}

# Load data ---------------------------------------------------------------
alignment <- read.FASTA('95_nt_supermatrix.fasta')
#alignment <- read.FASTA(opt$alignment)
sequence_ids <- names(alignment)
#write.csv(sequence_ids, file = "sequence_ids.csv", row.names = FALSE)

## The name in backbone tree is different from the database. To scynchronize them.
subset_backbone <- c()
# Traverse every value in sequence ids 遍历Values2中的每个对象
for (value in sequence_ids) {
  # Check if the value is in the names in backbone  检查当前对象是否存在于Values1中
  if (any(grepl(value, out_backbone))) {
    # extract the value
    subset_backbone <- c(subset_backbone, value)
  }
}
#write.csv(subset_backbone, file = "subset_backbone.csv", row.names = FALSE)



# Compute pairwise distances
distmat <- dist.dna(alignment, model = opt$model, pairwise.deletion = T)
distmat[is.nan(distmat)] <- ceiling(max(distmat, na.rm = T))

# Subset the backbone matrix 
# Transform it to matrix 转换dist对象为矩阵
dist_matrix <- as.matrix(distmat)

# Get the subset indices 获取子集序列的索引
subset_indices <- match(subset_backbone, rownames(dist_matrix))

# Extract the subset distance matrix 根据子集索引提取子集距离矩阵
subset_distmat <- dist_matrix[subset_indices, subset_indices]

# Transform the matrix to dist 将子集距离矩阵重新转换为dist对象
subset_dist <- as.dist(subset_distmat)

# Output in csv打印子集距离矩阵
#write.csv(as.matrix(subset_dist), file = "subset_dist.csv", row.names = FALSE)



## Calculating the threshold

# k for testing 期望的k值 ？如何决定大多少比较好？
k <- 200*1.5
# Make sure k is not larger than the data
if k >= nrow(dist_matrix){
  k <- max(optk*1.5, nrow(dist_matrix)*0.8 )
}

# Set the initial percentage 初始阈值百分比
threshold_percent <- 80

# Set the final threshold and k 初始化最终的阈值百分比和k值
final_threshold_percent <- 0
final_k <- 0

#write.csv(dist_matrix, file = "rowname.csv", row.names = FALSE)

stop_loop <- FALSE  #Initializes the variable that stops the loop 初始化停止循环的变量

while (threshold_percent >= 0 && !stop_loop) {
  threshold <- quantile(subset_dist, probs = threshold_percent / 100) #Calculate the threshold value
  
  distmat_thresholded <- distmat
  distmat_thresholded[distmat < threshold] <- 0 #Sets elements that are less than the threshold to 0
  distmat_threshold_matrix <- as.matrix(distmat_thresholded)
  
  kmeans_result <- tryCatch({ #Catch possible error situations, output error if errors occur and upadate the threshold_percent
    kmeans(distmat_thresholded, centers = k)
  }, error = function(e) {
    cat("Clustering failed with threshold percent:", threshold_percent, "%\n")
    threshold_percent <- threshold_percent - 1
    if (threshold_percent < 0) {
      cat("Threshold Percent reached 0%. Clustering failed.\n")
      stop_loop <- TRUE  # If the threshold_percent is less than 0, a message indicating that the threshold is 0% is displayed. Set stop_loop to TRUE to stop the loop.
    }
    return(NULL)
  })
  
  if (!is.null(kmeans_result)) {#If no errors occur
    cluster_labels <- kmeans_result$cluster #Gets the label vector of the clustering result
    zero_dist_indices <- which(distmat_threshold_matrix == 0) #The corresponding labels are extracted according to the index zero_dist_indices with a value of 0 in the distance matrix
    zero_dist_names <- names(cluster_labels[zero_dist_indices])
    merged_clusters <- cluster_labels
    merged_clusters[zero_dist_names] <- max(cluster_labels) + 1 #Update the label value in zero_dist_clusters to the maximum label value plus 1 to avoid conflicts with existing clustering labels
    non_na_indices_names <- zero_dist_names[!is.na(zero_dist_names)] #Get an index that does not contain NA values
    kmeans_result$cluster[non_na_indices_names] <- merged_clusters[non_na_indices_names]  #Reassign the clustering label 重新分配 cluster_labels
    
    current_k <- length(unique(kmeans_result$cluster)) # Calculate the number of different tags in the kmeans_result$cluster
    
    
    if (current_k >= 200) {
      cat("Threshold Percent:", threshold_percent, "%\n")
      cat("k:", current_k, "\n") #Output the values
      stop_loop <- TRUE  # Set the stopping condition to TRUE设置停止循环的条件为 TRUE
    } else {
      threshold_percent <- threshold_percent - 1
      kmeans_result <- NULL  # Reinitialize the kmeans result重新初始化 kmeans_result
    }
  }
}

if (!stop_loop) {
  cat("Threshold Percent reached 0%. Clustering failed.\n")
}

#Load the output to a dataframe
data <- data.frame(data_name= names(kmeans_result$cluster),cluster = kmeans_result$cluster)
# Groups according to the original discontinuous group label 按照原始的非连续分组标签进行分组
grouped_data <- data %>% group_by(cluster)

# Add a new contiguous group label and remove the original non-contiguous group label添加一个新的连续的组别标签，并删除原始的非连续分组标签
new_clusters <- grouped_data %>% 
  mutate(new_cluster = cur_group_id()) %>% 
  ungroup() %>% 
  select(-cluster)


# Prepare data for CSV output
output_data <- data.frame(Sequence = rowname, Cluster = kmeans_result$cluster)

# Write data to CSV file
write.csv(output_data, file = "output4.csv", row.names = FALSE)
#write.csv(distmat_threshold_matrix, file = "thresholded.csv", row.names = FALSE)

#write.csv(dist_matrix, file = "original dist.csv", row.names = FALSE)

#write.csv(zero_dist_indices, file = "zero_dist_indices.csv", row.names = FALSE)
#write.csv(zero_dist_clusters, file = "zero_dist_clusters.csv", row.names = FALSE)



