# Load libraries ----------------------------------------------------------
library(ape)
library(dplyr)
library(getopt)
library(cluster)
suppressMessages(require(getopt))
suppressMessages(require(ape))
suppressMessages(require(cluster))
suppressMessages(require(dplyr))
# Setup -------------------------------------------------------------------

rm(list = ls())
options(stringsAsFactors = F)

##Set up options------------------------
spec0 <- matrix(c(
  'help'     , 'h', 0, "logical",   "show this helpful message",
  'phy'    , 'p', 1, "character", "tree",
  'nodes'  , 'n', 2, "boolean", "nodes",
  'taxonomy' , 't', 1, "character", "path to a taxonomy csv with tree tip labels in the first column and taxonomic levels as other columns",
  'taxlevel' , 'l', 1, "character", "the taxonomic level for which to calculate an index",
  'threshold', 's', 1, "integer"  , "the desired starting threshold percent of backbone distmat",
  'ngroups'  , 'k', 1, 'integer'  , "the number of groups, typically equals to the backbone sequence count and the desired sequence count",
  'nsequences','ns' ,1, 'integer'  , "the desired number of sequence",
  'alignment', 'a', 1, "character","fasta",
  'model'    , 'm', 2, "character","model",
  'distmax'  , 'd', 2, "integer","distance matrix",
  'cores'    , 'c', 2, "integer","cores"
  ), byrow = T, ncol = 5)


#Read the options and set the defuul操作设置--------------------------
opt0 <- getopt(spec0)
if ( is.null(opt0) | !is.null(opt0$help) ){
  message(getopt(spec0, usage = T))
  q(status = 1)
}

if (!is.null(opt0$threshold)) {# 检查 threshold 参数是否存在
  threshold_percent <- as.integer(opt0$threshold)# 将 threshold 参数的值赋值给变量 k
} else {  # 如果 threshold 参数未提供，默认赋值为 1
  threshold_percent <- 80
}
if (!is.null(opt0$ngroups)) {# 检查 ngroups 参数是否存在
  k <- as.integer(opt0$ngroups)# 将 ngroups 参数的值赋值给变量 k
} else {  # 如果 ngroups 参数未提供，默认赋值为 2000
  k <- 300
}
if (!is.null(opt0$nsequences)) {# 检查 nsequences 参数是否存在
  n <- as.integer(opt0$nsequences)# 将 nsequences 参数的值赋值给变量 k
} else {  # 如果 threshold 参数未提供，默认赋值为 1500
  n <- 200
}
if (is.null(opt0$taxonomy)) {
  stop("Error: path to taxonomy file is required")
}
Tax_data <- read.csv(opt0$taxonomy)

###Extract the tip names of backbone tree------------------------------------------
# Set up options ----------------------------------------------------------
# col1: long flag name
# col2: short flag name
# col3: 0 = no argument, 1 = required, 2 = optional
# col3: logical, integer, double, complex, character
# col5: optional, brief description

#spec_backbone <- matrix(c('dummy'  , '', 0, "logical" ), byrow = T, ncol = 4)


# Read options and do help -----------------------------------------------

#opt_backbone <- getopt(spec_backbone)

#if ( !is.null(opt_backbone$help) ){
 # cat(getopt(spec_backbone, usage = T))
 # q(status = 1)
#}
if (is.null(opt0$phy)) {
  stop("Error: path to constraint tree file is required")
}
# Set defaults -----------------------------------------------------------

if ( is.null( opt0$nodes   ) ) { opt0$nodes  = F   }

#testing the file
#phy_backbone <- read.tree('9_backbone_v3_multifurcatingconsensus_2022-02-16_TJCauto_GBMID_RN (2).tre')

# Load in backbone data -----------------------------------------------------------
phy_backbone <- read.tree(opt0$phy)

if( opt0$nodes){
  out_backbone <- phy_backbone$node.label
} else {
  out_backbone <- phy_backbone$tip.label
}

#write.csv(out_backbone, file = "output_backbone.csv", row.names = FALSE)

###VARIANCE
##Load the sequence data-----------------------------------------------
#spec <- matrix(c(
 # 'help'     , 'h', 0, "logical",

#), byrow = TRUE, ncol = 4)

# Read options
#opt <- getopt(spec)

# Set defaults ------------------------------------------------------------

if (is.null(opt0$alignment)) {
  stop("Error: path to alignment is required")
}

if (is.null(opt0$model)) opt0$model <- "F84"
if (is.null(opt0$distmax)) opt0$distmax <- 65536
if (is.null(opt0$cores)) opt0$cores <- 1

if (opt0$distmax < 2 | opt0$distmax > 65536) {
  stop("Error: -d/--distmax must be greater than 1 and less than 65,537")
}

# Load data ---------------------------------------------------------------
#alignment <- read.FASTA('95_nt_supermatrix.fasta')
alignment <- read.FASTA(opt0$alignment)
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



# Compute pairwise distances---------------------------------
distmat <- dist.dna(alignment, model = opt0$model, pairwise.deletion = T)
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
#k <- 200*1.5
# Make sure k is not larger than the data
if (k>= nrow(dist_matrix)){
  k <- max(k, nrow(dist_matrix)*0.8 )
}
print(k)
print(n)
# Set the initial percentage 初始阈值百分比
#threshold_percent <- 80

# Set the final threshold and k 初始化最终的阈值百分比和k值
final_threshold_percent <- 0
final_k <- 0

#write.csv(dist_matrix, file = "rowname.csv", row.names = FALSE)

##Clustering-------------------------------------------
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
    stop_loop <- TRUE 
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
    
    
    if (current_k >= k) {
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

##Delete the cluster where the backbone sequences are in & output-----------------------
# Prepare data for CSV output
rowname<-rownames(dist_matrix)
output_data <- data.frame(Sequence = rowname, Cluster = kmeans_result$cluster)

# 凡是在out_backbone中存在的名称，都使它和在output中相同组别的名称都被删掉
#output_filtered <- output_data[!(output_data$Sequence %in% subset_backbone) | !(output_data$Cluster %in% output_data$Cluster[output_data$Sequence %in% subset_backbone]), ]

#Get the Sequence column in out_backbone 获取out_backbone中的Sequence列
out_backbone_df <-as.data.frame(subset_backbone)
sequences_to_remove <- out_backbone_df$subset_backbone

#Find the Cluster where the matching Sequence resides in output_data 找到output_data中匹配的Sequence所在的Cluster
clusters_to_remove <- output_data %>%
  filter(Sequence %in% sequences_to_remove) %>%
  distinct(Cluster)

#Delete the Cluster 从output_data中删除clusters_to_remove中的Cluster
output_filtered <- output_data %>%
  anti_join(clusters_to_remove, by = "Cluster")



# 输出过滤后的DataFrame
print(output_filtered)

# Write data to CSV file
write.csv(output_data, file = "output.csv", row.names = FALSE)
write.csv(output_filtered, file = "output_filtered.csv", row.names = FALSE)

#write.csv(distmat_threshold_matrix, file = "thresholded.csv", row.names = FALSE)

#write.csv(dist_matrix, file = "original dist.csv", row.names = FALSE)

#write.csv(zero_dist_indices, file = "zero_dist_indices.csv", row.names = FALSE)
#write.csv(zero_dist_clusters, file = "zero_dist_clusters.csv", row.names = FALSE)
###Taxonomy information--------------------------------------------------------
##Prepare the files ----------------------------------------
#Import taxonomy file and select the information 导入文件，选择分类信息
#Tax_data <- read.csv('SITE-100 Database - mitogenomes.csv')
Tax_selected <- select(Tax_data, db_id,species,subgenus, genus, subtribe,tribe, 
                       subfamily, family, superfamily, infraorder, suborder, order)
#Replace the space with NA 将空格替换成NA
Tax_selected <- Tax_selected %>%
  mutate(across(everything(), ~ case_when(. %in% c("", " ") ~ NA_character_, TRUE ~ .)))
#Summarize the number of valid cells in each column统计每列数目
tax_summary <- Tax_selected %>%
  summarise(across(.fns = ~sum(!is.na(.)), .cols = everything()))
#Summarize the taxa of each column统计分类数目
tax_counts <- Tax_selected %>%
  summarise(across(everything(), ~ n_distinct(.)))
print(tax_counts)
#Specify the taxonomic level you want ???要指定使用哪一列哦

#Change the "Sequence" column name to "db_id"  使用rename()函数将"Sequence"列名改为"db_id"
output_filtered <- output_filtered %>%
  rename(db_id = Sequence)
#Connect the two data frames according to the ID 使用inner_join()函数将两个数据框根据标本ID进行内连接
joined_data <- inner_join(Tax_selected, output_filtered, by = 'db_id')

##Compare the clustering result and taxonomy information------------------
#Compare the joined_data 对joined_data进行比对判断
comparison_result <- joined_data %>%
  group_by(subfamily) %>%
  summarise('if taxa in one cluster' = ifelse(length(unique(Cluster)) == 1, "Y", "N")) %>%
  left_join(joined_data, by = 'subfamily',multiple = "all") 

comparison_result<-select(comparison_result,'db_id', 'Cluster', 'subfamily','if taxa in one cluster')
#check the result 输出比对结果
print(comparison_result)

#Compare the joined_data reversely 反方向比对对joined_data进行分组和比对判断
comparison_result_reverse <- joined_data %>%
  group_by(Cluster) %>%
  summarise('if cluster in one taxa' = ifelse(n_distinct(subfamily) == 1, "Y", "N")) %>%
  left_join(joined_data, by = 'Cluster',multiple = "all") 
comparison_result_reverse<-select(comparison_result_reverse,'db_id','if cluster in one taxa')
#Match the two data frames into a new one 使用left_join()函数将两个数据框匹配到一个新的数据框中
merged_result <- left_join(comparison_result, comparison_result_reverse, by = "db_id")

##Group the clusters step by step--------------------------
#taxa perfectly match with clusters
#Extract the sequences where taxa perfectly match with clusters 使用filter()函数提取互相唯一的行
taxa_perfect <- merged_result %>%
  filter(`if taxa in one cluster` == "Y", `if cluster in one taxa` == "Y")

#ALL NAs
#Extract rows where a single cluster is all "NA" 使用group_by()函数和filter()函数提取满足条件的行
NAs_in_one_cluster <- merged_result %>%
  group_by(Cluster) %>%
  filter(all(is.na(subfamily)))
#Exclude rows in NAs_in_one_cluster and taxa_perfect 使用anti_join()函数排除NAs_in_one_cluster和taxa_perfect中的行，得到rest数据框
rest <- anti_join(merged_result, NAs_in_one_cluster, by = "db_id") %>%
  anti_join(taxa_perfect, by = c("db_id", "Cluster", "subfamily"))

#Taxa and NAs are in the same cluster  taxa与na在同一cluster里面的
#Extract rows where a single cluster contains only one unique taxa and "NA"  使用group_by()函数按照Cluster分组，然后进行条件筛选
taxa_nas_unique <- rest %>%
  group_by(Cluster) %>%
  filter(any(`if taxa in one cluster` == 'Y')& n_distinct(subfamily, na.rm = TRUE) == 1)
#Exclude them 删除这些
rest<-anti_join(rest,taxa_nas_unique,by='db_id')

#a single cluster contains multiple non-dispersed taxa
#Extract taxa where a single cluster contains multiple non-dispersed taxa taxa不分开，但同cluster里有多个taxa
ntaxa_inonecluster <- rest %>%
  filter(`if taxa in one cluster` == 'Y')
#Exclude the clusters instead of mere taxa 删除这些分组
rest<-anti_join(rest,ntaxa_inonecluster,by="Cluster")

#Dispersed taxa which contain a pure cluster (with NAs)
#Extract clusters where dispersed taxa contain a pure cluster (NAs can exist)
NNbutpurewithNAs <- rest %>%
  group_by(Cluster) %>%
  filter(any(n_distinct(subfamily, na.rm = TRUE) == 1))
#Exclude the corresponding taxa instead of mere clusters
rest<-anti_join(rest,NNbutpurewithNAs,by="subfamily")


##Select sequences from each group in turn -------------------
random_rows1 <- taxa_perfect %>%
  group_by(Cluster) %>%
  sample_n(size = 1) 
#length(unique(taxa_perfect$Cluster))
random_rows2 <- NAs_in_one_cluster %>%
  group_by(Cluster) %>%
  sample_n(size = 1) 
#length(unique(NAs_in_one_cluster$Cluster))
random_rows3 <- taxa_nas_unique %>%
  group_by(subfamily) %>%
  sample_n(size = 1) 
#length(unique(taxa_nas_unique$subfamily))
random_rows4 <- ntaxa_inonecluster %>%
  group_by(subfamily) %>%
  sample_n(size = 1) 
#length(unique(ntaxa_inonecluster$subfamily))
random_rows5 <- NNbutpurewithNAs %>%
  group_by(subfamily) %>%
  sample_n(size = 1) 
#length(unique(NNbutpurewithNAs$subfamily))
random_rows6 <- rest %>%
  group_by(subfamily) %>%
  sample_n(size = 1) 
#length(unique(rest$subfamily))

#Merge the id together
totalcandidates <- c(random_rows1$db_id, random_rows2$db_id, random_rows3$db_id, random_rows4$db_id, random_rows5$db_id,random_rows6$db_id)
if (length(totalcandidates) > n){
  #Specifies the number of elements to delete 指定要删除的最后几个元素个数
  num_to_remove <- length(totalcandidates) - n
  #Delete the last few elements of the list 使用切片操作来删除列表的最后几个元素
  totalcandidates <- totalcandidates[1:(length(totalcandidates) - num_to_remove)]
  
}
print(length(totalcandidates))
#Final output
write.csv(totalcandidates, file = "finallist.csv", row.names = FALSE)
