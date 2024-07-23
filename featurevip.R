# Set the working directory to where the module data is stored
setwd("D:/lab/AD/IBD/module/20240121")

# Define names for different result sets presumably from various analyses or models
# Merge and process the results from a specific model or analysis

###GWAS
name = "GWAS"
resmerge(name)

name = "ADIBD-KMEALL-DEG"
resmerge(name)

name = "ADIBD-KMEALL-DEG-MT"
resmerge(name)

name = "ADIBD-KMEALLDEG-MT+ALL"
resmerge(name)

name = "ADIBD-KMEALL-DEG"
resmerge(name)  


# Define a function to merge and process the results from random forest, rpart, and glm models
resmerge <- function(name){
  # For each dataset (e.g., GWAS, ADIBD-KMEALL-DEG, etc.), read the importance files for three different studies (GSE140829, GSE33000, GSE63063)
  # Merge the results based on the feature index, handling NA values, and calculate the sum of importance scores
  # Perform similar operations for rpart and glm models
  
  # For random forest importance
  rf14 = read.csv(...)  # Read the importance file for the first study and model
  # ... (other similar read.csv calls for the remaining studies and models)
  merge_res <- merge(merge(rf14, rf33, by = "X",all = TRUE), rf63, by = "X",all = TRUE)  # Merge the datasets
  merge_res [is.na(merge_res )] <- 0  # Replace NA values with 0
  # Calculate the sum of MeanDecreaseAccuracy and MeanDecreaseGini importance scores
  # Sort the results based on the summed importance scores in descending order
  
  # Repeat the process for rpart and glm models, calculating the sum of their respective importance scores
  # Sort the results based on the summed importance scores in descending order
  
  # Merge the sorted results from all models, calculate the overall sum of MeanDecreaseAccuracy scores
  # Sort the merged results based on the overall summed importance scores in descending order
  
  # Write the final sorted dataframe to a CSV file without row names
  write.csv(sorted_df4,file = paste(name,"-importanceall.csv"),row.names = F)
}