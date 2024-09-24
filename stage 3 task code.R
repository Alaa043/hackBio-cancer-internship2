# Loading necessary libraries
BiocManager::install("TCGAbiolinks")
library(TCGAbiolinks)
library(DESeq2)
library(dplyr)
library(biomaRt)

# Query and download TCGA BRCA data
query_TCGA <- GDCquery(project='TCGA-BRCA',
                       data.category = 'Transcriptome Profiling',
                       data.type = 'Gene Expression Quantification',
                       experimental.strategy = 'RNA-Seq',
                       workflow.type = 'STAR - Counts',
                       access = 'open')
GDCdownload(query_TCGA)
expression_data <- GDCprepare(query_TCGA)

# Convert the data into a data frame
data <- assay(expression_data)
data <- as.data.frame(data)

# Check for missing values
sum(is.na(data))

# Handling missing values: Replace NA values with 0
data[is.na(data)] <- 0

# Preprocessing: Normalize the gene expression data using DESeq2's variance stabilizing transformation
dds <- DESeqDataSetFromMatrix(countData = data, 
                              colData = colData(expression_data), 
                              design = ~1)  # No design as we are just transforming

# Apply variance stabilizing transformation (vst) or rlog transformation
vsd <- vst(dds, blind=TRUE)
normalized_data <- assay(vsd)  # Extract normalized data

# Reducing the dataset to 20 primary and 20 recurrent tumor samples
# Assuming that 'sample_type' column in colData defines sample classification, adjust based on your data
# Get sample information (metadata)
sample_info <- as.data.frame(colData(expression_data))  # Convert DataFrame to data.frame

# primary and recurrent samples
primary_samples <- sample_info %>% filter(sample_type == "Primary Tumor") %>% head(20)
recurrent_samples <- sample_info %>% filter(sample_type == "Recurrent Tumor") %>% head(20)

# Combine the selected primary and recurrent samples
selected_samples <- rbind(primary_samples, recurrent_samples)

# Filter expression data based on selected samples
reduced_data <- normalized_data[, rownames(selected_samples)]

# Output dimensions of the reduced data
dim(reduced_data)
