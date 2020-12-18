#### Install the libraries #####################################################
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("flowCore")
# BiocManager::install("ggplot2")
# BiocManager::install("ggpubr")
# BiocManager::install("pheatmap")
# BiocManager::install("tidyr")
# BiocManager::install("FlowRepositoryR")
# if (!requireNamespace("devtools", quietly = TRUE))
#   install.packages("devtools")
# devtools::install_github("saeyslab/FlowSOM", ref = "FlowSOM_v2")
# devtools::install_github("saeyslab/PeacoQC")


#### Download the data #########################################################
ds <- FlowRepositoryR::download(FlowRepositoryR::flowRep.get("FR-FCM-ZZQY"), 
                                "Data/Raw")
ds <- FlowRepositoryR::download(FlowRepositoryR::flowRep.get("FR-FCM-Z2TQ"), 
                                "Data/Raw")


#### Prepare data ##############################################################
#microbenchmark::microbenchmark({
# 1. Load the libraries
library(flowCore)
library(FlowSOM)
library(ggplot2)

# 2. Define the general and preprocessing variables
file_pattern <- "\\d.fcs" #digit at the end and fcs extension
reference_file <- read.FCS("Data/Raw/21-10-15_Tube_011.fcs")
reference_marker <- "PE-A" # Scatter values will be scaled to have the same range

markers_of_interest <- c("SSC-A", "MHCII", "CD49b", "CD11b", "CD64",
                         "FcERI", "CD161", "Ly-6G", "CD3", "CD19", "CD11c")

live_gate <- flowCore::polygonGate(filterId = "Live",
                                   .gate = matrix(data = c(60000, 100000, 150000, 
                                                           250000, 250000, 60000, 
                                                           60000, 1.6, 1.9, 2.5,
                                                           2.5, -0.3, -0.3, 1.6),
                                                  ncol = 2,
                                                  dimnames = list(c(), 
                                                                  c("FSC-A", 
                                                                    "APC-Cy7-A"))))

# 2. Define and create the directories
dir_prepr <- "Data/Preprocessed/" #where the preprocessed data will be stored
dir_QC <- "Data/Preprocessed/QC/" #where the data QC results will be stored
dir_RDS <- "RDS/" #where the R objects will be stored
dir_results <- "Results/" #where the results will be stored
dir_raw <- "Data/Raw/" #where the raw data is located
path_comp <- "Data/Raw/attachments/Compensation.csv" #where comp matrix is located

for (path in c(dir_prepr, dir_QC, dir_RDS, dir_results)){
  dir.create(path)
}

# 4. Prepare some additional information for preprocessing the files 
# given the variable choices of step 3.
files <- list.files(path = dir_raw,
                    pattern = file_pattern)
channels_of_interest <- GetChannels(object = reference_file,
                                    markers = markers_of_interest, 
                                    exact = FALSE)
compensation_matrix <- read.csv(path_comp, 
                                check.names = FALSE, row.names = 1)
colnames(compensation_matrix) <- sub(" :: .*", "",         
                                     colnames(compensation_matrix))

# Compute transformation list
ff_m <- PeacoQC::RemoveMargins(reference_file, channels_of_interest)
ff_c <- flowCore::compensate(ff_m, compensation_matrix)
translist <- estimateLogicle(ff_c, colnames(compensation_matrix))
ff_t <- flowCore::transform(ff_c, translist)
q5_goal <- quantile(exprs(ff_t)[,reference_marker], 0.05)
q95_goal <- quantile(exprs(ff_t)[,reference_marker], 0.95)
q5_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.05)
q95_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.95)
SSCA_a <- (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
SSCA_b <- q5_goal - q5_SSCA * (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
translist <- c(translist, 
               transformList("SSC-A", flowCore::linearTransform(a = SSCA_a,
                                                                b = SSCA_b)))
# 5. Execute the following for-loop, which will preprocess each fcs file.
for (file in files){
  # 6. Read the fcs file into a flowframe
  ff <- read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE)

  # 7. Remove margin events
  ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest)
  
  # 8. Compensate
  ff_c <- flowCore::compensate(ff_m, compensation_matrix)
  
  # 9. Transform, logicle for marker channels, linear for scatter channel
  ff_t <- flowCore::transform(ff_c, translist)
  
  # 10. Remove doublets and filter live cells
  ff_s <- PeacoQC::RemoveDoublets(ff_t, nmad = 6) #increase the accepted width
  
  selected_live <- filter(ff_s, live_gate)
  ff_l <- ff_s[selected_live@subSet, ]
  
  # 11. QC with PeacoQC
  PQC <- PeacoQC::PeacoQC(ff = ff_l,
                          channels = channels_of_interest,
                          plot = TRUE, save_fcs = FALSE,
                          output_directory = dir_QC)
  
  # 12. Save the preprocessed data
  write.FCS(PQC$FinalFF,
            file = paste0(dir_prepr, file))

  # 13. Visualize the preprocessing
  filter_plot <- function(ff_pre, ff_post, title, marker_x, marker_y){
    df <- data.frame(x = exprs(ff_pre)[,marker_x],
                     y = exprs(ff_pre)[,marker_y])
    i <- sample(nrow(df), 10000)
    if (!"Original_ID" %in% colnames(exprs(ff_pre))) {
      ff_pre@exprs <- cbind(ff_pre@exprs, 
                             Original_ID = seq_len(nrow(ff_pre@exprs)))
    }
    p <- ggplot(df[i,], aes(x = x, y = y)) +
      geom_point(size = 0.5,
                 color = ifelse(exprs(ff_pre)[i,"Original_ID"] %in% 
                                  exprs(ff_post)[,"Original_ID"], 'blue', 'red')) +
      xlab(marker_x) + ylab(marker_y) +
      theme_minimal() + theme(legend.position = "none") +
      ggtitle(title)
    return(p)
  }
  to_plot <- list(list(ff_pre = ff,
                       ff_post = ff_m,
                       title = "Removed margin events",
                       marker_x = "PerCP-Cy5-5-A",
                       marker_y = "BV605-A"),
                  list(ff_pre = ff_t,
                       ff_post = ff_s,
                       title = "Removed doublets",
                       marker_x = "FSC-A",
                       marker_y = "FSC-H"),
                  list(ff_pre = ff_s,
                       ff_post = ff_l,
                       title = "Removed debris and dead cells",
                       marker_x = "FSC-A",
                       marker_y = "APC-Cy7-A"),
                  list(ff_pre = ff_l,
                       ff_post = PQC$FinalFF,
                       title = "Removed low quality events",
                       marker_x = "Time",
                       marker_y = "PerCP-Cy5-5-A"))
  
  plot_list <- list()
  for (plot in to_plot) {
    plot_list[[length(plot_list) + 1]] <- filter_plot(ff_pre = plot$ff_pre, 
                                                      ff_post = plot$ff_post,
                                                      title = plot$title,
                                                      marker_x = plot$marker_x,
                                                      marker_y = plot$marker_y)
  }
  
  png(paste0(dir_QC, sub("fcs", "png", file)), width = 1920)
  print(ggpubr::ggarrange(plotlist = plot_list, nrow = 1))
  dev.off()
}

# 14. Perform quality control between all files
# 14.(A) Plot the signal per channel and per file
# 14.(A)(i) Define the variables
file_names <- sub(".*15_(.*).fcs", "\\1", files)
file_groups <- rep(c("KO", "WT"), times = c(3, 4))

# 14.(A)(ii) Make the overview plot
PlotFileScatters(input = paste0(dir_prepr, files),
                 channels = channels_of_interest,
                 names = file_names, legend = TRUE,
                 groups = file_groups, nrow = 2,
                 plotFile = paste0(dir_QC, "file_scatters.png"))

# 14.(B) Perform principal commponent analysis (PCA)
# 14.(B)(i) Retrieve the median marker expression values per file
medians <- matrix(data = NA,
                  nrow = length(files), ncol = length(channels_of_interest),
                  dimnames = list(files, channels_of_interest))

for (file in files){
  ff <- read.FCS(paste0(dir_prepr, file))
  medians[file,] <- apply(exprs(ff)[,channels_of_interest], 2, median)
}

# 14.(B)(ii) Calculate the PCs
pc <- prcomp(medians, scale. = TRUE)

# 14.(B)(iii) Visualize the PCs
ggplot(data.frame(pc$x[,1:2], file_groups)) + 
  geom_point(aes(x= PC1, y = PC2, col = file_groups)) +
  theme_minimal()
ggsave(paste0(dir_QC, "file_PCA.png"), width = 5)
  
#}, times = 10)

#### Create an aggregate file ##################################################
#microbenchmark::microbenchmark({

# 15. Choose the number of cells to include in the aggregate file
n <- 700000

# 16. Make an aggregate file
set.seed(2020)
agg <- AggregateFlowFrames(paste0(dir_prepr, files),
                           cTotal = n,
                           writeOutput = TRUE,
                           outputFile = paste0(dir_prepr, "aggregate.fcs"))

#}, times = 10)

#### Train FlowSOM model #######################################################
#microbenchmark::microbenchmark({
  
# 17. Specify the FlowSOM variables
SOM_x <- 10
SOM_y <- 10
n_meta <- 8
seed <- 2020
scaling <- FALSE

# 18. Compute the FlowSOM object
fsom <- FlowSOM(input = agg,
                scale = scaling,
                colsToUse = markers_of_interest,
                seed = seed,
                nClus = n_meta,
                xdim = SOM_x, ydim = SOM_y)
saveRDS(fsom, paste(dir_RDS, "fsom.rds"))

# 19. Visualize the FlowSOM object
PlotStars(fsom = fsom,
          backgroundValues = fsom$metaclustering)
ggsave(paste0(dir_results, "fsom_tree.pdf"),height = 8.5, width = 11)

#}, times = 10)

#### Test quality ##############################################################
#microbenchmark::microbenchmark({
  
# 20. Check the FlowSOM quality
# 20.(A) Make 2D scatter plots
# 20.(A)(i) Specify the parameters
channel_pairs = list(c("CD19", "SSC-A"),
                     c("CD3", "CD161"),
                     c("CD64", "CD49b"),
                     c("CD11c", "MHCII"),
                     c("Ly-6G", "CD11b"))
metaclusters_of_interest <- seq_len(n_meta)
clusters_of_interest <- NULL

# 20.(A)(ii) Make the 2D scatter plots
Plot2DScatters(fsom = fsom,
               channelpairs = channel_pairs,
               metaclusters = metaclusters_of_interest,
               clusters = clusters_of_interest,
               plotFile = paste0(dir_results, "fsom_2D_scatters.png"))

# 20.(B) Check the consistency with manual labeling
# 20.(B)(i) Extract the gating information from the wsp file
gating <- GetFlowJoLabels(files = files,
                          wspFile = "Data/Raw/attachments/General_panel.wsp",
                          path = dir_raw)

# 20.(B)(ii) Get an overview of the gatenames and define the cell types of interest
print(levels(gating[[1]][["manual"]]))
cell_types_of_interest <- c("B cells", "NK cells", "T cells", "Macrophages", 
                            "DCs", "Neutrophils","Non neutrophils")

# 20.(B)(iii) Compile the labels of the aggregate file
aggregate_labels <- c()
for (file in unique(exprs(agg)[, "File"])) {
  aggregate_labels <- c(aggregate_labels, 
                        as.character(ManualVector(gating[[file]][["matrix"]],
                                                  cell_types_of_interest)
                                     [exprs(agg)[, "Original_ID"]
                                       [exprs(agg)[, "File"] == file]]))
}

# 20.(B)(iv) Show the manual labeling on the FlowSOM tree
PlotPies(fsom = fsom,
         cellTypes = factor(aggregate_labels, levels = c("Unlabeled",
                                                         cell_types_of_interest)))
ggsave("Results/fsom_manual.pdf")

# 19.(B)(v) Calculate the purity of the FlowSOM clustering
Purity(realClusters = aggregate_labels,
       predictedClusters = GetClusters(fsom))

# 20.(C) Inspect the file contribution per cluster
# 20.(C)(i) Specify a color vector (optional)
file_colors <- c("#990000", "#cc0000", "#ff0000", #Different shades within the groups
                 "#1d1d77", "#2b3b92", "#3859ac", "#4677c7") 

# 20.(C)(ii) Show the file contribution
p <- PlotPies(fsom = fsom,
              cellTypes = factor(files[fsom$data[,"File"]]),
              colorPalette = file_colors)
AddStarsPies(p = p, # Legend to show how it should be
             arcs = data.frame(
               x0 = rep(0, length(files)),
               y0 = rep(0, length(files)),
               start = seq(0, 2 * pi, length.out = 8)[-8],
               end = seq(0, 2 * pi, length.out = 8)[-1],
               value = rep(1, length(files)),
               Markers = files),
             colorPalette = file_colors)
ggsave(paste0(dir_results, "fsom_filecontribution.pdf"))

#}, times = 10)

#### Discovery and downstream analysis #########################################
#microbenchmark::microbenchmark({
  
# 21. Explore the FlowSOM result
# 21.(A) Create the FlowSOMmary
FlowSOMmary(fsom = fsom,
            plotFile = paste0(dir_results, "fsom_summary.pdf"))

# 21.(B) Look for nodes with a specific pattern
# 21.(B)(i) Specify the query
query <- list("B cells" = c("CD19" = "high", "CD3" = "low"),
              "NK cells" = c("CD19" = "low", "CD161" = "high", 
                             "MHCII" = "low"),
              "T cells" = c("CD3" = "high", "MHCII" = "low", "Ly-6G" = "low"),
              "Macrophages" = c("CD64" = "high", "FcERI" = "high", 
                                "MHCII" = "high", "CD49b" = "high", 
                                "Ly-6G" = "low"),
              "Dendritic cells" = c("CD11c" = "high", "MHCII" = "high", 
                                    "CD11b" = "high", "FcERI" = "low"),
              "Neutrophils" = c("Ly-6G" = "high", "CD11b" = "high", 
                                "CD3" = "low"))

# 21.(B)(ii) Retrieve the cluster labels based on the query
labels <- QueryMultiple(fsom = fsom,
                        cellTypes = query,
                        plotFile = paste0(dir_results, "fsom_QueryStarPlot.pdf"))

# 21.(B)(iii) Show the retrieved labels on the FlowSOM tree
PlotVariable(fsom = fsom,
             variable = labels)
ggsave(paste0(dir_results, "fsom_query.pdf"))

# 22. Get features per fcs file
# Specify the variables of interest
types <- c("counts", "percentages", "MFIs")
MFIs <- c("CD49b", "Ly-6G")

# Get the features
features <- GetFeatures(fsom = fsom,
                        files = paste0(dir_prepr, files),
                        filenames = file_names,
                        type = types,
                        MFI = MFIs)

# 23. Define the groups and feature you would want to compare.
feature <- "cluster_percentages"
grouplist <- list("KO" = file_names[1:3],
                  "WT" = file_names[4:7])
stat <- "fold changes"

# 24. Compare the 2 groups of interest
# Perform the statistical tests
stats <- GroupStats(features = features[[feature]],
                    groups = grouplist)

# Define the plotting variables
stat_levels <- c(paste0(names(grouplist)[2], " underrepresented compared to ",
                        names(grouplist)[1]),
                 paste0(names(grouplist)[1], " underrepresented compared to ",
                        names(grouplist)[2]),
                 "--")
colors <- c("blue", "red", "white")

# Show statistical findings on FlowSOM trees
cluster_stat <- stats[stat,]
cluster_stat <- factor(ifelse(cluster_stat < -2.5, stat_levels[1],
                              ifelse(cluster_stat > 2.5, stat_levels[2],
                                     stat_levels[3])), 
                       levels = stat_levels)
cluster_stat[is.na(cluster_stat)] <- stat_levels[3]
gr_1 <- PlotStars(fsom = fsom, title = names(grouplist)[1], 
                  nodeSizes = stats[paste0("medians ", names(grouplist)[1]),], 
                  backgroundValues = cluster_stat,
                  backgroundColors = colors, 
                  list_insteadof_ggarrange = TRUE)
gr_2 <- PlotStars(fsom = fsom, title = names(grouplist)[2], 
                  nodeSizes = stats[paste0("medians ", names(grouplist)[2]),],
                  backgroundValues = cluster_stat,
                  backgroundColors = colors,
                  list_insteadof_ggarrange = TRUE)
ggpubr::ggarrange(plotlist = list(gr_1$tree, gr_2$tree, gr_2$starLegend, 
                                  gr_2$backgroundLegend), 
                  heights = c(3,1))
ggsave(paste0(dir_results, "fsom_groups.pdf"), width = 10, height = 7.5)


# 25. Map new data on the FlowSOM object
for (file in files){
  ff_prepr <- read.FCS(paste0(dir_prepr, file))
  ff_raw <- read.FCS(paste0(dir_raw, file))
  fsom_tmp <- NewData(fsom = fsom,
                      input = ff_prepr)
  clustering <- GetClusters(fsom_tmp)
  clustering_raw <- matrix(data = rep(0, nrow(exprs(ff_raw))),
                           ncol = 1, dimnames = list(c(), "FlowSOM"))
  clustering_raw[exprs(ff_prepr)[,"Original_ID"]] <- clustering
  ff_tmp <- flowCore::fr_append_cols(ff_raw, clustering_raw)
  write.FCS(ff_tmp, paste0(dir_prepr, "FlowSOM_", file))
}

#}, times = 10)

#### Additional FlowSOM approaches #############################################
#microbenchmark::microbenchmark({
  
### Applying FlowSOM to files or groups separately and then meta-cluster on all ####
# Compute separate FlowSOM objects
fsom_KO <- FlowSOM(input = paste0(dir_prepr, files[1:3]),
                   scale = FALSE, colsToUse = channels_of_interest,
                   seed = 2020)

fsom_WT <- FlowSOM(input = paste0(dir_prepr, files[4:7]),
                   scale = FALSE, colsToUse = channels_of_interest,
                   seed = 2020)

# Extract the cluster median fluorescence intensity values (MFIs)
MFI_KO <- GetClusterMFIs(fsom = fsom_KO, prettyColnames = TRUE, colsUsed = TRUE)
rownames(MFI_KO) <- paste0("KO", rownames(MFI_KO))
MFI_WT <- GetClusterMFIs(fsom = fsom_WT, prettyColnames = TRUE, colsUsed = TRUE)
rownames(MFI_WT) <- paste0("WT", rownames(MFI_WT))

# Obtain the meta-clusters by hierarchical clustering
all_clusters <- rbind(MFI_KO, MFI_WT)
hclust <- hclust(dist(all_clusters))
metaclustering <- cutree(hclust, 15) #MC 14 corresponds to the NK cells

# Generate one clustering heatmap from all clusters
ann <- data.frame(cohort = rep(c("KO", "WT"), each = 100), 
                  row.names = rownames(all_clusters))
p <- pheatmap::pheatmap(t(all_clusters), cluster_rows = F, cutree_cols = 15, 
                        cellwidth = 5, fontsize_col = 3, annotation_col = ann,
                        cluster_cols = hclust)
ggsave(p, filename = paste0(dir_results, "Higher_level_clustering.pdf"), width = 17)

# Generate the meta-cluster percentages boxplots
fsom_KO$metaclustering <- factor(unname(metaclustering[1:100]), levels = 1:15)
fsom_WT$metaclustering <- factor(unname(metaclustering[101:200]), levels =  1:15)
perc_KO <- GetFeatures(fsom = fsom_KO, 
                       files = paste0(dir_prepr, files[1:3]),
                       population = "metaclusters", type = "percentages", 
                       filenames = files[1:3])
perc_WT <- GetFeatures(fsom = fsom_WT, 
                       files = paste0(dir_prepr, files[4:7]),
                       population = "metaclusters", type = "percentages", 
                       filenames = files[4:7])

df <- data.frame(rbind(perc_KO[[1]], perc_WT[[1]])*100, 
                 cohort = rep(c("KO", "WT"), c(3, 4)), check.names = FALSE)
df_g <- tidyr::gather(df, "MC", "percentage", -cohort)
ggplot(df_g, aes(x = cohort, y = percentage)) +
  geom_boxplot() +
  facet_wrap(~MC, scales = "free") +
  theme_minimal()
ggsave(filename = "Results/FlowSOM_boxplot.pdf", width = 10, height = 10)

#}, times = 10)

## Hierarchical approach #######################################################
#microbenchmark::microbenchmark({
  
# Read in preprocessed fcs file, lymphocyte panel
ff <- read.FCS(paste0(dir_raw, "lympho.fcs"))
manual_labels <- readRDS(paste0(dir_raw, "attachments/lympho_labels.rds"))

# Perform a first level clustering to isolate the lymphocytes
fsom_level1 <- FlowSOM(input = ff,
                       scale = FALSE,
                       colsToUse = c("CD11b", "CD3", "CD161", "CD19"),
                       seed = 2020)

# Inspect the 2D scatter plots to identify the meta-clusters of interest
Plot2DScatters(fsom = fsom_level1, 
               channelpairs = list(c("CD3", "CD161")),
               metaclusters = 1:10, 
               plotFile = paste0(dir_results, "hierarchy_level1.png")) 
#MC 1, 4 and 5 are the lymphocytes (CD3+, CD161-)

# Subset the original fcs file
fsom_tmp <- NewData(fsom_level1, ff)
clustering <- GetMetaclusters(fsom_tmp)
ff_tmp <- ff[clustering %in% c(1, 4, 5),]

# Perform a second level clustering to characterize the lymphocytes
fsom_level2 <- FlowSOM(input = ff_tmp,
                       scale = FALSE,
                       colsToUse = c("TCRyd", "CD44", "CD4", "CD62L", "CD8"),
                       seed = 2020)

# Plot the lymphocytes FlowSOM tree
PlotStars(fsom = fsom_level2,
          backgroundValues = fsom_level2$metaclustering)

# Show the manual labels on the FlowSOM trees
PlotPies(fsom = fsom_level1, cellTypes = manual_labels,
         title = "First level clustering")
PlotPies(fsom = fsom_level2, cellTypes = factor(manual_labels[clustering %in% c(1, 4, 5)]),
         backgroundValues = fsom_level2$metaclustering, title = "Second level clustering")

#}, times = 10)

#### Figures for paper #########################################################
### 2 
tiff(paste0(dir_results, "FIGURE_2.tiff"), width = 4000, height = 1000, res = 330)
print(ggpubr::ggarrange(plotlist = plot_list, nrow = 1)) +
  theme(plot.margin = unit(c(5.5,12,5.5,5.5),"pt"))
dev.off()

### 3 
ff <- AggregateFlowFrames(paste0(dir_prepr,files), 
                          cTotal = 50000, 
                          channels = channels_of_interest)
data <- ff@exprs
file_values <- data[, "File"]
subset <- sample(seq_len(nrow(data)), min(50000, nrow(data)))
data <- data[subset,]

file_values <- file_values[subset]
channels <- channels_of_interest
names = file_names
groups = file_groups
plots_list <- list()
for (channel in channels) {
  df <- data.frame(intensity = data[, channel], 
                   names = factor(names[file_values], levels = unique(names)), 
                   group = factor(groups[file_values], levels = unique(groups)))
  p <- ggplot(df, aes(.data$names, .data$intensity)) + 
    geom_jitter(position = position_jitter(width = 0.1), 
                alpha = 0.5, aes(colour = .data$group), shape = ".") + 
    ylab(GetMarkers(ff, channel)) + theme_classic() + 
    theme(axis.title.x = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, size = 12),
          axis.title.y = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12),
          legend.text = element_text(size = 12)) + 
    guides(colour = guide_legend(override.aes = list(size = 5, shape = 15, 
                                                     alpha = 1)))
  plots_list[[length(plots_list) + 1]] <- p
}

tiff(paste0(dir_results, "FIGURE_3.tiff"), width = 800, height = 500)
ggpubr::annotate_figure(ggpubr::ggarrange(plotlist = plots_list, 
                                          common.legend = TRUE, 
                                          ncol = 6, nrow = 2), 
                        bottom = ggpubr::text_grob("Files", size = 15), 
                        left = ggpubr::text_grob("Fluorescence intensity", size = 15,
                                                 rot = 90))
dev.off()


### 4
### 4A Scatter channel not scaled
dir.create("Data/Preprocessed_wrong")
dir.create("Data/Preprocessed_wrong/QC")
dir_wprepr <- "Data/Preprocessed_wrong"
dir_wQC <- "Data/Preprocessed_wrong/QC"

translist <- NULL
for (file in files){
  # Read the fcs file into a flowframe
  ff <- read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE)
  
  # Remove margin events
  ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest)
  
  # Compensate
  ff_c <- flowCore::compensate(ff_m, compensation_matrix)
  
  # Transform, logicle for marker channels, "forget" the SSC-A scaling
  if (is.null(translist)) { # Calculate transformation on first file and apply to all
    translist <- estimateLogicle(ff_c,  
                                 colnames(compensation_matrix))
  }
  ff_t <- flowCore::transform(ff_c, translist)
  
  # Remove doublets and filter live cells
  ff_s <- PeacoQC::RemoveDoublets(ff_t)
  
  selected_live <- filter(ff_s, live_gate)
  ff_l <- ff_s[selected_live@subSet, ]
  
  # QC with PeacoQC
  PQC <- PeacoQC::PeacoQC(ff = ff_l,
                          channels = channels_of_interest,
                          plot = FALSE, save_fcs = FALSE,
                          output_directory = dir_wQC)
  
  # Save the preprocessed data
  write.FCS(PQC$FinalFF,
            file = paste0(dir_wprepr, file))
}

set.seed(2020)
wagg1 <- AggregateFlowFrames(paste0(dir_wprepr, files),
                             cTotal = n,
                             writeOutput = TRUE,
                             outputFile = paste0(dir_wprepr, "waggregate1.fcs"))

wfsom1 <- FlowSOM(input = wagg1,
                  scale = scaling,
                  colsToUse = markers_of_interest[1:6],
                  seed = seed,
                  nClus = n_meta,
                  xdim = SOM_x, ydim = SOM_y)
saveRDS(wfsom1, paste(dir_RDS, "wfsom1.rds"))

p1 <- PlotStars(fsom = wfsom1,
                maxNodeSize = 20, list_insteadof_ggarrange = TRUE)
p1 <- ggpubr::ggarrange(p1$starLegend, p1$tree,
                        nrow = 2, heights = c(1,5)) +
  theme(plot.margin=margin(10,0,0,0))

### 4B Batch effect, here introduced by an additional transformation
translist <- NULL
translist_batch <- c(transformList("BV605-A", flowCore::linearTransform(a = 1.2,
                                                                        b = 0.5)),
                     transformList("BV786-A", flowCore::linearTransform(a = 0.7,
                                                                        b = -0.3)),
                     transformList("APC-A", flowCore::linearTransform(a = 1.3,
                                                                      b = 0)),
                     transformList("PE-A", flowCore::linearTransform(a = 1,
                                                                     b = 0.5)))
file_names_ex <- c("day1 1", "day1 2", "day1 3",
                   "day2 1", "day2 2", "day2 3")
for (file in files){
  # Read the fcs file into a flowframe
  ff <- read.FCS(paste0(dir_raw, file), truncate_max_range = FALSE)
  
  # Remove margin events
  ff_m <- PeacoQC::RemoveMargins(ff, channels_of_interest)
  
  # Compensate
  ff_c <- flowCore::compensate(ff_m, compensation_matrix)
  
  # Transform, logicle for marker channels, "forget" the SSC-A scaling
  if (is.null(translist)) { # Calculate transformation on first file and apply to all
    translist <- estimateLogicle(ff_c,  
                                 colnames(compensation_matrix))
    ff_t <- flowCore::transform(ff_c, translist)
    q5_goal <- quantile(exprs(ff_t)[,"PE-A"], 0.05)
    q95_goal <- quantile(exprs(ff_t)[,"PE-A"], 0.95)
    q5_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.05)
    q95_SSCA <- quantile(exprs(ff_t)[,"SSC-A"], 0.95)
    SSCA_a <- (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
    SSCA_b <- q5_goal - q5_SSCA * (q95_goal - q5_goal) / (q95_SSCA - q5_SSCA)
    translist <- c(translist, 
                   transformList("SSC-A", flowCore::linearTransform(a = SSCA_a,
                                                                    b = SSCA_b)))
  }
  ff_t <- flowCore::transform(ff_c, translist)
  if (file %in% files[4:7]){
    ff_t <- flowCore::transform(ff_t, translist_batch)
  }
  
  # Remove doublets and filter live cells
  ff_s <- PeacoQC::RemoveDoublets(ff_t)
  
  selected_live <- filter(ff_s, live_gate)
  ff_l <- ff_s[selected_live@subSet, ]
  
  # QC with PeacoQC
  PQC <- PeacoQC::PeacoQC(ff = ff_l,
                          channels = channels_of_interest,
                          plot = FALSE, save_fcs = FALSE,
                          output_directory = dir_wQC)
  
  # Save the preprocessed data
  write.FCS(PQC$FinalFF,
            file = paste0(dir_wprepr, file))
}

plots <- PlotFileScatters(input = paste0(dir_wprepr, files[1:6]),
                          channels = channels_of_interest,
                          names = file_names_ex, legend = TRUE,
                          groups = rep(c("day 1", "day 2"), each = 3),
                          plotFile = NULL, color = c("red", "navy"))
p2 <- ggpubr::annotate_figure(ggpubr::ggarrange(plotlist = plots[4:9], 
                                                ncol = 3, nrow = 2,
                                                common.legend = TRUE) +
                                theme(plot.margin = margin(10,10,10,10)),
                              bottom = ggpubr::text_grob("Files"),
                              left = ggpubr::text_grob("Fluorescence intensity", rot = 90))

set.seed(2020)
wagg2 <- AggregateFlowFrames(paste0(dir_wprepr, files[1:6]),
                             cTotal = n,
                             writeOutput = TRUE,
                             outputFile = paste0(dir_wprepr, "waggregate2.fcs"))

wfsom2 <- FlowSOM(input = wagg2,
                  scale = scaling,
                  colsToUse = markers_of_interest,
                  seed = seed,
                  nClus = n_meta,
                  xdim = SOM_x, ydim = SOM_y)
saveRDS(wfsom2, paste(dir_RDS, "wfsom2.rds"))

file_factor <- factor(files[wfsom2$data[,"File"]])
levels(file_factor) <- file_names_ex
p3 <- PlotPies(fsom = wfsom2,
               cellTypes = file_factor,
               colorPalette = file_colors,
               maxNodeSize = 1.2)


### 4C Manual vector not in correct order
mock_manual_types <- c("Unlabeled", paste0("Cell Type ", 1:6))
mock_manual_labels <- factor(rep(mock_manual_types,
                                 c(50000, 120000, 100000, 130000, 
                                   90000, 130000, 80000)),
                             levels = mock_manual_types)
p4 <- PlotPies(fsom = fsom,
               cellTypes = mock_manual_labels)

tiff(paste0(dir_results, "FIGURE_4.tiff"), width = 4000, height = 2000, res = 330)
ggpubr::ggarrange(p1, p2, p4, p3,
                  labels = c("a", "c", "b", "d"))
dev.off()


### 5
p <- PlotStars(fsom = fsom,
               backgroundValues = fsom$metaclustering,
               list_insteadof_ggarrange = TRUE)

tiff(paste0(dir_results, "FIGURE_5.tiff"), width = 4000, height = 2400, res = 330)
ggpubr::ggarrange(p$tree, 
                  ggpubr::ggarrange(p$starLegend, p$backgroundLegend, 
                                    nrow = 2),
                  ncol = 2, widths = c(1.5,1))
dev.off()

### 6
tiff(paste0(dir_results, "FIGURE_6.tiff"), width = 3000, height = 1800, res = 400)
PlotPies(fsom = fsom,
         cellTypes = factor(aggregate_labels, 
                            levels = c("Unlabeled", cell_types_of_interest)))
dev.off()

### 7
tiff(paste0(dir_results, "FIGURE_7.tiff"), width = 3000, height = 1800, res = 400)
ggpubr::ggarrange(plotlist = list(gr_1$tree, gr_2$tree, 
                                  gr_2$starLegend, gr_2$backgroundLegend), 
                  heights = c(3,1))
dev.off()


### Session info ###############################################################
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")

#### Timing ####################################################################
# microbenchmark::microbenchmark(PlotFileScatters(input = paste0(dir_prepr, files),
#                                                 channels = channels_of_interest,
#                                                 names = file_names, legend = TRUE,
#                                                 groups = file_groups, nrow = 2,
#                                                 plotFile = paste0(dir_QC, "file_scatters.png")),
#                                times = 10)
# 
# microbenchmark::microbenchmark(AggregateFlowFrames(paste0(dir_prepr, files),
#                                                    cTotal = n,
#                                                    writeOutput = TRUE,
#                                                    outputFile = paste0(dir_prepr, "aggregate.fcs")),
#                                times = 10)
# 
#  microbenchmark::microbenchmark(FlowSOM(input = agg,
#                                         scale = scaling,
#                                         colsToUse = markers_of_interest,
#                                         seed = seed,
#                                         nClus = n_meta,
#                                         xdim = SOM_x, ydim = SOM_y),
#                                 times = 10)
#  
# microbenchmark::microbenchmark(PlotStars(fsom = fsom,
#                                          backgroundValues = fsom$metaclustering),
#                                times = 10)
# 
# microbenchmark::microbenchmark(Plot2DScatters(fsom = fsom,
#                                               channelpairs = channel_pairs,
#                                               metaclusters = metaclusters_of_interest,
#                                               clusters = clusters_of_interest,
#                                               plotFile = paste0(dir_results, "fsom_2D_scatters.png")),
#                                times = 10)
# 
# microbenchmark::microbenchmark(GetFlowJoLabels(files = files,
#                                                wspFile = "Data/Raw/attachments/General_panel.wsp",
#                                                path = dir_raw),
#                                times = 10)
# 
# microbenchmark::microbenchmark(Purity(realClusters = aggregate_labels,
#                                       predictedClusters = GetClusters(fsom)),
#                                times = 10)
# 
# microbenchmark::microbenchmark(PlotPies(fsom = fsom,
#                                         cellTypes = factor(aggregate_labels, levels = c("Unlabeled",
#                                                                                         cell_types_of_interest))),
#                                times = 10)
# 
# microbenchmark::microbenchmark(FlowSOMmary(fsom = fsom,
#                                            plotFile = paste0(dir_results, "fsom_summary.pdf")),
#                                times = 10)
# 
# microbenchmark::microbenchmark(QueryMultiple(fsom = fsom,
#                                              cellTypes = query,
#                                              plotFile = paste0(dir_results, "fsom_QueryStarPlot.pdf")),
#                                times = 10)
# 
# microbenchmark::microbenchmark(PlotVariable(fsom = fsom,
#                                             variable = labels),
#                                times = 10)
# 
# microbenchmark::microbenchmark(GetFeatures(fsom = fsom,
#                                            files = paste0(dir_prepr, files),
#                                            filenames = file_names,
#                                            type = types,
#                                            MFI = MFIs),
#                                times = 10)
# 
# microbenchmark::microbenchmark(GroupStats(features = features[[feature]],
#                                           groups = grouplist),
#                                times = 10)
# 
# microbenchmark::microbenchmark(fsom_tmp <- NewData(fsom = fsom,
#                                                    input = ff_prepr),
#                                times = 10)







