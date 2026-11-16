
woodpecker_metadata<-read.csv('/media/ssd/Raw_data/Picoides/Nuclear_data_02_2025/woodpeckers_metadata.csv',
                              sep=',')

#subset to only keep BBWO
arcticus_metadata<-woodpecker_metadata[woodpecker_metadata$Species=='Picoides arcticus',]
p_arcticus_sf<-st_as_sf(arcticus_metadata,coords=c('longitude','latitude'),crs=4269)


##make a metadata fior population groups

meta <- arcticus_metadata[, c("strpID", "State_Province")]
meta$strpID <- trimws(meta$strpID)
meta$State_Province <- gsub("Ê", "", meta$State_Province)


#assign population
### clean state_province names
#p_arcticus_sf$State_Province <- sub("Ê$", "", p_arcticus_sf$State_Province)

###add new column ##this is if i watn separated like xavier in his memoire
#p_arcticus_sf <- p_arcticus_sf %>%
 # mutate(
  #  POP = case_when(
   #   State_Province %in% c("California", "Washington", "Oregon") ~ "West",
    #  State_Province %in% c("Michigan", "Minnesota") ~ "South_central",
     # State_Province %in% c("Manitoba", "Alberta") ~ "North_central",
    #  State_Province == "New York" ~ "NewYork",
    #  State_Province == "Quebec" ~ "Quebec",
    #  TRUE ~ NA_character_
    #)
#  )

#p_arcticus_sf$POP[p_arcticus_sf$strpID=='MB19073'] <- "South_central"



########
##PCA##
########

cov <- as.matrix(read.table("/media/ssd/Bioinformatics/downstream_analyses/PCAngsd/PCAngsd.cov", 
                            header = F))

bamlist <- readLines("/media/ssd/popglen_workflow/results/datasets/p_arcticus_02_2026/bamlists/p_arcticus_02_2026.dataset-ref_all.bamlist")

samples <- basename(bamlist)
samples <- sub("\\.dataset.*", "", samples)   # removes suffix after ID

cov_names <- as.matrix(read.table(
  "/media/ssd/Bioinformatics/downstream_analyses/PCAngsd/PCAngsd.cov",
  header = FALSE
))

meta_pca <- meta[match(samples, meta$strpID), ]

west_states <- c("California", "Oregon", "Washington")

subgroups_pca <- ifelse(
  meta_pca$State_Province %in% west_states,
  "West",
  meta_pca$State_Province
)

pca_nuDNA <- eigen(cov_names)

eigenvectors <- pca_nuDNA$vectors

pca_vectors <- as_tibble(
  cbind(
    sample = samples,
    pop = subgroups_pca,
    data.frame(eigenvectors)
  )
)

##look at PC variance explained

pca_eigenval_sum = sum(pca_nuDNA$values) #sum of eigenvalues
varPC1 <- (pca_nuDNA$values[1]/pca_eigenval_sum)*100 #Variance explained by PC1
varPC2 <- (pca_nuDNA$values[2]/pca_eigenval_sum)*100 #Variance explained by PC2
varPC3 <- (pca_nuDNA$values[3]/pca_eigenval_sum)*100 #Variance explained by PC3
varPC4 <- (pca_nuDNA$values[4]/pca_eigenval_sum)*100 #Variance explained by PC4

pop_colors_pca <- c(
  "West" = "#377eb8", 
  "Alberta" = "#a6cee3",
  "Manitoba" = "#b2df8a",
  "Minnesota" = "#33a02c",
  "Michigan" = "#fb9a99",
  "New York" = "#e31a1c",
  "Quebec" = "#fdbf6f"
)

ggplot(pca_vectors, aes(X1, X2, color = pop)) +
  geom_point(size = 3) +
  xlab("PC1 (4.57%)") +
  ylab("PC2 (2.46%)")+
  scale_color_manual(values = pop_colors_pca) +
  theme_grey()


#########
#NGSadmix
#########

#same groupings as pca 

subgroups_ngsadmix <- meta$State_Province[
  match(samples, meta$strpID)
]


west_states <- c("California", "Oregon", "Washington")

meta_pca$pop_group <- ifelse(
  meta_pca$State_Province %in% west_states,
  "West",
  meta_pca$State_Province
)

#order from west to east 
pop_levels <- c( "California","Oregon","Washington" ,"Alberta", "Manitoba", "Minnesota", "Michigan", "New York", "Quebec" ) 

pop_colors <- c(
  "California" = "#377eb8", 
  "Oregon" = "#377eb8", 
  "Washington" = "#377eb8", 
  "Alberta" = "#a6cee3",
  "Manitoba" = "#b2df8a",
  "Minnesota" = "#33a02c",
  "Michigan" = "#fb9a99",
  "New York" = "#e31a1c",
  "Quebec" = "#fdbf6f"
)


pop_group <- factor(
  subgroups_ngsadmix,
  levels = pop_levels
)

ord <- order(pop_group)

#all(meta_pca$strpID == samples) #check to see if it is well aligned 

##now load the data
ngsADMIX_files<-list.files(path="/media/ssd/Bioinformatics/downstream_analyses/NGSadmix",pattern='qopt',full.names=TRUE)

#for k=2
k2_files <- ngsADMIX_files[grepl("k2",ngsADMIX_files)]
data_k2 <- read.table(k2_files[1])

data_k2 <- data_k2[ord, ]
pop_group <- pop_group[ord]

###better plot
par(mar = c(10, 4, 2, 1))  # bottom margin bigger
bp <- barplot(
  t(data_k2),
  col = c("grey30", "grey80"),
  border = NA,
  space = 0,
  xaxt = "n",
  ylab = "Admixture proportion"
)
separators <- tapply(bp, pop_group, range)
separators <- sapply(separators, function(x) x[2])
separators <- separators[-length(separators)]

abline(v = separators + 0.5, col = "white", lwd = 1)

rect(
  xleft = bp - 0.5,
  xright = bp + 0.5,
  ybottom = -0.06,
  ytop = 0,
  col = pop_colors[as.character(pop_group)],
  border = NA,
  xpd = TRUE
)
group_centers <- tapply(bp, pop_group, mean)

axis(
  1,
  at = group_centers,
  labels = pop_levels,
  tick = FALSE,
  las = 2,
  cex.axis = 0.8,
  line = 1   # <-- THIS is what fixes overlap
)



#for k=3
k3_files <- files[grepl("k3",ngsADMIX_files)]
data_k3 <- read.table(k3_files[1])

data_k3 <- data_k3[ord, ]
pop_group <- pop_group[ord]

###better plot
par(mar = c(10, 4, 2, 1))  # bottom margin bigger
bp <- barplot(
  t(data_k3),
  col = c("grey30", "grey60","grey90"),
  border = NA,
  space = 0,
  xaxt = "n",
  ylab = "Admixture proportion"
)
separators <- tapply(bp, pop_group, range)
separators <- sapply(separators, function(x) x[2])
separators <- separators[-length(separators)]

abline(v = separators + 0.5, col = "white", lwd = 1)

rect(
  xleft = bp - 0.5,
  xright = bp + 0.5,
  ybottom = -0.06,
  ytop = 0,
  col = pop_colors[as.character(pop_group)],
  border = NA,
  xpd = TRUE
)
group_centers <- tapply(bp, pop_group, mean)

axis(
  1,
  at = group_centers,
  labels = pop_levels,
  tick = FALSE,
  las = 2,
  cex.axis = 0.8,
  line = 1   # <-- THIS is what fixes overlap
)


##evaluating the ngsadmix model with Eval ADMIX 

mat <- as.matrix(read.table("/media/ssd/Bioinformatics/downstream_analyses/EvalAdmix/ALL/evaladmix_LD_pruned_ngsADMIX_k1_rep1_BBWO"))


ord <- order(pop_province)
mat_ord <- mat[ord, ord]

# IMPORTANT FIX HERE
lim <- max(abs(mat_ord), na.rm = TRUE)

cols <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

pheatmap(
  mat_ord,
  color = cols,
  breaks = seq(-lim, lim, length.out = 101),
  
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  border_color = NA,
  na_col = "grey90",   # optional but nice for diagonal NA
  
  show_rownames = FALSE,
  show_colnames = FALSE,
  legend = TRUE
)


###Isolation by distance BBWO

## mantel test

cov <- as.matrix(read.table("/media/ssd/Bioinformatics/downstream_analyses/PCAngsd/PCAngsd.cov", 
                            header = F))

gen_dist<-1-cov

bamlist <- readLines("/media/ssd/popglen_workflow/results/datasets/p_arcticus_02_2026/bamlists/p_arcticus_02_2026.dataset-ref_all.bamlist")

samples <- basename(bamlist)
samples <- sub("\\.dataset.*", "", samples)   # removes suffix after ID

##Compute geogrpahic distance

library(geosphere) #in order to take into account the sphere 

woodpecker_metadata<-read.csv('/media/ssd/Raw_data/Picoides/Nuclear_data_02_2025/woodpeckers_metadata.csv',
                              sep=',')
arcticus_metadata<-woodpecker_metadata[woodpecker_metadata$Species=='Picoides arcticus',]
arcticus_metadata_sorted <- arcticus_metadata[order(match(arcticus_metadata$strpID, samples)), ]

geo_coords <- data.frame(arcticus_metadata_sorted$longitude, arcticus_metadata_sorted$latitude) 
geo_dist <- distm(geo_coords, fun = distHaversine)
geo_dist <- as.matrix(geo_dist)

##Compute environmental distances 

env_data<-readRDS('/media/ssd/Bioinformatics/nuclear_BBWO_fire_proportion.RDS')
env_data_sorted <- env_data[match(samples, env_data$ID), ]

fire_dist<-dist(env_data_sorted$burned_area_m2, method = "euclidean")
fire_matrix<-as.matrix(fire_dist)


###Mantel test
#spearman is used since its not assumed to be linear the relationship
library(vegan)
IBD<- mantel(gen_dist, geo_dist, method = "spearman", permutations = 9999, na.rm = TRUE)
IBE<- mantel(gen_dist, fire_matrix, method = "spearman", permutations = 9999, na.rm = TRUE)


IBD_partial<-mantel.partial(gen_dist, geo_dist, fire_matrix,method = "spearman", permutations = 9999, na.rm = TRUE)
IBE_partial<-mantel.partial(gen_dist, fire_matrix, geo_dist,method = "spearman", permutations = 9999, na.rm = TRUE)

#plot 
gen_vec <- as.vector(gen_dist[lower.tri(gen_dist)])
geo_vec <- as.vector(geo_dist[lower.tri(geo_dist)])
fire_vec <- as.vector(fire_matrix[lower.tri(fire_matrix)])

df_ibd <- data.frame(
  geo = geo_vec,
  gen = gen_vec
)

df_ibe <- data.frame(
  fire = fire_vec,
  gen = gen_vec
)

library(ggplot2)

ggplot(df_ibd, aes(x = geo, y = gen)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("r = ", round(IBD_partial$statistic, 3),
                          "\np = ", IBD_partial$signif),
           hjust = 1.1, vjust = 1.1)+
  labs(title = "Isolation by Distance",
       x = "Geographical distance (m)",
       y = "Genetic distance")

ggplot(df_ibe, aes(x = fire, y = gen)) +
  geom_point(alpha = 0.6) +
  geom_smooth(method = "lm", se = FALSE) +
  theme_classic() +
  annotate("text",
           x = Inf, y = Inf,
           label = paste0("r = ", round(IBE_partial$statistic, 3),
                          "\np = ", IBE_partial$signif),
           hjust = 1.1, vjust = 1.1)+
  labs(title = "Isolation by Environment (fire)",
       x = "Environmental distance",
       y = "Genetic distance")

#################################
###########FST PLOT##############
#################################

######Xavier's populations##########

library(tidyverse)

# Directory containing your FST results
fst_all_dir <- "/media/ssd/Bioinformatics/p_arcticus_02/downstream_analysis/SFS/outputs/xavier_pop_FST"

# Population order
pops <- c('West', 'South_central', 'North_central', 'NewYork', 'Quebec')


#Get only the .global.txt files
files_all <- list.files(
  path = fst_all_dir,
  pattern = "\\.global\\.txt$",
  full.names = TRUE
)

# Check that 45 files were found
length(files_all)
basename(files_all) #ok all files are there


# Read the pairwise FST results
fst_data_all <- map_dfr(files_all, function(files_all) {
  
  values <- scan(files_all, quiet = TRUE)
  
  populations <- str_remove(
    basename(files_all),
    "\\.global\\.txt$"
  ) %>%
    str_split("\\.") %>%
    unlist()
  
  tibble(
    Pop1 = populations[1],
    Pop2 = populations[2],
    FST_unweighted = values[1],
    FST_weighted = values[2]
  )
})

# Create the symmetric FST matrix
fst_matrix_all <- matrix(
  NA,
  nrow = length(pops),
  ncol = length(pops),
  dimnames = list(pops, pops)
)

# Set diagonal to zero
diag(fst_matrix_all) <- 0

# Fill matrix with weighted FST values
for (i in seq_len(nrow(fst_data_all))) {
  
  p1 <- fst_data_all$Pop1[i]
  p2 <- fst_data_all$Pop2[i]
  
  fst_matrix_all[p1, p2] <- fst_data_all$FST_weighted[i]
  fst_matrix_all[p2, p1] <- fst_data_all$FST_weighted[i]
}

# Inspect matrix
round(fst_matrix_all, 4)

#get only the upper matrix since up and down are the same values
fst_matrix_all[lower.tri(fst_matrix_all)] <- NA

#flatten the matrix
library(reshape2)

melt_fst_matrix_all<-melt(fst_matrix_all)

fst_all_plot <- ggplot(melt_fst_matrix_all, aes(x = Var1, y = Var2, fill = value)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(
    aes(label = ifelse(is.na(value), "", sprintf("%.4f", value))),
    size = 3
  ) +
  scale_fill_gradient(
    low = "#deebf7",
    high = "#08519c",
    na.value = "white",
    name = expression(F[ST])
  ) +
  coord_fixed() +
  theme_minimal(base_size = 13) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank()
  )

print(fst_all_plot)
