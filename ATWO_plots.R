library(dplyr)

##########
#metadata
##########

woodpecker_metadata<-read.csv('/media/ssd/Raw_data/Picoides/Nuclear_data_02_2025/woodpeckers_metadata.csv',
                              sep=',')

#subset to only keep BBWO
dorsalis_metadata<-woodpecker_metadata[woodpecker_metadata$Species=='Picoides dorsalis',]

##in this woodpecker metadata there are some missing coordinates, but in MArie-Pier's memoir, they're present. 
###add missing coordinates 
dorsalis_metadata[is.na(dorsalis_metadata$latitude), ] #look at which row is missing 

dorsalis_metadata[dorsalis_metadata$strpID == "UWBM125469", 9] <- 48.692 #latitude
dorsalis_metadata[dorsalis_metadata$strpID == "UWBM125469", 10] <- 117.519 #longitude

dorsalis_metadata[dorsalis_metadata$strpID == "UWBM123111", 9] <- 61.5 #latitude
dorsalis_metadata[dorsalis_metadata$strpID == "UWBM123111", 10] <- -147 #longitude

#make it a spatial object
p_dorsalis_sf<-st_as_sf(dorsalis_metadata,coords=c('longitude','latitude'),crs=4269)

##make a metadata fior population groups

meta <- dorsalis_metadata[, c("strpID", "State_Province")]
meta$strpID <- trimws(meta$strpID)
meta$State_Province <- gsub("Ê", "", meta$State_Province) #clean state and prvince names
p_dorsalis_sf$State_Province <- sub("Ê$", "", p_dorsalis_sf$State_Province)

#assign population
###population as Marie-Pier's memoir 
p_dorsalis_sf <- p_dorsalis_sf %>%
 mutate(
  POP = case_when(
   State_Province %in% c("New Mexico", "Arizona","Utah","Colorado") ~ "South",
   State_Province %in% c("Washington", "Oregon","Idaho","British Columbia") ~ "Pacific",
  State_Province %in% c("Alaska") ~ "Alaska",
 State_Province %in% c("Manitoba") ~ "Manitoba",
  State_Province == "Alberta" ~ "Alberta",
  State_Province == "Quebec" ~ "Quebec",
  TRUE ~ NA_character_
)
  )

###population as morphological subspecies
p_dorsalis_sf <- p_dorsalis_sf %>%
  mutate(
    subspecies = case_when(
      State_Province %in% c("New Mexico", "Arizona","Utah","Colorado") ~ "dorsalis",
      State_Province %in% c("Washington", "Oregon","Idaho","British Columbia","Alaska","Manitoba","Alberta") ~ "fasciatus",
      State_Province == "Quebec" ~ "bacatus",
      TRUE ~ NA_character_
    )
  )


########
#PCA
######

cov <- as.matrix(read.table("/media/ssd/Bioinformatics/p_dorsalis_07/downstream_analyses/PCAngsd/PCAngsd_all.cov", 
                            header = F))

bamlist <- readLines("/media/ssd/popglen_workflow/results/datasets/p_dorsalis_07_2026/bamlists/p_dorsalis_07_2026.dataset-ref_all.bamlist")

samples <- basename(bamlist)
samples <- sub("\\.dataset.*", "", samples)   # removes suffix after ID

cov_names <- as.matrix(read.table(
  "/media/ssd/Bioinformatics/p_dorsalis_07/downstream_analyses/PCAngsd/PCAngsd_all.cov",
  header = FALSE
))

meta_pca <- meta[match(samples, meta$strpID), ]


subgroups_pca_population <- p_dorsalis_sf$POP #marie-pier's population
subgroups_pca_subspecies<-p_dorsalis_sf$subspecies #by subspecies

pca_nuDNA <- eigen(cov_names)

eigenvectors <- pca_nuDNA$vectors

pca_vectors <- as_tibble(
  cbind(
    sample = samples,
    pop = subgroups_pca_population,
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
  "Pacific" = "#377eb8", 
  "Alberta" = "#a6cee3",
  "Manitoba" = "#b2df8a",
  "South" = "#33a02c",
  "Alaska" = "#fb9a99",
  "Quebec" = "#fdbf6f"
)

ggplot(pca_vectors, aes(X1, X2, color = pop)) +
  geom_point(size = 3) +
  xlab("PC1 (13.50%)") +
  ylab("PC2 (1.95%)")+
  scale_color_manual(values = pop_colors_pca) +
  theme_grey()
