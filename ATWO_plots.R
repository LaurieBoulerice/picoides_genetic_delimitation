library(dplyr)
library(sf)

########sfheaders##########
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
   State_Province %in% c("Washington", "Oregon","Idaho","British Columbia") ~ "Pacific", #marie-Pier wrote Idaho in her memoir but there is no Idaho samples
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


pca_nuDNA <- eigen(cov_names)

eigenvectors <- pca_nuDNA$vectors


##look at PC variance explained

pca_eigenval_sum = sum(pca_nuDNA$values) #sum of eigenvalues
varPC1 <- (pca_nuDNA$values[1]/pca_eigenval_sum)*100 #Variance explained by PC1
varPC2 <- (pca_nuDNA$values[2]/pca_eigenval_sum)*100 #Variance explained by PC2
varPC3 <- (pca_nuDNA$values[3]/pca_eigenval_sum)*100 #Variance explained by PC3
varPC4 <- (pca_nuDNA$values[4]/pca_eigenval_sum)*100 #Variance explained by PC4

########################################
#By Marie-Pier's population assignment
########################################
subgroups_pca_population<- p_dorsalis_sf$POP[ #natch to follow the same order 
  match(samples, p_dorsalis_sf$strpID)
]

pca_vectors <- as_tibble(
  cbind(
    sample = samples,
    pop = subgroups_pca_population,
    data.frame(eigenvectors)
  )
)

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



########################################
#By Msubspecies geographical range
########################################

subgroups_pca_subspecies<- p_dorsalis_sf$subspecies[ #natch to follow the same order 
  match(samples, p_dorsalis_sf$strpID)
]


pca_vectors_subspecies <- as_tibble(
  cbind(
    sample = samples,
    pop = subgroups_pca_subspecies,
    data.frame(eigenvectors)
  )
)

pop_colors_pca_subspecies <- c(
  "fasciatus" = "#377eb8", 
  "bacatus" = "#b2df8a",
  "dorsalis" = "#fb9a99"
)

ggplot(pca_vectors_subspecies, aes(X1, X2, color = pop)) +
  geom_point(size = 3) +
  xlab("PC1 (13.50%)") +
  ylab("PC2 (1.95%)")+
  scale_color_manual(values = pop_colors_pca_subspecies) +
  theme_grey()


########################################
#By state and province 
########################################

subgroups_pca_state<- p_dorsalis_sf$State_Province[ #natch to follow the same order 
  match(samples, p_dorsalis_sf$strpID)
]


pca_vectors_state <- as_tibble(
  cbind(
    sample = samples,
    pop = subgroups_pca_state,
    data.frame(eigenvectors)
  )
)


pop_colors_pca_state <- c(   "British Columbia" = "#377eb8",
                             "Idaho"            = "#4daf4a",
                             "Alaska"           = "#e41a1c",
                             "Washington"       = "#984ea3",
                             "Oregon"           = "#ff7f00",
                             "Manitoba"         = "#b2df8a",
                             "Alberta"          = "#a6cee3",
                             
                             # Southern states — same color
                             "New Mexico"       = "#fb9a99",
                             "Utah"             = "#fb9a99",
                             "Colorado"         = "#fb9a99",
                             "Arizona"          = "#fb9a99",
                             
                             "Quebec"           = "#fdbf6f"
)


ggplot(pca_vectors_state, aes(X1, X2, color = pop)) +
  geom_point(size = 3) +
  xlab("PC1 (13.50%)") +
  ylab("PC2 (1.95%)")+
  scale_color_manual(values = pop_colors_pca_state) +
  theme_grey()
geom_text(aes(label = sample), vjust = -0.7, size = 3)

##########
#NGSadmix
##########

###################################
#Marie-Pier's population assignment
###################################


subgroups_ngsadmix <- p_dorsalis_sf$State_Province[
  match(samples, p_dorsalis_sf$strpID)
]

# Order from west to east
pop_levels <- c(
  "New Mexico", "Arizona", "Colorado", "Utah",
  "Oregon", "Washington", "British Columbia", "Alaska",
  "Alberta", "Manitoba", "Quebec"
)

pop_colors_ngsadmix <- c(
  "Oregon" = "#377eb8",
  "Washington" = "#377eb8",
  "British Columbia" = "#377eb8",
  "Alberta" = "#a6cee3",
  "Manitoba" = "#b2df8a",
  "New Mexico" = "#33a02c",
  "Arizona" = "#33a02c",
  "Colorado" = "#33a02c",
  "Utah" = "#33a02c",
  "Quebec" = "#fdbf6f",
  "Alaska" = "#fb9a99"
)

pop_group <- factor(
  subgroups_ngsadmix,
  levels = pop_levels
)

ord <- order(pop_group)

######
#K=2
#####

#Load NGSadmix data 

ngsADMIX_files <- list.files(
  path = "/media/ssd/Bioinformatics/p_dorsalis_07/downstream_analyses/NGSadmix/all_samples/",
  pattern = "qopt",
  full.names = TRUE
)

k2_files <- ngsADMIX_files[grepl("k2", ngsADMIX_files)]

data_k2 <- read.table(k2_files[1])

#reorder both

data_k2 <- data_k2[ord, ]
pop_group <- pop_group[ord]

#plot
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
  col = pop_colors_ngsadmix[as.character(pop_group)],
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

#####
#K=3
#####

#Load NGSadmix data 

ngsADMIX_files <- list.files(
  path = "/media/ssd/Bioinformatics/p_dorsalis_07/downstream_analyses/NGSadmix/all_samples/",
  pattern = "qopt",
  full.names = TRUE
)

k3_files <- ngsADMIX_files[grepl("k3", ngsADMIX_files)]

data_k3 <- read.table(k3_files[1])


#reorder both

data_k3 <- data_k3[ord, ]
pop_group <- pop_group[ord]

###plot
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
  col = pop_colors_ngsadmix[as.character(pop_group)],
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
  line = 1   # fix overlap
)

###############################################
#Morphological subspecies population assignment
###############################################

subgroups_ngsadmix_subspecies <- p_dorsalis_sf$State_Province[
  match(samples, p_dorsalis_sf$strpID)
]

# Order from west to east
pop_levels <- c(
  "New Mexico", "Arizona", "Colorado", "Utah",
  "Oregon", "Washington", "British Columbia", "Alaska",
  "Alberta", "Manitoba", "Quebec"
)

pop_colors_ngsadmix_subspecies <- c(
  "Oregon" = "#377eb8",
  "Washington" = "#377eb8",
  "British Columbia" = "#377eb8",
  "Alberta" = "#377eb8",
  "Manitoba" = "#377eb8",
  "New Mexico" = "#fb9a99",
  "Arizona" = "#fb9a99",
  "Colorado" = "#fb9a99",
  "Utah" = "#fb9a99",
  "Quebec" = "#b2df8a",
  "Alaska" = "#377eb8"
)

pop_group <- factor(
  subgroups_ngsadmix_subspecies,
  levels = pop_levels
)

ord <- order(pop_group)

######
#K=2
#####

#Load NGSadmix data 

ngsADMIX_files <- list.files(
  path = "/media/ssd/Bioinformatics/p_dorsalis_07/downstream_analyses/NGSadmix/all_samples/",
  pattern = "qopt",
  full.names = TRUE
)

k2_files <- ngsADMIX_files[grepl("k2", ngsADMIX_files)]

data_k2 <- read.table(k2_files[1])

#reorder both

data_k2 <- data_k2[ord, ]
pop_group <- pop_group[ord]

#plot
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
  col = pop_colors_ngsadmix_subspecies[as.character(pop_group)],
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

#####
#K=3
#####

ngsADMIX_files <- list.files(
  path = "/media/ssd/Bioinformatics/p_dorsalis_07/downstream_analyses/NGSadmix/all_samples/",
  pattern = "qopt",
  full.names = TRUE
)

k3_files <- ngsADMIX_files[grepl("k3", ngsADMIX_files)]

data_k3 <- read.table(k3_files[1])


#reorder both

data_k3 <- data_k3[ord, ]
pop_group <- pop_group[ord]

###plot
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
  col = pop_colors_ngsadmix_subspecies[as.character(pop_group)],
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
  line = 1   # fix overlap
)

################################
##evaluating the ngsadmix model 
################################

###############
#Evanno method
###############

#get data


data_all<-list.files("/media/ssd/Bioinformatics/p_dorsalis_07/downstream_analyses/NGSadmix/all_samples", pattern = ".log", full.names = T)

bigData_all<-lapply(1:160, FUN = function(i) readLines(data_all[i]))
foundset<-sapply(1:160, FUN= function(x) bigData_all[[x]][which(str_sub(bigData_all[[x]], 1, 1) == 'b')])



as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) )

logs<-data.frame(K = rep(1:8, each=20))

#add to it our likelihood values

logs$like<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ))

#calculate
#Evanno method can't identify the best number of K if it's K=1, so it has its limitations. 

library(dplyr)


summary_df <- logs %>%
  group_by(K) %>%
  summarise(
    mean_ll = mean(like),
    sd_ll = sd(like),
    .groups = "drop"
  ) %>%
  arrange(K)

summary_df$deltaK <- NA_real_
for (k in 2:(nrow(summary_df)-1)) {
  L_second <- summary_df$mean_ll[k + 1] -
    2 * summary_df$mean_ll[k] +
    summary_df$mean_ll[k - 1]
  
  summary_df$deltaK[k] <- abs(L_second) / summary_df$sd_ll[k]
} #gives the same result as clumpak! 



###########
#EvalAdmix
###########

mat <- as.matrix(read.table("/media/ssd/Bioinformatics/p_dorsalis_07/downstream_analyses/NGSadmix/all_samples/evaladmix_LD_pruned_ngsADMIX_k1_rep1_BBWO"))


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


