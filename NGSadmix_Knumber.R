library(stringr)

data_all<-list.files("/media/ssd/Bioinformatics/p_dorsalis_07/downstream_analyses/NGSadmix/all_samples", pattern = ".log", full.names = T)

bigData_all<-lapply(1:160, FUN = function(i) readLines(data_all[i]))
foundset<-sapply(1:160, FUN= function(x) bigData_all[[x]][which(str_sub(bigData_all[[x]], 1, 1) == 'b')])



as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) )

logs<-data.frame(K = rep(1:8, each=20))

#add to it our likelihood values

logs$like<-as.vector(as.numeric( sub("\\D*(\\d+).*", "\\1", foundset) ))

#and now we can calculate our delta K and probability

tapply(logs$like, logs$K, FUN= function(x) mean(abs(x))/sd(abs(x)))

tapply(logs$like, logs$K, mean)


library(stringr)

files <- list.files(
  "/media/ssd/Bioinformatics/p_dorsalis_07/downstream_analyses/NGSadmix/all_samples",
  pattern = "\\.log$",
  full.names = TRUE
)

get_ll <- function(file) {
  x <- readLines(file, warn = FALSE)
  
  # adjust pattern depending on your log format
  line <- x[str_detect(x, "loglikelihood|Log likelihood|like")]
  
  # extract numeric value (handles scientific notation too)
  as.numeric(str_extract(line, "-?\\d+\\.?\\d*(e[+-]?\\d+)?"))
}

ll <- sapply(files, get_ll)

K <- rep(1:8, each = 20)

df <- data.frame(
  K = K,
  loglik = ll
)

library(dplyr)

summary_df <- df %>%
  group_by(K) %>%
  summarise(
    mean_ll = mean(loglik, na.rm = TRUE),
    sd_ll   = sd(loglik, na.rm = TRUE),
    .groups = "drop"
  )

summary_df

plot(summary_df$K, summary_df$mean_ll, type = "b",
     xlab = "K",
     ylab = "Mean log-likelihood")


##################
### Evanno method
##################

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


##################
##eval admix plot
##################

# read evalAdmix output
mat <- as.matrix(read.table("/media/ssd/Bioinformatics/downstream_analyses/EvalAdmix/evaladmix_LD_pruned_ngsADMIX_k2_rep2_BBWO"))


# optional: order individuals by population
ord <- order(pop_province_nowest)

mat_ord <- mat[ord, ord]

# colors
library(RColorBrewer)
cols <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# basic heatmap
heatmap(
  mat_ord,
  Rowv = NA,
  Colv = NA,
  scale = "none",
  col = cols
)




library(pheatmap)
library(RColorBrewer)

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



#NGS ADMOX PLOT 

### =========================
### 1. User-defined pop order
### =========================
##plotting
### custom pop leves: 
pop_levels <- c(
  "California",
  "Oregon",
  "Washington",
  "Alberta",
  "Manitoba",
  "Minnesota",
  "Michigan",
  "New York",
  "Quebec"
)

pop_levels_nowest <- c(
  "Alberta",
  "Manitoba",
  "Minnesota",
  "Michigan",
  "New York",
  "Quebec"
)

pop <- factor(pop_province, levels = pop_levels)
pop_nowest <- factor(pop_province_nowest, levels = pop_levels_nowest)

ord <- order(pop)
ord_nowest<-order(pop_nowest)

files <- list.files(
  path="/media/ssd/Bioinformatics/downstream_analyses/NGSadmix",
  pattern='qopt',
  full.names=TRUE
)

k2_files <- files[grepl("k4", files)]
data_k2 <- read.table(k2_files[1])

cols <- rainbow(4)

### population boundaries
pop_breaks <- cumsum(table(pop[ord]))

### population label positions (center of each block)
pop_mid <- pop_breaks - table(pop[ord])/2

barplot(
  t(as.matrix(data_k2))[ , ord],
  col = cols,
  names = pop[ord],
  las = 2,
  space = 0,
  ylab = "Admixture proportions",
  xlab = "Individuals",
  border = NA,
  cex.names = 0.75,
  main = "NGSAdmix K=2"
)


barplot(
  t(as.matrix(data_k2))[ , ord],
  col = cols,
  names = rep("", ncol(t(as.matrix(data_k2))[ , ord])),   # REMOVE repeated labels
  las = 2,
  space = 0,
  ylab = "Admixture proportions",
  xlab = ,
  border = NA,
  cex.names = 0.75,
  main = "NGSAdmix K=4"
)

### white separators
abline(v = pop_breaks, col = "white", lwd = 2)

### ONE label per population
axis(
  1,
  at = pop_mid,
  labels = names(table(pop[ord])),
  tick = FALSE,
  las = 2
)

###From PCAngsd

Q <- as.matrix(read.table("/media/ssd/Bioinformatics/downstream_analyses/PCAngsd/ADMIXTURE/ALL/K5_rep1.admix.5.Q"))
Q_ord <- t(Q)[, ord]
cols <- rainbow(5)

barplot(
  Q_ord,
  col = cols,
  border = NA,
  space = 0,
  ylab = "Admixture proportions",
  xlab = "Individuals",
  main = "PCAngsd Admixture K=3"
)
pop_breaks <- cumsum(table(pop[ord]))
pop_mid <- pop_breaks - table(pop[ord]) / 2

abline(v = pop_breaks, col = "white", lwd = 2)

axis(
  1,
  at = pop_mid,
  labels = names(table(pop[ord])),
  tick = FALSE,
  las = 2
)
