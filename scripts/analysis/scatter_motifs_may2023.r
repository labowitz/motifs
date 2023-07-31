# 1. The data frame contains 412 unique cell type IDs -- I remember we had 622 (or something like that) for TGFb. Why is the number different here? Is that because this is only tabula muris?

# 2. I see some cell types from developmental datasets, is there any reason we are including all. 
# Is the idea to subset from this data frame?

# 3. The plots I see with the basic scatterplot dont match the ones you posted on the google doc. I would assume this is realted to 1 and 2? 



# Code 


df_list <- readRDS('tgfb_distance_plt_data.RDS')


df <- df_list[[1]]

# add dataset column
df %>%  mutate(dataset2 = sapply(strsplit(cell_id_2, "_"), tail, n=1), dataset1 = sapply(strsplit(cell_id_1, "_"), tail, n=1)) -> df_dist 

# filter for adult datasets
adult_datasets <- c("18m 10x" ,"1m 10x" , "3m 10x" , "21m 10x" ,"30m 10x", "24m 10x")

df_dist %>% dplyr::filter(dataset1 %in% adult_datasets & dataset2 %in% adult_datasets) -> adult_dist 



# paired histograms
adult_dist %>% dplyr::filter(global_dist >40 & global_dist <90 ) %>%  ggplot(aes(x = pathway_dist)) + geom_density() + geom_density(aes(x = random_dist), color = 'red')

# non disperse
adult_dist %>% dplyr::filter(global_dist <40 ) %>%  ggplot(aes(x = pathway_dist)) + geom_density() + geom_density(aes(x = random_dist), color = 'red')



# Make it tidy 
adult_dist %>% pivot_longer(cols = c("pathway_dist", "random_dist"), names_to = "type", values_to = "distance") -> adult_tidy

adult_tidy %>% dplyr::filter(global_dist >40 & global_dist <90 ) %>% ggplot(aes(x = distance, color = type, fill = type)) + geom_histogram(position = 'dodge') 





# TEST: bins 
# Define the ranges and labels for the bins
breaks <- c(0,30, 60, 90, Inf)
labels <- c("low", "medium", "high", "very high")
# Create a new column called "bins" that bins the values of "values"
adult_tidy$bins <- cut(adult_tidy$global_dist, breaks = breaks, labels = labels, right = FALSE)



# TEST: 2-D densities 
# 1. estimate 2-d prob density with MASS
dens1 <- MASS::kde2d(adult_tidy$global_dist[adult_tidy$type == "pathway_dist"], adult_tidy$distance[adult_tidy$type == "pathway_dist"], n = 100)
dens2 <- MASS::kde2d(adult_tidy$global_dist[adult_tidy$type == "random_dist"], adult_tidy$distance[adult_tidy$type == "random_dist"], n = 100)


# 2. plot individual distributions 
# create a tidy data frame from the density matrices
dens_df <- expand.grid(x = dens1$x, y = dens1$y)
dens_df$z <- as.vector(dens1$z)

# create heatmap
ggplot(dens_df, aes(x = x, y = y, fill = z)) +
    geom_tile() +    
    theme_minimal()


# 3. Plot the difference 

# Compute the difference in the 2D densities between the two distributions
diff_dens <- dens1$z - dens2$z

# create a tidy data frame from the density matrices
dens_df <- expand.grid(x = dens1$x, y = dens1$y) # same axis for both
dens_df$z <- as.vector(diff_dens)

# create heatmap
ggplot(dens_df, aes(x = x, y = y, fill = z)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = median(dens_df$z)) +
    theme_minimal()
