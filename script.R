### PREPARE ###

#set the directory
setwd("M:/R/VA")

#load the R libraries required
library(kohonen) 
library(ggplot2) 
library(rgdal) 
library(gridExtra) 
library(grid) 
library(viridis) 
library(dplyr)
library(sf)
library(maptools)
library(corrplot) 
library(psych) 
library(rgeos)
library(reshape2)


###  LOAD COLOR PALLETTE  ###
coolBlueHotRed <- function(n, alpha = 1) {rainbow(n, end=4/6, alpha=alpha)[n:1]}


###  LOAD DATA  ###
#read the census data 
read.table(file = "Data.csv", sep = ",", header=TRUE)

#load shapefile
mapa <- readOGR("SG_SIMD_2016_EDINBURGH.shp",stringsAsFactors = FALSE)

#check the data
class(Data)
class(mapa)
names(Data)
names(mapa)

#plot the shapefile
proj4string(mapa)
plot(mapa)

#change the projection
mapa <- spTransform(mapa, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# fortify the spatial polygon of mapa
Edinburgh_map <- fortify(mapa, region="DataZone")
Edinburgh_map <- merge(Edinburgh_map, Data, by.x="id", by.y="Data_Zone")


###  VARIBLES SELECTION ### 

# create different training data sets

variable <- c("Income_count",
              "Employment_count", 
              "CIF", "ALCOHOL","DEPRESS", "EMERG", "Attendance", 
              "Noquals", "NEET", "PT_GP", "crime_count", 
              "overcrowded_count")

data_test <- select(Data, Income_count, Income_rate, 
                    Employment_count,Employment_rate,
                    CIF, ALCOHOL, DRUG, SMR, EMERG, LBWT, DEPRESS,
                    Attendance, Attainment, Noquals, NEET, HESA, 
                    drive_petrol, drive_GP, drive_PO, drive_primary, drive_retail, drive_secondary, PT_GP, PT_Post, PT_retail, 
                    crime_count,crime_rate, 
                    overcrowded_count, nocentralheat_count, overcrowded_rate, nocentralheat_rate)

data_train <- select(Data,Income_count,
                     Employment_count, 
                     CIF, ALCOHOL,DEPRESS, EMERG, Attendance, 
                     Noquals, NEET, PT_GP, crime_count, 
                     overcrowded_count)

data_a <- select(Data,Data_Zone,
                 Income_count,
                 Employment_count, 
                 CIF, ALCOHOL,DEPRESS, EMERG, Attendance, 
                 Noquals, NEET, PT_GP, crime_count, 
                 overcrowded_count)

# change the data frame with test data sets to a matrix
# also scale all variables to standardize the data using z-value
data_test_matrix <- as.matrix(scale(data_test))

# keep the column names
names(data_test_matrix) <- names(data_test)

# set the size and topology of the som grid
som_grid_test <- somgrid(xdim = 11, ydim = 11, topo="hexagonal")

som_model_test <- som(data_test_matrix, 
                 grid = som_grid_test, 
                 rlen = 500, 
                 alpha = c(0.05,0.01), 
                 keep.data = TRUE)

# plot all the variables to see their distribution for collinearity check
par(mfrow = c(4,8))
for (i in 1:31) {
  plot(som_model_test, type = "property", property = getCodes(som_model_test)[,i],
       main=colnames(getCodes(som_model_test))[i],shape="straight", palette.name=coolBlueHotRed)
} 

dev.off()

#select the variables in health domain from file called 'health' 
health_a <- cor (health[1:7])

#calculate the correlation and plot the result
corr.test(health_a, use = "complete", method = "pearson", adjust = "none")
corrplot(corr=health_a,method = "color",order = "hclust",tl.col="black",addrect=4,addCoef.col = "grey")

#filter the variables and plot the result
var_a <- cor (var[1:16])
corrplot(var_a, method = "color")  


###  TRAIN THE TRAINING DATASET of SOM ###

# change the data frame with test data sets to a matrix
# also scale all variables to standardize the data using z-value
data_train <- mutate_all(data_train, as.numeric)
data_train_matrix <- as.matrix(scale(data_train))

# keep the column names
names(data_train_matrix) <- names(data_train)

# set the size and topology of the som grid
som_grid <- somgrid(xdim = 11, ydim = 11, topo="hexagonal")

som_model <- som(data_train_matrix,
                       grid=som_grid,
                       rlen=800,
                       alpha=c(0.05,0.01),
                       keep.data = TRUE )

# set the format of plot
par(mfrow = c(1,1))

# plot of the training progress 
plot(som_model, type="changes")
# plot "counts"
plot(som_model, type="count", main = "Node Counts", shape="straight",palette.name=coolBlueHotRed)
# plot "quality"
plot(som_model, type="quality", shape="straight",main = "Node Quality",palette.name=coolBlueHotRed )
# plot "neighbour distances"
plot(som_model, type="dist.neighbours", main = "SOM neighbour distances", shape="straight", palette.name=coolBlueHotRed)
# plot prevalence of one variable over the other using pie chart
plot(som_model, type="codes", shape="straight", palette.name=coolBlueHotRed)
       
# plot selected variables ina grid
par(mfrow = c(4,3))
for (i in 1:12) {
  plot(som_model, type = "property", property = getCodes(som_model)[,i],
       main=colnames(getCodes(som_model))[i], shape="straight",palette.name=coolBlueHotRed)
} 

dev.off()

###   CLUSTERING   ### 
# determining the number of cluster by showing the WCSS metric for k-means
mydata <- getCodes(som_model)
wss <- (nrow(mydata)-1)*sum(apply(mydata,2,var))

for (i in 2:15) wss[i] <- sum(kmeans(mydata,
                                     centers=i)$withinss)
par(mfrow = c(1,1))
par(mar=c(5.1,4.1,4.1,2.1))
plot(1:15, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares", main="Within cluster sum of squares (WCSS)")

# form clusters on grid
# use hierarchical clustering to cluster the codebook vectors
som_cluster <- cutree(hclust(dist(getCodes(som_model))), 7)

cbPalette <- c("#CC79A7", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#800080")

# plot divided in clusters colors
plot(som_model, type="mapping", codeRendering="segments", bgcol = cbPalette[som_cluster], main = "Mapping Clusters",shape="straight")
# add separation lines
add.cluster.boundaries(som_model, som_cluster)


# plot with pie chart instead of just colours and add separation lines
plot(som_model, type="codes", bgcol = cbPalette[som_cluster], main = "Codes Clusters", codeRendering="segments", shape="straight")
add.cluster.boundaries(som_model, som_cluster)

# create data frame of the Data's DataZone conlume and of the cluster unit
cluster_details <- data.frame(id=Data$Data_Zone, cluster=som_cluster[som_model$unit.classif])

# combine the census and cluster data
map <- merge(Edinburgh_map,cluster_details,by="id")

# finally map the areas and colour by cluster
ggplot(data=map, aes(x=long, y=lat, group=group, fill=factor(cluster))) + geom_polygon() + scale_fill_manual(values = cbPalette) +geom_path(colour="white", alpha=0.5, size=0.05) # if you want an outline


###  EXPORT SHAPE FILE OF RESULT MAP   ### 
#merge the tables to get outcome
output <- merge(mapa, Data, by.x="DataZone", by.y="Data_Zone")
output <- merge(output, cluster_details, by.x = "DataZone", by.y="id")

# Export the shapefile
writeOGR(obj=output,
         dsn="shapefile",
         layer="results",
         driver="ESRI Shapefile")


###  CLUSTER RESULT ANALYSIS  ### 
data_a_train <- as.matrix(scale(select_if(data_a, is.numeric)))
names(data_a_train) <- names(data_a)

# merge the data in different files
clustera <- merge(data_a, cluster_details, by.x="Data_Zone", by.y="id")
clusterb <- clustera %>%
            mutate_at(variable, funs(scale))

# calculate the mean of Z-value for each cluster
Means <- clusterb %>%
        group_by(cluster) %>%
        summarise_at(vars(variable), funs(mean))

# calculate the percentage of each cluster
pro <- clusterb %>%
  group_by(cluster) %>%
  summarize(proportion = n() / nrow(clusterb))

#save the results
write.csv(Means, file = "mean.csv",row.names=FALSE)



###   PLOT THE STATISTIC ANALYSIS OF EACH CLUSTER   ###

mydata <- read.csv("mean.csv")
mydata$cluster <- as.factor(mydata$cluster)
melted_data <- melt(mydata, id.vars = "cluster")

# plot the result and adjust the range and color of the legend
ggplot(melted_data, aes(x = cluster, y = variable, fill = value)) + 
  geom_tile() + 
  scale_fill_gradient(low = "white", high = "blue", breaks = seq(-3, 15, 2))+ 
  labs(x = "", y = "Cluster", fill = "Value") + 
  theme_bw()

