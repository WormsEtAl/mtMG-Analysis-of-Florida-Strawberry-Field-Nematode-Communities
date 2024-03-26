# Required Libraries
library(zCompositions)
library(coda.base)
library(pheatmap)
library(compositions)
library(ape)
library(vegan)
library(corrplot)
library(fossil)
library(lmodel2)
library(Hmisc)
library(psych)

## BLAST Stats Summary ##
# Data Import
dat<-read.csv("annotatedcox1_blaststats.csv", header = T)
sampleids<-unique(dat$SampleID)

# Loop by Sample to Summarize Characteristics
for ( i in 1:length(sampleids)) {
sampleid<-sampleids[i]
sampledat<-subset(dat, dat$SampleID == sampleid)
print(mean(sampledat$LENGTH, na.rm = TRUE)) # Change column name for desired variable
}

# Heatmaps
dat<-read.csv("annotatedcox1_results.csv", header = T, row.names =1 )
rownames(dat) <- factor(rownames(dat), levels = rownames(dat))

pdf("cox1morpho_pp_compare.pdf")
pheat.m1<-pheatmap(dat, cluster_rows = FALSE, cluster_cols = FALSE, fontsize = 5, color = c("white","black"), cellwidth = 5, cellheight = 5)
dev.off()




## Community Analysis ##
#Importing Data
#Define community table and map file
# Tell R where to find you metadata file
map_name<-"compare_map.csv"
# Tell R where to find your community data file
table_name<-"treefilteredv2_annotatedCOX1_trophictable.csv"

#Import community table
otutablewt <- read.csv(file=table_name, comment.char="", header=T, row.names=1, stringsAsFactors=T, check.names=FALSE)
dim(otutablewt)
# Sort the columns by alphanumeric sorting
otutable<-otutablewt[,order(colnames(otutablewt))]


#Import map file
unorder.map<-read.csv(file=map_name, header=T, comment.char="", , row.names=1, stringsAsFactors=T, check.names=FALSE)
# Sort the rows by alphanumeric sorting
map<-unorder.map[order(rownames(unorder.map)),]
# subset map for only samples found in the community table
sub.map<-subset(map, rownames(map) %in% colnames(otutable))

# Subset community table
sub.otutable<-subset(otutable, select = colnames(otutable) %in% rownames(sub.map))

# Confirm symmetry of metadata and community tables
rownames(sub.map) == colnames(sub.otutable)
otutable<-sub.otutable


#### Betadiv calculations ####
### Jaccard distances ###
# Filter empty rows and columns
nonemptysamples<-which(colSums(otutable) != 0)
nonemptyspecies<-which(rowSums(otutable) != 0)
nonempty.otutable<-otutable[nonemptyspecies,nonemptysamples]
nonempty.map<-subset(sub.map, rownames(sub.map) %in% colnames(nonempty.otutable))
jac.dist<-vegdist(t(nonempty.otutable), method="jaccard", na.rm = TRUE)

# Plots

#### Plot dissmilartiy matrix as PCoA DPF ####
###		####	Jaccard 	###		###
## Define variables for plot code
plot.dm<-as.matrix(jac.dist) # Change the source of the plot data here
plot.map<-nonempty.map

# Make plot
pdf("annotatedcox1_ppcomm.pdf") # Name of file output, CHANGE ME as needed
# Plotting 
# FUNCTION:
    Plot_ConvexHull<-function(xcoord, ycoord, lcolor){
      hpts <- chull(x = xcoord, y = ycoord)
      hpts <- c(hpts, hpts[1])
      lines(xcoord[hpts], ycoord[hpts], col = lcolor)
    } 
    # END OF FUNCTION
# Now run the Principle Coordinate Analysis ("pcoa") in package ape:
    pcoa.pts<-pcoa(plot.dm,correction="cailliez")
    x.minimum<-min(pcoa.pts$vectors[,1])
    x.maximum<-max(pcoa.pts$vectors[,1])
    y.minimum<-min(pcoa.pts$vectors[,2])
    y.maximum<-max(pcoa.pts$vectors[,2])

    # And to be fancy, b/c Jack put the time into figuring this out,
    # calculate the % variation explained by each axis
    # getting percent explained by each PC vector
    vars_percent <- (pcoa.pts$values$Eigenvalues / sum(pcoa.pts$values$Eigenvalues)) * 100
# Step 1: make very simple plot, no colors
par(mar = c(5,5,5,1)) 
   
    plot(pcoa.pts$vectors[,2] ~ pcoa.pts$vectors[,1],
         xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
         ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
         cex=0.1, cex.lab=1.5, cex.axis=1.5, col="white"
         ,main="C. Annotated COX1 Plant Parasites", cex.main=2 # CHANGE this to be the title of the plot
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    # Step 2: now subset vectors for the positions in the alpha file that correspond to cohcryo
    map.F1D1 <- which(plot.map$Field == "Field1" & plot.map$Depth == "Depth1") ## CHANGE these depending on the number of categories in your variable of interest
    map.F1D2 <- which(plot.map$Field == "Field1" & plot.map$Depth == "Depth2") ## CHANGE these depending on the number of categories in your variable of interest
	map.F1D3 <- which(plot.map$Field == "Field1" & plot.map$Depth == "Depth3") ## CHANGE these depending on the number of categories in your variable of interest
    map.F1D4 <- which(plot.map$Field == "Field1" & plot.map$Depth == "Depth4") ## CHANGE these depending on the number of categories in your variable of interest

    map.F2D1 <- which(plot.map$Field == "Field2" & plot.map$Depth == "Depth1") ## CHANGE these depending on the number of categories in your variable of interest
    map.F2D2 <- which(plot.map$Field == "Field2" & plot.map$Depth == "Depth2") ## CHANGE these depending on the number of categories in your variable of interest
	map.F2D3 <- which(plot.map$Field == "Field2" & plot.map$Depth == "Depth3") ## CHANGE these depending on the number of categories in your variable of interest
    map.F2D4 <- which(plot.map$Field == "Field2" & plot.map$Depth == "Depth4") ## CHANGE these depending on the number of categories in your variable of interest

    map.F3D1 <- which(plot.map$Field == "Field3" & plot.map$Depth == "Depth1") ## CHANGE these depending on the number of categories in your variable of interest
    map.F3D2 <- which(plot.map$Field == "Field3" & plot.map$Depth == "Depth2") ## CHANGE these depending on the number of categories in your variable of interest
	map.F3D3 <- which(plot.map$Field == "Field3" & plot.map$Depth == "Depth3") ## CHANGE these depending on the number of categories in your variable of interest
    map.F3D4 <- which(plot.map$Field == "Field3" & plot.map$Depth == "Depth4") ## CHANGE these depending on the number of categories in your variable of interest

    map.F4D1 <- which(plot.map$Field == "Field4" & plot.map$Depth == "Depth1") ## CHANGE these depending on the number of categories in your variable of interest
    map.F4D2 <- which(plot.map$Field == "Field4" & plot.map$Depth == "Depth2") ## CHANGE these depending on the number of categories in your variable of interest
	map.F4D3 <- which(plot.map$Field == "Field4" & plot.map$Depth == "Depth3") ## CHANGE these depending on the number of categories in your variable of interest
    map.F4D4 <- which(plot.map$Field == "Field4" & plot.map$Depth == "Depth4") ## CHANGE these depending on the number of categories in your variable of interest

    # Step 3: and plot the points for cohcryo samples in red over the top using points()
    # Add a set of points and Plot_ConvexHull functions for each category in your variable
    # If may have to remove or add category specific plots as needed
   
   # Field 1
    points(pcoa.pts$vectors[map.F1D1,2] ~ pcoa.pts$vectors[map.F1D1,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=15, col="black", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F1D1,1], ycoord = pcoa.pts$vectors[map.F1D1,2], lcolor = "black")

    points(pcoa.pts$vectors[map.F1D2,2] ~ pcoa.pts$vectors[map.F1D2,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=15, col="#646060", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F1D2,1], ycoord = pcoa.pts$vectors[map.F1D2,2], lcolor = "#646060")

    points(pcoa.pts$vectors[map.F1D3,2] ~ pcoa.pts$vectors[map.F1D3,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=15, col="#817F7F", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F1D3,1], ycoord = pcoa.pts$vectors[map.F1D3,2], lcolor = "#817F7F")

    points(pcoa.pts$vectors[map.F1D4,2] ~ pcoa.pts$vectors[map.F1D4,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=15, col="#B5B0B0", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F1D4,1], ycoord = pcoa.pts$vectors[map.F1D4,2], lcolor = "#B5B0B0")

	#Field 2
	points(pcoa.pts$vectors[map.F2D1,2] ~ pcoa.pts$vectors[map.F2D1,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=16, col="#F20808", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F2D1,1], ycoord = pcoa.pts$vectors[map.F2D1,2], lcolor = "#F20808")

    points(pcoa.pts$vectors[map.F2D2,2] ~ pcoa.pts$vectors[map.F2D2,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=16, col="#EE5353", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F2D2,1], ycoord = pcoa.pts$vectors[map.F2D2,2], lcolor = "#EE5353")

    points(pcoa.pts$vectors[map.F2D3,2] ~ pcoa.pts$vectors[map.F2D3,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=16, col="#F58B8B", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F2D3,1], ycoord = pcoa.pts$vectors[map.F2D3,2], lcolor = "#F58B8B")

    points(pcoa.pts$vectors[map.F2D4,2] ~ pcoa.pts$vectors[map.F2D4,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=16, col="#F9C3C3", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F2D4,1], ycoord = pcoa.pts$vectors[map.F2D4,2], lcolor = "#F9C3C3")

	# Field 3
	points(pcoa.pts$vectors[map.F3D1,2] ~ pcoa.pts$vectors[map.F3D1,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=17, col="#065F09", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F3D1,1], ycoord = pcoa.pts$vectors[map.F3D1,2], lcolor = "#065F09")

    points(pcoa.pts$vectors[map.F3D2,2] ~ pcoa.pts$vectors[map.F3D2,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=17, col="#358438", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F3D2,1], ycoord = pcoa.pts$vectors[map.F3D2,2], lcolor = "#358438")

    points(pcoa.pts$vectors[map.F3D3,2] ~ pcoa.pts$vectors[map.F3D3,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=17, col="#6CAC6E", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F3D3,1], ycoord = pcoa.pts$vectors[map.F3D3,2], lcolor = "#6CAC6E")

    points(pcoa.pts$vectors[map.F3D4,2] ~ pcoa.pts$vectors[map.F3D4,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=17, col="#CEF6D0", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F3D4,1], ycoord = pcoa.pts$vectors[map.F3D4,2], lcolor = "#CEF6D0")

	# Field 4
    points(pcoa.pts$vectors[map.F4D1,2] ~ pcoa.pts$vectors[map.F4D1,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=18, col="#07185A", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F4D1,1], ycoord = pcoa.pts$vectors[map.F4D1,2], lcolor = "#07185A")

    points(pcoa.pts$vectors[map.F4D2,2] ~ pcoa.pts$vectors[map.F4D2,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=18, col="#31407C", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F4D2,1], ycoord = pcoa.pts$vectors[map.F4D2,2], lcolor = "#31407C")

    points(pcoa.pts$vectors[map.F4D3,2] ~ pcoa.pts$vectors[map.F4D3,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=18, col="#7482BD", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F4D3,1], ycoord = pcoa.pts$vectors[map.F4D3,2], lcolor = "#7482BD")

    points(pcoa.pts$vectors[map.F4D4,2] ~ pcoa.pts$vectors[map.F4D4,1],
           xlab=paste("PCoA1: ", round(vars_percent[1], 2), "%", sep=""),
           ylab=paste("PCoA2: ", round(vars_percent[2], 2), "%", sep=""),
           pch=18, col="#C4CFF9", cex.lab=1.5, cex.axis=1.5
         ,xlim=c(x.minimum-0.05, x.maximum+0.05), ylim=c(y.minimum-0.05, y.maximum+0.05)
)
    Plot_ConvexHull(xcoord = pcoa.pts$vectors[map.F4D4,1], ycoord = pcoa.pts$vectors[map.F4D4,2], lcolor = "#C4CFF9")




	legend(x.maximum-0.4,y.maximum, legend=c("Field1: 0-25cm","Field1: 25-50cm","Field1: 50-75cm", "Field1: 75-100cm","Field2: 0-25cm","Field2: 25-50cm","Field2: 50-75cm", "Field2: 75-100cm","Field3: 0-25cm","Field3: 25-50cm","Field3: 50-75cm", "Field3: 75-100cm","Field4: 0-25cm","Field4: 25-50cm","Field4: 50-75cm", "Field4: 75-100cm"), 
		col=c("black","#646060","#817F7F","#B5B0B0","#F20808","#EE5353","#F58B8B","#F9C3C3","#065F09","#358438","#6CAC6E","#CEF6D0","#07185A","#31407C","#7482BD","#C4CFF9"),
			pch=c(15,15,15,15,16,16,16,16,17,17,17,17,18,18,18,18), cex=1.2,
			box.lwd = 0, box.lty=0) # Add labels, colors, and point shape for each category you want to plot
dev.off()





##### Alphadiv calculations #####
# Here are examples of how to get diversity metrics in R
# Shannon Diversity Index
Shannon<-diversity(t(otutable), index = "shannon")
# Simpson Diversity index
Simpson<-diversity(t(otutable), index = "simpson")
# Inverse Simpson (simpson_e) or Evenness Index
Invert.Simpson<-diversity(t(otutable), index = "invsimpson")
# Richness
Richness<-apply(t(otutable)>0,1,sum)
# Chao1
Chao1 <- apply(t(otutable), 1, chao1)
# combine into a single data table
alternative.alphadiv<-as.data.frame(cbind(Shannon,Simpson,Invert.Simpson,Richness,Chao1))
sub.map<-subset(map, rownames(map) %in% rownames(alternative.alphadiv))
rownames(alternative.alphadiv) == rownames(sub.map)

# Loop by Sample to Summarize Richness
sampleids<-rownames(alternative.alphadiv)
for ( i in 1:length(sampleids)) {
sampleid<-sampleids[i]
sampledat<-subset(alternative.alphadiv, rownames(alternative.alphadiv) == sampleid)
print(sampledat$Richness) # Change column name for desired variable
}

# Plots
map<-read.csv("compare_map.csv", header = T, row.names = 1)
map<-map[order(rownames(map)),]
map$Metavariable <- factor(map$Metavariable, levels = c("1.0 – 25 cm.Annotated","1.0 – 25 cm.Morpho","1.25 – 50 cm.Annotated","1.25 – 50 cm.Morpho","1.50 -75 cm.Annotated","1.50 -75 cm.Morpho","1.75-100 cm.Annotated","1.75-100 cm.Morpho","2.0 – 25 cm.Annotated","2.0 – 25 cm.Morpho","2.25 – 50 cm.Annotated","2.25 – 50 cm.Morpho","2.50 -75 cm.Annotated","2.50 -75 cm.Morpho","2.75-100 cm.Annotated","2.75-100 cm.Morpho","4.0 – 25 cm.Annotated","4.0 – 25 cm.Morpho","4.25 – 50 cm.Annotated","4.25 – 50 cm.Morpho","4.50 -75 cm.Annotated","4.50 -75 cm.Morpho","4.75-100 cm.Annotated","4.75-100 cm.Morpho","3.0 – 25 cm.Annotated","3.0 – 25 cm.Morpho","3.25 – 50 cm.Annotated","3.25 – 50 cm.Morpho","3.50 -75 cm.Annotated","3.50 -75 cm.Morpho","3.75-100 cm.Annotated","3.75-100 cm.Morpho"))

myColors<-c("lightgray","lightblue")
pdf("richness_trophic.pdf")
boxplot(map$Treefiltered.Richness.Trophic ~ map$Metavariable, col = myColors, las = 2, cex.axis = 0.3, cex.lab = 0.1)
dev.off()

pdf("richness_plantparasites.pdf")
boxplot(map$Treefiltered.Richness.PP ~ map$Metavariable, las = 2, col = myColors, cex.axis = 0.3, cex.lab = 0.1)
dev.off()

sub.map<-subset(map, map$Method == "Annotated")
sub.map<-droplevels(sub.map)
sampleorder<-c("F111-025","F111-2550","F111-5075","F111-75100","F112-025","F112-2550","F112-5075","F112-75100","F121-025","F121-2550","F121-5075","F121-75100","F122-025","F122-2550","F122-5075","F122-75100","F211-025","F211-2550","F211-5075","F211-75100","F212-025","F212-2550","F212-5075","F212-75100","F221-025","F221-2550","F221-5075","F221-75100","F222-025","F222-2550","F222-5075","F222-75100","F411-025","F411-2550","F411-5075","F411-75100","F412-025","F412-2550","F412-5075","F412-75100","F421-025","F421-2550","F421-5075","F421-75100","F422-025","F422-2550","F422-5075","F422-75100","F311-025","F311-2550","F311-5075","F311-75100","F312-025","F312-2550","F312-5075","F312-75100","F321-025","F321-2550","F321-5075","F321-75100","F322-025","F322-2550","F322-5075","F322-75100")
sub.map2<-sub.map[sampleorder,]
sub.map2$Metavariable <- factor(sub.map2$Metavariable, levels = c("1.0 – 25 cm.Annotated","1.25 – 50 cm.Annotated","1.50 -75 cm.Annotated","1.75-100 cm.Annotated","2.0 – 25 cm.Annotated","2.25 – 50 cm.Annotated","2.50 -75 cm.Annotated","2.75-100 cm.Annotated","4.0 – 25 cm.Annotated","4.25 – 50 cm.Annotated","4.50 -75 cm.Annotated","4.75-100 cm.Annotated","3.0 – 25 cm.Annotated","3.25 – 50 cm.Annotated","3.50 -75 cm.Annotated","3.75-100 cm.Annotated"))

myColors = c("black","#646060","#817F7F","#B5B0B0","#F20808","#EE5353","#F58B8B","#F9C3C3","#065F09","#358438","#6CAC6E","#CEF6D0","#07185A","#31407C","#7482BD","#C4CFF9")

pdf("richness_annotatedcox1_v3.pdf")
boxplot(sub.map2$Richness ~ as.factor(sub.map2$Metavariable), col = myColors, las = 2, cex.axis = 0.3, cex.lab = 0.1)
dev.off()

# Alphadiv Stats
# Generalized linear model
summary(glm(alternative.alphadiv$Richness ~ Treatment, data = plot.map))


#### Metadata Analysis ####
# Metadata Correlations
sub.sub.map<-sub.map[,4:ncol(sub.map)]
pval <- psych::corr.test(sub.sub.map, adjust="none")$p
corrplot(cor(sub.sub.map, use = "complete.obs", method = "pearson"), type="lower", p.mat=pval, tl.pos="n")
corrplot(cor(sub.sub.map, use = "complete.obs", method = "pearson"), type="lower", p.mat=pval, sig.level=0.05,
add=T, tl.pos="ld", cl.pos="n", cl.cex = 0.5, pch.cex = 0.5, tl.cex = 0.5, number.cex = 0.1)

