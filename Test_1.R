##################################################
######## Testing the genomic potato of Endelman ##
##################################################
setwd("G:/Meu Drive/Doutorado/Tese/BEPE/MResende")
### testing to git
### Using the GWASpoly
library(knitr)

install.packages("devtools")
devtools::install_github("jendelman/GWASpoly", build_vignettes=FALSE)
library(GWASpoly)

#obtaining the data set from GWASpoly package
genofile <- system.file("extdata", "new_potato_geno.csv", package = "GWASpoly")
phenofile <- system.file("extdata", "new_potato_pheno.csv", package = "GWASpoly")

library(GWASpoly)
data <- read.GWASpoly(ploidy=4, pheno.file=phenofile, geno.file=genofile,
                      format="numeric", n.traits=1, delim=",")
##

datageno <- data@geno #just the genomic data
datageno <- as.data.frame(datageno)

## marker curation
N <- 957
params <- set.params(geno.freq = 1 - 5/N, fixed = "env", fixed.type = "factor", MAF = 0.1)
data.loco <- set.K(data,LOCO=TRUE,n.core=2)
data.loco.scan <- GWASpoly(data=data.loco,models=c("additive","1-dom"),
                           traits=c("vine.maturity"),params=params,n.core=2)
library(ggplot2)
qq.plot(data.loco.scan,trait="vine.maturity") + ggtitle(label="Original")


### trying the loop for the classes in genomic data

molecular_data <- data.frame(
  Marker1 = c(1, 0, NA, 1, 0),
  Marker2 = c(0, NA, 1, 1, 0),
  Marker3 = c(1, 1, 1, 0, NA),
  Marker4 = c(0, 1, NA, NA, NA)
)

print(molecular_data)
# Calculate the number of missing values (NAs) for each marker
na_counts <- colSums(is.na(molecular_data))

# Create a histogram of NA counts
hist(na_counts, 
     breaks = seq(min(na_counts), max(na_counts)),
     xlab = "Number of NAs",
     ylab = "Frequency",
     main = "Distribution of NAs across Markers")

###

install.packages("AGHmatrix")
library(AGHmatrix)


#Because of the dimension of data, the best way could be plotting the G matrix 

#Computing the additive relationship matrix based on VanRaden 2008
# adapted by Ashraf 2016
G_VanRaden <- Gmatrix(datageno, method="VanRaden", ploidy=4)

heatmap(G_VanRaden)

###### Testing the Histogram plot
levels <- unique(df$category)

# Create a list to store the histograms
histograms <- list()

# Loop over each level and create histograms
for (level in levels) {
  # Subset data for the current level
  subset_data <- df[df$category == level, ]
  
  # Create histogram for the subset
  histograms[[level]] <- hist(subset_data$values, main = paste("Histogram for", level))
}

# Plot all histograms
par(mfrow = c(2, ceil(length(levels)/2)))  # Adjust layout based on the number of levels
for (i in seq_along(levels)) {
  plot(histograms[[i]])
}
