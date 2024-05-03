#########################################
## Converting in the FastStructure format
#########################################
# STEP 1

#Should be a numerical or genmat matrix
matrix <- read.csv("/blue/mresende/share/Givanildo/M_85K_matrix.csv", row.names = 1, header = T)

#Creating the function:
numeric2structure <- function(genmat, 
                              file,
                              indNames = dimnames(genmat)[[1]],
                              addtlColumns = NULL, ploidy = 4, #if diploid, change it
                              exportMarkerNames = TRUE){
  nInd <- dim(genmat)[1] # number of individuals
  if(length(indNames) != nInd){
    stop("Number of individuals does not match between indNames and genmat.")
  }
  if(!is.null(addtlColumns) && dim(addtlColumns)[1] != nInd){
    stop("Number of individuals does not match between addtlColumns and genmat.")
  }
  genmat <- as.matrix(genmat)
  if(!all(genmat %in% c(0:ploidy,NA))){
    stop("genmat must only contain 0, 1, 2... ploidy and NA")
  }
  if(length(file) != 1 || !is.character(file)){
    stop("file must be a single character string.")
  }
  if(length(ploidy) != 1 || !is.numeric(ploidy)){
    stop("ploidy must be a single number")
  }
  if(!exportMarkerNames %in% c(TRUE, FALSE)){
    stop("exportMarkerNames must be TRUE or FALSE")
  }
  
  # make sets of possible genotypes
  G <- list()
  for(i in 0:ploidy){
    G[[i + 1]] <- c(rep(1, ploidy - i), rep(2, i))
  }
  G[[ploidy + 2]] <- rep(-9, ploidy) # for missing data
  
  # set up data frame for Structure
  StructTab <- data.frame(ind = rep(indNames, each = ploidy))
  # add any additional columns
  if(!is.null(addtlColumns)){
    for(i in 1:dim(addtlColumns)[2]){
      StructTab <- data.frame(StructTab, rep(addtlColumns[,i], each = ploidy))
      if(!is.null(dimnames(addtlColumns)[[2]])){
        names(StructTab)[i + 1] <- dimnames(addtlColumns)[[2]][i]
      } else {
        names(StructTab)[i + 1] <- paste("X", i, sep = "")
      }
    }
  }
  
  # add genetic data
  for(i in 1:dim(genmat)[2]){
    thesegen <- genmat[,i] + 1
    thesegen[is.na(thesegen)] <- ploidy + 2
    StructTab[[dimnames(genmat)[[2]][i]]] <- unlist(G[thesegen])
  }
  
  # add marker name header
  if(exportMarkerNames){
    cat(paste(dimnames(genmat)[[2]], collapse = "\t"), sep = "\n", file = file)
  }
  
  # export all data
  write.table(StructTab, row.names = FALSE, col.names = FALSE, append = TRUE,
              sep = "\t", file = file, quote = FALSE)
}

## STEP 2 - Converting the genomic matrix to FastStructure/STRUCTURE format
numeric2structure(genmat = matrix, 
                  file = "/blue/mresende/share/Givanildo/Structure_files/structure_output.txt", 
                  indNames = row.names(matrix), 
                  addtlColumns = NULL, 
                  ploidy = 4, 
                  exportMarkerNames = TRUE)

#Verifying the file:
table <- read.table("/blue/mresende/share/Givanildo/Structure_files/structure_output.txt", row.names = NULL, header = T)
