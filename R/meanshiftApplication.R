############################################################
                # MEANSHIFT APPLICATION #
############################################################

# made by Casey Pei, June 2022

# adapted from algorithm to be run on FullData.csv to combine islands by distance between position via the meanshift algorithm
# relies on the methylation package made by Jiang Yuan in Zhang Lab > jiangyuan > Shared Code > methylttest
# only works on linux or mac machines (for the parallel package to work)


# INPUT
#
# file - file (path of csv file) - EDIT ON LINE 116
# [Chromosome, Start Position, End Position, X1235 Methylated, X1235 Unmethylated, ...]
# expects the columns, "Chromosome," "Start Position", "End Position"
#
# setpsize - size of window for the meanshift algorithm (int) - EDIT ON LINE 133
# default value of 50
#

# OUTPUT
#
# file - meanshift_result.csv (file in directory)
# [Chromosome, Start Position, End Position, Sets, X1235 Methylated, X1235 Unmethylated, ...]
# Sets is the number of sets combined into each row
#


# Rcpp::sourceCpp("src/mean_shift0324.cpp")


############################################################
                      # FUNCTIONS #
############################################################

######## sum the columns ########
computeplus <- function(x){
  x <- na.omit(x)
  if (length(x)==0) {return(NA)}
  else {return(sum(x))}
}

######## get the start of the column ########
getstart <- function(x){
  return(x[1])
}

######## get the end of the column ########
getend <- function(x){
  return(x[length(x)])
}


#' Apply Meanshift
#'
#' This function generates a list of lists of islands (as determined by a meanshift
#' algorithm) for each chromosome from a data.frame or data.table of CpG sites. Expects 
#' `Chromosome`, Start Position`, and `End Position` columns.
#'
#' @param x -------- data.frame / data.table
#'                   format: [Chromosome, Start Position, End Position, X1235 Methylated, X1235 Unmethylated, ...]
#' @param stepsize - integer (window size that defines range single sets
#'                   can differ from the island center to be in the island)
#' @return a list of lists of indices that make up islands for each island
#' @export
applyMeanshift <- function(x, stepsize=50){
  tmppos <- x$`Start Position`

  # Rcpp::sourceCpp("src/mean_shift0324.cpp")
  outpos = mean_shift(tmppos, stepsize) #.Call("src/mean_shift0324",tmppos,package="RcppArmadillo")
}


#' Get Mode
#'
#' This function takes a sequence of numerical values and returns the statistical mode
#'
#' @param x - list of numeric values
#' @return the statistical mode of the numerical values
#' @export
getMode <- function(x) {
  u <- unique(x)
  tab <- tabulate(match(x, u))
  u[tab == max(tab)]
}


#' Format Islands
#'
#' This function takes in a list with data.frames of the original CpG sites and 
#' a list of the island indices, separated by chromosome, formatting it into a
#' single final dataframe with the island indices combined into one row.
#' Expects original CpG site data.table to have `Chromosome`, `Start Position`,
#' `End Position`, and methylation read columns -- but no summarized data columns.
#'
#' @param x - list of dataframe for each chromosome, and list of island indices
#'            columns: [Chromosome1, [20, 21]]
#' @return a list of data.tables for each chromosome with combined islands and 
#'         unchanged single sets
#' @export
formatIslands <- function(x) {
  chrom <- x[[1]]
  islands <- x[[2]]

  islands_unlisted <- unlist(islands)
  if (length(islands_unlisted) != 0) {
    finalobj <- chrom[-islands_unlisted, ] # all rows but those in islands
    finalobj <- finalobj %>%
                mutate(Sites=1, .after=`End Position`)

    summed_rows <- lapply(islands, function(x) c(chrom[getstart(x),"Chromosome"], chrom[getstart(x), "Start Position"], chrom[getend(x), "End Position"], length(x), dplyr::summarise_all(chrom[getstart(x):getend(x),4:ncol(chrom)],computeplus))) # list of summed rows
    summed_rows <- lapply(summed_rows, function(x) setNames(x, colnames(finalobj))) # change it so the column names match
    summed_rows <- rbindlist(summed_rows, use.names=TRUE, fill=TRUE, idcol=NULL) # put the summed_rows into DF

    finalobj <- rbind(finalobj, summed_rows) # add to final DF
    finalobj <- finalobj[order(finalobj$`Start Position`),] # put them in the correct order

    finalobj <- data.table(finalobj)
  }
  else {
    finalobj <- chrom
  }

  return(finalobj)
}


#' Get Islands
#'
#' This function takes in a list with data.frames of the original CpG sites and 
#' a list of the island indices, separated by chromosome, formatting it into a
#' single final dataframe with the island indices combined into one row.
#' Expects original CpG site data.table to have `Chromosome`, `Start Position`,
#' `End Position`, and methylation read columns -- but no summarized data columns.
#'
#' @param file ----- string that is pathfile to a table of CpG sites
#'                   format: [Chromosome, Start Position, End Position, X1235 Methylated, X1235 Unmethylated, ...]
#' @param stepsize - integer of window size for meanshift algorithm
#'                   default value = 50
#' @return a data.table which is the original data.table but with the islands (as
#'         determined by the meanshift algorithm using the given stepsize) added into one row
#' @export
getIslands <- function(file, stepsize=50) {

  message("\nload data")
  data <- data.table::fread(file)
  data <- data[with(data, order(Chromosome, `Start Position`)), ]

  numCores <- parallel::detectCores()

  message("\nseparate the data by chromosome")
  ori.list <- list()
  for(i in c(1:22,"X","Y")){
    ori.list[[i]] <- data[which(data$Chromosome==paste(i)),]
  }

  message("\nget islands")
  # list of indices in islands
  island.list <- parallel::mclapply(ori.list, applyMeanshift, stepsize, mc.cores = numCores)


  # list of dataframe for each chromosome and indices in the island
  format.list <- list()
  for(i in c(1:22,"X","Y")){

    format.list[[i]] <- list()
    format.list[[i]][[1]]  <-  ori.list[[i]]
    format.list[[i]][[2]] <- island.list[[i]]
  }

  message("\nstart formatting")

  out.list <- parallel::mclapply(format.list, formatIslands, mc.cores = numCores)

  outframe <- data.frame(matrix(0, nrow = 0, ncol = ncol(out.list[[1]]), dimnames = list( NULL, paste0(colnames(out.list[[1]])) ) ) )
  class(outframe$Chromosome) = "character"
  result <- data.table::setDT(outframe)
  result <- rbindlist(out.list, use.names=TRUE, fill=TRUE, idcol=NULL)

  return(result)
}

############################################################
                   # Statistics #
############################################################

# message("\n-----#Site for Row (Single Sites and Islands)----")
# nrow(result)
# mean(result$Sites)
# median(result$Sites)
# getMode(result$Sites)
# max(result$Sites)
# min(result$Sites)
# hist(result$Sites, main=paste("#Sites per Row for step=", stepsize))

# message("\n-----Size for Single Sites and Islands----")
# nrow(result)
# mean(result$`End Position` - result$`Start Position`)
# median(result$`End Position` - result$`Start Position`)
# find_mode(result$`End Position` - result$`Start Position`)
# max(result$`End Position` - result$`Start Position`)
# min(result$`End Position` - result$`Start Position`)
# hist(result$`End Position` - result$`Start Position`, main=paste("Size ('End Position' - 'Start Position') per Row for step=", stepsize))

# result <- result[result$Sites > 1, ]

# message("\n-----#Site for Island----")
# nrow(result)
# mean(result$Sites)
# median(result$Sites)
# find_mode(result$Sites)
# max(result$Sites)
# min(result$Sites)
# hist(result$Sites, main=paste("#Sites per Island for step=", stepsize))

# message("\n-----Size for Island----")
# nrow(result)
# mean(result$`End Position` - result$`Start Position`)
# median(result$`End Position` - result$`Start Position`)
# find_mode(result$`End Position` - result$`Start Position`)
# max(result$`End Position` - result$`Start Position`)
# min(result$`End Position` - result$`Start Position`)
# hist(result$`End Position` - result$`Start Position`, main=paste("Size ('End Position' - 'Start Position') per Island for step=", stepsize))