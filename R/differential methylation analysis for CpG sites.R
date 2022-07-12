#Differential methylation analysis for CpG sites

#NOTE:This package is only compatible with MethylDackel output files, which are .bedgraph files. 

#1. Read in .bedgraphs and join them together
#This function requires "name", which is the sample name chosen under the user's discretion.

read.bedgraph <- function(x,name){
  library("tibble")
  meth.table <- read.table(x, header = TRUE)
  meth.table <- meth.table [,c(1,2,3,5,6)]
  colnames(meth.table) <- c("Chromosome","Start Position","End Position")
  colnames(meth.table)[4] <- paste(name,"Methylated")
  colnames(meth.table)[5] <- paste(name,"Nonmethylated")
  return(meth.table)
}

#Theoretically, each sample has its own .bedgraph file. 
#Thus, this function joins the methylation results for all the samples into one table.
join.bedgraph <- function(x,...){
  library(plyr)
  full=list(...)
  FullDat=join_all(full,type="full")
  return(FullDat)
}

#2.Calculate total reads for each CpG site and add to the full table
#(the output of this function will be used as the input for the mean shift application)
total.reads <- function(x){
  tot=c()
  tot=apply(x[,-c(1:3)],1, sum, na.rm = TRUE)
  full <- add_column(x, tot, .after = 3)     
  colnames(full)[4]<-"Total Reads"
  return(full)
}

#3.Remove NAs if there are NAs for all samples in one group
#The user needs to define groups (Ex. group1 = normal, group2 = cancer)
#df is the table obtained from total.reads()
#x: methylated columns for group 1
#y: methylated columns for group 2
filter.dataset <- function(x,y,df){
  count_func <- function(x) sum(!is.na(x)) 
  df$group1.count <<- apply(x, 1, count_func) 
  df$group2.count <<- apply(y, 1, count_func) 
  df$group1.count <<- as.numeric(z$group1.count)
  df$group2.count <<- as.numeric(z$group2.count)
  filtered <- subset(df,group1.count != 0 & group2.count != 0) #subset rows that has at least one set of values (meth + unmeth) for both groups
  return(filtered)
}

#4.Calculate methylation ratio for each sample 
meth.ratio <- function(x){
  for(i in 1:length(x)){
    ratio <- x[4][n+2]/
  }
}

#5.Calculate methylation difference between groups
meth.difference <- 