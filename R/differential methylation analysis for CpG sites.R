#Differential methylation analysis for CpG sites

#NOTE:This portion of the package is only compatible with MethylDackel output files, which are .bedgraph files. 
#Recommended to follow the step-by-step workflow of this R script.

#1. Read in .bedgraphs and join them together
#The input x corresponds to .bedgraph files.
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
#Note: this function will introduce NAs since not all samples have the same CpG sites. 
#NAs will be filtered out with the function filter.dataset() in the later steps.
join.bedgraph <- function(x,...){
  library(plyr)
  full=list(x,...)
  FullDat=join_all(full,type="full")
  return(FullDat)
}

#2.Calculate total reads for each CpG site and add to the full table
#The output of this function will be used as the input for the mean shift application.
total.reads <- function(x){
  tot=c()
  tot=apply(x[,-c(1:3)],1, sum, na.rm = TRUE)
  full <- add_column(x, tot, .after = 3)     
  colnames(full)[4]<-"Total Reads"
  return(full)
}

#3. Remove NAs if there are NAs for all samples in one group
#The user needs to define groups (Ex. group1 = normal, group2 = cancer)
#df is the table obtained from total.reads()
#The user has to define which columns correspond to group 1 or group 2.
#x: Columns of df corresponding to group 1 (Ex. 5:8)
#y: Columns of df corresponding to group 2 (Ex. 9:12)
#Columns 1:4 of df should be disregarded since it should be the same for both groups.
filter.dataset <- function(x,y,df){
  count_func <- function(x) {
    sum(!is.na(x))
  }
  x <- df[,x]
  y <- df[,y]
  x <- apply(x[seq(1,ncol(x),2)], 1, count_func) #apply the function on only the column "Methylated" to avoid double counting samples
  df$group1.count <- as.numeric(x)
  y <- apply(y[seq(1,ncol(y),2)], 1, count_func) 
  df$group2.count <- as.numeric(y)
  filtered <- subset(df,group1.count != 0 & group2.count != 0) #subset rows that has at least one set of values (meth + unmeth) for both groups
  filtered <- filtered[-rev(seq_len(ncol(filtered)))[1:2]]
  return(filtered)
}

#4.Calculate methylation ratio for each sample 
#The methylation ratio is derived by dividing the number of methylated reads from the total number of reads (methylated reads + unmethylated reads).
#which gives the proportion of methylated reads. 
meth.ratio <- function(df){
  out=list() 
  temp = df[,-c(1:4)]
  for(i in seq(1,ncol(temp),2)){
    ratio <- temp[,i]/(temp[,i] + temp[,i+1])
    out <- c(out,list(ratio))
  }
  out <- do.call(cbind.data.frame, out)
  for(i in 1:ncol(out)){
    colnames(out)[i] <- i
  }
  return(out)
}

#5.Calculate methylation difference between groups
#Input: x = range for group 1, y = range for group 2, df = table obtained from meth.ratio()
meth.diff <- function(x,y,df){
  library("parallel")
  numCores <- parallel::detectCores()
  cl <- parallel::makeCluster(numCores)
  
  #subset df into different groups
  g1 <- df[,x]
  g2 <- df[,y]
  #conditional statement for different sample sizes per group
  if (ncol(g1) == 1 && ncol(g2) > 1) {
    g1methyl <- g1
    g2methyl <- parallel::parApply(cl,g2,1,mean,na.rm = TRUE)
    dif=round((g2methyl-g1methyl)*100)
    return(dif)
  } else if (ncol(g1) > 1 && ncol(g2) == 1) {
    g1methyl <- parallel::parApply(cl,g1,1,mean,na.rm = TRUE)
    g2methyl <- g2
    dif=round((g2methyl-g1methyl)*100)
    return(dif)
  } else if (ncol(g1) == 1 && ncol(g2) == 1) {
    g1methyl <- g1
    g2methyl <- g2
    dif=round((g2methyl-g1methyl)*100)
    return(dif)
  } else {
    g1methyl <- parallel::parApply(cl,g1,1,mean,na.rm = TRUE)
    g2methyl <- parallel::parApply(cl,g2,1,mean,na.rm = TRUE)
    dif=round((g2methyl-g1methyl)*100)
    return(dif)
  }
}

meth.diff <- function(x,y,df){
  library("parallel")
  numCores <- parallel::detectCores()
  cl <- parallel::makeCluster(numCores)
  
  #subset df into different groups
  g1 <- df[,x]
  g2 <- df[,y]
  
  if (isTRUE(ncol(g1) == 1 && ncol(g2) > 1) == TRUE) {
    g1methyl <- g1
    g2methyl <- parallel::parApply(cl,g2,1,mean,na.rm = TRUE)
    dif=round((g2methyl-g1methyl)*100)
    return(dif)
  } else if (isTRUE(ncol(g1) > 1 && ncol(g2) == 1) == TRUE) {
    g1methyl <- parallel::parApply(cl,g1,1,mean,na.rm = TRUE)
    g2methyl <- g2
    dif=round((g2methyl-g1methyl)*100)
    return(dif)
  } else if (isTRUE(ncol(g1) == 1 && ncol(g2) == 1) == TRUE) {
    g1methyl <- g1
    g2methyl <- g2
    dif=round((g2methyl-g1methyl)*100)
    return(dif)
  } else {
    g1methyl <- parallel::parApply(cl,g1,1,mean,na.rm = TRUE)
    g2methyl <- parallel::parApply(cl,g2,1,mean,na.rm = TRUE)
    dif=round((g2methyl-g1methyl)*100)
    return(dif)
  }
}

#6. Visualizing frequency of methylation patterns
#The user can set the criteria. For example, methylation difference >= 10 can be considered as hypermethylation.
#Input: meth.diff()
#Pie chart for hyper, hypo, no change
meth.pattern <- function(x,df){
  hyper <- df[df >= x,]
  hypo <- df[df <= -x,]
  no_change <- df[df < x,]
  library(ggplot2)
  data <- data.frame(
    group=c("Hypermethylated", "Hypomethylated", "No change")
    value=c(hyper,hypo,no_change)
    ggplot(data, aes(x="", y=value, fill=group)) +
    geom_bar(stat="identity", width=1, color="white") +
    coord_polar("y", start=0) +
    theme_void() # remove background, grid, numeric labels
  )
  #pie(c(hyper,hypo,no_change))
}

#7. Fisher exact test
#8. T-test
#9. Finding singificant CpG sites