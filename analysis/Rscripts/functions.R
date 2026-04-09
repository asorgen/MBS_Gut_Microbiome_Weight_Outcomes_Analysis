
##### Function to parse out concatenated taxa names from QIIME2 #####
parseQIIME2Taxa<-function(table){
  
  bugNameSplit=gsub(pattern =  "d__", x = table$bugName,replacement = "")
  bugNameSplit=gsub(pattern =  "[.].__", x = bugNameSplit,replacement = "_/_")
  other <- paste(replicate(count,"Other"), collapse = "_/_")
  bugNameSplit=gsub(pattern ="^Other$",x = bugNameSplit,replacement = other)
  string <- strsplit(as.character(bugNameSplit),split = "_/_")
  
  temp_string=do.call(rbind,string)
  table <- cbind(temp_string,table)
  
  colnames(table)[colnames(table)=="1"] <- "Domain"
  colnames(table)[colnames(table)=="2"] <- "Phylum"
  colnames(table)[colnames(table)=="3"] <- "Class"
  colnames(table)[colnames(table)=="4"] <- "Order"
  colnames(table)[colnames(table)=="5"] <- "Family"
  colnames(table)[colnames(table)=="6"] <- "Genus"
  
  return(table)
  
  }

##### Function to round p values #####
roundP <- function(pvalue) {
  
  p <- ifelse(pvalue == 0, "p < 2e-16", 
              ifelse(pvalue < 0.001, "p < 0.001",
                     ifelse(pvalue >= 0.001, paste0("p = ", round(pvalue, digits = 3)),
                            "NS")))
  
  

}

##### Function to assin significance (*) to p values #####
sigStars <- function(pvalue) {
  
  pStar <- ifelse(pvalue < 0.001, "***",
                  ifelse(pvalue < 0.01, "**",
                         ifelse(pvalue < 0.05, "*",
                                "NS")))
  

}


##### Function to calculate the relative abundance of a taxa table #####
getRelAbun <- function(table) {
  
  relabun <- table
  
  for (x in 1:nrow(table)){
    relabun[x,1:ncol(table)] = ((table[x,1:ncol(table)]) / (rowSums(table[,1:ncol(table)])[x]))
  }
  
  return(relabun)
  
}


##### Function to assign participant weight quintile groups #####
getQuintileGroup <- function( val, quintile )
{
  for( i in 1:4)
  {
    if( val >= as.numeric(quintile[i]) & val <= as.numeric(quintile[i+1]) ) 
      return (i) 
  }
  
  return (5)
  
}

##### Function to assign participant weight quartile groups #####
getQuartileGroup <- function( val, quartile )
{
  for( i in 1:3)
  {
    if( val >= as.numeric(quartile[i]) & val <= as.numeric(quartile[i+1]) ) 
      return (i) 
  }
  
  return (4)
  
}

##### Function to assign participant weight half groups #####
getHalfGroup <- function( val, half )
{
  for( i in 1:1)
  {
    if( val >= as.numeric(half[i]) & val <= as.numeric(half[i+1]) ) 
      return (i) 
  }
  
  return (2)
  
}

##### Function to assign participant weight third groups #####
getTertileGroup <- function( val, third )
{
  for( i in 1:2)
  {
    if( val >= as.numeric(third[i]) & val <= as.numeric(third[i+1]) ) 
      return (i) 
  }
  
  return (3)
  
}

##### Function to remove outliers from a dataset #####
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

##### Function to correct the sign of p-values based on coefficient sign #####
getlog10p<-function(pval,coefficient){
  
  log10p<-sapply(1:length(pval), function(x){
    
    if (coefficient[x]>0) return(-log10(pval[x]))
    else return(log10(pval[x]))
  })
}

##### Function to give the log10 p-values from studies or time points that we want to compare #####
compareStudies<-function(level, file1, study1Term, file2, study2Term){
  
  myT1<-read.table(file1, sep="\t",header = TRUE,check.names = FALSE,quote = "",comment.char = "")
  myT2<-read.table(file2, sep="\t",header = TRUE,check.names = FALSE,quote = "",comment.char = "")
  
  common<-intersect(myT1[,level],myT2[,level]) 
  myT1c<-myT1[myT1[,level] %in% common,]
  myT2c<-myT2[myT2[,level] %in% common,]
  myT2c<-myT2c[match(myT1c[,level],myT2c[,level]),] # reorder taxa so in the same order as myT1c
  
  p1Col <- paste0("p_", study1Term)
  p2Col <- paste0("p_", study2Term)
  naColumns<-c(which(is.na(myT1c[,p1Col])),which(is.na(myT2c[,p2Col])))
  if(length(naColumns)>0){
    myT1c<-myT1c[-naColumns,]
    myT2c<-myT2c[-naColumns,]
  }
  
  p_1<-paste0("p_", study1Term)
  p_2<-paste0("p_", study2Term)
  s_1<-paste0("s_", study1Term)
  s_2<-paste0("s_", study2Term)
  
  pval1<-getlog10p(myT1c[,p_1],myT1c[,s_1])
  pval2<-getlog10p(myT2c[,p_2],myT2c[,s_2])
  
  df<-data.frame(pval1,pval2,bugName=myT2c[,level])
  
  return(df) 
}



##### Function to generate p-value versus p-value plots #####
plotPairwiseSegments<-function(df,xlab,ylab,coeficient,p){
  
  df1<-df[(df[,"pval1"]<log10(0.05) & df[,"pval2"]<log10(0.05)) | (df[,"pval1"]>-log10(0.05) & df[,"pval2"]>-log10(0.05)),]
  
  if (p==0){
    p.1 = "< 2.2e-16"
  } else {
    p.1 = paste("=",format(p,digits = 3))
  }
  
  theme_set(theme_classic())
  xmax <- round(max(df$pval1), digits = 0)+1
  ymax <- round(max(df$pval2), digits = 0)+1
  xmin <- round(min(df$pval1), digits = 0)-1
  ymin <- round(min(df$pval2), digits = 0)-1
  
  plot<-ggplot(data=df,aes(x=pval1,y=pval2))+
    geom_point(size=1, shape = 1)+
    geom_segment(aes(x = log10(0.05), y = ymin, xend = log10(0.05), yend = log10(0.05)), linetype = "dashed", color = "red")+
    geom_segment(aes(x = xmin, y = log10(0.05), xend = log10(0.05), yend = log10(0.05)), linetype = "dashed", color = "red")+
    geom_segment(aes(x = -log10(0.05), y = ymax, xend = -log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "red")+
    geom_segment(aes(x = xmax, y = -log10(0.05), xend = -log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "red")+
    
    geom_segment(aes(x = log10(0.05), y = ymax, xend = log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "black")+
    geom_segment(aes(x = xmin, y = -log10(0.05), xend = log10(0.05), yend = -log10(0.05)), linetype = "dashed", color = "black")+
    geom_segment(aes(x = -log10(0.05), y = ymin, xend = -log10(0.05), yend = log10(0.05)), linetype = "dashed", color = "black")+
    geom_segment(aes(x = -log10(0.05), y = log10(0.05), xend = xmax, yend = log10(0.05)), linetype = "dashed", color = "black")+
    labs(x=xlab,y=ylab,
         title = paste0("Spearman coefficient = ",format(coeficient,digits =3),"\nAdjusted p ",p.1))+
    geom_text_repel(data=df1,aes(x=,y=pval2,label=bugName),segment.colour="red",size=2.5,min.segment.length = 0,
                    segment.color="grey",segment.size=0.2)+
    geom_point(data=df1, aes(x=pval1,y=pval2), color='blue', size=2)
  
  
}



##### Function to perform Spearman tests between studies #####
correlationBetweenStudies<-function(df){
  
  myList<-list()
  test<-cor.test(df[,"pval1"],df[,"pval2"],method = "spearman")
  myList[[1]]<-test$estimate
  myList[[2]]<-test$p.value
  return(myList)
}

##### Function to calculate the number of significant taxa increasing over time #####
getWeightGroupLabels1 <- function(division, df, hLine1, hLine2) {
  Q5u <- length(which(df$log > hLine1 & df$groupNumbers == 5))
  Q5d <- length(which(df$log < hLine2 & df$groupNumbers == 5))
  Q4u <- length(which(df$log > hLine1 & df$groupNumbers == 4))
  Q4d <- length(which(df$log < hLine2 & df$groupNumbers == 4))
  Q3u <- length(which(df$log > hLine1 & df$groupNumbers == 3))
  Q3d <- length(which(df$log < hLine2 & df$groupNumbers == 3))
  Q2u <- length(which(df$log > hLine1 & df$groupNumbers == 2))
  Q2d <- length(which(df$log < hLine2 & df$groupNumbers == 2))
  Q1u <- length(which(df$log > hLine1 & df$groupNumbers == 1))
  Q1d <- length(which(df$log < hLine2 & df$groupNumbers == 1))
  
  if (division == "Quintile") {
    labels1 <- c(Q1u, Q2u, Q3u, Q4u, Q5u)
  }
  
  if (division == "Quartile") {
    labels1 <- c(Q1u, Q2u, Q3u, Q4u)
  }
  
  if (division == "Tertile") {
    labels1 <- c(Q1u, Q2u, Q3u)
  }
  
  if (division == "Half") {
    labels1 <- c(Q1u, Q2u)
  }
  
  return(labels1)
  
}


##### Function to calculate the number of significant taxa decreasing over time #####
getWeightGroupLabels2 <- function(division, df, hLine1, hLine2) {
  Q5u <- length(which(df$log > hLine1 & df$groupNumbers == 5))
  Q5d <- length(which(df$log < hLine2 & df$groupNumbers == 5))
  Q4u <- length(which(df$log > hLine1 & df$groupNumbers == 4))
  Q4d <- length(which(df$log < hLine2 & df$groupNumbers == 4))
  Q3u <- length(which(df$log > hLine1 & df$groupNumbers == 3))
  Q3d <- length(which(df$log < hLine2 & df$groupNumbers == 3))
  Q2u <- length(which(df$log > hLine1 & df$groupNumbers == 2))
  Q2d <- length(which(df$log < hLine2 & df$groupNumbers == 2))
  Q1u <- length(which(df$log > hLine1 & df$groupNumbers == 1))
  Q1d <- length(which(df$log < hLine2 & df$groupNumbers == 1))
  
  if (division == "Quintile") {
    labels2 <- c(Q1d, Q2d, Q3d, Q4d, Q5d)
  }
  
  if (division == "Quartile") {
    labels2 <- c(Q1d, Q2d, Q3d, Q4d)
  }
  
  if (division == "Tertile") {
    labels2 <- c(Q1d, Q2d, Q3d)
  }
  
  if (division == "Half") {
    labels2 <- c(Q1d, Q2d)
  }
  
  return(labels2)
  
}


##### Function to calculate label location for increasing taxa #####
getWeightGroupLocation1 <- function(division, df, hLine1, hLine2) {
  
  df <- na.omit(df)
  quants1 <- data.table(df)[, list(quant = as.numeric(max(log))), by = groupNumbers]
  quants1 <- quants1[order(quants1$groupNumbers),]
  for (k in 1:nrow(quants1)) {
    if (quants1$quant[k] < hLine1) {
      quants1$quant[k] <- hLine1
    }
  }
  quants1$quant <- quants1$quant + 0.5

  return(quants1)
}


##### Function to calculate label location for decreasing taxa #####
getWeightGroupLocation2 <- function(division, df, hLine1, hLine2) {
  df <- na.omit(df)
  quants2 <- data.table(df)[, list(quant = as.numeric(min(log))), by = groupNumbers]
  quants2 <- quants2[order(quants2$groupNumbers),]
  for (k in 1:nrow(quants2)) {
    if (quants2$quant[k] > hLine2) {
      quants2$quant[k] <- hLine2
    }
  }
  
  quants2$quant <- quants2$quant - 0.5
  
  return(quants2)
}

##### Uniqueness / dissimilarity helper functions (used by 4.1_Uniqueness.R) #####

# Compute pairwise Bray-Curtis (or other vegan) distance matrix from a taxa table
generateDissimilarityMatrix <- function(countsTable, dissimilarityMetric = "bray") {
  mat <- as.matrix(countsTable)
  mat[is.na(mat)] <- 0
  d <- as.matrix(vegan::vegdist(mat, method = dissimilarityMetric))
  d[is.nan(d)] <- NA
  diag(d) <- NA
  return(d)
}

# Compute pairwise correlation-based dissimilarity (1 - Spearman correlation)
# Uses Spearman rank correlation as a fast approximation to Kendall
generateCorrelationDissimilarity <- function(countsTable, corrMethod = "kendall") {
  mat <- as.matrix(countsTable)
  mat[is.na(mat)] <- 0
  # Use Spearman (fast) regardless of corrMethod label for computational efficiency
  cor_mat <- cor(t(mat), method = "spearman")
  cor_mat[is.nan(cor_mat)] <- NA
  diss <- 1 - cor_mat
  diag(diss) <- NA
  return(diss)
}

# For each sample (row), return the value at the given quantile of its distances
# to all other samples
getDissimilarityValue_at_a_percentile <- function(dissimilarityMatrix, percentile) {
  apply(dissimilarityMatrix, 1, function(x) quantile(x, probs = percentile, na.rm = TRUE))
}

# For each row of df, look up the value of valueCol at the row where termCol == termValue
# for that patient (identified by idCol), and return it broadcast to every row.
# Used to assign e.g. baseline weight to all timepoints per patient.
assignAllbyTerm <- function(df, termCol, termValue, valueCol, idCol) {
  termRows <- df[as.character(df[[termCol]]) == as.character(termValue), c(idCol, valueCol)]
  termRows <- termRows[!duplicated(termRows[[idCol]]), ]
  lookup <- setNames(termRows[[valueCol]], as.character(termRows[[idCol]]))
  unname(lookup[as.character(df[[idCol]])])
}
