library(survey)


`%notin%` <- Negate(`%in%`)


MergeSubjectTime <- function(df1, df2, mergecol, timecol1, timecol2) {
  
  
  m_data <- lapply(intersect(df1[[mergecol]], df2[[mergecol]]),function(id) {
    d1   <- subset(df1,df1[[mergecol]] == id)
    d2   <- subset(df2,df2[[mergecol]] == id)
    
    d1$indices <- sapply(d1[[timecol1]],function(d) which.min(abs(d2[[timecol2]] - d)))
    d2$indices <- 1:nrow(d2)
    base::merge(d1, d2, by = c(mergecol, 'indices'), all.x=FALSE, all.y=FALSE)
    
  })
  mergeddata <- do.call(rbind, m_data)
  mergeddata$indices <- NULL
  return(mergeddata)
}

TimeSinceBaseline <- function(data, timecol) {
  subjlist <- split(data, data$RID)
  return.list <- list()
  for(i in 1:length(subjlist)) {
    subj <- subjlist[[i]]
    min.time <- min(subj[[timecol]])
    subj$new_time <- subj[[timecol]] - min.time
    return.list[[i]] <- subj
  }
  return.list <- do.call(rbind, return.list)
  return(return.list)
}





