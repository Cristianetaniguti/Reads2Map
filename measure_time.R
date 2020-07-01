#' Measuring cromwell tasks time
#' 
#' @param log.file character with log file name in cromwell-workflow-logs directory after 
#' running WDL workflow
#' 
#' @return ggplot boxplot with time spent by each task
#' 
workflow_times <- function(log.file){
  log.file <- "cromwell-workflow-logs/workflow.ee9b225d-0aba-4f99-ab6a-050757c53f96.log"
  system(paste("grep 'job id'",log.file, "> temp.starting"))
  system(paste("grep 'to Done'",log.file, "> temp.done"))
  
  start <- read.table("temp.starting")
  dates.start <- as.POSIXct(paste(start$V1, start$V2))
  start.df <- data.frame(dates.start, id=factor(start$V6, levels = unique(as.character(start$V6))))
  
  done <- read.table("temp.done")
  dates.done <- as.POSIXct(paste(done$V1, done$V2))
  done.df <- data.frame(dates.done, id= factor(done$V6, levels = unique(as.character(done$V6))))
  
  file.remove(c("temp.starting", "temp.done"))
  
  tot.df <- merge(done.df, start.df)
  
  tot.df <- cbind(tot.df, diff=tot.df[,2] - tot.df[,3])
  
  tot.df <- cbind(task=sapply(strsplit(as.character(tot.df$id), ":"), "[", 1), tot.df)
  
  p <- ggplot(tot.df, aes(x=task, y=diff)) + geom_boxplot() +
    guides(fill=FALSE) + coord_flip() + ylab("Time in seconds") + scale_y_continuous()

  # fig <- ggplotly(p)
  # fig
  
  return(p)
}