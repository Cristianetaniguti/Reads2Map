
# Packages
library(lubridate)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)

# Functions
transform_day <- function(x){
  day <- as.numeric(x[1])*24
  split.t <- strsplit(x[2], ":")
  new.h <- as.numeric(sapply(split.t, "[[", 1)) + day
  paste0(new.h, ":",sapply(split.t, "[",2),":",sapply(split.t, "[",3))
}

efficiency_table <- function(infos_seff, jobs){
  
  df <- read.table(infos_seff, sep = "\n")
  
  cpu_used <- sapply(strsplit(df$V1[grep("CPU Utilized", df$V1)],": "), "[[",2)
  time <- sapply(strsplit(df$V1[grep("Job Wall-clock time", df$V1)],": "), "[[",2)
  
  cpu_used[grep(cpu_used, pattern = "-")] <- sapply(strsplit(cpu_used[grep(cpu_used, pattern = "-")], "-"), transform_day)
  time[grep(cpu_used, pattern = "-")] <- sapply(strsplit(time[grep(time, pattern = "-")], "-"), transform_day)
  
  df.infos <- data.frame(
    cores = sapply(strsplit(df$V1[grep("Cores", df$V1)], ":"), "[[",2),
    cpu_used = period_to_seconds(hms(cpu_used)),
    cpu_eff = sapply(strsplit(df$V1[grep("CPU Efficiency", df$V1)],": "), "[[",2),
    time = period_to_seconds(hms(time)),
    mem_used_GB = sapply(strsplit(df$V1[grep("Memory Utilized", df$V1)],": "), "[[",2),
    mem_eff = sapply(strsplit(df$V1[grep("Memory Efficiency", df$V1)],": "), "[[",2)
  )
  
  df.infos <- cbind(jobs, df.infos)
  
  df.infos$mem_GB <- sapply(strsplit(df.infos$mem_eff, "of"),"[[",2)
  
  idx.mb <- grep(df.infos$mem_GB, pattern = "MB")
  idx.kb <- grep(df.infos$mem_GB, pattern = "KB")
  
  df.infos$mem_GB <- as.numeric(sapply(strsplit(df.infos$mem_GB, " "), "[[",2))
  df.infos$mem_GB[idx.mb] <- df.infos$mem_GB[idx.mb]/1000
  df.infos$mem_GB[idx.kb] <- df.infos$mem_GB[idx.kb]/1000000
  
  idx.mb <- grep(df.infos$mem_used, pattern = "MB")
  idx.kb <- grep(df.infos$mem_used, pattern = "KB")
  
  df.infos$mem_used_GB <- as.numeric(sapply(strsplit(df.infos$mem_used_GB, " "), "[[",1))
  df.infos$mem_used_GB[idx.mb] <- df.infos$mem_used_GB[idx.mb]/1000
  df.infos$mem_used_GB[idx.kb] <- df.infos$mem_used[idx.kb]/1000000
  
  df.infos$cpu <- sapply(strsplit(df.infos$cpu_eff, "of"),"[[",2)
  df.infos$cpu <- gsub(df.infos$cpu, pattern = " core-walltime", replacement = "")
  
  df.infos$cpu <- sapply(strsplit(df.infos$cpu, " "), "[[",2)
  df.infos$cpu[grep(df.infos$cpu, pattern = "-")] <- sapply(strsplit(df.infos$cpu[grep(df.infos$cpu, pattern = "-")], "-"), transform_day)
  
  df.infos$cpu <- period_to_seconds(hms(df.infos$cpu))
  
  df.infos <- df.infos[,-c(5,8)]
  colnames(df.infos)[2] <- "job_id"
  
  return(df.infos)
}

efficiency_graphics <- function(eff_table, out_name){
  perc <- eff_table %>% 
    mutate(mem = (mem_used_GB/mem_GB)*100 , cpu = (cpu_used/cpu)*100) %>% 
    select(tasks, mem, cpu) %>% pivot_longer(cols = c(2:3), names_to="Usage (%)", values_to = "percentage") %>%
    ggplot(aes(x=tasks, y=percentage, fill=`Usage (%)`, color = `Usage (%)`)) +   coord_flip() +
    geom_boxplot()  + theme_bw()
  
  mem <- eff_table %>%  
    select(tasks, mem_used_GB) %>% pivot_longer(cols = c(2), values_to="usage (GB)") %>%
    ggplot(aes(x=tasks, y=`usage (GB)`)) +   coord_flip() +
    geom_boxplot()  + ylab("usage mem (GB)")  + theme_bw() 
  
  cpu <- eff_table %>%  
    select(tasks, cpu_used) %>% pivot_longer(cols = c(2), values_to = "usage (min)") %>%
    ggplot(aes(x=tasks, y = `usage (min)`/60)) +   coord_flip() +
    geom_boxplot() + ylab("usage cpus (min)") + theme_bw() +
    theme(axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) 
  
  p <- ggarrange(mem, cpu, nrow = 1)
  (p2 <- ggarrange(p, perc, nrow = 2))
  
  ggsave(p2, filename = out_name, width = 10, height = 7)
}

# grep -rw "Submitted batch job" * > jobs
# df <- read.table("/scratch/user/chtaniguti/cromwell-executions/SNPCalling/7cefff16-8d18-46d7-adeb-b22e4f35b663/jobs")

df <- read.table("jobs_acca")
tasks <- sapply(strsplit(df$V1, "/"), function(x) paste0(x[grep("call-", x)], collapse = "/"))
tasks <- gsub("call-", "", tasks)
jobs_acca <- data.frame(tasks, df$V4)
write.table(df$V4, quote = F, file = "jobs_acca_id", col.names = F, row.names = F)


df <- read.table("jobs_euc") 
tasks <- sapply(strsplit(df$V1, "/"), function(x) paste0(x[grep("call-", x)], collapse = "/"))
tasks <- gsub("call-", "", tasks)
jobs_euc <- data.frame(tasks, df$V4)
write.table(df$V4, quote = F, file = "jobs_euc_id", col.names = F, row.names = F)

# for i in $(cat jobs_euc_id);do
#   seff $i >> infos.temp.euc
# done

df.acca <- efficiency_table(infos_seff = "infos.temp.acca", jobs_acca)
efficiency_graphics(df.acca, out_name = "acca_efficiency_wf.png")

df.euc <- efficiency_table(infos_seff = "infos.temp.euc", jobs_euc)
efficiency_graphics(df.euc, out_name = "euc_efficiency_wf.png")

