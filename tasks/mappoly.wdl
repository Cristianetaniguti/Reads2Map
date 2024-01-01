version 1.0

task MappolyReport {
  input {
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String parent1
    String parent2
    Float prob_thres = 0.8
    Int repetitions = 100
    Int sample_size = 30
    Int max_cores
    Int ploidy
    String filt_segr = "TRUE"
    Array[String] global_errors = ["0.05"]
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 2)
  Int memory_size = ceil(size(vcf_file, "MiB") * max_cores + 12000*max_cores)

  command <<<
    R --vanilla --no-save <<RSCRIPT
      library(mappoly)
      library(parallel)
      
      run_tests <- function(seed, s_all, sample_size, dat){
        set.seed(seed[1])
        library(mappoly)
        map.err.p1.idx <- map.err.p2.idx <- map.prob.p1.idx <- map.prob.p2.idx <- 1
        res.p1 <- res.p2 <- res.prob.p1 <- res.prob.p2 <- NULL
        map.err.p1 <- map.err.p2 <- map.prob.p1 <- map.prob.p2 <- vector("list", length(seed))
        while(map.err.p1.idx <= length(seed) | map.err.p2.idx <= length(seed) | map.prob.p1.idx <= length(seed) | map.prob.p2.idx <= length(seed)){  
          tryCatch({
            mrk.id<-sample(s_all[["seq.mrk.names"]], sample_size)
            s<-get_genomic_order(make_seq_mappoly(s_all,mrk.id))
            
            ## Parent 1
            s1<-make_seq_mappoly(s, info.parent = "p1")
            tpt <- est_pairwise_rf(s1)
            map <- est_rf_hmm_sequential(input.seq = s1,
                                         start.set = 5,
                                         thres.twopt = 10,
                                         thres.hmm = 10,
                                         extend.tail = 30,
                                         info.tail = TRUE,
                                         twopt = tpt,
                                         phase.number.limit = 10,
                                         reestimate.single.ph.configuration = TRUE,
                                         tol = 10e-2,
                                         tol.final = 10e-4, 
                                         verbose = FALSE)
            map <- filter_map_at_hmm_thres(map, thres.hmm = 0.0001)
            
            # Global error
            if(map.err.p1.idx <= length(seed)){
               map2 <- est_full_hmm_with_global_error(map, error = 0.05, tol = 10e-3)
               map3 <- split_and_rephase(map2, gap.threshold = 20, size.rem.cluster = 3, twopt = tpt)
               map.err.p1[[map.err.p1.idx]] <- est_full_hmm_with_global_error(map3, error = 0.05, tol = 10e-4)
               x<-summary_maps(list(map.err.p1[[map.err.p1.idx]]))
               res.p1 <-rbind(res.p1,x[1,])
               map.err.p1.idx <- map.err.p1.idx + 1
            }
            
            # Prob error
            if(!is.null(dat[["geno"]]) & map.prob.p1.idx <= length(seed)){
              map2 <- est_full_hmm_with_prior_prob(map, tol = 10e-3)
              map3 <- split_and_rephase(map2, gap.threshold = 20, size.rem.cluster = 3, twopt = tpt)
              map.prob.p1[[map.prob.p1.idx]] <- est_full_hmm_with_prior_prob(map3, tol = 10e-4)
              x<-summary_maps(list(map.prob.p1[[map.prob.p1.idx]]))
              res.prob.p1 <-rbind(res.prob.p1,x[1,])
              map.prob.p1.idx <- map.prob.p1.idx + 1          
            } else if (is.null(dat[["geno"]])) map.prob.p1.idx <- map.prob.p1.idx + 1
            
            # Parent 2
            s2<-make_seq_mappoly(s, info.parent = "p2")
            tpt <- est_pairwise_rf(s2)
            map <- est_rf_hmm_sequential(input.seq = s2,
                                         start.set = 5,
                                         thres.twopt = 10,
                                         thres.hmm = 10,
                                         extend.tail = 30,
                                         info.tail = TRUE,
                                         twopt = tpt,
                                         phase.number.limit = 10,
                                         reestimate.single.ph.configuration = TRUE,
                                         tol = 10e-2,
                                         tol.final = 10e-4, 
                                         verbose = FALSE)
            map <- filter_map_at_hmm_thres(map, thres.hmm = 0.0001)
            
            # Global error
            if(map.err.p2.idx <= length(seed)){
               map2 <- est_full_hmm_with_global_error(map, error = 0.05, tol = 10e-3)
               map3 <- split_and_rephase(map2, gap.threshold = 20, size.rem.cluster = 3, twopt = tpt)
               map.err.p2[[map.err.p2.idx]] <- est_full_hmm_with_global_error(map3, error = 0.05, tol = 10e-4)
               x<-summary_maps(list(map.err.p2[[map.err.p2.idx]]))
               res.p2<-rbind(res.p2,x[1,])
               map.err.p2.idx <- map.err.p2.idx + 1            
            }
            
            # Prob error
            if(!is.null(dat[["geno"]]) & map.prob.p2.idx <= length(seed)){
              map2 <- est_full_hmm_with_prior_prob(map, tol = 10e-3)
              map3 <- split_and_rephase(map2, gap.threshold = 20, size.rem.cluster = 3, twopt = tpt)
              map.prob.p2[[map.prob.p2.idx]] <- est_full_hmm_with_prior_prob(map3, tol = 10e-4)
              x<-summary_maps(list(map.prob.p2[[map.prob.p2.idx]]))
              res.prob.p2 <-rbind(res.prob.p2,x[1,])
              map.prob.p2.idx <- map.prob.p2.idx + 1          
            } else if (is.null(dat[["geno"]])) map.prob.p2.idx <- map.prob.p2.idx + 1
          }, error=function(e){})
        }
        return(list(map.err.p1, map.err.p2, map.prob.p1, map.prob.p2, 
                    res.p1, res.p2, res.prob.p1, res.prob.p2))
      }
      
      if("~{GenotypeCall_program}" == "supermassa") 
        prob.thres = ~{prob_thres} - ~{prob_thres}*0.3 else prob.thres = ~{prob_thres}
      
      dat <- read_vcf(file = "~{vcf_file}", 
                      parent.1 = "~{parent1}", 
                      parent.2 = "~{parent2}", 
                      verbose = FALSE, 
                      read.geno.prob = TRUE, 
                      prob.thres = prob.thres, 
                      ploidy = ~{ploidy})
      
      info <- data.frame(dat = paste0("~{SNPCall_program}", "_", "~{GenotypeCall_program}", "_", "~{CountsFrom}"), 
                         step = "raw", 
                         n.markers = dat[["n.mrk"]], 
                         mis.perc = round((sum(dat[["geno.dose"]] == dat[["ploidy"]]+1)/(nrow(dat[["geno.dose"]])*ncol(dat[["geno.dose"]])))*100,2),
                         n.redundant = round(100*(nrow(dat[["elim.correspondence"]])/(length(dat[["kept"]])+nrow(dat[["elim.correspondence"]]))),2), 
                         n.ind = dat[["n.ind"]])
      
      dat <- filter_missing(input.data = dat, type = "marker", 
                            filter.thres = 0.25, inter = FALSE)
      
      info_temp <- data.frame(dat = paste0("~{SNPCall_program}", "_", "~{GenotypeCall_program}", "_", "~{CountsFrom}"), 
                              step = "miss filtered", 
                              n.markers = dat[["n.mrk"]], 
                              mis.perc = round((sum(dat[["geno.dose"]] == dat[["ploidy"]]+1)/(nrow(dat[["geno.dose"]])*ncol(dat[["geno.dose"]])))*100,2),
                              n.redundant = round(100*(nrow(dat[["elim.correspondence"]])/(length(dat[["kept"]])+nrow(dat[["elim.correspondence"]]))),2), 
                              n.ind = dat[["n.ind"]])
      
      info <- rbind(info, info_temp)
      
      
      if("TRUE"){
        pval.bonf <- 0.05/dat[["n.mrk"]]
        mrks.chi.filt <- filter_segregation(dat, 
                                            chisq.pval.thres =  pval.bonf, 
                                            inter = FALSE)
        
        seq.init <- make_seq_mappoly(mrks.chi.filt)
      } else {
        seq.init <- make_seq_mappoly(dat, "all")
      }
      
      info_temp <- data.frame(dat = paste0("~{SNPCall_program}", "_", "~{GenotypeCall_program}", "_", "~{CountsFrom}"), 
                              step = "segr filtered", 
                              n.markers = length(seq.init[["seq.mrk.names"]]), 
                              mis.perc = round((sum(dat[["geno.dose"]][which(rownames(dat[["geno.dose"]]) %in% seq.init[["seq.mrk.names"]]),] == dat[["ploidy"]]+1)/(ncol(dat[["geno.dose"]])*length(seq.init[["seq.mrk.names"]])))*100,2), 
                              n.redundant = round(100*(nrow(dat[["elim.correspondence"]])/(length(dat[["kept"]])+nrow(dat[["elim.correspondence"]]))),2), 
                              n.ind = dat[["n.ind"]])
      
      info <- rbind(info, info_temp)
      
      sp1<-make_seq_mappoly(seq.init, info.parent = "p1")
      
      info_temp <- data.frame(dat = paste0("~{SNPCall_program}", "_", "~{GenotypeCall_program}", "_", "~{CountsFrom}"), 
                              step = "p1", 
                              n.markers = length(sp1[["seq.mrk.names"]]), 
                              mis.perc = round((sum(dat[["geno.dose"]][which(rownames(dat[["geno.dose"]]) %in% sp1[["seq.mrk.names"]]),] == dat[["ploidy"]]+1)/(ncol(dat[["geno.dose"]])*length(sp1[["seq.mrk.names"]])))*100,2), 
                              n.redundant = round(100*(nrow(dat[["elim.correspondence"]])/(length(dat[["kept"]])+nrow(dat[["elim.correspondence"]]))),2), 
                              n.ind = dat[["n.ind"]])
      info <- rbind(info, info_temp)
      
      sp2<-make_seq_mappoly(seq.init, info.parent = "p2")
      
      info_temp <- data.frame(dat = paste0("~{SNPCall_program}", "_", "~{GenotypeCall_program}", "_", "~{CountsFrom}"), 
                              step = "p2", 
                              n.markers = length(sp2[["seq.mrk.names"]]), 
                              mis.perc = round((sum(dat[["geno.dose"]][which(rownames(dat[["geno.dose"]]) %in% sp2[["seq.mrk.names"]]),] == dat[["ploidy"]]+1)/(ncol(dat[["geno.dose"]])*length(sp2[["seq.mrk.names"]])))*100,2), 
                              n.redundant = round(100*(nrow(dat[["elim.correspondence"]])/(length(dat[["kept"]])+nrow(dat[["elim.correspondence"]]))),2), 
                              n.ind = dat[["n.ind"]])
      info <- rbind(info, info_temp)
      
      # Estimate two-point recombination fraction
      tpt <- est_pairwise_rf(input.seq = seq.init, ncpus = ~{max_cores})
      mat <- rf_list_to_matrix(input.twopt = tpt)
      
      sample_size <- if(length(seq.init[["seq.mrk.names"]]) <= ~{sample_size}) length(seq.init[["seq.mrk.names"]]) else ~{sample_size}
      seeds <- split(1:~{repetitions}, rep(1:~{max_cores}, each=~{repetitions}/~{max_cores}))
      
      clust <- makeCluster(~{max_cores})
      clusterExport(clust, c("seq.init", "sample_size", "dat", "run_tests"))
      results <- parLapply(clust, seeds, function(x) run_tests(seed = x, s_all = seq.init, 
                                                               sample_size = sample_size, dat = dat))
      stopCluster(clust)
      
      map.err.p1 <- sapply(results, "[[", 1)
      idx <- sapply(map.err.p1, is.null)
      map.err.p1 <- map.err.p1[which(!idx)]
      
      map.err.p2 <- sapply(results, "[[", 2)
      idx <- sapply(map.err.p2, is.null)
      map.err.p2 <- map.err.p2[which(!idx)]
      
      if(!is.null(dat[["geno"]])){
        map.prob.p1 <- sapply(results, "[[", 3)
        idx <- sapply(map.prob.p1, is.null)
        map.prob.p1 <- map.prob.p1[which(!idx)]
        
        map.prob.p2 <- sapply(results, "[[", 4)
        idx <- sapply(map.prob.p2, is.null)
        map.prob.p2 <- map.prob.p2[which(!idx)]
        
        maps <- list(map.err.p1 = map.err.p1, 
                     map.err.p2 = map.err.p2, 
                     map.prob.p1 = map.prob.p1, 
                     map.prob.p2 = map.prob.p2)
      } else {
        maps <- list(map.err.p1 = map.err.p1, 
                     map.err.p2 = map.err.p2) 
      }
      
      res.p1 <- lapply(results, "[[", 5)
      res.p1 <- do.call(rbind, res.p1)
      res.p1[["map"]] <- "error.p1"
      
      res.p2 <- lapply(results, "[[", 6)
      res.p2 <- do.call(rbind, res.p2)
      res.p2[["map"]] <- "error.p2"
      
      if(!is.null(dat[["geno"]])){
        res.prob.p1 <- lapply(results, "[[", 7)
        res.prob.p1 <- do.call(rbind, res.prob.p1)
        res.prob.p1[["map"]] <- "prob.p1"
        
        res.prob.p2 <- lapply(results, "[[", 8)
        res.prob.p2 <- do.call(rbind, res.prob.p2)
        res.prob.p2[["map"]] <- "prob.p2"
        summaries <- rbind(res.p1, res.p2, res.prob.p1, res.prob.p2)
      } else summaries <- rbind(res.p1, res.p2)
      
      summaries[["data"]] <- "~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}"
      
      saveRDS(summaries, file= "~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_summaries.rds")
      saveRDS(info, file="~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_info.rds")
      saveRDS(dat, file= "~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_dat.rds")
      saveRDS(mat, file="~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_mat2.rds")
      saveRDS(maps, file="~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_maps.rds")
      
      system("mkdir results")
      system("mv *.rds  results")
      
      system(paste0("tar -czvf ", "~{SNPCall_program}", "_", "~{GenotypeCall_program}", "_", "~{CountsFrom}","_poly_results.tar.gz results"))

    RSCRIPT
  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.1.0"
    singularity: "docker://cristaniguti/reads2map:0.1.0"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "MappolyReport"
    mem:"~{memory_size}G"
    time: 24
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Build linkage map using MAPpoly"
  }

  output {
    File results = "~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_poly_results.tar.gz"
  }
}