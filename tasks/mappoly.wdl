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
    Int max_cores
    Int ploidy
    String filt_segr = "TRUE"
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 2)
  Int memory_size = ceil(size(vcf_file, "MiB") + 2000)

  command <<<
    R --vanilla --no-save <<RSCRIPT

        library(mappoly)

        if("~{GenotypeCall_program}" == "supermassa") prob.thres = ~{prob_thres} - 0.3 else prob.thres = ~{prob_thres}
        
        dat <- read_vcf(file = "~{vcf_file}", 
                        parent.1 = "~{parent1}", 
                        parent.2 = "~{parent2}", 
                        verbose = FALSE, 
                        read.geno.prob = TRUE, 
                        prob.thres = prob.thres, 
                        ploidy = ~{ploidy})

        dat <- filter_missing(input.data = dat, type = "marker", 
                            filter.thres = 0.25, inter = FALSE)

        dat <- filter_missing(dat, type = 'individual', filter.thres = 0.1, inter = FALSE)

        if("~{filt_segr}"){
          pval.bonf <- 0.05/dat[[3]]
          mrks.chi.filt <- filter_segregation(dat, 
                                              chisq.pval.thres =  pval.bonf, 
                                              inter = FALSE)

          seq.init <- make_seq_mappoly(mrks.chi.filt)
        } else {
          seq.init <- make_seq_mappoly(dat, "all")
        }
        
        # Estimate two-point recombination fraction
        tpt <- est_pairwise_rf(input.seq = seq.init, ncpus = ~{max_cores})
        mat <- rf_list_to_matrix(input.twopt = tpt)

        # Filter markers by recombination fraction values
        seq.filt <- rf_snp_filter(input.twopt = tpt, diagnostic.plot = FALSE, probs = c(0.05, 0.95))
        mat2 <- make_mat_mappoly(mat, seq.filt)

        # Run MDS for ordering according to recombination fractions
        seq_test_mds <- mds_mappoly(mat2)
        seq_mds <- make_seq_mappoly(seq_test_mds)

        # Sequence with genomic order
        geno_order <- get_genomic_order(seq_mds)
        seq_geno_order <- make_seq_mappoly(geno_order)

        init.map.list <- framework_map(input.seq = seq_geno_order,
                                       twopt = tpt,
                                       start.set = 5,
                                       inflation.lim.p1 = 10,
                                       inflation.lim.p2 = 10,
                                       verbose = FALSE)

        res <- update_framework_map(input.map.list = init.map.list,
                                    input.seq = seq_geno_order,
                                    twopt = tpt,
                                    thres.twopt = 5,
                                    init.LOD = 100,
                                    max.rounds = 3,
                                    size.rem.cluster = 3,
                                    gap.threshold = 3,
                                    verbose = FALSE)

        # Get last interaction
        iter <- length(res[[2]][[1]])
        map_error <- est_full_hmm_with_global_error(res[[2]][[1]][[iter]], error = 0.05, verbose = FALSE)
        map_prob <- est_full_hmm_with_prior_prob(res[[2]][[1]][[iter]], dat.prob = dat, verbose = FALSE)

        # Diagnostic graphics - overall from plot(dat)
        # Heatmap of mds with plot(mat2, ord=seq_mds)
        # Heatmap of genomic order order plot(mat2, ord = map_error$info$mk.names)
        # Relation between mds and genome plot(seq_mds$genome.pos)

        saveRDS(dat, file= "~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_dat.rds")
        saveRDS(mat2, file="~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_mat2.rds")
        saveRDS(seq_mds, file="~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_seq_mds.rds")
        saveRDS(map_error, file="~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_map_error.rds")
        saveRDS(map_prob, file= "~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_map_prob.rds")

        system("mkdir results")
        system("mv *.rds  results")

        system(paste0("tar -czvf ", "~{SNPCall_program}", "_", "~{GenotypeCall_program}", "_", "~{CountsFrom}","_results.tar.gz results"))

    RSCRIPT

  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.0.9"
    singularity: "docker://cristaniguti/reads2map:0.0.9"
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
    File results = "~{SNPCall_program}_~{GenotypeCall_program}Poly_~{CountsFrom}_results.tar.gz"
  }
}