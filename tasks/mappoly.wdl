version 1.0

task MappolyReport {
  input {
    File vcf_file
    String SNPCall_program
    String GenotypeCall_program
    String CountsFrom
    String parent1
    String parent2
    Int max_cores
    Int ploidy
  }

  Int disk_size = ceil(size(vcf_file, "GiB") * 2)
  Int memory_size = ceil(size(vcf_file, "MiB") + 2000)

  command <<<
    R --vanilla --no-save <<RSCRIPT

        library(mappoly)

        dat <- read_vcf(file = "~{vcf_file}", 
                        parent.1 = "~{parent1}", 
                        parent.2 = "~{parent2}", 
                        verbose = FALSE, 
                        read.geno.prob = TRUE, 
                        prob.thres = 0.8, ploidy = ~{ploidy})

        png(paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"_raw_data.png"))
        plot(dat)
        dev.off()

        dat <- filter_missing(input.data = dat, type = "marker", 
                            filter.thres = 0.25, inter = FALSE)

        pval.bonf <- 0.05/dat[[3]]
        mrks.chi.filt <- filter_segregation(dat, 
                                            chisq.pval.thres =  pval.bonf, 
                                            inter = FALSE)

        seq.init <- make_seq_mappoly(mrks.chi.filt)

        png(paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"_","filters.png"))
        plot(seq.init)
        dev.off()

        all.rf.pairwise <- est_pairwise_rf(input.seq = seq.init, ncpus = ~{max_cores})
        mat <- rf_list_to_matrix(input.twopt = all.rf.pairwise)

        id<-get_genomic_order(seq.init)
        s.o <- make_seq_mappoly(id)

        png(paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"_","rf.png"))
        plot(mat, ord = s.o[[3]])
        dev.off()

        tpt <- make_pairs_mappoly(all.rf.pairwise, input.seq = s.o)
        temp2 <- rf_snp_filter(input.twopt = tpt, diagnostic.plot = FALSE)
        lgtemp <- get_genomic_order(temp2)
        s.o <- make_seq_mappoly(lgtemp)

        png(paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"_","rf.filt.png"))
        plot(mat, ord = s.o[[3]])
        dev.off()

        est.map <- est_rf_hmm_sequential(input.seq = s.o,
                                        start.set = 5,
                                        thres.twopt = 10,
                                        thres.hmm = 50,
                                        extend.tail = 30,
                                        twopt =  all.rf.pairwise,
                                        verbose = F,
                                        phase.number.limit = 20,
                                        sub.map.size.diff.limit = 5)

        map.err <- est_full_hmm_with_global_error(input.map = est.map, error = 0.05)
        map.prob <- est_full_hmm_with_prior_prob(input.map = est.map, dat.prob = dat)

        png(paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"_no_error_cMbyMb.png"))
        plot_genome_vs_map(est.map, same.ch.lg = TRUE)
        dev.off()

        png(paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"global_error_cMbyMb.png"))
        plot_genome_vs_map(map.err, same.ch.lg = TRUE)
        dev.off()

        png(paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"_probs_cMbyMb.png"))
        plot_genome_vs_map(map.prob, same.ch.lg = TRUE)
        dev.off()

        summary <- summary_maps(list(est.map, map.err, map.prob))
        summary <- cbind(method = c("no_error", "global_error", "probs", "-"), summary)

        write.csv(summary, file = paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"_map_summary.csv"))

        export_map_list(est.map, file = paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"_no_error_","map_file.csv"))
        export_map_list(map.err, file = paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"_global_error_","map_file.csv"))
        export_map_list(map.prob, file = paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"_probs_","map_file.csv"))

        png(paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"_map_draw.png"))
        plot_map_list(list(default = est.map, 
                        global = map.err,
                        probs = map.prob), col = "ggstyle")
        dev.off()

        genoprob <- calc_genoprob_error(input.map = est.map, error = 0)
        genoprob.err <- calc_genoprob_error(input.map = map.err, error = 0.05)
        genoprob.prob <- calc_genoprob_dist(input.map = map.prob, dat.prob = dat)

        homoprobs = calc_homologprob(genoprob)
        homoprobs.err = calc_homologprob(genoprob.err)
        homoprobs.prob = calc_homologprob(genoprob.prob)

        save(homoprobs, file = paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"homoprobs.RData"))
        save(homoprobs.err, file = paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"homoprobs.err.RData"))
        save(homoprobs.prob, file = paste0("~{SNPCall_program}", "_","~{GenotypeCall_program}", "_", "~{CountsFrom}" ,"homoprobs.prob.RData"))

        system("mkdir results")
        system("mv *.png *.RData *csv results")
        system(paste0("tar -czvf ", "~{SNPCall_program}", "_", "~{GenotypeCall_program}", "_", "~{CountsFrom}","_results.tar.gz results"))

    RSCRIPT

  >>>

  runtime {
    docker:"cristaniguti/reads2map:0.0.5"
    cpu: max_cores
    # Cloud
    memory:"~{memory_size} MiB"
    disks:"local-disk " + disk_size + " HDD"
    # Slurm
    job_name: "MappolyReport"
    mem:"~{memory_size}G"
    time:"24:00:00"
  }

  meta {
        author: "Cristiane Taniguti"
        email: "chtaniguti@tamu.edu"
        description: "Build linkage map using MAPpoly"
  }

  output {
    File results = "~{SNPCall_program}_~{GenotypeCall_program}_~{CountsFrom}_results.tar.gz"
  }
}