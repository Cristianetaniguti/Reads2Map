def test_genotyping_with_r(workflow_data, workflow_runner):
    inputs = {
        "analysis_bam": workflow_data["vcf_bi_bam_counts"],
        "analysis_vcf": workflow_data["vcf_bi"],
        "analysis_multi_vcf": workflow_data["vcf_multi"],
        "true_vcf": workflow_data["true_vcf"],
        "ref_alt_alleles": workflow_data["ref_alt_alleles"],
        "simulated_phases": workflow_data["simulated_phases"],
        "method": "gatk",
        "parent1": "P1",
        "parent2": "P2",
        "cross": "F1",
        "seed": 22,
        "depth": 10,
        "max_cores": 4,
    }

    expected = {"gusmap_out": [workflow_data["gusmap_rdata_1"], workflow_data["gusmap_rdata_2"]]}
    workflow_runner(
        "tests/workflows/TestGenotypingR.wdl",
        inputs,
        expected
    )
