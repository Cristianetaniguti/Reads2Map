def test_read_simu(workflow_data, workflow_runner):
    inputs = {
        "parent1": "P1",
        "parent2": "P2",
        "sampleNames": ["P1", "P2", "F1_01", "F1_02"],
        "program": "gatk",
        "references": {
            "ref_fasta": workflow_data["ref_fasta"],
            "ref_dict": workflow_data["ref_dict"],
            "ref_ann": workflow_data["ref_ann"],
            "ref_sa": workflow_data["ref_sa"],
            "ref_amb": workflow_data["ref_amb"],
            "ref_bwt": workflow_data["ref_bwt"],
            "ref_fasta_index": workflow_data["ref_fasta_index"],
            "ref_pac": workflow_data["ref_pac"]
        },
        "bams": [workflow_data["F1_01_bam"], workflow_data["F1_02_bam"], workflow_data["P1_bam"], workflow_data["P2_bam"]],
        "bais": [workflow_data["F1_01_bai"], workflow_data["F1_02_bai"], workflow_data["P1_bai"], workflow_data["P2_bai"]]
    }

    expected = {"vcf_bi": workflow_data["smallest_vcf_bi"]}
    workflow_runner(
        "tasks/gatk_genotyping.wdl",
        inputs,
        expected
    )
