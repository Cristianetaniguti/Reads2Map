def test_read_simulation_workflow(workflow_data, workflow_runner):
    inputs = {
        "family": {
            "cmBymb": None,
            "cross": "F1",
            "doses": None,
            "ploidy": 2,
            "popsize": 20,
            "seed": 22
        },
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
        "sequencing":{
            "emp_vcf": workflow_data["emp_vcf_no_indels"],
            "enzyme1": "HinDIII",
            "enzyme2": "NlaIII",
            "library_type": "ddRAD",
            "chromosome": "Chr10",
            "pcr_cycles": 0,
            "read_length": 300,
            "insert_size": 500,
            "insert_size_dev": 30,
            "depth": 10,
            "depth_parents": 10,
            "ref_map": workflow_data["ref_map"],
            "multiallelics": "TRUE",
            "vcf_parent1": "PT_F",
            "vcf_parent2": "PT_M"
        },
        "max_cores": 4
    }

    expected = {
        "simu_haplo": workflow_data["simu_haplo"],
        "simulated_phases": workflow_data["simulated_phases"],
        "true_vcf": workflow_data["true_vcf"],
        "ref_alt_alleles": workflow_data["ref_alt_alleles"],
    }
    workflow_runner(
        "tasks/create_alignment_from_read_simulations.wdl",
        inputs,
        expected
    )
