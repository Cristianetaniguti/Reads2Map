from pathlib import Path

import pytest

# Parei aqui! Adaptar os caminhos
@pytest.mark.workflow("SNP and genotype calling with GATK")
def test_ensure_vcf_1_are_equal(workflow_dir, compare_files):
    vcf1 = Path(workflow_dir, "_LAST/out/vcf_norm/vcf_norm.vcf.gz")
    vcf2 = Path("tests/data/Ptremula_PRJNA395596_subset/snpcalling_vcfs/gatk/vcf_norm.vcf.gz")
    vcf3 = Path(workflow_dir, "_LAST/out/vcf_norm_bamcounts/gatk_bam_vcf.vcf.gz")
    vcf4 = Path("tests/data/Ptremula_PRJNA395596_subset/snpcalling_vcfs/gatk/gatk_bam_vcf.vcf.gz")

    # Asserts
    assert compare_files.assert_vcf_equal(vcf1, vcf2)
    assert compare_files.assert_vcf_equal(vcf3, vcf4)
