from pathlib import Path

import pytest


@pytest.mark.workflow("Freebayes SNP calling")
def test_ensure_vcf_1_are_equal(workflow_dir, compare_files):
    vcf1 = Path(workflow_dir, "_LAST/out/vcfs/0/vcf_norm.vcf.gz")
    vcf2 = Path("tests/data/Ptremula_PRJNA395596_subset/snpcalling_vcfs/freebayes/vcf_norm.vcf.gz")
    vcf3 = Path(workflow_dir, "_LAST/out/vcfs/1/freebayes_bam_vcf.vcf.gz")
    vcf4 = Path("tests/data/Ptremula_PRJNA395596_subset/snpcalling_vcfs/freebayes/freebayes_bam_vcf.vcf.gz")


    # Asserts
    assert compare_files.assert_vcf_equal(vcf1, vcf2)
    assert compare_files.assert_vcf_equal(vcf3, vcf4)
