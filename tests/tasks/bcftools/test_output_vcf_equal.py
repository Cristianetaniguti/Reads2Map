from pathlib import Path

import pytest


@pytest.mark.workflow("Use bcftools to normalize variant representation in VCF")
def test_ensure_vcf_1_are_equal(workflow_dir, compare_files):
    vcf1 = Path(workflow_dir, "_LAST/out/vcf_norm/vcf_norm.vcf.gz")
    vcf2 = Path("tests/data/ecoli/ecoli_first_10k_bases.normalized.vcf.gz")

    # Asserts
    assert compare_files.assert_vcf_equal(vcf1, vcf2)
