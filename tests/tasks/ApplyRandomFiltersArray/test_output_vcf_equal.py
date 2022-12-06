from pathlib import Path

import pytest


@pytest.mark.workflow("Use bcftools to filter variants in VCF file")
def test_ensure_vcf_1_are_equal(workflow_dir, compare_files):
    vcf1 = Path(workflow_dir, "_LAST/out/vcfs_filt/0/vcf_filt_gatk_vcf.vcf.gz")
    vcf2 = Path("tests/data/polyploid/ApplyRandomFiltersArray/vcf_filt_gatk_vcf.vcf.gz")

    # Asserts
    assert compare_files.assert_vcf_equal(vcf1, vcf2)
