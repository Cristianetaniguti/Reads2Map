from pathlib import Path

import pytest


@pytest.mark.workflow("Filter non-informative in VCF file ")
def test_ensure_vcf_1_are_equal(workflow_dir, compare_files):
    vcf1 = Path(workflow_dir, "_LAST/out/vcf_filtered/filtered.vcf.gz")
    vcf2 = Path("tests/data/Ptremula_PRJNA395596_subset/RemoveNonInformative_results/filtered_gatk_bam.vcf.gz")

    # Asserts
    assert compare_files.assert_vcf_equal(vcf1, vcf2)
