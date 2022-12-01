from pathlib import Path

import pytest

# Parei aqui! Adaptar os caminhos
@pytest.mark.workflow("Use updog polyRAD or SuperMASSA to call dosages")
def test_ensure_vcf_1_are_equal(workflow_dir, compare_files):
    vcf1 = Path(workflow_dir, "_LAST/out/regeno_vcf/regeno.vcf")
    vcf2 = Path("tests/data/Ptremula_PRJNA395596_subset/genotyping_vcfs/freebayes/updog/regeno.vcf")

    # Asserts
    assert compare_files.assert_vcf_equal(vcf1, vcf2)
