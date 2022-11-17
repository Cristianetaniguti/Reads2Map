# from pathlib import Path

# import pytest

# from utils.vcf import assert_vcf_equal


# @pytest.mark.workflow("Use bcftools view to create bgzip and tabixed vcfs from the provided VCF inputs")
# def test_ensure_vcf_1_are_equal(workflow_dir):
#     vcf1 = Path(workflow_dir, "_LAST/out/vcf_files/0/TEST-ONLY.another-test.conversion.vcf.gz")
#     vcf2 = Path("cached-data/circleci-test-data/utils_tests/expected_left_align_vcf.vcf.gz")

#     # Asserts
#     assert assert_vcf_equal(vcf1, vcf2)


# @pytest.mark.workflow("Use bcftools view to create bgzip and tabixed vcfs from the provided VCF inputs")
# def test_ensure_vcf_2_are_equal(workflow_dir):
#     vcf1 = Path(workflow_dir, "_LAST/out/vcf_files/1/TEST-ONLY.testing.conversion.vcf.gz")
#     vcf2 = Path("cached-data/circleci-test-data/utils_tests/mito_variants_other.vcf.gz")

#     # Asserts
#     assert assert_vcf_equal(vcf1, vcf2)
