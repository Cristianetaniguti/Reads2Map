from pathlib import Path

import pandas as pd
from pytest import fixture


class CompareFiles:
    @staticmethod
    def assert_vcf_equal(vcf_a: Path, vcf_b: Path):
        df_a = pd.read_csv(vcf_a, comment="#", sep="\t", header=None)
        df_b = pd.read_csv(vcf_b, comment="#", sep="\t", header=None)

        if df_a.equals(df_b):
            return True

        print(df_a)
        print(df_b)
        return False


@fixture
def compare_files() -> CompareFiles:
    return CompareFiles()
