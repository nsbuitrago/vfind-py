from Bio.Seq import Seq
from Bio import pairwise2
import polars as pl
import pyfastx
from collections import defaultdict
import time


def revcomp(seq: str) -> str:
    """
    Take reverse complement of a DNA sequence

    :param seq: DNA sequence (str)
    :return: Reverse complement of DNA sequence (str)
    """

    return seq[::-1].translate(str.maketrans("ATCG", "TAGC"))


class Adapters:
    """
    Adapter sequences flanking variable region

    :attr left: Left adapter sequence
    :attr right: Right adapter sequence
    """

    def __init__(self, left: str, right: str):
        self.left = left.upper()
        self.right = right.upper()


def get_variable_region(
    seq: str, adapters: Adapters,
    max_left_align_score: float, max_right_align_score: float,
    match_score: int, mismatch_score: int,
    gap_open_penalty: int, gap_extend_penalty: int
) -> str:

    start = seq.find(adapters.left)
    end = seq.find(adapters.right)

    if start != -1 and end != -1:
        return seq[start + len(adapters.left):end]

    else:
        left_alignment = pairwise2.align.localms(seq, adapters.left, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty, one_alignment_only=True)
        right_alignment = pairwise2.align.localms(seq, adapters.right, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty, one_alignment_only=True)

        if left_alignment[0][2] > max_left_align_score and right_alignment[0][2] > max_right_align_score:
            variable_region_start = left_alignment[0][4]
            variable_region_end = right_alignment[0][3]
            if variable_region_start < variable_region_end:
                return seq[variable_region_start:variable_region_end]

def variant_qc(variant: str, length: int, scheme: str) -> bool:
    """
    Quality control for variant sequences

    :param variant: Variant sequence
    :return: True if variant passes QC, False otherwise
    """
    if variant:
        if len(variant) == length:
            if scheme == "NNK":
                for i in range(2, len(variant)-2, 3):
                    last_base = variant[i]
                    if last_base in ["G", "T"]:
                        return True
    return False

def find_variants(
    adapters: Adapters,
    fq_1: str,
    fq_2: str = None,
    save_path: str = None,
    match_score: int = 3,
    mismatch_score: int = -2,
    gap_open_penalty: int = -5,
    gap_extend_penalty: int = -2,
    align_score_threshold: float = 0.75,
    length: int = 24,
    scheme: str = "NNK",
) -> pl.DataFrame:
    """
    Find variants flanked by adapter sequences

    :param adapters: Adapters object containing left and right flanking sequences
    :param fq_1: Path to forward fastq file
    :param fq_2: Path to reverse fastq file in the case of paired-end reads (Default: None)

    :return: None
    """
    rc_adapters = Adapters(revcomp(adapters.left), revcomp(adapters.right))
    peptides = defaultdict(int)
    max_left_align_score = align_score_threshold * len(adapters.left) * match_score
    max_right_align_score = align_score_threshold * len(adapters.right) * match_score

    for name, seq, qual in pyfastx.Fastq(fq_1, build_index=False):
        variable_region = get_variable_region(seq, adapters, max_left_align_score, max_right_align_score, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty) 
        if variant_qc(variable_region, length, scheme):
            peptide = str(Seq(variable_region).translate())
            print(peptide)
            peptides[peptide] += 1

    if fq_2:
        for name, seq, qual in pyfastx.Fastq(fq_2, build_index=False):
            variable_region_rc = get_variable_region(seq, rc_adapters, max_left_align_score, max_right_align_score, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty) 
            if variant_qc(variable_region_rc, length, scheme):
                variable_region = revcomp(variable_region_rc)
                peptide = str(Seq(variable_region).translate())
                print(peptide)
                peptides[peptide] += 1

    peptides_df = pl.from_dict({"peptide": list(peptides.keys()), "count": list(peptides.values())})

    if save_path:
        peptides_df.write_csv(save_path)

    return peptides_df

if __name__ == "__main__":
    adapters = Adapters("GATTATGCTGGGGCCCAGCCGGCCGGATCC", "GGAGGCGGAGGTTCAGGAGGAGGGGGATCG" )
    fq_1 = "30-982598147/17L_R1_clean.fq.gz"
    fq_2 = "30-982598147/17L_R2_clean.fq.gz"
    
    start = time.perf_counter()
    find_variants(
        adapters,
        fq_1,
        fq_2,
        save_path="30-982598147/17L_variants.csv"
    )
    end = time.perf_counter()
    print(f"Time: {end - start}")

