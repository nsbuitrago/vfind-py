from Bio.Seq import Seq
from Bio import pairwise2
import polars as pl
import pyfastx
from collections import defaultdict
import time

class Adapters:
    """
    Adapter sequences flanking variable region

    :attr left (str): left adapter sequence
    :attr right (str):  right adapter sequence
    :attr id (str, optional): identifiers for adapter pair
    """

    def __init__(self, left: str = None, right: str = None, id: str = None):
        self.adapters = [(left, right)] if left and right else []
        self.id = [id] if id else []

    def add(self, left: str, right: str, id: str = None):
        """
        Add adapter pairs to Adapters object

        :param left (str): Left adapter sequence
        :param right (str): Right adapter sequence
        :param id (str, optional): Identifier for adapter pair
        """
        self.adapters.append((left, right))
        self.id.append(id)


class AlignParams:
    """
    Parameters for pairwise alignment

    :attr match_score (int): Score for match
    :attr mismatch_score (int): Score for mismatch
    :attr gap_open_penalty (int): Penalty for gap opening
    :attr gap_extend_penalty (int): Penalty for gap extension
    :attr left_align_threshold (float): Threshold for accepting left alignments
    :attr right_align_threshold (float): Threshold for accepting right alignments 
    """

    def __init__(
        self,
        match_score: int = 3,
        mismatch_score: int = -2,
        gap_open_penalty: int = -5,
        gap_extend_penalty: int = -2,
        left_align_threshold: float = 0.75,
        right_align_threshold: float = 0.75
    ):
        self.match_score = match_score
        self.mismatch_score = mismatch_score
        self.gap_open_penalty = gap_open_penalty
        self.gap_extend_penalty = gap_extend_penalty
        self.left_align_threshold = left_align_threshold
        self.right_align_threshold = right_align_threshold

    def __post_init__(self):
        if self.match_score <= 0:
            raise ValueError("Match score must be positive")
        if self.mismatch_score >= 0:
            raise ValueError("Mismatch score must be negative")
        if self.gap_open_penalty >= 0:
            raise ValueError("Gap open penalty must be negative")
        if self.gap_extend_penalty >= 0:
            raise ValueError("Gap extend penalty must be negative")
        if self.left_align_threshold < 0 or self.left_align_threshold > 1:
            raise ValueError("Left alignment threshold must be between 0 and 1")
        if self.right_align_threshold < 0 or self.right_align_threshold > 1:
            raise ValueError("Right alignment threshold must be between 0 and 1")


def get_variable_region(
        seq: str, left_adapter: str, right_adapter: str, max_left_align_score: float, max_right_align_score: float, align_params: AlignParams) -> str:

    start = seq.find(left_adapter)
    end = seq.find(right_adapter)

    if start != -1 and end != -1:
        return seq[start + len(left_adapter):end]
    elif start == -1:
        left_alignment = pairwise2.align.localms(seq, left_adapter, align_params.match_score, align_params.mismatch_score, align_params.gap_open_penalty, align_params.gap_extend_penalty, one_alignment_only=True)
        if left_alignment[0][2] > max_left_align_score:
            start = left_alignment[0][4]
    elif end == -1:
        right_alignment = pairwise2.align.localms(seq, right_adapter, align_params.match_score, align_params.mismatch_score, align_params.gap_open_penalty, align_params.gap_extend_penalty, one_alignment_only=True)
        if right_alignment[0][2] > max_right_align_score:
            end = right_alignment[0][3]

    if start != -1 and end != -1 and start < end:
        return seq[start:end]


def find_variants(
    adapters: Adapters,
    fq_path: str,
    save_path: str = None,
    align_params: AlignParams = AlignParams(),
    skip_translation: bool = False
) -> pl.DataFrame:
    """
    Find variants flanked by adapter sequences

    :param adapters: Adapters object containing left and right flanking sequences
    :param fq_1: Path to forward fastq file
    :param fq_2: Path to reverse fastq file in the case of paired-end reads (Default: None)

    :return: None
    """
    variants = defaultdict(int)

    for name, seq, qual in pyfastx.Fastq(fq_path, build_index=False):
        for left, right in adapters.adapters:
            max_left_align_score = align_params.left_align_threshold * len(left) * align_params.match_score
            max_right_align_score = align_params.right_align_threshold * len(right) * align_params.match_score

            variant = get_variable_region(seq.upper(), left, right, max_left_align_score, max_right_align_score, align_params)
            if variant:
                if skip_translation:
                    variants[variant] += 1
                else:
                    variants[str(Seq(variant).translate())] += 1

    variants_df = pl.from_dict({"seq": list(variants.keys()), "count": list(variants.values())})

    if save_path:
        variants_df.write_csv(save_path)

    return variants_df

if __name__ == "__main__":
    adapters = Adapters("GATTATGCTGGGGCCCAGCCGGCCGGATCC", "GGAGGCGGAGGTTCAGGAGGAGGGGGATCG", id="17L")
    fq_1 = "tests/toy_R1.fq.gz"
    
    start = time.perf_counter()
    find_variants(
        adapters,
        fq_1,
        save_path="new_test.csv"
    )
    end = time.perf_counter()
    print(f"Time: {end - start}")

