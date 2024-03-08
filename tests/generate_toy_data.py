from Bio import pairwise2
from Bio.Seq import Seq
from collections import defaultdict
import polars as pl
import random

class Adapters:
    """
    Adapter sequences flanking variable region

    :attr left: Left adapter sequence
    :attr right: Right adapter sequence
    """

    def __init__(self, left: str, right: str):
        self.left = left.upper()
        self.right = right.upper()

def generate_toy_data(
    adapters: Adapters,
    lib_size: int = 1000,
    match_score: int = 3,
    mismatch_score: int = -2,
    gap_open_penalty: int = -5,
    gap_extend_penalty: int = -3,
    accept_alignment: float = 0.75,
    mut_freq: int = 1,
    adapter_id: str = "fwd",
    seq_len: int = 24,
    reverse: bool = False
) -> (list, dict):

    fastq = []
    recovered_peptides = defaultdict(int)

    for i in range(lib_size):
        if adapter_id == "fwd":
            adapter = adapters.left
        else:
            adapter = adapters.right

        for m in range(mut_freq):
            idx = random.randint(0, len(adapter)-1)
            mut_adapter = adapter[:idx] + random.choice("ACGT") + adapter[idx+1:]

        rand_seq = "".join(random.choices("ACGT", k=seq_len))
        peptide = Seq(rand_seq).translate()

        if adapter_id == "fwd":
            rand_seq = mut_adapter + rand_seq + adapters.right
        else:
            rand_seq = adapters.left + rand_seq + mut_adapter

        alignment = pairwise2.align.localms(
            rand_seq, adapter, match_score, mismatch_score, gap_open_penalty, gap_extend_penalty, one_alignment_only=True
        )

        max_alignment_score = match_score * len(adapter) * accept_alignment
        if alignment[0][2] > max_alignment_score:
            recovered_peptides[str(peptide)] += 1

        if reverse:
            rand_seq = Seq(rand_seq).reverse_complement()

        fastq.append(rand_seq)
    
    return fastq, recovered_peptides

if __name__ == "__main__":
    variants_path = "variants_ground_truth.csv",
    fastq_fwd_path = "toy_R1.fq",
    fastq_rev_path = "toy_R2.fq"

    fastq_fwd = open("toy_R1.fq", "w+")
    fastq_rev = open("toy_R2.fq", "w+")

    adapters = Adapters("GATTATGCTGGGGCCCAGCCGGCCGGATCC", "GGAGGCGGAGGTTCAGGAGGAGGGGGATCG")
    seq_len = 24
    mut_freqs = range(1, seq_len, 4)
    #peptides_df = pl.DataFrame(columns=["peptide", "count"])
    peptides_df = pl.DataFrame()

    for mut_freq in mut_freqs:
        fwd_seqs, peptides_R1 = generate_toy_data(adapters, lib_size=1000, mut_freq=mut_freq, seq_len=seq_len)
        rev_seqs, peptides_R2 = generate_toy_data(adapters, lib_size=1000, mut_freq=mut_freq, seq_len=seq_len, reverse=True)
        for i, seq in enumerate(fwd_seqs):
            fastq_fwd.write(f"@seq_{i}\n{seq}\n+\n{'F'*len(seq)}\n")

        for i, seq in enumerate(rev_seqs):
            fastq_rev.write(f"@seq_{i}\n{seq}\n+\n{'F'*len(seq)}\n")

        peptides_df = peptides_df.vstack(pl.DataFrame({"peptide": list(peptides_R1.keys()), "count": list(peptides_R1.values())}))
        peptides_df = peptides_df.vstack(pl.DataFrame({"peptide": list(peptides_R2.keys()), "count": list(peptides_R2.values())}))

    fastq_fwd.close()
    fastq_rev.close()

    breakpoint()
    peptides_df.write_csv(variants_path)

