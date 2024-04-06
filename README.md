# vfind-py

A simple variant finder. 

## Examples

```python
from vfind import find_variants, Adapters, AlignParams
import polars as pl # variants are returned in a polars dataframe

# we want to find all variable regions that are flanked by constant adapters

# define the Adapters
adapters = Adapters(
    left="GGG",
    right="CCC"
    id="lib_1" # id is optional
)

# you can also set additional adapters with the `add` method
adapters.add("AAA", "TTT", id="lib_2")

# now, let's find some variants
fq_path = "./path/to/your/fastq/file.fq.gz" # we support gzipped fastq files
variants = find_variants(fq_path, adapters)

# variants is a polars dataframe with 'seq' and 'count' columns.
# 'seq' contains the amino acid sequence of the variable region
# 'count' contains the number of times the sequence was found in the fastq file

# if you only want the DNA sequence, you can set `skip_translation` to True
variants = find_variants(fq_path, adapters, skip_translation=True)

# if you want to adjust the alignment parameters, use the AlignParams class
align_params = AlignParams(
    match=1,
    mismatch=-1,
    gap_open=-1,
    gap_extend=-1
)

variants = find_variants(fq_path, adapters, align_params=align_params)

# you can also set how stringent you want to be when accepting alignments.
# the AlignParams class takes a left_align_threshold and right_align_threshold.
# these are floats between 0 and 1 that represent the threshold for accepting an
# alignment. The default is 0.8 for both. If you set them to 1, you will only
# accept perfect alignments (i.e., exact matches to the adapter sequence).
# these thresholds are actually a percentage of the maximum alignment score.
# So, if you want to take sequences that are 90% an exact match based on your
# alignment parameters, then you would set the threshold to 0.9.

align_params = AlignParams(
    left_align_threshold=0.9,
    right_align_threshold=0.9
)

variants = find_variants(fq_path, adapters, align_params=align_params)
```

## Installation

I have not packaged this, so you will need to clone the repo or just the `vfind.py`
file into your project directory to use it. We have a new version that is significantly
faster, but is not quite ready for release. This new version will be packaged
when ready for ease of use.

## Contribution

Contributions are more than welcome. We are currently working on a new version
that will likely go into a new repo. If you would like to contribute to this
version, please open an issue or PR.

## License

vfind-py is licensed under the MIT license.
