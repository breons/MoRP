# MoRP
Python implementation of Median of Ratio's normalisation method.

## Props to:
Anders, S. and Huber, W. (2010) ‘Differential expression analysis for sequence count data’, GENOME BIOLOGY, 11(10). doi: 10.1186/gb-2010-11-10-r106.

Josh Starmer - https://www.youtube.com/watch?v=UFB993xufUU

## Use it!

1. Clone it somewhere, navigate to root directory.
2. Install it `pip install .`
3. `from morp import morp`
4. `normalised = morp.MORP(your_counts_matrix)` Note: Rows = Genes, Columns = Samples.
5. `normalised.mor` for normalised counts. `normalised.gm` for geometric means used in calculation.
