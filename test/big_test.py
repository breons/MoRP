#=============================================================================
#
#   MOR Big Test
#
#	Input: Big matrix of randomly generated numbers m
#	Test: Profile algorithm
#
#=============================================================================

import pandas as pd
import numpy as np
import cProfile

from context import morp 

''' --------------------------------------------------------------------------
Run Test
---------------------------------------------------------------------------'''

print("\nLoading Counts\n------------------------")
counts = pd.read_csv("randomcounts/test_big.csv", index_col=0).transpose()

print("Rows: ", counts.shape[0], "Columns: ", counts.shape[1], "\n")

print("Start Profiling\n------------------------")
profile = cProfile.run('morp.MORP(counts)')



