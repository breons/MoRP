#=============================================================================
#
#   MOR Small Test
#
#	Input: Small matrix of randomly generated numbers m
#	Output: Small matrix mor normalised. 
#	Test: Matrix must equal hand calculations
#
#=============================================================================

import pandas as pd
import numpy as np

from context import morp 

''' --------------------------------------------------------------------------
Run Test
---------------------------------------------------------------------------'''

counts = pd.read_csv("randomcounts/test_sm.csv", index_col=0).transpose()
answers = pd.read_csv("randomcounts/answer_sm.csv", index_col=0)

normalised = morp.MORP(counts)

print("\nRaw Counts\n------------------------")
print(counts.round(2))

print("\nMOR attempt\n------------------------")
print(normalised.mor.round(2))

print("\nMOR manual calc\n------------------------")
print(answers.round(2))

print("\nOutcome\n------------------------")
if answers.round(2).equals(normalised.mor.round(2)):
	print("Test passed.\n")
else:
	print("Test failed.\n")