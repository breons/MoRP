#=============================================================================
#
#   Median of Ratios - Python | RNA-Seq normalisation
#
#   Anders, S. and Huber, W. (2010) 
#.  ‘Differential expression analysis for sequence count data’, 
#   GENOME BIOLOGY, 11(10)
#
#   License: GNU LGPLv3
#
#   Input: Pandas Dataframe of counts (samples x genes)
#   Output: Normalised pandas dataframe and Geometric means
#
#=============================================================================

import numpy as np
import pandas as pd


class MORP():

    def __init__(self, counts):
        self.mor, self.gm = self._mor(counts)


    def _mor(self, counts):
        
        # Step 1 - Take the log of every count, note np.log is the natural log
        step_1_logs = counts.apply(np.log)
        
        # Step 2 - Take the average of each gene, that is, each row.
        step_2_geometric_mean = step_1_logs.mean(axis=1)
        
        # Step 3 - Remove genes with a -infinity geometric mean
        step_2_geometric_mean = step_2_geometric_mean.replace(-np.inf, np.nan)
        step_2_geometric_mean = step_2_geometric_mean.dropna()
        
        geometric_means = {}

        for gene, g_mean in step_2_geometric_mean.iteritems():
            geometric_means[gene] = g_mean  
            
        step_3_logs = step_1_logs[step_1_logs.index.isin(geometric_means.keys())]       
            
        # Transform 
        counts_normalised, success = self._applyMor(counts, 
                                                    step_3_logs, 
                                                    geometric_means)
        
        # Finish up, report to user.
        if success:
            print("Normalisation successful.")
        else:
            print("Error.")

        return counts_normalised, geometric_means 


    def _applyMor(self, counts, step_3_logs, geometric_means):

        # Step 4 - For each of these remaining genes (We assume housekeeping), 
        # subtract them from log matrix (step 1)

        step_3_logs_t = step_3_logs.transpose() # Because of pandas stupid slice copy thing

        for sample, gene in step_3_logs.iterrows():
            try:
                step_3_logs_t[sample] -= geometric_means[sample] 
            except KeyError:
                continue 

        step_3_logs = step_3_logs_t.transpose()

        # Step 5 - Perform median of ratios
        sample_scaler_log = step_3_logs.median(axis=0)
        
        # Step 6 - Convert the scalers to counts
        sample_scaler = sample_scaler_log.apply(np.exp)
        gene_names = counts.index
        sample_names = counts.columns
        sample_scaler_matrix = pd.np.tile(sample_scaler, (len(counts.index), 1))
        
        # Step 7, divide each read count per sample in the raw matrix!
        counts_normalised = counts.copy().values
        counts_normalised = counts_normalised.astype("float")
        np.divide(counts_normalised, 
                  sample_scaler_matrix, 
                  counts_normalised, 
                  casting="unsafe")

        # Format output
        sample_scaler_matrix = pd.DataFrame(sample_scaler_matrix, 
                                            index=gene_names, 
                                            columns=sample_names)

        counts_normalised = pd.DataFrame(counts_normalised, 
                                         index=gene_names, 
                                         columns=sample_names)

        # Sanity check
        sample_check = step_3_logs.columns[0]
        gene_check = step_3_logs.index[0]

        original = counts.loc[gene_check, sample_check]
        scaler = sample_scaler_matrix.loc[gene_check, sample_check]
        normalised = counts_normalised.loc[gene_check, sample_check]  
        success = False
        
        if (original/scaler) == normalised:
            success = True

        return counts_normalised, success


    def transformMor(self, counts):

        # Step 1 - Take the log of every count, note np.log is the natural log
        step_1_logs = counts.apply(np.log)
        step_3_logs = step_1_logs[step_1_logs.index.isin(self.gm.keys())]    
        
        # Transform 
        counts_normalised, success = self._applyMor(counts, step_3_logs, self.gm)
        
        return counts_normalised
