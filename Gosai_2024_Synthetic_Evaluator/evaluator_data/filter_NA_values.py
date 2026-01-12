#script to make downsampled data from original
import json
import pandas as pd

#Read in measurement file from paper
measured_values = pd.read_csv("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/game_dev/Evaluators/Gosai_2024_Synthetic_Evaluator/evaluator_data/41586_2024_8070_MOESM14_ESM.txt", sep = '\t')
print(measured_values)
#[77157 rows x 16 columns]
cols = ["K562_l2fc", "HepG2_l2fc", "SKNSH_l2fc"]
#Drop rows where the measurments are NA for all the cell types

measured_values_clean = measured_values.dropna(subset=cols, how="all")
print(measured_values_clean)
measured_values_clean.to_csv("/scratch/st-cdeboer-1/iluthra/game_apis/RestAPI/game_dev/Evaluators/Gosai_2024_Synthetic_Evaluator/evaluator_data/41586_2024_8070_MOESM14_ESM_clean.txt", sep = "\t")
#[76193 rows x 16 columns]