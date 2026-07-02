'''Calculate and save the final evaluation metrics.

This module computes Pearson r between the Predictor's expression predictions and
the measured log2(RNA/DNA) values from the Gosai 2024 synthetic MPRA dataset, for each
of the three cell types (K562, HepG2, SK-N-SH).

Cell-type specificity is NOT computed for Gosai — the dataset's measurements are not
made across all three cell types in a way that supports paired-difference correlation.

NOTE: Every evaluator will do this slightly differently depending on how the data is presented.
'''

import os
import sys
import json
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
from datetime import datetime, timezone

from config import EVALUATOR_NAME, EVALUATOR_INPUT_PATH

def _calculate_task_correlation(
    measured_df: pd.DataFrame,
    single_task_data: dict,
    measured_value_column: str,
    seq_id_column: str,
    chromosome_column: str = None,
    chromosomes_to_filter: list = None
    ):
    
    """
    Calculates Pearson r and extracts metadata for single prediction task,
    using pre-loaded measured_df and a single task data dictionary.
    
    Args:
        measured_df (pd.DataFrame): The dataframe of measured values loaded once by the caller (Evaluator __main__ block).
        single_task_data (dict): The dictionary with predictions and metadata for a single task from the `prediction_tasks` 
                                 list of a predictions JSON.
        measured_value_column (str): Column name in measured_df for correlation.
        id_column (str): Common identifier for sequences column (IDs or sequences).
        chromosome_column (str, optional): Chromosome column name in measured_df. Defaults to 'None'.
        chromosomes_to_filter (list, optional): List of chromosomes for filtering. Defaults to 'None'. 

    Returns:
        correlation_details (dict): Dictionary containing 'pearson_r' (float or None) and task metadata ("task_name", "task_type", "cell_type_actual")
    """
    
    # Extract metadata from single task data
    print(f"\n--- Extracting prediction_task metadata ---")
    task_name = single_task_data.get("name")
    task_type_actual = single_task_data.get("type_actual")
    cell_type_actual = single_task_data.get("cell_type_actual")
    predictions_dict =  single_task_data.get("predictions")
    scale_prediction_actual = single_task_data.get("scale_prediction_actual")
    scale_prediction_requested = single_task_data.get("scale_prediction_requested")
    
    pearson_r_value = None # If there's an error, default to None
    
    if "error" in predictions_dict:
        print("No predictions were returned for this task -> Skipping evaluation calculation")
        correlation_details = {
        'task_name': task_name, 
        'task_type': task_type_actual,
        'cell_type_actual': cell_type_actual,
        'pearson_r': pearson_r_value
        }
        return correlation_details

    if scale_prediction_actual != scale_prediction_requested:
        print("Predictions scale does not match requested scale (log): Skipping evaluation calculation")
        correlation_details = {
        'task_name': task_name, 
        'task_type': task_type_actual,
        'cell_type_actual': cell_type_actual,
        'pearson_r': pearson_r_value
        }
        return correlation_details
    
    # --- Data Validation and processing ---
    print("--- Validating data ---")
    if not isinstance(predictions_dict, dict):
        print(f"WARNING: 'predictions' in task: '{task_name}'\
            \nCell type: {cell_type_actual}\
            \nType: {task_type_actual}\
            \nis not a valid dictionary or is missing.")
    elif not predictions_dict:
        print(f"WARNING: 'predictions' dictionary is empty in task: '{task_name}'\
            \nCell type: {cell_type_actual}\
            \nType: {task_type_actual}")
    elif seq_id_column not in measured_df.columns:
        print(f"ERROR: Sequence ID column '{seq_id_column}' not found in measured_df.\
            \nCannot merge for task: {task_name}.")
    elif measured_value_column not in measured_df.columns:
         print(f"ERROR: Measured value column '{measured_value_column}'\
             \nfor cell type '{cell_type_actual}' not found in measured_df.\
             \nCannot correlate task: '{task_name}'.") 
         # NOTE: More checks can be added.
    else:
        # Proceed with calculation if checks pass
        # Create DataFrame from Predictions
        print("--- Creating predictions_df ---")
        predictions_df = pd.DataFrame(list(predictions_dict.items()), columns=[seq_id_column, 'Predicted_Value'])
        predictions_df['Predicted_Value'] = predictions_df['Predicted_Value'].apply(
            lambda x: x[0] if isinstance(x, list) and len(x) > 0 else x
            )

        #check here is there is NA is any of the prediction values
        na_rows = predictions_df[predictions_df['Predicted_Value'].isna()]
        if not na_rows.empty:
            print("NA values were found in the predictions, skipping evaluation")
            print(na_rows)
            correlation_details = {
            'task_name': task_name, 
            'task_type': task_type_actual,
            'cell_type_actual': cell_type_actual,
            'pearson_r': pearson_r_value
            }
            return correlation_details
    
        # Now select only the necessary columns of measured_df
        columns_to_keep = [seq_id_column, measured_value_column]
        if (chromosome_column and chromosomes_to_filter and (chromosome_column in measured_df)):
            if chromosome_column not in columns_to_keep:
                columns_to_keep.append(chromosome_column)
                
        measured_df_subset = measured_df[columns_to_keep].copy()
        merged_df = pd.merge(
            measured_df_subset, predictions_df,
            on=seq_id_column, how="inner"
        )

        # Filter by chromosome (if needed)
        if (chromosomes_to_filter and chromosome_column and (chromosome_column in merged_df.columns)):
            print(f"Filtering chromosomes: {chromosomes_to_filter}...")
            merged_df[chromosome_column] = merged_df[chromosome_column].astype(str)
            # In order to handle None
            str_chromosomes_to_filter = [str(c) for c in chromosomes_to_filter]
            filtered_df = merged_df[merged_df[chromosome_column].isin(str_chromosomes_to_filter)]
        else:
            print("No chromosomes to filter.")
            filtered_df = merged_df # No chromosome filter
        
        # Columns for correlation calculation
        # Drop any rows that have NaNs for either column
        correlation_columns = [measured_value_column]
        print("Original size of the measurement file is:")
        print(filtered_df.shape)
        na_rows = filtered_df[filtered_df[correlation_columns].isna().any(axis=1)]
        if not na_rows.empty:
            print("Rows with NaN values in any of measured value column (will be dropped):")
            print(na_rows)
            
        final_df = filtered_df.dropna(subset=correlation_columns)
        print("Size of data after dropping rows with NaN measured values")
        print(final_df.shape)

        # Sanitize the final_df in case values are non-numeric
        if not final_df.empty:
            print("Sanitizing final_df in case values are non-numeric for correlation...")
            final_df.loc[:, 'Predicted_Value'] = pd.to_numeric(final_df['Predicted_Value'], errors='coerce')
            final_df.loc[:, measured_value_column] = pd.to_numeric(final_df[measured_value_column], errors='coerce')
            
            # Drop rows that became NaN after numeric coercion, with a count for visibility
            rows_before = len(final_df)
            final_df = final_df.dropna(subset=['Predicted_Value', measured_value_column])
            rows_dropped = rows_before - len(final_df)
            if rows_dropped > 0:
                print(f"Warning: Dropped {rows_dropped} rows with non-numeric values after coercion.")
            print("Size of data that will be used to calculate Pearson r")
            print(final_df.shape)
            
            # Check for 0 variance in either column (all values identical)
            std_predicted = final_df['Predicted_Value'].std()
            std_measured = final_df[measured_value_column].std()

            if std_predicted == 0 or std_measured == 0:
                print(f"WARNING: Zero variance detected for task '{task_name}'. Setting Pearson r to 0 to reflect no correlation.")
                pearson_r_value = 0.0
            
            else:
                # Calculate pearson r
                try:
                    r, _ = pearsonr(final_df['Predicted_Value'], final_df[measured_value_column])
                    pearson_r_value = 0.0 if np.isnan(r) else float(r)
                    print(f"Calculated Pearson r for {task_name}: {pearson_r_value}")
                except Exception as e:
                    print(f"Error during Pearson correlation calculation for task '{task_name}': {e}")
            
        else:
            print(f"DataFrame is empty after numeric conversion and NaN drop for task: '{task_name}'")
            
    correlation_details = {
        'task_name': task_name, 
        'task_type': task_type_actual,
        'cell_type_actual': cell_type_actual,
        'pearson_r': pearson_r_value
    }
    return correlation_details

def calculate_and_save_metrics(saved_predictions_path, output_dir):
    """
    Calculates custom evaluation metrics and saves them to CSV files.
    This is the primary function to customize for a new evaluator.
    """
    print("----- Starting Evaluation Calculation and Saving as CSV -----")
    MEASURED_DATA_PATH = EVALUATOR_INPUT_PATH # NOTE: This may not be the same for other evaluators
    print(f"Using measured data from: {MEASURED_DATA_PATH}")
    print(f"Using predictions from: {saved_predictions_path}")
    print(f"Correlation metadata will be saved in {output_dir}")

    seq_column = "ID" # This can change depending on data
    measured_value_columns_map = {
        "K562 (erythroid precursors)": "K562_l2fc",
        "HepG2 (hepatocytes)": "HepG2_l2fc", 
        "SK-N-SH (neuroblastoma)": "SKNSH_l2fc"
    }
    # Define output paths
    evaluation_metrics_filename = f"evaluation_summary_{EVALUATOR_NAME}.csv"
    evaluation_metrics_filepath = os.path.join(output_dir, evaluation_metrics_filename)
    
    # Initialize an empty list to get summary for all tasks
    all_task_correlation_results = []
    
    try:
        try:
            # Load measured data file and predictions file ONCE (not with every function call).
            # NOTE: Evaluator builders: If measured_file_path is not a tab-separated file,
            # this line (pd.read_csv) will need to be adjusted or replaced with the
            # appropriate pandas read function (e.g., pd.read_excel, pd.read_csv with different sep)
            # or custom loading logic (e.g., for .npy files).
            measured_df = pd.read_csv(MEASURED_DATA_PATH, sep='\t', header=0)
            print(measured_df)
            
            print("\nFiltering to synthetic sequences only (matching sequences sent to Predictor)...")
            print(f"Original measured_df shape: {measured_df.shape}")
            measured_df = measured_df[measured_df['origin'].isin(["Simulated_Annealing", "FastSeqProp", "AdaLead"])]
            print(f"Filtered measured_df shape: {measured_df.shape}")
            print("This ensures we only evaluate on sequences that were sent for prediction.\n")

            # Now load predictions
            with open(saved_predictions_path, 'r') as f:
                predictions_file_content = json.load(f)
            
            # Extract Predictor Name
            predictor_name_base = predictions_file_content.get("predictor_name", "UnknownPredictor")
            predictor_name = predictor_name_base.replace(" ", "_").replace("/", "_")
        
        except Exception as e:
            print(f"Error loading data files: {e}", file=sys.stderr)
            raise
            
        try:
            if (
                "prediction_tasks" not in predictions_file_content or
                # Also flag cases in case prediction_tasks key is returned empty
                not predictions_file_content["prediction_tasks"] or
                # And flag if any 'predictions' keys are empty
                any(not key.get("predictions") for key in predictions_file_content["prediction_tasks"])
            ):
                print("WARNING: 'prediction_tasks' key missing, empty, or one of the tasks has empty predictions.")
            else:
                # Loop through each prediction_task from Predictor
                # Calculate the correlation for each task separately
                for task_index, single_task_data_dict in enumerate(predictions_file_content["prediction_tasks"]):
                    if not isinstance(single_task_data_dict, dict):
                        print(f"WARNING: Task item at index {task_index} is not a dictionary. Skipping!")
                        continue
                    
                    # Extract metadata from this task
                    task_type_actual = single_task_data_dict.get("type_actual")
                    predicted_cell_type = single_task_data_dict.get("cell_type_actual")
                    # We also want to extract the cell_type_requested to map it to measured_value_columns_map
                    requested_cell_type = single_task_data_dict.get("cell_type_requested")
                    
                    # Find the corresponding measured data column from the map
                    measured_col_for_task = measured_value_columns_map.get(requested_cell_type)
                    
                    print(f"\nProcessing task {task_index+1} (Cell Type: {predicted_cell_type}). Correlating against measured column '{measured_col_for_task}'\
                    for requested cell type '{requested_cell_type}'.")
                    prediction_task_data_nopredictions = [{k: v for k, v in single_task_data_dict.items() if k != "predictions"}]
                    # Call the correlation calculation function
                    task_correlation_dict = _calculate_task_correlation(
                        measured_df=measured_df,
                        single_task_data=single_task_data_dict,
                        measured_value_column=measured_col_for_task,
                        seq_id_column=seq_column # chromosome_column and chromosomes_to_filter_list can be added as arguments
                    )
                    
                    if task_correlation_dict:
                        pearson_r_value = task_correlation_dict.get('pearson_r')
                        val_str = "NaN" if pearson_r_value is None else str(pearson_r_value)
                        
                        # Get UTC timestamp for predictor_name
                        timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S.%f")

                        description = f"Gosai Synthetic MPRA ({requested_cell_type})"
                        all_task_correlation_results.append({
                            'evaluator_name': EVALUATOR_NAME,
                            'description': description,
                            'predictor_name': predictor_name,
                            'time_stamp': timestamp,
                            'metric': 'pearson_r',
                            'value': val_str,
                            'prediction_task(s)_data': prediction_task_data_nopredictions,
                        })

        except Exception as e:
            print(f"An error occurred during correlation calculation: {e}")
            
       # Once all the data is received, save them all into a summary CSV
        # print(all_task_correlation_results)
        if all_task_correlation_results:
            summary_df = pd.DataFrame(all_task_correlation_results)
            csv_file_exists: bool = os.path.isfile(evaluation_metrics_filepath)
            try:
                summary_df.to_csv(evaluation_metrics_filepath, mode='a',
                                  sep='\t', header=(not csv_file_exists), index=False)
                if csv_file_exists:
                    print("Appended to existing summary CSV file")
                else:
                    print("Created a new summary CSV file")
                print(f"Saved correlation summary to {evaluation_metrics_filepath}!")
            except IOError as e:
                print("\nNo correlation results were saved!")

    except Exception as e:
        print(f"An unexpected error occurred during evaluation calculations: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()

