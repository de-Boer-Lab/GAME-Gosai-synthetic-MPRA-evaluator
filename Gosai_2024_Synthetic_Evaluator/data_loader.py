'''Handle Loading and Validating Evaluator Input/Request Data'''

import os
import json
from collections import Counter
import functools
import pandas as pd

from config import EVALUATOR_INPUT_PATH, PLASMID_BACKBONE_INPUT_PATH

class DuplicateKeysError(ValueError):
    """Raised when duplicate keys are found in a JSON object."""
    pass

# Internal helper function to detect duplicates during JSON parsing
def _detect_duplicates(pairs, duplicate_keys_state):

    """
    Detects duplicate keys during JSON parsing and counts occurrences of each key.

    This function intercepts the key-value pairs provided by `json.loads` and ensures that
    duplicate keys are flagged. It constructs the dictionary normally but counts how often
    each key appears, recording any keys that occur more than once.

    Args:
        pairs (list of tuple): A list of key-value pairs at the current level of the JSON.
        duplicate_keys_state (dict): The dictionary to update with any duplicates found.

    Returns:
        result_dict: A dictionary created from the key-value pairs.
    """

    # Use a local Counter to count occurrences of keys at this level
    local_counts = Counter()
    result_dict = {}
    for key, value in pairs:
        # Increment the count for each key
        local_counts[key] += 1
        # If the key is a duplicate, record it in the duplicate_keys dictionary
        if local_counts[key] > 1:
            duplicate_keys_state[key] = local_counts[key]
        # Add the key-value pair to the resulting dictionary
        result_dict[key] = value
    return result_dict

def _process_results(data, duplicate_keys):
    """
    Checks the duplicate_keys dictionary and prints a report.

    Args:
        data (dict): The dictionary of parsed data.
        duplicate_keys (dict): The dictionary of duplicates.

    Raises:
        DuplicateKeysError: If duplicate keys are found in the JSON structure.

    Returns:
        data: The parsed data if no errors or duplicates are found.
    """
    # Report duplicates if any were found
    if duplicate_keys:
        print("Duplicate keys found:")
        error_messages = [f"Key: '{key}', Count: {count}" for key, count in duplicate_keys.items()]
        raise DuplicateKeysError(f"Duplicate keys found:\n" + "\n".join(error_messages))
    else:
        print("No duplicates found.")
        return data # Return the parsed data if no duplicates.


# Function to check for duplicate keys in JSON object

def check_duplicates_from_string(json_string):

    """
    Parses a JSON string to detect and report any duplicate keys at the same level in the same object.
    This function ensures that no keys are silently overwritten in dictionaries.

    The function uses a helper to track the number of times each key appears during parsing,
    leveraging the `object_pairs_hook` parameter of `json.loads()` to intercept key-value pairs
    before they are processed into a dictionary. If duplicates are detected at any level, they
    are reported with their counts. Keys reused in separate objects within arrays (e.g. lists) 
    are not considered duplicates.

    Args:
        json_string (str): The JSON content as a string to parse and check for duplicates.

    Raises:
        json.JSONDecodeError: If the string is not valid JSON.
        DuplicateKeysError: If duplicate keys are found in the JSON structure.

    Returns:
        dict or list: The parsed data if no errors or duplicates are found.
    """

    # Initialize a dictionary to track duplicate keys and their counts
    duplicate_keys = {}
    
    # Create a 1-argument hook callable by "freezing" the duplicate_keys dict
    # as the second argument to the helper.
    hook = functools.partial(_detect_duplicates, duplicate_keys_state=duplicate_keys)

    # Parse the JSON string using the helper to track duplicates
    data = json.loads(json_string, object_pairs_hook=hook)
    
    return _process_results(data, duplicate_keys)
    
# Function for check for duplicate keys if input file is in JSON format

def check_duplicates_from_json(json_file_path):
    """
    Parses a JSON file to detect and report any duplicate keys at the same level in the same object.
    This function ensures that no keys are silently overwritten in dictionaries.

    The function uses a helper to track the number of times each key appears during parsing,
    leveraging the `object_pairs_hook` parameter of `json.load()` to intercept key-value pairs 
    before they are processed into a dictionary. If duplicates are detected at any level, they
    are reported with their counts and paths. Keys reused in separate objects within arrays 
    (e.g. lists) are not considered duplicates.
    
    Args:
        json_file_path (str): The path to the JSON file to parse and check for duplicates.

    Raises:
        FileNotFoundError: If the specified file does not exist.
        json.JSONDecodeError: If the file content is not valid JSON.
        DuplicateKeysError: If duplicate keys are found in the JSON structure.

    Returns:
        dict or list: The parsed data if no errors or duplicates are found.
    """

    # Initialize a dictionary to track duplicate keys and their counts
    duplicate_keys = {}
    
    # Create a 1-argument hook callable by "freezing" the duplicate_keys dict
    # as the second argument to the helper.
    hook = functools.partial(_detect_duplicates, duplicate_keys_state=duplicate_keys)

    # Open and parse the JSON file, using the helper to track duplicates
    with open(json_file_path, 'r') as file:
        data = json.load(file, object_pairs_hook=hook)
        
    return _process_results(data, duplicate_keys)


def create_json():
    """
    Loads the input file specified in `config.py`, filters to synthetic sequences,
    loads the plasmid backbone, and constructs the JSON object to be sent to a Predictor.

    Returns:
        evaluator_dict (dict): The constructed request payload.

    Raises:
        FileNotFoundError: If the input file does not exist.
        ValueError: If duplicate sequence IDs are found, or if the data is malformed.
    """

    # Validate evaluator input file exists
    if not os.path.exists(EVALUATOR_INPUT_PATH):
        print(f"ERROR: Evaluator input file '{EVALUATOR_INPUT_PATH}' not found.")
        raise FileNotFoundError(f"Evaluator input file not found: {EVALUATOR_INPUT_PATH}")

    try:
        input_data = pd.read_csv(EVALUATOR_INPUT_PATH, delimiter='\t', header=0)
        # Filter sequences to only synthetic sequences that would have no train test leakage (51k total)
        input_data_synthetic = input_data[input_data['origin'].isin([
            "Simulated_Annealing", "FastSeqProp" , "AdaLead"
            ])]
        print("Data loaded shape:")
        print(input_data_synthetic.shape)
        
        # Check for duplicates in the 'ID' column
        if input_data_synthetic['ID'].duplicated().any():
            duplicates = input_data_synthetic['ID'][input_data_synthetic['ID'].duplicated()].tolist()
            raise DuplicateKeysError(
                f"Duplicate values found in 'ID' column: {duplicates[:5]}"
                f"{'...' if len(duplicates) > 5 else ''}. "
                "These would be silently overwritten in the request payload."
            )
        else:
            print("No duplicates found in 'ID' column.")

        #Read in plasmid backbone
        backbone = pd.read_csv(PLASMID_BACKBONE_INPUT_PATH, header=0, sep= '\t', index_col=0)
        print(backbone)
        # Get backbone sequences and gene coordinates from the backbone file
        upstream_seq = backbone.loc["upstream", "sequence"]
        downstream_seq = backbone.loc["downstream_padded", "sequence"]
        promoter_coordinates = json.loads(backbone.loc["gene_coordinates", "sequence"])
        
        # These parameters are decided based on the sequence dataset.
        # Define prediction_tasks as a raw JSON string and validate BEFORE parsing,
        # so duplicate keys inside task definitions are actually caught (v0 ran the check
        # after json.loads, by which point any duplicates had already been silently dropped).
        prediction_tasks_str = """
        [
            {
                "name": "gosai_synthetic_sequences_k562",
                "type": "expression",
                "cell_type": "K562 (erythroid precursors)",
                "scale": "log",
                "species": "homo_sapiens"
            },
            {
                "name": "gosai_synthetic_sequences_hepg2",
                "type": "expression",
                "cell_type": "HepG2 (hepatocytes)",
                "scale": "log",
                "species": "homo_sapiens"
            },
            {
                "name": "gosai_synthetic_sequences_sknsh",
                "type": "expression",
                "cell_type": "SK-N-SH (neuroblastoma)",
                "scale": "log",
                "species": "homo_sapiens"
            }
        ]
        """
        prediction_tasks = check_duplicates_from_string(prediction_tasks_str)
    
        # Build the sequences dictionary from the DataFrame
        sequence_dict = dict(zip(input_data_synthetic.ID, input_data_synthetic.sequence))
        
        # gene span is read from backbone, not hardcoded
        # Map every sequence to the same gene coordinates (read from the backbone file)
        prediction_ranges = {name: promoter_coordinates for name in sequence_dict.keys()}
        
        evaluator_dict = {
            "readout": "point",
            "prediction_tasks": prediction_tasks,
            "upstream_seq": upstream_seq,
            "downstream_seq": downstream_seq,
            "sequences": sequence_dict,
            "prediction_ranges": prediction_ranges
        }
        
        # Final safety-net check on the fully assembled payload
        json_string = json.dumps(evaluator_dict, indent=4)
        check_duplicates_from_string(json_string)

        print("Input data loaded and validated successfully.")
        # Return the dict directly, not the result of re-parsing the JSON
        return evaluator_dict
    
    except (json.JSONDecodeError, DuplicateKeysError, KeyError) as e:
        # Raise a general ValueError that the main script's handler
        # will catch and report cleanly
        raise ValueError(f"Input data is invalid.\nDetails: {e}") from e