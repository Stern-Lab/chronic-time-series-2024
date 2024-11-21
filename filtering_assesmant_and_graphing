## IMPORTS
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import numpy as np

##Initial loading of the data

# Define base directory
base_directory = "" # Enter path to directory that includes nested sub-directories for patients and timepoint 

# Initialize an empty list to collect dataframes and track missing files
dataframes = []
missing_data = []

# Traverse through each patient folder in the base directory
for patient_folder in os.listdir(base_directory):
    patient_path = os.path.join(base_directory, patient_folder)
    if os.path.isdir(patient_path):  # Ensure it’s a directory
        # Traverse through each timepoint folder within the patient folder
        for timepoint_folder in os.listdir(patient_path):
            timepoint_path = os.path.join(patient_path, timepoint_folder)
            if os.path.isdir(timepoint_path):  # Ensure it’s a directory
                # Search for files containing "phase1.csv" and identify those with "freq1" or "freq2"
                phase_files = [f for f in os.listdir(timepoint_path) if "phase1.csv" in f]
                
                # Check for required files 'freq1' and 'freq2' as part of larger strings
                required_files = {"rep1": None, "rep2": None}
                for file_name in phase_files:
                    if "freq1" in file_name:
                        required_files["rep1"] = file_name
                    elif "freq2" in file_name:
                        required_files["rep2"] = file_name
                
                # Process each required file if found, else log missing data
                for replicate, file_name in required_files.items():
                    if file_name:
                        file_path = os.path.join(timepoint_path, file_name)
                        
                        # Read the CSV file
                        try:
                            df = pd.read_csv(file_path)
                            
                            # Add the 'Patient', 'timepoint', and 'replicate' columns
                            df['Patient'] = patient_folder
                            df['timepoint'] = timepoint_folder
                            df['replicate'] = replicate
                            
                            # Append the dataframe to the list
                            dataframes.append(df)
                        except Exception as e:
                            print(f"    Error reading {patient_folder}, {timepoint_folder}, {file_name}: {e}")
                    else:
                        # Log missing files if they aren't found
                        missing_data.append((patient_folder, timepoint_folder, replicate))

# Concatenate all dataframes into a single large dataframe if any dataframes are found
if dataframes:
    combined_df = pd.concat(dataframes, ignore_index=True)
    print("Final DataFrame created successfully.")
else:
    combined_df = pd.DataFrame()  # Empty DataFrame if no data found
    print("No data found to create the DataFrame.")

# Display missing files information if any
if missing_data:
    print("Missing files:")
    for patient, day, replicate in missing_data:
        print(f"  Could not find data for Patient '{patient}' on day '{day}' for replicate '{replicate}'")

##Removing unwanted samples
#Filteration is done using either by removing specific patients or by loading a df of relevant patients:

def remove_patients_timepoints(df, patients_to_remove):
    """
    Removes specified timepoints for given patients from the dataframe.
    If a patient has an empty list as their value, all data for that patient is removed.
    
    Parameters:
        df (pd.DataFrame): The original dataframe with columns 'Patient' and 'timepoint'.
        patients_to_remove (dict): A dictionary where keys are patient names and values are lists of timepoints to remove.
    
    Returns:
        pd.DataFrame: The modified dataframe with specified patient-timepoint pairs removed.
    """
    # Filter out specified patient-timepoint pairs
    for patient, timepoints in patients_to_remove.items():
        if not timepoints:  # Check if timepoints list is empty
            # Remove all data for the patient
            df = df[df['Patient'] != patient]
        else:
            # Remove only specified timepoints for the patient
            df = df[~((df['Patient'] == patient) & (df['timepoint'].isin(timepoints)))]
    
    return df
		
def filter_patients_timepoints(df, patients_timepoint_dict):
    # Ensure 'Patient' and 'timepoint' columns are strings
    df['Patient'] = df['Patient'].astype(str)
    df['timepoint'] = df['timepoint'].astype(str)
    
    # Create a set of valid (Patient, timepoint) pairs from the dictionary
    valid_pairs = {(str(patient), str(timepoint)) for patient, timepoints in patients_timepoint_dict.items() for timepoint in timepoints}
    
    # Print valid pairs to debug
    print("Valid pairs:", valid_pairs)

    # Filter the dataframe by checking if each row's (Patient, timepoint) pair is in valid_pairs
    filtered_df = df[df.apply(lambda row: (row['Patient'], row['timepoint']) in valid_pairs, axis=1)]
    
    # Debugging: Show the number of rows before and after filtering
    print(f"Original row count: {len(df)}, Filtered row count: {len(filtered_df)}")
    
    return filtered_df

# Example usage
patients_to_remove = {
    #'P5': ['27', '44','75'],
    'P6': []
    #'P4':['13'],
    #'P3':['9']
}

#Use one of the methods (or a combination of both) :
#filtered_df = filter_patients_timepoints(combined_df, patients_timepoint_dict) #in a comment since I used only the second one
filtered_df = remove_patients_timepoints(combined_df, patients_to_remove)


###Processing of Samples:
def find_incomplete_samples(df):
    """
    Identifies samples (Patient + timepoint) in the dataframe that do not have exactly two replicates.
    Parameters:
        df (pd.DataFrame): The dataframe with columns 'Patient', 'timepoint', and 'replicate'.
    Returns:
        List[Tuple[str, str]]: A list of tuples, each containing the patient and timepoint of samples without two replicates.
    """
    # Group by Patient and timepoint and count unique replicates for each group
    replicate_counts = df.groupby(['Patient', 'timepoint'])['replicate'].nunique().reset_index()
    # Filter for groups that do not have exactly 2 replicates
    incomplete_samples = replicate_counts[replicate_counts['replicate'] != 2]
    # Convert to a list of (Patient, timepoint) tuples
    missing_replicates = list(zip(incomplete_samples['Patient'], incomplete_samples['timepoint']))
    
    return missing_replicates

# Example usage
# Assuming combined_df is the dataframe created earlier
missing_replicates = find_incomplete_samples(filtered_df)
# Filter `filtered_df` to keep only (Patient, timepoint) combinations with exactly two replicates
filtered_df = filtered_df[filtered_df.groupby(['Patient', 'timepoint'])['replicate'].transform('nunique') == 2]
# Assuming combined_df is the dataframe created earlier
filtered_df = remove_patients_timepoints(filtered_df, patients_to_remove)
# Display the result
print("Removed Samples without two replicates:", missing_replicates)


### Plot rep1 vs rep 2:

def plot_patient_timepoints(df, output=None):
    """
    Plots patient-timepoint subplots and optionally saves as a high-resolution PNG.
    Parameters:
        df (pd.DataFrame): Dataframe containing 'Patient', 'timepoint', 'mutation', 'replicate', and 'new_freq' columns.
        output (str): Optional. File path to save the plot as a PNG. If None, the plot will only display.
    """
    # Filter out frequencies of -1
    df = df[df['new_freq'] != -1]
    # Convert 'timepoint' to numeric if it's not already
    df['timepoint'] = pd.to_numeric(df['timepoint'], errors='coerce')
    # Sort the dataframe by 'Patient' and 'timepoint' to ensure desired plotting order
    df = df.sort_values(by=['Patient', 'timepoint'], ascending=[True, True])
    # Get unique patients and timepoints
    patients = df['Patient'].unique()
    timepoints = df['timepoint'].unique()
    # Determine number of rows and columns for the subplots grid
    num_subplots = len(df.groupby(['Patient', 'timepoint']))
    num_cols = 6  # Set the number of columns to 6
    num_rows = int(np.ceil(num_subplots / num_cols))
    # Create the subplots figure
    fig, axes = plt.subplots(num_rows, num_cols, figsize=(20, num_rows * 4))
    axes = axes.flatten()  # Flatten axes array for easy indexing
    # Plot each patient-timepoint combination
    subplot_idx = 0
    for (patient, timepoint), group in df.groupby(['Patient', 'timepoint']):
        ax = axes[subplot_idx]
        # Find mutations that are present in at least one replicate and get their frequencies in each
        mut_freqs = group.pivot_table(index='mutation', columns='replicate', values='new_freq', fill_value=0)
        # Scatter plot for frequencies in rep1 vs rep2, using get with default value of 0 to avoid KeyErrors
        rep1_freqs = mut_freqs.get('rep1', pd.Series(0, index=mut_freqs.index))
        rep2_freqs = mut_freqs.get('rep2', pd.Series(0, index=mut_freqs.index))
        ax.scatter(rep1_freqs, rep2_freqs)
        # Add y=x line as a regression line
        ax.plot([0, 1], [0, 1], color='black', linestyle='--', linewidth=1)
        # Set x and y limits
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        # Set title for each subplot
        ax.set_title(f"Patient {patient} - Timepoint {timepoint}")
        # Count outliers (frequency > 0.2 in one replicate and 0 in the other)
        outliers = ((rep1_freqs > 0.2) & (rep2_freqs == 0)) | \
                   ((rep2_freqs > 0.2) & (rep1_freqs == 0))
        outlier_count = outliers.sum()
        # Display outlier count
        ax.text(0.05, 0.9, f"Outliers: {outlier_count}", transform=ax.transAxes, fontsize=10,
                bbox=dict(facecolor='white', alpha=0.8))
        # Move to the next subplot
        subplot_idx += 1
    # Hide any unused subplots if the grid has extra spaces
    for i in range(subplot_idx, len(axes)):
        axes[i].axis('off')
    # Add a single set of x and y labels for the entire figure, with increased font size and padding
    fig.text(0.5, 0.02, 'Frequency in Replicate 1', ha='center', fontsize=16, fontweight='bold')
    fig.text(0.02, 0.5, 'Frequency in Replicate 2', va='center', rotation='vertical', fontsize=16, fontweight='bold')
    # Adjust layout to add space for the labels outside the subplots
    plt.subplots_adjust(left=0.07, bottom=0.07, right=0.95, top=0.95, wspace=0.3, hspace=0.4)
    # Save as high-resolution PNG if output path is provided
    if output:
        plt.savefig(output, dpi=300, format='png', bbox_inches='tight')
    # Show the plot
    plt.show()

# Example usage assuming `filtered_df` is the filtered dataframe
plot_patient_timepoints(filtered_df, output="/sternadi/home/volume3/adi_bz/Low_Freqs_24/Graphs/comparing_freqs_graphs/patients_timepoints_freqs_plot.png")


