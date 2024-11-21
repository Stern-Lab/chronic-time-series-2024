import Filter_Usecase.filter_replicate_script.filter_functions as ff
import Filter_Usecase.filter_replicate_script.protein_data as prot
import time
from datetime import datetime
import pandas as pd
import numpy as np
import os

def get_freq_file(sample_id, rep_dirs, rep_path):
    if sample_id in rep_dirs:
        freq_file = rep_path + "/" + sample_id + "/freqs.tsv"
        return freq_file, True
    else:
        print(f"Directory {sample_id} wasn't found in {rep_path}. Skipping...")
        return "", False

def main(ui=False):
    log_txt = ""
    print("Data filtering and Usecase table creation script is starting...")
    date_time_str = datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
    print(date_time_str)
    log_txt += f"Data filtering and Usecase table creation script is starting...\n{date_time_str}\n"

    start = time.time()
    try:
        while True:
            if ui:
                user_input = "1"
            else:
                user_input = input("Enter 1 for cluster and 2 for local run: ")
                
            # Prepartion actions
            if (user_input == "1"):
                # REP_PATH = r"/sternadi/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/V3+V4"
                # REP_PATH = r"/sternadi/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/first_timepoint_as_reference"
                REP_PATH = r"/sternadi/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/AccuNGS/1st_timepoint_as_ref_w-manual_trim_of_trimmomatic_samples"
                PATIENTS = r"/sternadi/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_analysis/replicates_2023/all_patients_global_content_initials_V4.csv"
                break
            elif (user_input == "2"):
                REP_PATH = r"Z:/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_after_pipeline/replicates_2023/AccuNGS/1st_timepoint_as_ref_w-manual_trim_of_trimmomatic_samples"
                PATIENTS = r"Z:/nobackup/volume1/natalie/ichilov_chronic_omicron/libraries_analysis/replicates_2023/all_patients_global_content_initials_V4.csv"
                break
            else:
                user_input = input("Wrong input please try again\nEnter 1 for cluster and 2 for local run: ")
        
        log_txt += f"Replicates path: {REP_PATH}\n"
        log_txt += f"Patients info table path: {PATIENTS}\n"

        protein_dict = prot.create_protein_dict() # Protein dictionary
        all_patients_df = pd.read_csv(PATIENTS) # Load patients data
        sample_size = all_patients_df.shape[0]
        res_cols = ["Patient", "Timepoint", "sample_ID", "Date","Ct", "total_merged_mutations", "merged_mutations_with_f", "merged_mutations_NA", "merged_mutations_0"]
        results_df = pd.DataFrame(columns=res_cols) # Create results data frame
        
        # Get list of all directories
        replicate_dirs = os.listdir(REP_PATH)

        # Change default params
        while True:
            if ui:
                change_filter = "n"
            else:
                change_filter = input("Defaults filtering paramters are FREQ = 0.01, COVERAGE = 100, BASECOUNT = 50.\nDo you want to change filtering parameters (y/n)? ")
            
            if (change_filter == "n"):
                FREQ = 0.01
                COVERAGE = 100
                BASECOUNT = 50
                break

            elif (change_filter == "y"):
                FREQ = float(input("Minimal Frequnecy: "))
                COVERAGE = int(input("Minimal Coverage: "))
                BASECOUNT = int(input("Minimal Basecount: "))
                break

            else:
               print("Wrong input! Please enter (y/n)") 
        log_txt += f"Defaults filtering paramters are FREQ = {FREQ}, COVERAGE = {COVERAGE}, BASECOUNT = {BASECOUNT}.\n"
        
        # Keep or filter indels
        while True:
            if ui:
                indel_input = "n"
            else:
                indel_input = input("Do you want to filter (remove) indels mutations (y/n)? ")
            
            if (indel_input == "n"):
                filter_indels = False
                break
            elif (indel_input == "y"):
                filter_indels = True
                break
        
        log_txt += f"Filter (removing) indels: {filter_indels}\n"

        # Creates results directory
        run_dir = f"results_({FREQ}_{COVERAGE}_{BASECOUNT})_{date_time_str}"
        res_dir = r"./Filter_Usecase/results/" + run_dir
        if not os.path.exists(res_dir):
            os.makedirs(res_dir)

        log_txt += f"Results directory: {res_dir}\n"
        ind = 0
        for i, curr_row in all_patients_df.iterrows():

            print(f"Progress: {(i/sample_size*100):.2f}%")

            # Check if sample is not nan
            if (np.isnan(curr_row["sample"])):
                print("Sample id is empty. Skipping iteration...")
                continue

            # Get patient info
            curr_sample_id = str(int(curr_row["sample"]))
            if (pd.isna(curr_row["patient_ID"])):
                print("Patient id is empty. Skipping iteration...")
                continue
            curr_patient_id = curr_row["patient_ID"][:2]

            if (np.isnan(curr_row["time_since_first_sampling"])):
                print("Sample time is empty. Skipping iteration...")
                continue
            curr_timepoint = int(curr_row["time_since_first_sampling"])

            print(f"Filtering patient {curr_patient_id} timepoint {curr_timepoint}.")
            
            log_txt += f"Filtering patient {curr_patient_id} timepoint {curr_timepoint}.\n"
            
            # Find replicate1 file
            s1_rep1, found = get_freq_file(curr_sample_id, replicate_dirs, REP_PATH)
            if not (found):
                print("Replicate1 wasn't found. Skipping iteration...")
                log_txt += "Replicate1 wasn't found. Skipping iteration...\n"
                continue
            
            # Find replicate2 file
            s1_rep2, found = get_freq_file(curr_sample_id + "_L001", replicate_dirs, REP_PATH)
            if not (found):
                print("Replicate2 wasn't found. Skipping iteration...")
                log_txt += "Replicate2 wasn't found. Skipping iteration...\n"
                continue
            
            # Update index
            ind += 1
            
            # Update patients data to results_df
            results_df.loc[ind, "Timepoint"] = curr_timepoint
            results_df.loc[ind, "Patient"] = curr_patient_id
            results_df.loc[ind, "sample_ID"] = curr_sample_id
            results_df.loc[ind, "Date"] = curr_row["sampling date"]
            results_df.loc[ind, "Ct"] = curr_row["mean_ct"]
            
            # Create specific result dur for patient/timepoint
            specific_res_dir = f"{res_dir}/{curr_patient_id}/{curr_timepoint}"
            os.makedirs(specific_res_dir)
            
            # Filter timepoint and update data to results
            total_merged_mutations, merged_mutations_NA, merged_mutations_0, merged_mutations_with_f = ff.filter(s1_rep1, s1_rep2, FREQ, COVERAGE, BASECOUNT, protein_dict, specific_res_dir, curr_patient_id, curr_timepoint, filter_indels)
            results_df.loc[ind, "total_merged_mutations"] = total_merged_mutations
            results_df.loc[ind, "merged_mutations_with_f"] = merged_mutations_with_f
            results_df.loc[ind, "merged_mutations_NA"] = merged_mutations_NA
            results_df.loc[ind, "merged_mutations_0"] = merged_mutations_0

            if (total_merged_mutations != merged_mutations_with_f+merged_mutations_NA+merged_mutations_0):
                print("***Error***\nIn decision tree phase 2, there is a case that's not covered!!!")
        
        # Save data frame as a file
        results_df.to_csv(res_dir + "/Results.csv", index=False)
        
        # Save log data as a file
        tot_time = (time.time() - start)
        with open(res_dir + "/log.txt", 'w') as log_file:
            log_txt += f"Script elapsed time: {tot_time} sec"
            log_file.write(log_txt)
        
        print("***Filter Script finished successfully!***")
        print(f"Script elapsed time: {tot_time} sec")
    
    except Exception as e:
        # Save log data as a file
        tot_time = (time.time() - start)
        with open(res_dir + "/log.txt", 'w') as log_file:
            log_txt += f"Script elapsed time: {tot_time} sec"
            log_file.write(log_txt)
        print("An error has occured!\nTerminating script...")
        print(f"Filter script elapsed time: {(time.time() - start)} sec")
        print("Exception:")
        print(e)
        exit(1)

if __name__ == "__main__":
    main(True)