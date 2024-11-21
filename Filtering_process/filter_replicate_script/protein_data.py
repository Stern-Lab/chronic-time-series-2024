import pandas as pd
import re

PATH = r"./Filter_Usecase/filter_replicate_script/mutation_type_full_mutation.csv"

def create_protein_dict():
    print("Creating protein dictionary...")
    protein_df = pd.read_csv(PATH)
    protein_dict = dict()
    pattern = r"\['(.*?)'\]"

    for ind, row in protein_df.iterrows():
        match = re.search(pattern, row['mutation_type'])
        mut_type = ""
        if match: 
            mut_type = match.group(1)
            
        protein_dict[row['full_mutation']] = (mut_type, row['protein'])


    print("Dictionary created...")
    return protein_dict
    