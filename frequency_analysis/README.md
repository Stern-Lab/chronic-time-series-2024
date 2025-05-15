# 🧬 SARS-CoV-2 Chronic Patients Mutation's Frequency Analysis

This Jupyter notebook processes mutation frequency data from chronically infected SARS-CoV-2 patients. 
It includes data loading, manual corrections, filtering, annotation, table export, plotting.

---

## 📁 Input Files & Requirements

| File Path                                           | Description & Required Columns                                                                                                                                       |
|-----------------------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| 📋 `chronic_metadata.csv`                           | Patient metadata. **Required columns:**<br>• `sample_id`<br>• `patient`<br>• `timepoint`<br>• other sample info (e.g., collection date)                               |
| 🗂️ `*/mutation_analysis/*/*/*merged.csv`            | Per-sample mutation outputs. **Required columns:**<br>• `mutation` (e.g. `A23403G`)<br>• `final_freq` (frequency)<br>• other pipeline-specific metrics              |
| 🛠️ `FINAL_SUMMARY_hard_mutation_drops_w_tech_seq_p5.csv` | Manual corrections summary. **Required columns:**<br>• `patient_id` / `patient`<br>• `timepoint`<br>• `mutation`<br>• `fixed_freq`                                    |
| 🚨 `Recurrent_mutations_unreplicated_data.tsv`       | Known problematic mutations. **Required columns:**<br>• `mutation`<br>• any flags or counts for recurrence                                                            |
| 🔎 `patients_timepoints_to_analyze.csv`              | Whitelist of valid timepoints per patient. **Required columns:**<br>• `patient`<br>• `timepoints` (e.g., `0,45,90`)                                                    |
| 📖 `wuhan_ref.fasta`                                 | Reference genome FASTA. Used for annotation.                                                                                                                          |

> 💡 **Tip:** Ensure all paths are updated in the **Paths & Global Constants** section before running.

---

## 🛠️ Notebook Sections

1. **Imports & Dependencies**  
2. **Paths & Global Constants**  
3. **Load & Merge Data**  
4. **Manual Corrections**  
5. **Filter Relevant Rows**  
6. **Load Reference & Annotation Helpers**  
7. **Annotate & Prepare Final Table**  
8. **Export Timepoint Tables**  
9. **Plotting Functions & Main Loop**

---

## ⚙️ Dependencies

Install required packages:

pandas
numpy
biopython
fpdf
plotly
tqdm

---

## 🚀 Usage

1. **Clone/download** this notebook.  
2. **Update** file paths in **Section 2**.  
3. **Run** cells sequentially.  
4. **Outputs**:  
   - Per-patient Excel tables in `Frequency_Tables` 📊  
   - Publication-ready PNG & SVG plots in `Frequency_Graphs` 📈  
