### IMPORTS ###
import pandas as pd
import numpy as np
import re
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
import plotly.io as pio
from fpdf import FPDF
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy import stats

### Functions ###

def remove_patients_timepoints(df, patients_to_remove):
    # Iterate over the dictionary to remove specified timepoints for each patient
    for patient, timepoints in patients_to_remove.items():
        for timepoint in timepoints:
            # Remove rows where both the patient and the timepoint match
            df = df[~((df['patient'] == patient) & (df['days_since_first'] == timepoint))]
    return df
  
def remove_days_greater_than(df, max_days):
    # Filter out rows where 'days_since_first' is greater than the specified max_days
    df = df[df['days_since_first'] <= max_days]
    return df
  
def get_position(mutation):
    # Define a regular expression pattern to extract the numeric position
    position_pattern = re.compile(r'[ACGT](\d+)[ACGT+-]*')
    # Define a helper function to extract the position using the regex
    match = position_pattern.search(mutation)
    return int(match.group(1)) if match else None

def remove_non_coding_muts(df, orf_boundaries=orf_boundaries):
    # Helper function to extract the position from mutation
    df['position'] = df['mutation'].apply(get_position)
    # Filter mutations within any of the ORF boundaries
    coding_df = df[df['position'].apply(lambda pos: any(start <= pos <= end for start, end in orf_boundaries))]
    # Drop the temporary 'position' column
    coding_df = coding_df.drop(columns=['position'])
    return coding_df

# Define a function to check if the mutation length after + or - is divisible by 3
def is_triplet_indel(mutation_type):
    # Search for patterns starting with + or - followed by letters
    match = re.search(r'([+-])([A-Za-z]+)', mutation_type)
    if match:
        # Check if the length of the sequence is divisible by 3
        sequence = match.group(2)
        return len(sequence) % 3 == 0
    return True  # Keep rows that don't match the pattern as they are not indels

def remove_non_triplet_indels(df):
    # Apply the helper function and filter the dataframe
    filtered_df = df[df['mutation_type'].apply(is_triplet_indel)]
    return filtered_df

def Calc_divergence(df, mutation_kinds=["Syn"], mutation_percentage=0.22, genome_length=29264):
    # Calculate the actual genome size affected by mutation percentage
    actual_genome_size = genome_length * mutation_percentage
    # Capture all unique time points per patient before filtering by mutation type
    unique_timepoints_per_patient = df[['days_since_first', 'patient']].drop_duplicates()
    # Filter by the specified mutation kinds
    df_for_divergence_calc = df[df['mut_kind'].isin(mutation_kinds)]
    # Sum frequencies by 'days_since_first' and 'patient'
    summed_frequencies = df_for_divergence_calc.groupby(['days_since_first', 'patient'])["final_after_manual"].sum()
    # Convert summed frequencies to a dataframe
    summed_frequencies_df = summed_frequencies.reset_index()
    # Merge with unique time points to fill missing time points with a frequency of 0
    merged_df = unique_timepoints_per_patient.merge(summed_frequencies_df, on=['days_since_first', 'patient'], how='left').fillna(0)
    # Normalize the summed frequencies for divergence calculation
    merged_df['final_after_manual'] = merged_df['final_after_manual'] / actual_genome_size
    # Return the final dataframe with columns 'days_since_first', 'patient', and 'final_after_manual'
    divergence_df = merged_df[['days_since_first', 'patient', 'final_after_manual']]
    return divergence_df

def regress_divergence(df):
    # Prepare the data
    X = df['days_since_first'].values.reshape(-1, 1)
    y = df['final_after_manual'].values
    # Initialize the linear regression model, with intercept forced to 0
    model = LinearRegression(fit_intercept=False)
    model.fit(X, y)
    slope = model.coef_[0]
    intercept = model.intercept_
    y_pred = model.predict(X)
    residuals = y - y_pred
    s_err = np.sqrt(np.sum(residuals**2) / (len(X) - 1)) / np.sqrt(np.sum(X.flatten()**2))
    t_stat = slope / s_err
    p_value = 2 * (1 - stats.t.cdf(np.abs(t_stat), df=len(X) - 1))
    return model, slope, intercept, p_value, residuals.std()


def plot_regression_with_divergence(df, slope, intercept, p_value, model, mut_type, residual_std, confidence_level=0.95, is_uniqe=""):
    # Prepare the data and sort it for a smooth plot
    df = df.sort_values('days_since_first')
    X = df['days_since_first'].values
    y = df['final_after_manual'].values
    trend_line_y = slope * X  # Calculate the trend line
    r2 = r2_score(y, trend_line_y)  # Calculate R²
    equation = f'y = {slope:.4e}x + {intercept:.4e} (R² = {r2:.4e}, p-value = {p_value:.4e})'  # Equation with scientific notation
    # Calculate the confidence interval for sorted data
    se_line = np.sqrt(np.sum((y - model.predict(X.reshape(-1, 1)))**2) / (len(X) - 2)) * \
              np.sqrt(1 / len(X) + (X - np.mean(X))**2 / np.sum((X - np.mean(X))**2))
    conf_interval = stats.t.ppf(1 - (1 - confidence_level) / 2, df=len(X) - 2) * se_line
    upper_bound = trend_line_y + conf_interval
    lower_bound = trend_line_y - conf_interval
    # Create the figure
    fig = go.Figure()
    # Plot data points for each patient
    for patient in df['patient'].unique():
        patient_df = df[df['patient'] == patient]
        fig.add_trace(go.Scatter(
            x=patient_df['days_since_first'],
            y=patient_df['final_after_manual'],
            mode='markers',
            name=patient,
            marker=dict(size=15),
            legendgroup='Patients',
            showlegend=True
        ))
    # Add an asterisk annotation for patient N1 on day 241
    bal_sample = df[(df['patient'] == 'N1') & (df['days_since_first'] == 241)]
    if not bal_sample.empty:
        fig.add_annotation(
            x=bal_sample['days_since_first'].values[0],
            y=bal_sample['final_after_manual'].values[0],
            text="*",
            showarrow=True,
            arrowhead=2,
            ax=15,
            ay=-15,
            font=dict(size=16, color="red")  # Customize the asterisk's appearance
        )
    # Plot the confidence interval band
    fig.add_trace(go.Scatter(
        x=np.concatenate((X, X[::-1])),
        y=np.concatenate((upper_bound, lower_bound[::-1])),
        fill='toself',
        fillcolor='rgba(150,150,150,0.3)',  # Gray color for the confidence interval shading
        line=dict(color='rgba(150,150,150,0)'),  # No line for edges
        showlegend=False,
        name=f'Confidence Interval ({confidence_level*100:.0f}%)'
    ))
    # Plot the regression line
    fig.add_trace(go.Scatter(
        x=X,
        y=trend_line_y,
        mode='lines',
        name=equation,
        line=dict(color='black', width=3),
        legendgroup='Trendline',
        showlegend=True
    ))
    # Add ±1 standard deviation bands for 'Not-Syn' and 'Indels'
    if mut_type in ['Not-Syn', 'Indels']:
        y_upper_sd = trend_line_y + residual_std
        y_lower_sd = trend_line_y - residual_std
        # Plot the shaded region between y_upper_sd and y_lower_sd
        fig.add_trace(go.Scatter(
            x=np.concatenate((X, X[::-1])),
            y=np.concatenate((y_upper_sd, y_lower_sd[::-1])),
            fill='toself',
            fillcolor='rgba(0,100,250,0.2)',  # Light blue color for shading
            line=dict(color='rgba(0,100,250,0)'),  # No line for edges
            showlegend=False,
            name=f'{mut_type} ±1 SD'
        ))
    # Add an additional legend entry for the BAL sample
    fig.add_trace(go.Scatter(
        x=[None], y=[None],  # Invisible entry for legend
        mode='markers',
        marker=dict(size=10, color='rgba(0,0,0,0)'),  # Transparent marker
        showlegend=False,
        name="* - BAL sample"
    ))
    # Layout settings
    fig.update_layout(
        title=f'Divergence per Day with Regression Trend Line for {mut_type}',
        title_x=0.5,
        xaxis_title='Days',
        yaxis_title='Synonymous divergence (r<sub>Syn</sub>)',
        showlegend=True,
        plot_bgcolor='white',
        paper_bgcolor='white',
        xaxis=dict(showgrid=True, gridcolor='lightgrey', linecolor='lightgrey', automargin=True),
        yaxis=dict(showgrid=True, gridcolor='lightgrey', linecolor='lightgrey', automargin=True, exponentformat="e"),
        width=1200,
        height=800,
        font=dict(size=23),
        legend=dict(
            yanchor="top",
            y=-0.2,
            xanchor="center",
            x=0.5,
            orientation="h",
            traceorder="grouped",
            font=dict(size=20),
            itemsizing='constant'
        ),
        margin=dict(l=150, r=100, t=100, b=150)
    )
    # Save the figure as a high-quality PNG file
    fig.write_image(f"/sternadi/home/volume3/adi_bz/Low_Freqs_24/Graphs/divergence_pngs/divergence_plot_{mut_type}{is_uniqe}.png", scale=3) #change to fit export dir
    # Save the plot as an HTML file and open in a new tab
    html_file_path = f"/sternadi/home/volume3/adi_bz/Low_Freqs_24/Graphs/divergence_pngs/divergence_plot_{mut_type}{is_uniqe}.html" #change to fit export dir
    pio.write_html(fig, file=html_file_path, auto_open=False)
    # Show the plot
    fig.show()

def regress_and_table_non_syn(rel_df,kinds):
    patient_values= {}
    for patient in np.sort(rel_df.patient.unique()).tolist():
        final_df_per_p = rel_df[rel_df['patient']==patient]
        non_syn_df=Calc_divergence(final_df_per_p,kinds,mutation_percentage_non_synonymous)
        model_ns, divergence_rate_ns, intercept_ns, p_value_ns, residual_std_ns = regress_divergence(non_syn_df)
        patient_values[patient]=divergence_rate_ns
        if(p_value_ns>0.05):
            print(f"P-value for {patient} is not significant at {p_value_ns} - divergence at that point is {divergence_rate_ns}")
        print(f"P-value for {patient} is {p_value_ns} - divergence at that point is {divergence_rate_ns}")
    # Prepare data for the table
    patients = list(patient_values.keys())
    divergences = [f"{value:.4e}" for value in patient_values.values()]
    r_syn_chosen = r_dict['del_N1_241'][0]
    syn_divergences = [f"{r_syn_chosen:.4e}"] * len(patients)
    ratios = [value / r_syn_chosen for value in patient_values.values()]  # Keep as floats
    # Create the DataFrme
    syn_non_syn_df = pd.DataFrame({
        'Patient': patients,
        'r_ns': divergences,
        'r_s': syn_divergences,
        'ratios': [f"{ratio:.4f}" for ratio in ratios]
    })
    # Normalize the ratio values for color scaling
    min_val = 0
    max_val = 3.2
    norm_ratios = [(float(val) - min_val) / (max_val - min_val) for val in syn_non_syn_df['ratios']]
    scale = px.colors.sequential.Sunset
    colors = [scale[int(r * (len(scale) - 1))] for r in norm_ratios]
    # Define cell-specific font colors
    num_rows = len(patients)
    font_colors = [['black'] * num_rows, ['black'] * num_rows, ['black'] * num_rows]  # Default all black
    # Set specific cell color, e.g., cell at row i=1, col j=2 to white
    i, j = 8, 2  # Define the specific cell position
    font_colors[j][i] = 'white'  # Change color of cell (i=1, j=2) to white
    # Create the Plotly table with custom cell-specific font color
    fig = go.Figure(data=[go.Table(
        header=dict(values=['Patient', 'r<sub>Non_Syn</sub>', 'r<sub>Non_Syn</sub> / r<sub>Syn</sub>'],
                    fill_color='#4b494d',
                    font=dict(size=22, color='white'),
                    align='center',
                    height=20),
        cells=dict(values=[syn_non_syn_df['Patient'].tolist(),
                           syn_non_syn_df['r_ns'].tolist(),
                           syn_non_syn_df['ratios'].tolist()],
                   fill_color=[['#f4f5f6'] * num_rows,  
                               ['#f4f5f6'] * num_rows,
                               colors],                
                   font=dict(size=17,color=font_colors),      
                   align='center',
                   height=30),
        columnwidth=[40, 50, 50]
    )])
    # Set the overall size of the table
    fig.update_layout(width=700, height=1000)
    # Save the figure as a high-quality PNG file
    fig.write_image(f"/sternadi/home/volume3/adi_bz/Low_Freqs_24/Graphs/divergence_pngs/divergence_table_{kinds}.png", scale=3)
    # Display the table
    fig.show()
    return patient_values



### Parameters ###
# Define the non-coding regions of the SARS-CoV-2 genome
non_coding_regions = [(1,265+1),(29675, 29903+1)]

sars_cov_2_proteins = {
    'ORF1ab': [266, 21555],
    'S (Spike)': [21563, 25384],
    'ORF3a': [25393, 26220],
    'E (Envelope)': [26245, 26472],
    'M (Membrane)': [26523, 27191],
    'ORF6': [27202, 27387],
    'ORF7a': [27394, 27759],
    'ORF7b': [27756, 27887],
    'ORF8': [27894, 28259],
    'N (Nucleocapsid)': [28274, 29533],
    'ORF10': [29558, 29674]
}

# Flattened list of all ORF boundaries from the dictionary
orf_boundaries = [(start, end) for start, end in sars_cov_2_proteins.values()]

# Define the coding genes and their lengths
coding_genes_and_lengths = [
    ("ORF1ab", 21290),
    ("S", 3822),
    ("ORF3a", 828),
    ("E", 228),
    ("M", 669),
    ("ORF6", 186),
    ("ORF7a", 366),
    ("ORF7b", 132),
    ("ORF8", 366),
    ("N", 1260),
    ("ORF10", 117)
]
# OPTIONAL: Sum the lengths (should be 29264) (I just calculated the numbers by hand).
total_length = sum(length for gene, length in coding_genes_and_lengths)

# Define the number of possible non-synonymous and synonymous mutations
genome_length_without_non_coding = 29264  # SARS-CoV-2 genome length without non-coding regions
mutation_percentage_non_synonymous = 0.72
mutation_percentage_synonymous = 0.22

# Patients and timepoints to remove from analysis
patients_to_remove = {
    'P5': [27, 44,75],
    'P6': [0,22,30,35,37],
    'P4':[13],
    'P3':[9]
}

### Preprocessing ###

# Files Imports
final_df = pd.read_pickle("/sternadi/home/volume3/adi_bz/Low_Freqs_24/Datasets/freqs_for_graphing_with_zeros.pkl")
final_df = final_df[final_df['final_after_manual']>=0]

# Remove irelevant patients
final_df = remove_patients_timepoints(final_df,patients_to_remove)
    
# Remove values that are not-whithin the coding regions
final_df = remove_non_coding_muts(final_df)

# Remove Indels that are not in 3 Nucleotide
final_df =remove_non_triplet_indels(final_df)

# Initiallize the r_values dictionary
r_dict = {}

### Run Linear Models and Graph###

#Synonymous - All
syn_df=Calc_divergence(final_df,["Syn"],mutation_percentage_synonymous)
model_s,divergence_rate_s,intercept_s, p_val_s, res_s=regress_divergence(syn_df)
plot_regression_with_divergence(syn_df, divergence_rate_s,intercept_s, p_val_s, model_s,"Syn", res_s)
r_dict['no_filter']=[divergence_rate_s,intercept_s,p_val_s]

#Synonymous - w/o N1 T241 (BAL sample)
max_days = 230 
final_df_241 = remove_days_greater_than(final_df, max_days)
syn_df=Calc_divergence(final_df_241,["Syn"],mutation_percentage_synonymous)
model_s,divergence_rate_s,intercept_s, p_val_s, res_s=regress_divergence(syn_df)
plot_regression_with_divergence(syn_df, divergence_rate_s,intercept_s, p_val_s, model_s,"Syn", res_s,0.95, "_No_N1_241_")
r_dict['del_N1_241']=[divergence_rate_s,intercept_s,p_val_s]

#Synonymous - day limit - 100 days
max_days = 100 
final_df_100 = remove_days_greater_than(final_df, max_days)
syn_df=Calc_divergence(final_df_100,["Syn"],mutation_percentage_synonymous)
model_s,divergence_rate_s,intercept_s, p_val_s, res_s=regress_divergence(syn_df)
plot_regression_with_divergence(syn_df, divergence_rate_s,intercept_s, p_val_s, model_s,"Syn", res_s,0.95, "_max_days_100_")
r_dict['max_100']=[divergence_rate_s,intercept_s,p_val_s]

#Synonymous - day limit - 50 days
max_days = 50 
final_df_50 = remove_days_greater_than(final_df, max_days)
syn_df=Calc_divergence(final_df_50,["Syn"],mutation_percentage_synonymous)
model_s,divergence_rate_s,intercept_s, p_val_s, res_s=regress_divergence(syn_df)
plot_regression_with_divergence(syn_df, divergence_rate_s,intercept_s, p_val_s, model_s,"Syn", res_s,0.95, "_max_days_50_")
r_dict['max_50']=[divergence_rate_s,intercept_s,p_val_s]

#Synonymous - w/o N1 T241 (BAL sample) and w/o P5
max_days =230 
final_df_241 = remove_days_greater_than(final_df, max_days)
final_df_no_p5 = final_df_241[final_df_241['patient']!='P5']
syn_df=Calc_divergence(final_df_no_p5,["Syn"],mutation_percentage_synonymous)
model_s,divergence_rate_s,intercept_s, p_val_s, res_s=regress_divergence(syn_df)
plot_regression_with_divergence(syn_df, divergence_rate_s,intercept_s, p_val_s, model_s,"Syn", res_s,0.95, "_No_P5_")
r_dict['no_p5']=[divergence_rate_s,intercept_s,p_val_s]

#Non-Synonymous - w/o N1 T241 (BAL sample)
ns_df_rs = regress_and_table_non_syn(final_df_241,['Not-Syn']) # will also print P-values 

#Non-Synonymous and Indels - w/o N1 T241 (BAL sample)
ns_and_indels_df_rs = regress_and_table_non_syn(final_df_241,['Not-Syn','Indel']) # will also print P-values 


#Graph the different Syn divergence rates for later use
r_syns_df = pd.DataFrame.from_dict(r_dict, orient='index', columns=['Divergence Rate', 'Intercept', 'P-Value'])
divergence_rate_formatted = [f"{x:.5e}" for x in r_syns_df['Divergence Rate']]
fig = go.Figure(data=[go.Table(
    header=dict(
        values=['', 'Divergence Rate', 'Intercept', 'P-Value'],
        fill_color='#4b494d',
        font=dict(size=22, color='white'),
        align='center',
        height=40
    ),
    cells=dict(
        values=[
            r_syns_df.index.tolist(),
            divergence_rate_formatted,
            r_syns_df['Intercept'].tolist(),
            r_syns_df['P-Value'].tolist()
        ],
        fill_color='#f4f5f6',
        font=dict(size=17, color='black'),
        align='center',
        height=30
    ),
    columnwidth=[30, 70, 50, 70]
)])
# Remove margins
fig.update_layout(
    width=1000, height=500,
    margin=dict(l=0, r=0, t=0, b=0))  # Set all margins to zero
# Save the figure as a high-quality PNG file
fig.write_image(f"/sternadi/home/volume3/adi_bz/Low_Freqs_24/Graphs/syn_divergence_rates.png", scale=5)
# Display the table
fig.show()
r_syns_df.to_pickle("/sternadi/home/volume3/adi_bz/Low_Freqs_24/Datasets/syn_divergence_rates.pkl")
