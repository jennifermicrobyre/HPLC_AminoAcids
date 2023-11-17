import pandas as pd
import pathlib
import yaml
import numpy as np
import matplotlib.pyplot as plt
import os

def inputs():
    experimentID = 'exp-741748223893'
    date_time_stamp = '2023-05-11 12-18-24'

    # select the kind of plots you want. You need to answer 'Yes' to at least one
    plot_final = 'Yes'     #Answer 'Yes' or 'No'
    plot_final_minus_initial = 'Yes'   #Answer 'Yes' or 'No'

    # parameter extracted from Gaussian Process analysis to plot against HPLC
    # e.g. 'Max Specific Growth Rate (1/Hour)' or 'ln(Max Amplitude)'
    plot_x_label = 'Max Specific Growth Rate (1/Hour)'

    # This is for sources of organic acids other than those detailed in the definitions yaml and the fillstock yamls
    # Units must be micromolar (uM) !!!!  Also - looking for concentration of cystine (not cysteine) !!!!!!!
    other_sources = [['alanine', 0], ['arginine', 0], ['asparagine', 0], ['aspartate', 0], ['citrulline', 0], ['cystine', 0],
                     ['glutamate', 0], ['glutamine', 0], ['glycine', 0], ['histidine', 0], ['isoleucine', 0], ['leucine', 0],
                     ['lysine', 0], ['methionine', 0], ['norvaline', 0], ['ornithine', 0], ['phenylalanine', 0], ['proline', 0],
                     ['sarcosine', 0], ['serine', 0], ['threonine', 0], ['tryptophan', 0], ['tyrosine', 0], ['valine', 0]
                     ]

    # This should not need to be changed, but take a look at it and make sure nothing is missing
    # This is to match chemicals that are in the definitions and media yamls to the appropriate molecules detected by HPLC
    # The number provides the number of moles of the ion of interest. e.g. for 1 mole of calciumLactatePentahydrate there are 2 moles of lactate
    match_chemical_w_ion = [['ornithineMonohydrochloride', 'ornithine', 1], ['L-norvaline', 'norvaline', 1],
                                     ['histidineMonohydrochlorideMonohydrate', 'histidine', 1]
                                     ]

    ####################################################################################
    # Stuff that should not need to be changed
    ####################################################################################
    last_characters = 'uM'  # this helps identify the columns that have the amount of the molecule in uM
    num_last_characters = 2  # number of characters in last_characters

    filename = 'hplc_clean_data.csv'

    # 'root' location of HPLC data
    hplc_directory_path = pathlib.Path('/Volumes/data/prod/Experiments/RobotData/HPLC/Results/amino_acid/')
    hplc_data_path = hplc_directory_path / experimentID / date_time_stamp / filename

    # location of the yamls
    yaml_layouts_path = pathlib.Path('/Volumes/data/prod/Experiments/Layouts/')
    yaml_definitions_path = pathlib.Path('/Volumes/data/prod/Experiments/Definitions/')
    yaml_made_path = pathlib.Path('/Volumes/data/prod/Experiments/Made/')

    # biolector analysis location
    bl_location = pathlib.Path('/Volumes/data/prod/Experiments/Analytics/')

    # Dictionary with molecular weights
    conversion_dict = dict(
        [('alanine', 89.09), ('arginine', 174.2), ('asparagine', 132.12), ('aspartate', 133.1), ('citrulline', 175.2), ('cystine', 121.16),
         ('glutamate', 147.13), ('glutamine', 146.14), ('glycine', 75.07), ('histidine', 155.15), ('isoleucine', 131.17), ('leucine', 131.17),
         ('lysine', 146.19), ('methionine', 149.21), ('norvaline', 117.15), ('ornithine', 132.16), ('phenylalanine', 165.19), ('proline', 115.13),
         ('sarcosine', 89.093), ('serine', 105.09), ('threonine', 119.12), ('tryptophan', 204.23), ('tyrosine', 181.19), ('valine', 117.151)
         ])
    #######################################

    return plot_final, plot_final_minus_initial, experimentID, yaml_definitions_path, yaml_layouts_path, yaml_made_path, hplc_data_path, last_characters, num_last_characters, conversion_dict, other_sources, match_chemical_w_ion, bl_location, plot_x_label

def get_single_analysis(experiment, analysis_loc):
    files = os.listdir(analysis_loc / experiment)

    for i in range(len(files)):
        if files[i][:8] == 'analysis':
            analysis_filename = files[i]
            break

    analysis = pd.read_csv(analysis_loc / experiment / analysis_filename, keep_default_na=False)

    return analysis

def open_hplc_dataframe(hplc_data_path):
    hplc_df = pd.read_csv(hplc_data_path, keep_default_na=False)

    return hplc_df

def remove_calibration_samples(hplc_df):
    indices_to_drop = []
    for i in range(len(hplc_df)):
        plate = hplc_df['plate_name'].iloc[i]
        if plate[:5] != 'plate':
            indices_to_drop.append(i)
    calibrations_removed_df = hplc_df.drop(indices_to_drop)

    return calibrations_removed_df

def open_yaml_expcond(file_loc, file, type):
    # open the yaml files, if we are opening the layout yaml file '_layout.yaml' needs to be appended to the experiment ID
    if type == 'layout':
        filename = file + '_layout.yaml'
        file_to_open = file_loc / filename
    elif type == 'made':
        filename = file + '_made.yaml'
        file_to_open = file_loc / filename
    else:
        filename = file + '.yaml'
        file_to_open = file_loc / filename

    with open(file_to_open) as file:
        try:
            data = yaml.safe_load(file)
        except yaml.YAMLError as exception:
            print(exception)

    return data

def stack_columns(hplc_df, last_characters, num_last_characters):
    stacked_hplc_df = pd.DataFrame({'plate_name': [], 'well_id': [], 'Inoculant': [], 'Condition': [], 'Amino Acid': [], 'HPLC Measured Concentration (uM)': []})
    headers = list(hplc_df.columns)

    for i in range(len(headers)):
        if headers[i][-num_last_characters:] == last_characters:
            aa = []
            for j in range(len(hplc_df)):
                aa.append(headers[i].split('_')[0])
            df = pd.DataFrame({'plate_name': hplc_df['plate_name'], 'well_id': hplc_df['well_id'], 'Inoculant': hplc_df['Inoculant'], 'Condition': hplc_df['Condition'], 'Amino Acid': aa, 'HPLC Measured Concentration (uM)': hplc_df[headers[i]]})
            stacked_hplc_df = pd.concat([stacked_hplc_df, df])

    for i in range(len(stacked_hplc_df)):
        if stacked_hplc_df['HPLC Measured Concentration (uM)'].iloc[i] == 'none':
            stacked_hplc_df['HPLC Measured Concentration (uM)'].iloc[i] = 0
    stacked_hplc_df['HPLC Measured Concentration (uM)'] = stacked_hplc_df['HPLC Measured Concentration (uM)'].astype(float)

    return stacked_hplc_df

def get_concentration(value, aa, conversion_dict, multiplication):
    value = str(value)
    if value[-2:] == ' M':  # need to leave space before M otherwise it catches mM, uM and nM
        c = float(value[:-2]) * 1e6
    elif value[-2:] == 'mM':
        c = float(value[:-3]) * 1e3
    elif value[-2:] == 'uM':
        c = float(value[:-3])
    elif value[-2:] == 'nM':
        c = float(value[:-3]) / 1e3
    elif value[-3:] == 'g/L':
        c = float(value[:-4])
        c /= conversion_dict[aa]
        c *= 1e6
    elif value[-3:] == 'g/l':
        c = float(value[:-4])
        c /= conversion_dict[aa]
        c *= 1e6
    else:
        print('units on organic acid in made yaml had either (1) no units or (2) unusual units, value was set to zero: ', str)
        c = 0

    c *= multiplication

    return c

def extract_from_yaml(experimentID, yaml_layouts_path):
    layout = open_yaml_expcond(yaml_layouts_path, experimentID, 'layout')

    plates = []
    well_id = []
    inoculant = []
    cond = []
    keys1 = list(layout['Plates'].keys())
    for i in range(len(keys1)):
        keys2 = list(layout['Plates'][keys1[i]]['Wells'].keys())
        for j in range(len(keys2)):
            plates.append(keys1[i].split('.')[1])
            well_id.append(keys2[j])
            inoculant.append(layout['Plates'][keys1[i]]['Wells'][keys2[j]]['inoculant'])
            cond.append(layout['Plates'][keys1[i]]['Wells'][keys2[j]]['variant'])

    df = pd.DataFrame({'plate_name': plates, 'well_id': well_id, 'Inoculant': inoculant, 'Condition': cond})

    return df

def add_condition(hplc_df, experimentID, yaml_layouts_path):
    layout_df = extract_from_yaml(experimentID, yaml_layouts_path)
    hplc_df = pd.merge(hplc_df, layout_df, on=['well_id', 'plate_name'])

    return hplc_df

def expcond(data):
    dict = {}

    variants = data['Variants']
    keys = list(variants.keys())
    for i in range(len(keys)):
        keys2 = list(variants[keys[i]].keys())
        for j in range(len(keys2)):
            keys3 = list(variants[keys[i]][keys2[j]].keys())
            values = list(variants[keys[i]][keys2[j]].values())
            for t in range(len(values)):
                values[t] = str(values[t])
            for k in range(len(keys3)):
                new_value = ''
                for m in range(len(values[k])):
                    new_value += values[k][m]
                    if values[k][m] == ' ':
                        temp = float(new_value)
                        temp = np.round(temp, 3)
                        new_value = str(temp)
                        new_value += values[k][m]
                dict[keys3[k]] = keys[i] + " " + new_value

    return dict

def add_xaxis_info(yaml_definitions_path, experimentID, hplc_df):
    data = open_yaml_expcond(yaml_definitions_path, experimentID, 'definition')
    experiment_type = data['Design']

    if experiment_type == 'Dilution':
        expcond_dict = expcond(data)
        conditions = []
        for i in range(len(hplc_df)):
            conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + expcond_dict[hplc_df['Condition'].iloc[i]])

    elif experiment_type == 'Desolation':
        expcond_dict = expcond(data)
        conditions = []
        for i in range(len(hplc_df)):
            if hplc_df['Condition'].iloc[i] != 'cond000':
                conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + expcond_dict[hplc_df['Condition'].iloc[i]])
            else:
                conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + 'Control')

    else:
        conditions = hplc_df['Condition']

    hplc_df.insert(4, 'Experimental Conditions', conditions)

    return hplc_df

def make_dataframe_using_made(experimentID, yaml_made_path):
    made = open_yaml_expcond(yaml_made_path, experimentID, 'made')

    conditions = list(made['Actual'].keys())
    df = pd.DataFrame({'Condition': [], 'Chemical': [], 'Amount': [], 'Ion': []})
    for i in range(len(conditions)):
        chemicals = list(made['Actual'][conditions[i]].keys())
        ions = list(made['Actual'][conditions[i]].keys())
        amount = []
        for j in range(len(chemicals)):
            amount.append(made['Actual'][conditions[i]][chemicals[j]])
        df1 = pd.DataFrame({'Chemical': chemicals, 'Amount': amount, 'Ion': ions})
        df1.insert(0, 'Condition', conditions[i])
        df = pd.concat([df, df1])

    return df

def convert_to_ion(df, match_chemical_w_ion, amino_acids):
    df.insert(4, 'Multiplication', np.ones(len(df)))

    for i in range(len(df)):
        for j in range(len(match_chemical_w_ion)):
            if df['Chemical'].iloc[i] == match_chemical_w_ion[j][0]:
                df.iat[i, 3] = match_chemical_w_ion[j][1]
                df.iat[i, 4] = match_chemical_w_ion[j][2]
                break

    for i in range(len(amino_acids)):
        if amino_acids[i] == 'cystine':
            amino_acids[i] = 'cysteine'
    df_filtered = df[df['Ion'].isin(amino_acids)]

    return df_filtered

def convert_to_uM(df, conversion_dict):
    df.insert(5, 'Concentration (uM)', np.zeros(len(df)))
    for i in range(len(df)):
        df.iat[i, 5] = get_concentration(df['Amount'].iloc[i], df['Ion'].iloc[i], conversion_dict, df['Multiplication'].iloc[i])

    return df

def add_made_to_hplc_df(hplc_df, made_df):
    hplc_df.insert(7, 'Concentration (uM) from Made', np.zeros(len(hplc_df)))
    for i in range(len(hplc_df)):
        for j in range(len(made_df)):
            if hplc_df['Amino Acid'].iloc[i] != 'cystine':
                if hplc_df['Condition'].iloc[i] == made_df['Condition'].iloc[j] and hplc_df['Amino Acid'].iloc[i] == \
                        made_df['Ion'].iloc[j]:
                    hplc_df.iat[i, 7] = made_df['Concentration (uM)'].iloc[j]
                    break
            else:
                if hplc_df['Condition'].iloc[i] == made_df['Condition'].iloc[j] and made_df['Ion'].iloc[j] == 'cysteine':
                    hplc_df.iat[i, 7] = made_df['Concentration (uM)'].iloc[j] / 2
                    break

    return hplc_df

def get_initial_from_made(hplc_df, experimentID, yaml_made_path, conversion_dict, match_chemical_w_ion):
    made_df = make_dataframe_using_made(experimentID, yaml_made_path)
    made_df = convert_to_ion(made_df, match_chemical_w_ion, hplc_df['Amino Acid'].unique())
    made_df = convert_to_uM(made_df, conversion_dict)
    hplc_df = add_made_to_hplc_df(hplc_df, made_df)

    return hplc_df

def get_cond_from_other_sources(hplc_df, other_sources):
    conc_other_sources = np.zeros(len(hplc_df))
    for i in range(len(hplc_df)):
        for j in range(len(other_sources)):
            if hplc_df['Amino Acid'].iloc[i] == other_sources[j][0]:
                conc_other_sources[i] = other_sources[j][1]
                break

    hplc_df.insert(8, 'Concentration (uM) from Other Sources', conc_other_sources)

    return hplc_df

def get_initial_concentration(hplc_df):
    hplc_df['Initial Concentration (uM)'] = hplc_df['Concentration (uM) from Made'] + hplc_df['Concentration (uM) from Other Sources']

    return hplc_df

def update_HPLC_with_BL(hplc_df, bl_df):

    hplc_df.rename(columns={'plate_name': 'PlateNumber', 'well_id': 'Well Number'}, inplace=True)
    hplc_bl_df = pd.merge(hplc_df, bl_df, on=['PlateNumber', 'Well Number'])

    return hplc_bl_df

def colors_symbols():
    colors = ['blue', 'orange', 'skyblue', 'grey', 'black', 'red', 'green', 'purple', 'yellow', 'peru',
              'blue', 'orange', 'skyblue', 'grey', 'black', 'red', 'green', 'purple', 'yellow', 'peru',
              'blue', 'orange', 'skyblue', 'grey', 'black', 'red', 'green', 'purple', 'yellow', 'peru',
              'blue', 'orange', 'skyblue', 'grey', 'black', 'red', 'green', 'purple', 'yellow', 'peru',
              'blue', 'orange', 'skyblue', 'grey', 'black', 'red', 'green', 'purple', 'yellow', 'peru']
    symbols = ['x', 'v', '^', 'o', '>', '<', 'P', 'h', 'D', 'p', 'd', '*', '8', 'h', '1', '2', '3', '4',
               'x', 'v', '^', 'o', '>', '<', 'P', 'h', 'D', 'p', 'd', '*', '8', 'h', '1', '2', '3', '4',
               'x', 'v', '^', 'o', '>', '<', 'P', 'h', 'D', 'p', 'd', '*', '8', 'h', '1', '2', '3', '4']

    return colors, symbols

def individual_plots_w_BL(hplc_bl_df, plot_x_label):
    colors, symbols = colors_symbols()

    # get list of amino acids
    aa_tested = pd.unique(hplc_bl_df['Amino Acid'])
    aa_tested.sort()

    for i in range(len(aa_tested)):
        # set up the figure
        fig, axs = plt.subplots()
        fig.set_figwidth(8.8)
        fig.set_figheight(5.7)

        subset = hplc_bl_df.loc[(hplc_bl_df['Amino Acid'] == aa_tested[i])]

        subset_mean = subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
        subset_std = subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)

        if subset_mean.empty == False:
            for j in range(len(subset_mean)):
                leg = subset_mean['Experimental Conditions'][j]
                if pd.isna(subset_std['HPLC Measured Concentration (uM)'].iloc[j]) == False:
                    axs.scatter(subset_mean[plot_x_label][j], subset_mean['HPLC Measured Concentration (uM)'][j], linestyle='None', marker=symbols[j], color=colors[j], label=leg)
                    axs.errorbar(subset_mean[plot_x_label][j], subset_mean['HPLC Measured Concentration (uM)'][j], xerr=subset_std[plot_x_label][j], yerr=subset_std['HPLC Measured Concentration (uM)'][j], linestyle='None', color=colors[j])
                else:
                    axs.scatter(subset_mean[plot_x_label][j], subset_mean['HPLC Measured Concentration (uM)'][j], linestyle='None', marker=symbols[j], color=colors[j], label=leg)

        axs.set_ylabel('HPLC: Concentration (uM)')
        axs.set_xlabel(plot_x_label)
        axs.legend(bbox_to_anchor=(1, 1), loc='upper left')

        plt.subplots_adjust(top=0.95, bottom=0.25, left=0.09, right=0.6)
        plt.title(aa_tested[i])
        plt.show()

def individual_plots_final_minus_initial_w_BL(hplc_bl_df, plot_x_label):
    colors, symbols = colors_symbols()

    # get list of amino acids
    aa_tested = pd.unique(hplc_bl_df['Amino Acid'])
    aa_tested.sort()

    for i in range(len(aa_tested)):
        # set up the figure
        fig, axs = plt.subplots()
        fig.set_figwidth(8.8)
        fig.set_figheight(5.7)

        subset = hplc_bl_df.loc[(hplc_bl_df['Amino Acid'] == aa_tested[i])]

        subset_mean = subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
        subset_std = subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)

        if subset_mean.empty == False:
            for j in range(len(subset_mean)):
                leg = subset_mean['Experimental Conditions'][j]
                if pd.isna(subset_std['Final Minus Initial (uM)'].iloc[j]) == False:
                    axs.scatter(subset_mean[plot_x_label][j], subset_mean['Final Minus Initial (uM)'][j], linestyle='None', marker=symbols[j], color=colors[j], label=leg)
                    axs.errorbar(subset_mean[plot_x_label][j], subset_mean['Final Minus Initial (uM)'][j], xerr=subset_std[plot_x_label][j], yerr=subset_std['HPLC Measured Concentration (uM)'][j], linestyle='None', color=colors[j])
                else:
                    axs.scatter(subset_mean[plot_x_label][j], subset_mean['Final Minus Initial (uM)'][j], linestyle='None', marker=symbols[j], color=colors[j], label=leg)

        axs.set_ylabel('Concentration (uM): Final - Initial')
        axs.set_xlabel(plot_x_label)
        axs.legend(bbox_to_anchor=(1, 1), loc='upper left')
        axs.axhline(0, color='black', linestyle='dotted')

        plt.subplots_adjust(top=0.95, bottom=0.25, left=0.09, right=0.6)
        plt.title(aa_tested[i])
        plt.show()

def get_final_minus_initial(hplc_df):
    hplc_df['Final Minus Initial (uM)'] = hplc_df['HPLC Measured Concentration (uM)'] - hplc_df['Initial Concentration (uM)']

    return hplc_df

def save(hplc_df, experimentID):
    hplc_df.to_csv(experimentID + '.csv')

def main():
    # read in inputs
    plot_final, plot_final_minus_initial, experimentID, yaml_definitions_path, yaml_layouts_path, yaml_made_path, hplc_data_path, last_characters, num_last_characters, conversion_dict, other_sources, match_chemical_w_ion, bl_location, plot_x_label = inputs()
    # get HPLC data
    hplc_df = open_hplc_dataframe(hplc_data_path)
    # remove calibration samples
    hplc_df = remove_calibration_samples(hplc_df)
    # add condition to data frame
    hplc_df = add_condition(hplc_df, experimentID, yaml_layouts_path)
    # get amino acids & stack columns
    hplc_df = stack_columns(hplc_df, last_characters, num_last_characters)
    # add experimental conditions - details, e.g. serine 500 uM not cond001
    hplc_df = add_xaxis_info(yaml_definitions_path, experimentID, hplc_df)
    # get concentrations from made
    hplc_df = get_initial_from_made(hplc_df, experimentID, yaml_made_path, conversion_dict, match_chemical_w_ion)
    # get concentrations from other sources
    hplc_df = get_cond_from_other_sources(hplc_df, other_sources)
    # get initial concentrations of amino acids
    hplc_df = get_initial_concentration(hplc_df)
    # read BioLector analysis
    bl_df = get_single_analysis(experimentID, bl_location)
    # update HPLC data with BioLector data
    hplc_bl_df = update_HPLC_with_BL(hplc_df, bl_df)
    #plot stuff
    if plot_final == 'Yes':
        # individual plots
        individual_plots_w_BL(hplc_bl_df, plot_x_label)
    # final minus initial
    if plot_final_minus_initial == 'Yes':
        # get final minus initial concentration
        hplc_bl_df = get_final_minus_initial(hplc_bl_df)
        # individual plots final minus initial
        individual_plots_final_minus_initial_w_BL(hplc_bl_df, plot_x_label)
    # save data frame
    save(hplc_bl_df, experimentID)

main()

