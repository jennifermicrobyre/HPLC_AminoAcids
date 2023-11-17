import pandas as pd
import pathlib
import yaml
import numpy as np
import matplotlib.pyplot as plt
import math
import os

def inputs():
    experimentID = 'exp-251765846464'

    # HPLC
    hplc_filename_columnnames = 'REPORT00.csv'
    hplc_filename_data = 'REPORT01_recoded2.csv'
    hplc_directory_location = pathlib.Path('/Users/jennifer/Documents/Analysis/HPLC_AminoAcids/exp-251765846464/')

    ##########################################
    # Stuff that should not need to be changed
    # location of the yamls, !!! these need to end with a /
    yaml_definitions_path = pathlib.Path('/Volumes/data/prod/Experiments/Definitions/')
    yaml_fillstock_path = pathlib.Path('/Volumes/data/prod/Chemistry/Combinations/Media/')
    reads_path = pathlib.Path('/Volumes/data/prod/Experiments/RobotData/Biolector/Results/')

    return experimentID, hplc_filename_columnnames, hplc_filename_data, hplc_directory_location, yaml_definitions_path, yaml_fillstock_path, reads_path

def make_hplc_dataframe(hplc_filename_columnnames, hplc_filename_data, hplc_directory_location):
    header_filepath =  hplc_directory_location / hplc_filename_columnnames
    data_filepath =  hplc_directory_location / hplc_filename_data

    # merge file header and data
    header_df = pd.read_csv(header_filepath, encoding='UTF-16')
    headers = header_df.iloc[:, 2].to_list()  # brittle encoding to a specific column
    headers = [header.replace('|', '_').lower() for header in headers]  # clean and standardize
    #hplc_df = pd.read_csv(data_filepath, encoding='UTF-16', header=None)
    hplc_df = pd.read_csv(data_filepath, header=None)
    hplc_df.columns = headers  # set headers from other file

    hplc_df = hplc_df.rename(columns={'sample': 'Plate & Well Number'})

    return headers, hplc_df

def remove_calibration_samples(hplc_df):
    indices_to_drop = []
    for i in range(len(hplc_df)):
        if hplc_df['Plate & Well Number'].iloc[i][:5] != 'plate':
            indices_to_drop.append(i)
    calibrations_removed_df = hplc_df.drop(indices_to_drop)

    return calibrations_removed_df

def split_plate_and_well(hplc_df):
    plate_array = []
    well_number_array = []
    for i in range(len(hplc_df)):
        temp = hplc_df['Plate & Well Number'].iloc[i]
        well_number_array.append(temp.split("_")[1])
        plate_array.append(temp.split("_")[0])
    hplc_df.insert(0, 'Well Number', well_number_array)
    hplc_df.insert(0, 'PlateID', plate_array)

    return hplc_df

def open_yaml_expcond(file_loc, file, type):
    # open the yaml files, if we are opening the layout yaml file '_layout.yaml' needs to be appended to the experiment ID
    if type == 'layout':
        filename = file + '_layout.yaml'
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

def stack_columns(headers, hplc_df):
    stacked_hplc_df = pd.DataFrame({'PlateID': [], 'Well Number': [], 'Inoculant': [], 'Condition': [], 'Amino Acid': [], 'HPLC Measured Concentration (uM)': []})

    for i in range(len(headers)):
        if headers[i][-7:] == '_amount':
            aa = []
            for j in range(len(hplc_df)):
                aa.append(headers[i][:-7])
            df = pd.DataFrame({'PlateID': hplc_df['PlateID'], 'Well Number': hplc_df['Well Number'], 'Inoculant': hplc_df['Inoculant'], 'Condition': hplc_df['Condition'], 'Amino Acid': aa, 'HPLC Measured Concentration (uM)': hplc_df[headers[i]]})
            stacked_hplc_df = pd.concat([stacked_hplc_df, df])

    return stacked_hplc_df

def convert_to_uM(str):
    if str[-2:] == 'mM':
        c = float(str[:-3]) * 1000
    elif str[-2:] == 'uM':
        c = float(str[:-3])
    elif str[-2:] == 'nM':
        c = float(str[:-3]) / 1000
    else:
        print('Concentration in Media yaml file does not have attached unit - check that initial concentrations are correct')
        c = -1

    return c

def get_initial_from_fillstock(hplc_df, experimentID, yaml_definitions_path, yaml_fillstock_path):
    definitions_yaml = open_yaml_expcond(yaml_definitions_path, experimentID, 'definitions')
    fillstock_yaml = open_yaml_expcond(yaml_fillstock_path, definitions_yaml['Fillstock'], 'fillstock')

    dilution = definitions_yaml['Fillstock Dilution']

    initial_conc = np.zeros(len(hplc_df))

    chemicals_in_fillstock = list(fillstock_yaml['components'].keys())
    for i in range(len(hplc_df)):
        if hplc_df['Amino Acid'].iloc[i] != 'cystine':
            for j in range(len(chemicals_in_fillstock)):
                if hplc_df['Amino Acid'].iloc[i] == chemicals_in_fillstock[j]:
                    c = convert_to_uM(fillstock_yaml['components'][chemicals_in_fillstock[j]]['solvation'])
                    initial_conc[i] = c * dilution
        else:
            try:
                c = convert_to_uM(fillstock_yaml['components']['cysteine']['solvation'])
                initial_conc[i] = c * dilution / 2
            except:
                c = 0

    hplc_df.insert(6, 'Concentration (uM) from Fillstock', initial_conc)

    return hplc_df

def get_data(experimentID, reads_path):

    directories = os.listdir(path=reads_path / experimentID); directories.sort()

    if directories[0] == '.DS_Store':
        directories = directories[1:]

    data = []
    for i in range(len(directories)):
        d = pd.read_csv(reads_path / experimentID / directories[i] / 'reads.csv', sep="\t")
        d['PlateID'] = directories[i]
        d.sort_values(by=['time (hours)'], inplace=True)
        d.reset_index(drop=True, inplace=True)
        data.append(d)

    return data

def add_condition(hplc_df, experimentID, reads_path):

    d = get_data(experimentID, reads_path)

    inoculant = []
    condition = []
    for i in range(len(hplc_df)):
        keep_going = 'Yes'
        for j in range(len(d)):
            if keep_going == 'No':
                break
            df = d[j]
            for k in range(len(df)):
                if hplc_df['Well Number'].iloc[i] == df['wellnum'].iloc[k] and hplc_df['PlateID'].iloc[i] == df['PlateID'].iloc[k]:
                    inoculant.append(df['inoculant'].iloc[k])
                    condition.append(df['condition'].iloc[k])
                    keep_going = 'No'
                    break
    hplc_df.insert(2, 'Inoculant', inoculant)
    hplc_df.insert(3, 'Condition', condition)

    return hplc_df

def expcond(file_loc, experiment):
    data = open_yaml_expcond(file_loc, experiment, 'definition')
    dict = {}

    variants = data['Variants']
    keys = list(variants.keys())
    for i in range(len(keys)):
        keys2 = list(variants[keys[i]].keys())
        for j in range(len(keys2)):
            keys3 = list(variants[keys[i]][keys2[j]].keys())
            values = list(variants[keys[i]][keys2[j]].values())
            for k in range(len(keys3)):
                new_value = ''
                for m in range(len(values[k])):
                    if values[k][m] != '/':
                        new_value += values[k][m]
                    else:
                        new_value += ' per '

                    if values[k][m] == ' ':
                        temp = float(new_value)
                        temp = np.round(temp, 2)
                        new_value = str(temp)
                        new_value += values[k][m]

                dict[keys3[k]] = keys[i] + " " + new_value

    return dict

def create_condition_array(condition, length):
    array = []
    for i in range(length):
        array.append(condition)

    return array

def get_variant_concentrations(aa, amino_acids, concentration_variant, c):
    for i in range(len(amino_acids)):
        if amino_acids[i] == aa:
            concentration_variant[i] = c
            break

    return concentration_variant

def get_conc_from_definitions(hplc_df, experimentID, yaml_definitions_path):
    definitions_yaml = open_yaml_expcond(yaml_definitions_path, experimentID, 'definitions')

    amino_acids = list(definitions_yaml['Variants'].keys())
    concentration_control = np.zeros(len(amino_acids))
    condition_array = []

    df = pd.DataFrame({'Condition': [], 'Amino Acid': [], 'Concentration (uM) From Definitions': []})

    # get concentrations for cond000
    for i in range(len(amino_acids)):
        condition_array.append(list(definitions_yaml['Variants'][amino_acids[i]]['Control'].keys())[0])
        c = definitions_yaml['Variants'][amino_acids[i]]['Control'][condition_array[i]]
        concentration_control[i] = convert_to_uM(c)
        if amino_acids[i] == 'cysteine':  # new code
            concentration_control[i] /= 2  # new code
    df1 = pd.DataFrame({'Condition': condition_array, 'Amino Acid': amino_acids, 'Concentration (uM) From Definitions': concentration_control})
    df = pd.concat([df, df1])

    # get concentrations for the rest of the conditions
    for i in range(len(amino_acids)):
        conditions = list(definitions_yaml['Variants'][amino_acids[i]]['Variants'].keys())
        for j in range(len(conditions)):
            c = convert_to_uM(definitions_yaml['Variants'][amino_acids[i]]['Variants'][conditions[j]])
            if amino_acids[i] == 'cysteine': # new code
                c /= 2                       # new code
            condition_array = create_condition_array(conditions[j], len(amino_acids))
            concentration_variant = get_variant_concentrations(amino_acids[i], amino_acids, np.copy(concentration_control), c)
            df1 = pd.DataFrame({'Condition': condition_array, 'Amino Acid': amino_acids, 'Concentration (uM) From Definitions': concentration_variant})
            df = pd.concat([df, df1])

    # Add the concentrations from the definitions yaml to the data frame
    conc_from_def = np.zeros(len(hplc_df))
    for i in range(len(hplc_df)):
        for j in range(len(df)):
            if hplc_df['Condition'].iloc[i] == df['Condition'].iloc[j] and hplc_df['Amino Acid'].iloc[i] == df['Amino Acid'].iloc[j]:
                conc_from_def[i] = df['Concentration (uM) From Definitions'].iloc[j]
            elif hplc_df['Condition'].iloc[i] == df['Condition'].iloc[j] and hplc_df['Amino Acid'].iloc[i] == 'cystine' and df['Amino Acid'].iloc[j] == 'cysteine':  #new code
                conc_from_def[i] = df['Concentration (uM) From Definitions'].iloc[j]  # new code

    hplc_df.insert(7, 'Concentration (uM) from Definitions', conc_from_def)

    #hplc_df.to_csv('/Users/jennifer/Desktop/hplc_testing.csv')
    return hplc_df

def add_concentrations(hplc_df):
    hplc_df['Concentration of Added AAs (uM)'] = hplc_df['Concentration (uM) from Fillstock'] + hplc_df['Concentration (uM) from Definitions']

    return hplc_df

def add_yaxis_info(yaml_definitions_path, experimentID, hplc_df):
    expcond_dict = expcond(yaml_definitions_path, experimentID)
    conditions = []
    for i in range(len(hplc_df)):
        if hplc_df['Condition'].iloc[i] != 'cond000':
            conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + expcond_dict[hplc_df['Condition'].iloc[i]])
        else:
            conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + 'Control')

    hplc_df.insert(4, 'Experimental Conditions', conditions)

    return hplc_df

def grid_plot(hplc_df):
    # get list of amino acids
    aa_tested = pd.unique(hplc_df['Amino Acid'])
    aa_tested.sort()

    # find number of rows for plot
    aa_number = len(aa_tested)
    num_y = math.floor(aa_number/5)
    remainder = aa_number % 5
    if remainder != 0 :
        num_y += 1

    # set up the figure
    fig, axs = plt.subplots(num_y, 5, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0.6)
    fig.set_figwidth(15)
    fig.set_figheight(num_y * 2)

    x = 0; y= 0
    colors = np.random.rand(len(pd.unique(hplc_df['Condition'])), 3)
    for i in range(aa_number):
        subset = hplc_df.loc[(hplc_df['Amino Acid'] == aa_tested[i])]
        cond = pd.unique(subset['Condition'])
        for j in range(len(cond)):
            subset_of_subset = subset.loc[(subset['Condition'] == cond[j])]
            leg = subset_of_subset['Experimental Conditions'].iloc[0]
            axs[y, x].scatter(subset_of_subset['Experimental Conditions'], subset_of_subset['HPLC Measured Concentration (uM)'], marker='x', color=colors[j, :], label=leg)
            axs[y, x].scatter(subset_of_subset['Experimental Conditions'], subset_of_subset['Concentration of Added AAs (uM)'], marker='+', color=colors[j, :], label=leg)
        axs[y, x].set_title(aa_tested[i])
        axs[y, x].tick_params(axis='x', labelrotation = 90)

        x += 1
        if x > 4:
            x = 0
            y += 1

    if x > 0:
        for i in range(5 - x):
            axs[y, x].tick_params(axis='x', labelrotation=90)
            x += 1

    fig.supylabel('Concentration (uM)')
    plt.subplots_adjust(top=0.9, bottom=0.25, left=0.09, right=0.75)
    axs[0, 4].legend(bbox_to_anchor=(1, 1), loc='upper left')
    fig.suptitle('+ Initial Concentration (uM) \n x Concentration Measured by HPLC (uM) \n ' )
    plt.show()
    #plt.savefig('/Users/jennifer/Desktop/grid_plot')

    return colors

def grid_plot_final_minus_initial(hplc_df, colors):
    # get list of amino acids
    aa_tested = pd.unique(hplc_df['Amino Acid'])
    aa_tested.sort()

    # find number of rows for plot
    aa_number = len(aa_tested)
    num_y = math.floor(aa_number/5)
    remainder = aa_number % 5
    if remainder != 0 :
        num_y += 1

    # set up the figure
    fig, axs = plt.subplots(num_y, 5, sharex=True, sharey=True)
    fig.subplots_adjust(wspace=0, hspace=0.6)
    fig.set_figwidth(15)
    fig.set_figheight(num_y * 2)

    x = 0; y= 0
    for i in range(aa_number):
        subset = hplc_df.loc[(hplc_df['Amino Acid'] == aa_tested[i])]
        cond = pd.unique(subset['Condition'])
        for j in range(len(cond)):
            subset_of_subset = subset.loc[(subset['Condition'] == cond[j])]
            leg = subset_of_subset['Experimental Conditions'].iloc[0]
            axs[y, x].scatter(subset_of_subset['Experimental Conditions'], subset_of_subset['Final Minus Initial (uM)'], marker='o', color=colors[j, :], label=leg)
        axs[y, x].set_title(aa_tested[i])
        axs[y, x].tick_params(axis='x', labelrotation = 90)
        axs[y, x].hlines(0, -1, len(subset), linestyles='dotted')

        x += 1
        if x > 4:
            x = 0
            y += 1

    if x > 0:
        for i in range(5 - x):
            axs[y, x].tick_params(axis='x', labelrotation=90)
            x += 1

    fig.supylabel('Final Minus Initial Concentration (uM)')
    plt.subplots_adjust(top=0.9, bottom=0.25, left=0.09, right=0.75)
    axs[0, 4].legend(bbox_to_anchor=(1, 1), loc='upper left')
    plt.show()
    #plt.savefig('/Users/jennifer/Desktop/finalminusinitial_grid_plot')

    return colors

def individual_plots(hplc_df, colors):
    # get list of amino acids
    aa_tested = pd.unique(hplc_df['Amino Acid'])
    aa_tested.sort()

    for i in range(len(aa_tested)):
        # set up the figure
        fig, axs = plt.subplots(1, 1)
        fig.set_figwidth(8)
        fig.set_figheight(8)

        subset = hplc_df.loc[(hplc_df['Amino Acid'] == aa_tested[i])]
        cond = pd.unique(subset['Condition'])
        for j in range(len(cond)):
            subset_of_subset = subset.loc[(subset['Condition'] == cond[j])]
            #leg = subset_of_subset['Experimental Conditions'].iloc[0]
            axs.scatter(subset_of_subset['Experimental Conditions'], subset_of_subset['HPLC Measured Concentration (uM)'], marker='x', color=colors[j, :], s=64)
            axs.scatter(subset_of_subset['Experimental Conditions'], subset_of_subset['Concentration of Added AAs (uM)'], marker='+', color=colors[j, :], s=64)
        axs.tick_params(axis='x', labelrotation = 90)
        axs.set_ylabel('Concentration (uM)')

        plt.subplots_adjust(top=0.85, bottom=0.35, left=0.1, right=0.9)
        #axs.legend(bbox_to_anchor=(1, 1), loc='upper left')
        plt.title('+ Initial Concentration (uM) \n x Concentration Measured by HPLC (uM) \n \n' + aa_tested[i])
        plt.show()
        #plt.savefig('/Users/jennifer/Desktop/' + str(i))

def individual_plots_final_minus_initial(hplc_df, colors):
    # get list of amino acids
    aa_tested = pd.unique(hplc_df['Amino Acid'])
    aa_tested.sort()

    for i in range(len(aa_tested)):
        # set up the figure
        fig, axs = plt.subplots(1, 1)
        fig.set_figwidth(8)
        fig.set_figheight(8)

        subset = hplc_df.loc[(hplc_df['Amino Acid'] == aa_tested[i])]
        cond = pd.unique(subset['Condition'])
        for j in range(len(cond)):
            subset_of_subset = subset.loc[(subset['Condition'] == cond[j])]
            #leg = subset_of_subset['Experimental Conditions'].iloc[0]
            axs.scatter(subset_of_subset['Experimental Conditions'], subset_of_subset['Final Minus Initial (uM)'], marker='o', color=colors[j, :], s=64)
        axs.tick_params(axis='x', labelrotation = 90)
        axs.set_ylabel('Final Minus Initial Concentration (uM)')
        axs.hlines(0, -1, len(subset), linestyles='dotted')

        plt.subplots_adjust(top=0.85, bottom=0.35, left=0.15, right=0.9)
        #axs.legend(bbox_to_anchor=(1, 1), loc='upper left')
        plt.title(aa_tested[i])
        plt.show()
        #plt.savefig('/Users/jennifer/Desktop/finalminuinitial_' + str(i))

def final_minus_initial(hplc_df):
    hplc_df['Final Minus Initial (uM)'] = hplc_df['HPLC Measured Concentration (uM)'] - hplc_df['Concentration of Added AAs (uM)']

    return hplc_df

def main():
    # read in inputs
    experimentID, hplc_filename_columnnames, hplc_filename_data, hplc_directory_location, yaml_definitions_path, yaml_fillstock_path, reads_path = inputs()
    # get HPLC data
    headers, hplc_df = make_hplc_dataframe(hplc_filename_columnnames, hplc_filename_data, hplc_directory_location)
    # remove calibration samples
    hplc_df = remove_calibration_samples(hplc_df)
    # split plate & well
    hplc_df = split_plate_and_well(hplc_df)
    # add condition to data frame
    hplc_df = add_condition(hplc_df, experimentID, reads_path)
    # get amino acids & stack columns
    hplc_df = stack_columns(headers, hplc_df)
    # get concentrations from fillstock
    hplc_df = get_initial_from_fillstock(hplc_df, experimentID, yaml_definitions_path, yaml_fillstock_path)
    # get concentrations from definitions yaml
    hplc_df = get_conc_from_definitions(hplc_df, experimentID, yaml_definitions_path)
    # add concentrations from fillstock to concentrations from yaml
    hplc_df = add_concentrations(hplc_df)
    # add experimental conditions - details, e.g. serine 500 uM not cond001
    hplc_df = add_yaxis_info(yaml_definitions_path, experimentID, hplc_df)
    # grid plot
    colors = grid_plot(hplc_df)
    # individual plots
    individual_plots(hplc_df, colors)
    # final minus initial
    hplc_df = final_minus_initial(hplc_df)
    # grid plot final minus initial
    grid_plot_final_minus_initial(hplc_df, colors)
    # individual plots final minus initial
    individual_plots_final_minus_initial(hplc_df, colors)


main()

