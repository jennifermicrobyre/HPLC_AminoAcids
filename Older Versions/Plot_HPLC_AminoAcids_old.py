import pandas as pd
import pathlib
import yaml
import numpy as np
import matplotlib.pyplot as plt
import math
import os

def inputs():
    experimentID = 'exp-741748223893'
    time_date_stamp = '2023-05-11 12-18-24'

    # select the kind of plots you want. You need to answer 'Yes' to at least one
    plot_final_and_initial = 'Yes'     #Answer 'Yes' or 'No'
    plot_final_minus_initial = 'Yes'   #Answer 'Yes' or 'No'

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

    ##########################################
    # Stuff that should not need to be changed
    last_characters = 'uM'  # this helps identify the columns that have the amount of the molecule in uM
    num_last_characters = 2  # number of characters in last_characters

    file_name = 'hplc_clean_data.csv'

    # location of things
    yaml_layouts_path = pathlib.Path('/Volumes/data/prod/Experiments/Layouts/')
    yaml_definitions_path = pathlib.Path('/Volumes/data/prod/Experiments/Definitions/')
    yaml_fillstock_path = pathlib.Path('/Volumes/data/prod/Chemistry/Combinations/')
    hplc_directory_location = pathlib.Path('/Volumes/data/prod/Experiments/RobotData/HPLC/Results/amino_acid/')

    # Dictionary to convert M to g/L
    conversion_dict = dict(
        [('alanine', 89.09), ('arginine', 174.2), ('asparagine', 132.12), ('aspartate', 133.1), ('citrulline', 175.2), ('cystine', 121.16),
         ('glutamate', 147.13), ('glutamine', 146.14), ('glycine', 75.07), ('histidine', 155.15), ('isoleucine', 131.17), ('leucine', 131.17),
         ('lysine', 146.19), ('methionine', 149.21), ('norvaline', 117.15), ('ornithine', 132.16), ('phenylalanine', 165.19), ('proline', 115.13),
         ('sarcosine', 89.093), ('serine', 105.09), ('threonine', 119.12), ('tryptophan', 204.23), ('tyrosine', 181.19), ('valine', 117.151)
         ])
    #######################################

    return experimentID, time_date_stamp, plot_final_and_initial, plot_final_minus_initial, other_sources, match_chemical_w_ion, last_characters, num_last_characters, file_name, hplc_directory_location, yaml_layouts_path, yaml_definitions_path, yaml_fillstock_path, conversion_dict

def open_hplc_dataframe(hplc_directory_location, experimentID, time_date_stamp, file_name):
    data_filepath =  hplc_directory_location / experimentID / time_date_stamp / file_name
    hplc_df = pd.read_csv(data_filepath, keep_default_na=False)

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

    return stacked_hplc_df

def convert_to_uM(value, aa, conversion_dict):
    value = str(value)
    if value[-2:] == ' M':  # need to leave space before M otherwise it catched mM, uM and nM
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
        print('units on molecule in fillstock or definition yaml did not have any units attached, value was set to zero: ', str)
        c = 0

    return c

def convert_to_ion(chemicals_from_yaml, match_chemical_w_ion):
    chemicals_converted = chemicals_from_yaml.copy()
    multiplication = np.ones(len(chemicals_from_yaml))
    for i in range(len(chemicals_from_yaml)):
        for j in range(len(match_chemical_w_ion)):
            if chemicals_from_yaml[i] == match_chemical_w_ion[j][0]:
                chemicals_converted[i] = match_chemical_w_ion[j][1]
                multiplication[i] = match_chemical_w_ion[j][2]
                break

    return chemicals_converted, multiplication

def get_initial_from_fillstock(hplc_df, experimentID, yaml_definitions_path, yaml_fillstock_path, match_chemical_w_ion, conversion_dict):
    definitions_yaml = open_yaml_expcond(yaml_definitions_path, experimentID, 'definitions')

    if definitions_yaml['Fillstock'][:3] == 'Med':
        fillstock_yaml = open_yaml_expcond(yaml_fillstock_path / 'Media', definitions_yaml['Fillstock'], 'fillstock')
    elif definitions_yaml['Fillstock'][:3] == 'Fee':
        fillstock_yaml = open_yaml_expcond(yaml_fillstock_path / 'Feedstock', definitions_yaml['Fillstock'], 'fillstock')
    elif definitions_yaml['Fillstock'][:3] == 'Oth':
        fillstock_yaml = open_yaml_expcond(yaml_fillstock_path / 'Other', definitions_yaml['Fillstock'], 'fillstock')
    else:
        print('fillstock file cannot be found (not in Media, Feedstock or Other folder)')

    dilution = definitions_yaml['Fillstock Dilution']
    initial_conc = np.zeros(len(hplc_df))

    chemicals_in_fillstock = list(fillstock_yaml['components'].keys())
    chemicals_in_fillstock_converted_to_ion, multiplication = convert_to_ion(chemicals_in_fillstock, match_chemical_w_ion)
    for i in range(len(hplc_df)):
        if hplc_df['Amino Acid'].iloc[i] != 'cystine':
            for j in range(len(chemicals_in_fillstock)):
                if hplc_df['Amino Acid'].iloc[i] == chemicals_in_fillstock_converted_to_ion[j]:
                    c = convert_to_uM(fillstock_yaml['components'][chemicals_in_fillstock[j]]['solvation'], hplc_df['Amino Acid'].iloc[i], conversion_dict)
                    initial_conc[i] += c * dilution * multiplication[j]
        else:
            for j in range(len(chemicals_in_fillstock)):
                if chemicals_in_fillstock_converted_to_ion[j] == 'cysteine':
                    c = convert_to_uM(fillstock_yaml['components'][chemicals_in_fillstock[j]]['solvation'], hplc_df['Amino Acid'].iloc[i], conversion_dict)
                    initial_conc[i] += c * dilution * multiplication[j] / 2

    hplc_df.insert(6, 'Concentration (uM) from Fillstock', initial_conc)

    return hplc_df

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

def expcond(file_loc, experiment):
    data = open_yaml_expcond(file_loc, experiment, 'definition')
    experiment_type = data['Design']

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

    return dict, experiment_type

def create_condition_array(condition, length):
    array = []
    for i in range(length):
        array.append(condition)

    return array

def replace_control_concentration(df, chemical, cond, variant_concentration):
    for i in range(len(df)):
        if df['Condition'].iloc[i] == cond and df['Chemicals'].iloc[i] == chemical:
            df['Concentration'].iloc[i] = variant_concentration
            break

    return df

def get_mediafeedstockother_yaml(mediafeedstockother, yaml_path):
    yaml = []

    if mediafeedstockother[:3] == 'Med':
        yaml = open_yaml_expcond(yaml_path / 'Media', mediafeedstockother, 'fillstock')
    elif mediafeedstockother[:3] == 'Fee':
        yaml = open_yaml_expcond(yaml_path / 'Feedstock', mediafeedstockother, 'fillstock')
    elif mediafeedstockother[:3] == 'Oth':
        yaml = open_yaml_expcond(yaml_path / 'Other', mediafeedstockother, 'fillstock')
    else:
        print('File did not start with Med, Fed or Oth')

    return yaml


def get_conc_from_definitions(hplc_df, experimentID, yaml_definitions_path, match_chemical_w_ion, conversion_dict, yaml_fillstock_path):
    definitions_yaml = open_yaml_expcond(yaml_definitions_path, experimentID, 'definitions')
    if definitions_yaml['Design'] == 'Desolation':
        get_conc_from_definitions_desolation(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict)
    elif definitions_yaml['Design'] == 'Dilution':
        get_conc_from_definitions_dilution(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict, yaml_fillstock_path)
    else:
        print('Definition yaml is not a Desolation or Dilution - no initial values have been entered')
        hplc_df.insert(11, 'Concentration (g/L) from Definitions', np.zeros(len(hplc_df)))

    return hplc_df

def get_conc_from_definitions_desolation(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict):
    # get list of conditions
    conditions = hplc_df['Condition'].unique()
    conditions.sort()

    # get list of chemicals
    chemicals = list(definitions_yaml['Variants'].keys())
    chemicals.sort()

    # create data frame
    conditions_array = []
    chemicals_array = []
    for i in range(len(chemicals)):
        for j in range(len(conditions)):
            conditions_array.append(conditions[j])
            chemicals_array.append(chemicals[i])
    df = pd.DataFrame({'Condition': conditions_array, 'Chemicals': chemicals_array})

    # get control concentrations - going to put these in every space, then will overwrite specific ones with 'variant' concentrations below
    control_concentrations = []
    for i in range(len(df)):
        control_concentrations.append(definitions_yaml['Variants'][df['Chemicals'].iloc[i]]['Control']['cond000'])
    df.insert(2, 'Concentration', control_concentrations)

    # get variant concentrations and overwrite control concentrations where appropriate
    for i in range(len(chemicals)):
        cond = list(definitions_yaml['Variants'][chemicals[i]]['Variants'].keys())
        for j in range(len(cond)):
            variant_concentration = definitions_yaml['Variants'][chemicals[i]]['Variants'][cond[j]]
            df = replace_control_concentration(df, chemicals[i], cond[j], variant_concentration)

    # convert chemicals to ions
    chemicals_converted_to_ion, multiplication = convert_to_ion(chemicals_array, match_chemical_w_ion)
    df.insert(3, 'Ions', chemicals_converted_to_ion)
    df.insert(4, 'Multiplication', multiplication)

    # convert concentrations to uM
    concentrations_uM = np.zeros(len(df))
    for i in range(len(df)):
        concentrations_uM[i] = convert_to_uM(df['Concentration'].iloc[i], df['Ions'].iloc[i], conversion_dict)
        concentrations_uM[i] *= multiplication[i]
        if df['Ions'].iloc[i] == 'cysteine':
            concentrations_uM[i] /= 2
    df.insert(5, 'Concentration (uM) From Definitions', concentrations_uM)

    # Add the concentrations from the definitions yaml to the data frame
    conc_from_def = np.zeros(len(hplc_df))
    for i in range(len(hplc_df)):
        for j in range(len(df)):
            if hplc_df['Condition'].iloc[i] == df['Condition'].iloc[j] and hplc_df['Amino Acid'].iloc[i] == df['Ions'].iloc[j]:
                conc_from_def[i] += df['Concentration (uM) From Definitions'].iloc[j]
                break
            elif hplc_df['Condition'].iloc[i] == df['Condition'].iloc[j] and hplc_df['Amino Acid'].iloc[i] == 'cystine' and df['Ions'].iloc[j] == 'cysteine':
                conc_from_def[i] += df['Concentration (uM) From Definitions'].iloc[j]
                break

    hplc_df.insert(8, 'Concentration (uM) from Definitions', conc_from_def)

    return hplc_df

def get_conc_from_definitions_dilution(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict, yaml_path):
    variants = list(definitions_yaml['Variants'].keys())
    num_variants = len(variants)
    count = 0
    for i in range(num_variants):
        if variants[i][:3] == 'Med' or  variants[i][:3] == 'Fee' or variants[i][:3] == 'Oth':
            count += 1

    if count == num_variants:
        hplc_df = get_conc_from_definitions_dilution_mediafeedstocksother(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict, yaml_path)
    elif count == 0:
        hplc_df = get_conc_from_definitions_dilution_chemicals(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict)
    else:
        hplc_df.insert(11, 'Concentration (g/L) from Definitions', np.zeros(len(hplc_df)))
        print('!!! the definitions yaml is of a mixed type and could not be handled, values have all been set to zero')

    return hplc_df

def get_conc_from_definitions_dilution_mediafeedstocksother(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict, yaml_path):
    amino_acids = hplc_df['Amino Acid'].unique()
    for i in range(len(amino_acids)):
        if amino_acids[i] == 'cystine':
            amino_acids[i] = 'cysteine'

    mediafeedstockother = list(definitions_yaml['Variants'].keys())
    condition_array_bigloop = []
    chemcial_array_bigloop = []
    concentration_array_bigloop = []
    for i in range(len(mediafeedstockother)):
        mediafeedstockother_yaml = get_mediafeedstockother_yaml(mediafeedstockother[i], yaml_path)
        chemicals_in_mediafeedstockother = list(mediafeedstockother_yaml['components'].keys())
        chemicals_in_mediafeedstockother_converted_to_ion, multiplication_converted_to_ion = convert_to_ion(chemicals_in_mediafeedstockother, match_chemical_w_ion)
        solvation_array_small_loop = []
        for j in range(len(chemicals_in_mediafeedstockother)):
            solvation_array_small_loop.append(mediafeedstockother_yaml['components'][chemicals_in_mediafeedstockother[j]]['solvation'])

        for j in range(len(amino_acids)):
            # get the controls
            concentration = 0
            cond = list(definitions_yaml['Variants'][mediafeedstockother[i]]['Control'].keys())[0]
            for k in range(len(chemicals_in_mediafeedstockother_converted_to_ion)):
                if amino_acids[j] == chemicals_in_mediafeedstockother_converted_to_ion[k]:
                    temp = convert_to_uM(solvation_array_small_loop[k], amino_acids[j], conversion_dict)
                    temp *= multiplication_converted_to_ion[k]
                    temp *= float(definitions_yaml['Variants'][mediafeedstockother[i]]['Control'][cond])
                    concentration += temp
                    # no break b/c there might be more than one chemical corresponding to organic acid, e.g. sodiumAcetate & potassiumAcetate both become acetate
            condition_array_bigloop.append(cond)
            chemcial_array_bigloop.append(amino_acids[j])
            concentration_array_bigloop.append(concentration)

            # get the variants
            cond = list(definitions_yaml['Variants'][mediafeedstockother[i]]['Variants'].keys())
            for m in range(len(cond)):
                concentration = 0
                for k in range(len(chemicals_in_mediafeedstockother_converted_to_ion)):
                    if amino_acids[j] == chemicals_in_mediafeedstockother_converted_to_ion[k]:
                        temp = convert_to_uM(solvation_array_small_loop[k], amino_acids[j], conversion_dict)
                        temp *= multiplication_converted_to_ion[k]
                        temp *= float(definitions_yaml['Variants'][mediafeedstockother[i]]['Variants'][cond[m]])
                        concentration += temp
                        # no break b/c there might be more than one chemical corresponding to organic acid, e.g. sodiumAcetate & potassiumAcetate both become acetate
                condition_array_bigloop.append(cond[m])
                chemcial_array_bigloop.append(amino_acids[j])
                concentration_array_bigloop.append(concentration)

    df = pd.DataFrame({'Condition': condition_array_bigloop, 'Chemicals': chemcial_array_bigloop,
                       'Concentration (uM) From Definitions': concentration_array_bigloop})

    # Add the concentrations from the definitions yaml to the data frame
    conc_from_def = np.zeros(len(hplc_df))
    for i in range(len(hplc_df)):
        for j in range(len(df)):
            if hplc_df['Condition'].iloc[i] == df['Condition'].iloc[j] and hplc_df['Amino Acid'].iloc[i] == df['Chemicals'].iloc[j]:
                conc_from_def[i] += df['Concentration (uM) From Definitions'].iloc[j]
                break
            elif hplc_df['Condition'].iloc[i] == df['Condition'].iloc[j] and hplc_df['Amino Acid'].iloc[i] == 'cystine' and df['Chemicals'].iloc[j] == 'cysteine':
                conc_from_def[i] += df['Concentration (uM) From Definitions'].iloc[j] / 2
                break

    hplc_df.insert(8, 'Concentration (uM) from Definitions', conc_from_def)

    return hplc_df

def get_conc_from_definitions_dilution_chemicals(hplc_df, definitions_yaml, match_chemical_w_ion, conversion_dict):
    chemicals = list(definitions_yaml['Variants'].keys())
    chemicals_converted_to_ion, multiplication = convert_to_ion(chemicals, match_chemical_w_ion)

    amino_acids = hplc_df['Amino Acid'].unique()
    for i in range(len(amino_acids)):
        if amino_acids[i] == 'cystine':
            amino_acids[i] = 'cysteine'

    chemicals_array = []
    concentration_control = []
    condition_array = []
    for i in range(len(chemicals)):
        for j in range(len(amino_acids)):
            if chemicals_converted_to_ion[i] == amino_acids[j]:
                cond = list(definitions_yaml['Variants'][chemicals[i]]['Control'].keys())[0]
                c = definitions_yaml['Variants'][chemicals[i]]['Control'][cond]
                c = convert_to_uM(c, amino_acids[j], conversion_dict)
                c *= multiplication[i]

                condition_array.append(cond)
                concentration_control.append(c)
                chemicals_array.append(chemicals_converted_to_ion[i])

                cond = list(definitions_yaml['Variants'][chemicals[i]]['Variants'].keys())
                for k in range(len(cond)):
                    c = definitions_yaml['Variants'][chemicals[i]]['Variants'][cond[k]]
                    c = convert_to_uM(c, amino_acids[j], conversion_dict)
                    c *= multiplication[i]

                    condition_array.append(cond[k])
                    concentration_control.append(c)
                    chemicals_array.append(chemicals_converted_to_ion[i])

                break

    df = pd.DataFrame({'Condition': condition_array, 'Chemicals': chemicals_array,
                        'Concentration (uM) From Definitions': concentration_control})

    # Add the concentrations from the definitions yaml to the data frame
    conc_from_def = np.zeros(len(hplc_df))
    for i in range(len(hplc_df)):
        for j in range(len(df)):
            if hplc_df['Condition'].iloc[i] == df['Condition'].iloc[j] and hplc_df['Amino Acid'].iloc[i] == df['Chemicals'].iloc[j]:
                conc_from_def[i] += df['Concentration (uM) From Definitions'].iloc[j]
                break
            elif hplc_df['Condition'].iloc[i] == df['Condition'].iloc[j] and hplc_df['Amino Acid'].iloc[i] == 'cystine' and df['Chemicals'].iloc[j] == 'cysteine':
                conc_from_def[i] += df['Concentration (uM) From Definitions'].iloc[j] / 2
                break

    hplc_df.insert(8, 'Concentration (uM) from Definitions', conc_from_def)

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
    hplc_df['Initial Concentration (uM)'] = hplc_df['Concentration (uM) from Fillstock'] + hplc_df['Concentration (uM) from Definitions'] + hplc_df['Concentration (uM) from Other Sources']

    return hplc_df

def add_xaxis_info(yaml_definitions_path, experimentID, hplc_df):
    expcond_dict, experiment_type = expcond(yaml_definitions_path, experimentID)
    conditions = []

    if experiment_type == 'Dilution':
        for i in range(len(hplc_df)):
            conditions.append(hplc_df['Inoculant'].iloc[i] + ' ' + expcond_dict[hplc_df['Condition'].iloc[i]])
    else:
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
    fig.subplots_adjust(wspace=0, hspace=0.4)
    fig.set_figwidth(12)
    fig.set_figheight(num_y * 2)

    x = 0; y= 0
    for i in range(aa_number):
        subset = hplc_df.loc[(hplc_df['Amino Acid'] == aa_tested[i])]

        subset_mean = subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
        subset_mean = add_condition_for_sorting(hplc_df, subset_mean)
        subset_mean.sort_values(by=['Condition'], inplace=True)

        subset_std = subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)
        subset_std = add_condition_for_sorting(hplc_df, subset_std)
        subset_std.sort_values(by=['Condition'], inplace=True)

        if subset_mean.empty == False:
            if pd.isna(subset_std['HPLC Measured Concentration (uM)'].iloc[0]) == False:
                axs[y, x].scatter(subset_mean['Experimental Conditions'], subset_mean['HPLC Measured Concentration (uM)'], marker='s', color='blue')
                axs[y, x].scatter(subset_mean['Experimental Conditions'], subset_mean['Initial Concentration (uM)'], marker='s', color='orange')
                axs[y, x].errorbar(subset_mean['Experimental Conditions'], subset_mean['HPLC Measured Concentration (uM)'], yerr=subset_std['HPLC Measured Concentration (uM)'], color='blue', linestyle = 'None')
                axs[y, x].errorbar(subset_mean['Experimental Conditions'], subset_mean['Initial Concentration (uM)'], yerr=subset_std['Initial Concentration (uM)'], color='orange', linestyle = 'None')
            else:
                axs[y, x].scatter(subset_mean['Experimental Conditions'], subset_mean['HPLC Measured Concentration (uM)'], marker='s', color='blue')
                axs[y, x].scatter(subset_mean['Experimental Conditions'], subset_mean['Initial Concentration (uM)'], marker='s', color='orange')

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
    plt.subplots_adjust(top=0.9, bottom=0.25, left=0.09, right=0.91)
    fig.suptitle('Orange: Initial Concentration (uM) \n Blue: Concentration Measured by HPLC (uM) \n ')
    plt.show()
    #plt.savefig('/Users/jennifer/Desktop/grid_plot_final_minus_initial.png')

def add_condition_for_sorting(hplc_df, df):
    condition = []
    for i in range(len(df)):
        for j in range(len(hplc_df)):
            if df['Experimental Conditions'].iloc[i] == hplc_df['Experimental Conditions'].iloc[j]:
                condition.append(hplc_df['Condition'].iloc[j])
                break

    df.insert(0, 'Condition', condition)
    return df

def grid_plot_final_minus_initial(hplc_df):
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
    fig.subplots_adjust(wspace=0, hspace=0.4)
    fig.set_figwidth(12)
    fig.set_figheight(num_y * 2)

    x = 0; y= 0
    for i in range(aa_number):
        subset = hplc_df.loc[(hplc_df['Amino Acid'] == aa_tested[i])]

        subset_mean = subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
        subset_mean = add_condition_for_sorting(hplc_df, subset_mean)
        subset_mean.sort_values(by=['Condition'], inplace=True)

        subset_std = subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)
        subset_std = add_condition_for_sorting(hplc_df, subset_std)
        subset_std.sort_values(by=['Condition'], inplace=True)

        width = 0.25

        if subset_mean.empty == False:
            if pd.isna(subset_std['Final Minus Initial (uM)'].iloc[0]) == False:
                axs[y, x].bar(subset_mean['Experimental Conditions'], subset_mean['Final Minus Initial (uM)'], width, yerr=subset_std['Final Minus Initial (uM)'], color='blue')
            else:
                axs[y, x].bar(subset_mean['Experimental Conditions'], subset_mean['Final Minus Initial (uM)'], width, color='blue')

        axs[y, x].set_title(aa_tested[i])
        axs[y, x].tick_params(axis='x', labelrotation = 90)
        axs[y, x].axhline(0, color='black', linestyle='dotted')

        x += 1
        if x > 4:
            x = 0
            y += 1

    if x > 0:
        for i in range(5 - x):
            axs[y, x].tick_params(axis='x', labelrotation=90)

            x += 1

    fig.supylabel('Concentration (uM): Final - Initial')
    plt.subplots_adjust(top=0.95, bottom=0.32, left=0.09, right=0.95)
    plt.show()
    #plt.savefig('/Users/jennifer/Desktop/grid_plot_final_minus_initial.png')

def individual_plots(hplc_df):
    # get list of amino acids
    aa_tested = pd.unique(hplc_df['Amino Acid'])
    aa_tested.sort()

    for i in range(len(aa_tested)):
        # set up the figure
        fig, axs = plt.subplots(1, 1)
        fig.set_figwidth(8)
        fig.set_figheight(8)

        subset = hplc_df.loc[(hplc_df['Amino Acid'] == aa_tested[i])]

        subset_mean = subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
        subset_mean = add_condition_for_sorting(hplc_df, subset_mean)
        subset_mean.sort_values(by=['Condition'], inplace=True)

        subset_std = subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)
        subset_std = add_condition_for_sorting(hplc_df, subset_std)
        subset_std.sort_values(by=['Condition'], inplace=True)

        if subset_mean.empty == False:
            if pd.isna(subset_std['HPLC Measured Concentration (uM)'].iloc[0]) == False:
                axs.scatter(subset_mean['Experimental Conditions'], subset_mean['HPLC Measured Concentration (uM)'], marker='s', color='blue')
                axs.scatter(subset_mean['Experimental Conditions'], subset_mean['Initial Concentration (uM)'], marker='s', color='orange')
                axs.errorbar(subset_mean['Experimental Conditions'], subset_mean['HPLC Measured Concentration (uM)'], yerr=subset_std['HPLC Measured Concentration (uM)'], color='blue', linestyle = 'None')
                axs.errorbar(subset_mean['Experimental Conditions'], subset_mean['Initial Concentration (uM)'], yerr=subset_std['Initial Concentration (uM)'], color='orange', linestyle = 'None')
            else:
                axs.scatter(subset_mean['Experimental Conditions'], subset_mean['HPLC Measured Concentration (uM)'], marker='s', color='blue')
                axs.scatter(subset_mean['Experimental Conditions'], subset_mean['Initial Concentration (uM)'], marker='s', color='orange')

        axs.tick_params(axis='x', labelrotation = 90)
        axs.set_ylabel('Concentration (uM)')

        plt.subplots_adjust(top=0.85, bottom=0.35, left=0.1, right=0.9)
        plt.title('Orange: Initial Concentration (uM) \n Blue: Concentration Measured by HPLC (uM) \n \n' + aa_tested[i])
        plt.show()

def individual_plots_final_minus_initial(hplc_df):
    # get list of amino acids
    aa_tested = pd.unique(hplc_df['Amino Acid'])
    aa_tested.sort()

    for i in range(len(aa_tested)):
        # set up the figure
        fig, axs = plt.subplots(1, 1)
        fig.set_figwidth(9)
        fig.set_figheight(6)

        subset = hplc_df.loc[(hplc_df['Amino Acid'] == aa_tested[i])]

        subset_mean = subset.groupby(['Experimental Conditions'], as_index=False).mean(numeric_only=True)
        subset_mean = add_condition_for_sorting(hplc_df, subset_mean)
        subset_mean.sort_values(by=['Condition'], inplace=True)

        subset_std = subset.groupby(['Experimental Conditions'], as_index=False).std(numeric_only=True)
        subset_std = add_condition_for_sorting(hplc_df, subset_std)
        subset_std.sort_values(by=['Condition'], inplace=True)

        width = 0.25

        leg = 'Concentration (uM): Final - Initial'
        if subset_mean.empty == False:
            if pd.isna(subset_std['Final Minus Initial (uM)'].iloc[0]) == False:
                axs.bar(subset_mean['Experimental Conditions'], subset_mean['Final Minus Initial (uM)'], width, yerr=subset_std['Final Minus Initial (uM)'],
                        color='blue', label=leg)
            else:
                axs.bar(subset_mean['Experimental Conditions'], subset_mean['Final Minus Initial (uM)'], width, color='blue', label=leg)

        leg = 'Initial Concentration (uM)'
        axs.scatter(subset_mean['Experimental Conditions'], subset_mean['Initial Concentration (uM)'], marker='_', color='red', s=64, label=leg)

        axs.tick_params(axis='x', labelrotation=90)
        axs.set_ylabel('Concentration (uM): Final - Initial')
        axs.axhline(0, color='black', linestyle='dotted')
        axs.legend(bbox_to_anchor=(1, 1), loc='upper left')

        plt.subplots_adjust(top=0.95, bottom=0.35, left=0.1, right=0.68)
        plt.title(aa_tested[i])
        plt.show()
        # plt.savefig('/Users/jennifer/Desktop/finalminuinitial_' + str(i))

def get_final_minus_initial(hplc_df):
    hplc_df['Final Minus Initial (uM)'] = hplc_df['HPLC Measured Concentration (uM)'] - hplc_df['Initial Concentration (uM)']

    return hplc_df

def save(hplc_df, experimentID):
    hplc_df.to_csv(experimentID + '.csv')

def main():
    # read in inputs
    experimentID, date_time_stamp, plot_final_and_initial, plot_final_minus_initial, other_sources, match_chemical_w_ion, last_characters, num_last_characters, file_name, hplc_directory_location, yaml_layouts_path, yaml_definitions_path, yaml_fillstock_path, conversion_dict = inputs()
    # get HPLC data
    hplc_df = open_hplc_dataframe(hplc_directory_location, experimentID, date_time_stamp, file_name)
    # remove calibration samples
    hplc_df = remove_calibration_samples(hplc_df)
    # add condition to data frame
    hplc_df = add_condition(hplc_df, experimentID, yaml_layouts_path)
    # get amino acids & stack columns
    hplc_df = stack_columns(hplc_df, last_characters, num_last_characters)
    # add experimental conditions - details, e.g. serine 500 uM not cond001
    hplc_df = add_xaxis_info(yaml_definitions_path, experimentID, hplc_df)
    # get concentrations from fillstock
    hplc_df = get_initial_from_fillstock(hplc_df, experimentID, yaml_definitions_path, yaml_fillstock_path, match_chemical_w_ion, conversion_dict)
    # get concentrations from definitions yaml
    hplc_df = get_conc_from_definitions(hplc_df, experimentID, yaml_definitions_path, match_chemical_w_ion, conversion_dict, yaml_fillstock_path)
    # get concentrations from other sources
    hplc_df = get_cond_from_other_sources(hplc_df, other_sources)
    # get initial concentrations of amino acids
    hplc_df = get_initial_concentration(hplc_df)
    #plot stuff
    if plot_final_and_initial == 'Yes':
        # grid plot
        grid_plot(hplc_df)
        # individual plots
        individual_plots(hplc_df)
        # final minus initial
    if plot_final_minus_initial == 'Yes':
        # get final minus initial concentration
        hplc_df = get_final_minus_initial(hplc_df)
        # grid plot final minus initial
        grid_plot_final_minus_initial(hplc_df)
        # individual plots final minus initial
        individual_plots_final_minus_initial(hplc_df)
    # save data frame
    save(hplc_df, experimentID)
    hplc_df.to_csv('/Users/jennifer/Desktop/testing_aa_old.csv')

main()