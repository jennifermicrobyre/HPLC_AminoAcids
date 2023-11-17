import pathlib
import pandas as pd

#location of file
data_dir = pathlib.Path('/Users/jennifer/Documents/Analysis/HPLC_OrganicAcids/exp-238539175652/')
header_filepath = data_dir / 'REOT00.CSV'
result_filepath = data_dir / 'REOT01.CSV'
cleaned_filepath = data_dir / 'hplc_raw_data.csv'
#header_filepath = data_dir / 'HPLC' / 'REPORT00.CSV'
#result_filepath = data_dir / 'HPLC' / 'REPORT01.CSV'
#cleaned_filepath = data_dir / 'HPLC' / 'MERGEDREPORT.CSV'

#merge file header and data
header_df = pd.read_csv(header_filepath, encoding='UTF-16')
headers = header_df.iloc[:,2].to_list() #brittle encoding to a specific column
headers = [header.replace('|','_').lower() for header in headers] #clean and standardize
result_df = pd.read_csv(result_filepath, encoding='UTF-16', header=None)
result_df.columns = headers #set headers from other file

#export merged file
result_df.to_csv(cleaned_filepath, sep = ',', index=False, encoding='utf-8')