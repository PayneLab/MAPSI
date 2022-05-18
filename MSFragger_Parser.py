import pandas as pd
import numpy as np
import os
import pyteomics.mzml
import spectrum_utils.spectrum as sus
from pathlib import Path
import json

# loading functions
def load_mzml_df(mzml_file_path):
    # use pyteomics.mzml.read() to generate an iterator over the dicts with spectrum properties
    mzml_dicts = pyteomics.mzml.read(source=mzml_file_path)

    # load dataframe from the list of mzml dictionaires
    # drop the extra index column
    mzml_df = pd.DataFrame(mzml_dicts).drop(columns='index')

    # create a new dataframe containing only the ms/ms scans
    relevant_info = mzml_df.loc[(mzml_df['ms level'] == 2)]

    # reset the index to make up for the ms scans that were not included in this database
    relevant_info = relevant_info.reset_index(drop=True)

    # drop irrelevent columns (Note: We can change this if needed.)
    columns_to_drop = ["spectrum title", "count", "positive scan", "centroid spectrum", "defaultArrayLength", "MSn spectrum", "dataProcessingRef", "scanList", "MS1 spectrum", "ms level"]
    relevant_info = relevant_info.drop(columns=columns_to_drop)

    # create a new "Scan Number" column
    # the scan number info is contained within the "id" column so we will pull out the scan number and then delete the "id" column
    relevant_info["temp_split_column"] = (relevant_info["id"].str.split(" "))
    relevant_info["Scan Number"] = relevant_info["temp_split_column"].map(lambda x:x[2]).str.replace("scan=", "")
    relevant_info["Scan Number"] = relevant_info["Scan Number"].apply(pd.to_numeric)
    relevant_info = relevant_info.drop(columns=['temp_split_column', 'id'])

    # next, we'll want to pull out some info about the precursor in the "precursorList" columm
    # we will store the info we want under the "precursor info" column
    # then we'll drop the "precursorList" column
    relevant_info["precursor info"] = relevant_info["precursorList"].map(lambda x:x['precursor'][0]['selectedIonList']['selectedIon'][0]).astype(str)
    relevant_info = relevant_info.drop(columns=["precursorList"])

    # the precursor info is stored as a string in the format of a dictionary
    # json.loads() requires " instead of ' so we will fix that then convert the string into a dictionary
    dict_list = relevant_info["precursor info"].tolist()

    for index, dictionary in enumerate(dict_list):
        dictionary = dictionary.replace("'", '"')
        dict_list[index] = json.loads(dictionary)

    # we'll then load the dictionary data into a temporary dataframe
    three_column_df = pd.DataFrame.from_dict(dict_list)

    # next, we'll concatenate these two dataframes along the columns based on the index
    complete_mzml_df = pd.concat([relevant_info, three_column_df], axis="columns")

    # drop the "precursor info" column because we don't need it anymore
    complete_mzml_df = complete_mzml_df.drop(columns=['precursor info'])

    # as the scan number is the index we care about, we will use it as our index in the dataframe
    complete_mzml_df= complete_mzml_df.set_index("Scan Number")

    return complete_mzml_df 
def load_psm_df_msfragger(psm_file_path):
    # read in the psm file as a dataframe
    psm_df = pd.read_table(psm_file_path)

    # split the "Spectrum" column into a list at each period and store it under 
    # the "temp_split_column"
    psm_df["temp_split_column"] = psm_df["Spectrum"].str.split(".")
    # store the element located at index 1 of the "temp_split_column" in a 
    # "Scan Number" column
    psm_df["Scan Number"] = psm_df["temp_split_column"].map(lambda x:x[1]).apply(pd.to_numeric)
    # drop unneeded columns
    columns_to_drop = ['Spectrum', 'temp_split_column','Mapped Genes', 'Mapped Proteins']
    psm_df = psm_df.drop(columns=columns_to_drop)

    # rename 2 headers to match MM file formats
    psm_df = psm_df.rename({'Protein ID': 'Protein Accession', 'Spectrum File': 'File Name'}, axis=1)

    # drop duplicates
    psm_df = psm_df.drop_duplicates(subset=["Scan Number"], keep="first")
    
    return psm_df
def load_msfragger_peptideQ(peptideQ_file_path):
    # read the peptideQ file into a pandas dataframe
    peptideQ_dataframe = pd.read_table(peptideQ_file_path, delimiter='\t')

    # rename the "Protein Groups" header so we can use the df.merge function later
    peptideQ_dataframe = peptideQ_dataframe.rename({'Protein ID': 'Protein Accession', 'Peptide Sequence': 'Peptide'}, axis=1)
    
    return peptideQ_dataframe
def load_msfragger_protein(protein_file_path):
    # load the protein file into a pandas dataframe
    protein_dataframe = pd.read_table(protein_file_path)

    # Rename the "Protein ID" column to faciliate merging
    protein_dataframe = protein_dataframe.rename({'Protein ID': 'Protein Accession'}, axis=1)
    
    return protein_dataframe

# joining functions
def join_peptideQ_and_protein_dataframes(protein_df, peptideQ_df):
    # join based on the "Protein Accession"
    joined_dataframe = peptideQ_df.merge(right=protein_df, on="Protein Accession", how='inner', suffixes=('_protein', '_peptide'))
    
    return joined_dataframe
def join_psm_and_peptideQ_dataframes(psm_df, peptideQ_df):
    # find all the duplicate columns that are not the 'Peptide'
    duplicate_columns = []
    for column in psm_df.columns:
        if column in peptideQ_df.columns and column != 'Peptide':
            duplicate_columns.append(column)
        
    psm_df = psm_df.drop(columns=duplicate_columns)

    # join based on the "Base Sequence"
    joined_dataframe = psm_df.merge(right=peptideQ_df, on="Peptide", how='inner', )

    # generate multiIndex
    #joined_dataframe = joined_dataframe.set_index(['File Name','Protein Accession','Peptide', 'Scan Number']).drop(columns=["Protein Groups"])

    return joined_dataframe

# simple controller functions
def psm_and_peptideQ_controller(psm_file_path, peptideQ_file_path, columns_to_keep=None):
    
    ''' Joins a psm and a Peptide Quantification file. Files are joined into a pandas dataframe and saved as a tsv.
    
    Required Parameters:
        * psm_file_path: File path to the psm file
        * peptideQ_file_path: File path to the Peptide Quantification file
        * psm_and_peptideQ_file_path: Output file path
        
    Optional Parameters: 
        * columns_to_keep: List of columns to include in the dataframe. Note that column names may vary based on whether your files were generated with MetaMorpheus or MSFragger.'''
        
    # load dataframes
    peptideQ_df = load_msfragger_peptideQ(peptideQ_file_path=peptideQ_file_path)
    psm_df = load_psm_df_msfragger(psm_file_path)
    
    # join dataframes
    joined_df = join_psm_and_peptideQ_dataframes(psm_df, peptideQ_df)

    # select all columns to keep, if this parameter was not passed in, return dataframe with all columns
    if columns_to_keep != None:
        joined_df = joined_df[columns_to_keep]
    
    return joined_df
def mzml_and_psm_controller(mzml_file_path,psm_file_path, columns_to_keep=None):

    ''' Joins an mzml and psm file. Files are joined into a pandas dataframe and saved as a tsv.
    
    Required Parameters:
        * mzml_file_path: File path to the mzML file
        * psm_file_path: File path to the psm file
        * mzml_and_psm_file_path: Output file path
        
    Optional Parameters: 
        * columns_to_keep: List of columns to include in the dataframe. Note that column names may vary based on whether your files were generated with MetaMorpheus or MSFragger.'''

    # load psm dataframe based on psm file type
    psm_dataframe = load_psm_df_msfragger(psm_file_path=psm_file_path)

    # load mzML dataframe
    mzml_dataframe = load_mzml_df(mzml_file_path)

    # merge datafames based on "Scan Number"
    joined_dataframe = mzml_dataframe.join(other=psm_dataframe, on='Scan Number', how='inner')

    # select all columns to keep, if this parameter was not passed in, return dataframe with all columns
    if columns_to_keep != None:
        joined_dataframe = joined_dataframe[columns_to_keep]
  
    return joined_dataframe
def peptideQ_and_protein_controller (peptideQ_file_path, protein_file_path, columns_to_keep=None):
    
    ''' Joins a Peptide and Protein Quantification files. Files are joined into a pandas dataframe and saved as a tsv.
    
    Required Parameters:
        * peptideQ_file_path: File path to the Peptide Quantification file
        * protein_file_path: File path to the Protein Quantification file
        * peptideQ_and_protein_file_path: Output file path
        
    Optional Parameters: 
        * columns_to_keep: List of columns to include in the dataframe. Note that column names may vary based on whether your files were generated with MetaMorpheus or MSFragger. '''

    # load dataframes
    peptideQ_df = load_msfragger_peptideQ(peptideQ_file_path=peptideQ_file_path)
    protein_df = load_msfragger_protein(protein_file_path=protein_file_path)
    
    # create joined dataframe and save as csv
    joined_df = join_peptideQ_and_protein_dataframes(protein_df=protein_df, peptideQ_df=peptideQ_df)

    # select all columns to keep, if this parameter was not passed in, return dataframe with all columns
    if columns_to_keep != None:
        joined_df = joined_df[columns_to_keep]

    return joined_df
def mzml_psm_and_peptideQ_controller(mzml_file_path, psm_file_path, peptideQ_file_path, columns_to_keep=None):
    # merge mzml and psm dataframes
    mzml_and_psm_df = mzml_and_psm_controller(mzml_file_path=mzml_file_path, psm_file_path=psm_file_path, columns_to_keep=None)
    
    # load peptideQ dataframe
    peptideQ_df = load_msfragger_peptideQ(peptideQ_file_path=peptideQ_file_path)
    
    # merge dataframes
    merged_df = join_psm_and_peptideQ_dataframes(psm_df=mzml_and_psm_df, peptideQ_df=peptideQ_df)
    # get rid of duplicate entries and columns
    if 'Peptide_x' in merged_df.columns and 'Peptide_y' in merged_df.columns and 'Protein Accession_x' in merged_df.columns and 'Protein Accession_y' in merged_df.columns:
        merged_df = merged_df.drop_duplicates(subset=["Scan Number"], keep="first")
        merged_df = merged_df.rename({'Peptide_x' : 'Peptide', 'Protein Accession_x': 'Protein Accession'}, axis=1)
        merged_df = merged_df.drop(columns=['Peptide_y', 'Protein Accession_y'])

    print (merged_df.columns)
    merged_df = merged_df.set_index(['File Name', 'Protein Accession', 'Peptide', 'Scan Number'])
    
    # select all columns to keep, if this parameter was not passed in, return dataframe with all columns
    if columns_to_keep != None:
        merged_df = merged_df[columns_to_keep]

    return merged_df
def psm_peptideQ_and_protein_controller(psm_file_path, peptideQ_file_path, protein_file_path,columns_to_keep=None):
    # merge psm and peptideQ dataframes
    psm_and_peptideQ_df = psm_and_peptideQ_controller(psm_file_path=psm_file_path, peptideQ_file_path=peptideQ_file_path)

    # load protein dataframe
    protein_df = load_msfragger_protein(protein_file_path=protein_file_path)

    # merge all dataframes
    merged_df = join_peptideQ_and_protein_dataframes(protein_df=protein_df, peptideQ_df=psm_and_peptideQ_df)

    # select all columns to keep, if this parameter was not passed in, return dataframe with all columns
    if columns_to_keep != None:
        merged_df = merged_df[columns_to_keep]

    return merged_df
def master_df_controller(mzml_file_path, psm_file_path, peptideQ_file_path, protein_file_path, columns_to_keep=None):
    # merge mzml and psm dataframes
    mzml_and_psm_df = mzml_and_psm_controller(mzml_file_path=mzml_file_path, psm_file_path=psm_file_path)

    # merge peptideQ and protein dataframes
    peptideQ_and_protein_df = peptideQ_and_protein_controller(peptideQ_file_path=peptideQ_file_path, protein_file_path=protein_file_path)

    # merge all data
    merged_df = join_psm_and_peptideQ_dataframes(psm_df=mzml_and_psm_df, peptideQ_df=peptideQ_and_protein_df)
    
    return merged_df

# general parsing helper functions
def assign_file_types(input_files):

    file_path_list = [None, None, None, None]

    for index, file_path in enumerate(input_files):
        file_name = os.path.basename(file_path).lower()
        if 'protein' in file_name or 'prot' in file_name:
            file_path_list[3] = file_path
        elif 'peptide' in file_name or 'pep' in file_name:
            file_path_list[2] = file_path
        elif 'psm' in file_name:
            file_path_list[1] = file_path
        elif 'mzml' in file_path:
            file_path_list[0] = file_path
        else:
            # most of the mzML files do not contain the word 'mzml'
            file_path_list[0] = file_path
    
    return file_path_list
def generate_bool_file_list(interpreted_file_list):
    bool_file_list = []
    for file_path in interpreted_file_list:
        if file_path:
            bool_file_list.append(1)
        else:
            bool_file_list.append(0)

    bool_file_list = str(bool_file_list).replace(" ", "")
    return bool_file_list
def load_call_dictionary(interpreted_file_list, columns_to_keep):
    call_dictionary = {}

    # load one dataframe
    call_dictionary['[1,0,0,0]'] =  [load_mzml_df, interpreted_file_list[0], ['Scan Number']]
    call_dictionary['[0,1,0,0]'] = [load_psm_df_msfragger, interpreted_file_list[1], ['Protein Accession','Peptide', 'Scan Number']]
    call_dictionary['[0,0,1,0]'] = [load_msfragger_peptideQ, interpreted_file_list[2], ['Protein Accession', 'Peptide']]
    call_dictionary['[0,0,0,1]'] = [load_msfragger_protein, interpreted_file_list[3], ['Protein Accession']]
    
    # load and merge 2 dataframes
    call_dictionary['[1,1,0,0]'] = [mzml_and_psm_controller, [interpreted_file_list[0], interpreted_file_list[1], columns_to_keep],
    ['Protein Accession','Peptide', 'Scan Number']]
    call_dictionary['[0,1,1,0]'] = [psm_and_peptideQ_controller, [interpreted_file_list[1], interpreted_file_list[2], columns_to_keep],
    ['Protein Accession','Peptide', 'Scan Number']]
    call_dictionary['[0,0,1,1]'] = [peptideQ_and_protein_controller, [interpreted_file_list[2], interpreted_file_list[3], columns_to_keep],
    ['Protein Accession', 'Peptide']]

    # load and merge 3 dataframes
    call_dictionary['[1,1,1,0]'] = [mzml_psm_and_peptideQ_controller, [interpreted_file_list[0],interpreted_file_list[1], interpreted_file_list[2], columns_to_keep],
    ['Protein Accession','Peptide', 'Scan Number']]
    call_dictionary['[0,1,1,1]'] = [psm_peptideQ_and_protein_controller, [interpreted_file_list[1],interpreted_file_list[2], interpreted_file_list[3], columns_to_keep],
    ['Protein Accession','Peptide', 'Scan Number']]

    # load and merge all 3 dataframes
    call_dictionary['1,1,1,1'] = [master_df_controller, [interpreted_file_list[0],interpreted_file_list[1], interpreted_file_list[2], interpreted_file_list[3],columns_to_keep],
    ['Protein Accession','Peptide', 'Scan Number']]

    return call_dictionary
def select_columns_to_keep(user_dataframe, columns_to_keep):
    if columns_to_keep != None:
        # check that all of the columns listed exist in the user database
        mistake = False
        for column in columns_to_keep:
            if column not in user_dataframe.columns:
                mistake = True
                print(f"{column} column does not exist in this dataframe.")
        if mistake:
            print("To ensure that you are given all the infomation needed, the entire database will be returned.")
        else:
            user_dataframe = user_dataframe[columns_to_keep]
    return user_dataframe
def select_multiIndex(user_dataframe, multiIndex, default_multiIndex):
    if multiIndex != None:
        # check that the columns exist
        adjusted_default_multiIndex = multiIndex.copy()
        for column in multiIndex:
            if column not in user_dataframe.columns:
                adjusted_default_multiIndex.remove(column)
                print(f"{column} column does not exist in this dataframe and cannot be used as an index")
        user_dataframe = user_dataframe.set_index(adjusted_default_multiIndex)
    else:
        adjusted_default_multiIndex = default_multiIndex.copy()
        for column in default_multiIndex:
            if column not in user_dataframe.columns:
                adjusted_default_multiIndex.remove(column)
        user_dataframe = user_dataframe.set_index(adjusted_default_multiIndex)
    return user_dataframe
def save_df(joined_dataframe, file_path):
    joined_dataframe.to_csv(file_path, sep="\t", index=False)
    print(f'Dataframe saved.')

# msfragger parser
def parse_files(input_files, output_file_path, columns_to_keep=None, multiIndex=None, proteins_to_keep=None, peptides_to_keep=None):
    # check rows_to_keep parameter
    
    # interpret the file list
    interpreted_file_list = assign_file_types(input_files)

    # check that the input_files contained a valid list of files
    bool_file_list = generate_bool_file_list(interpreted_file_list)
    call_dictionary = load_call_dictionary(interpreted_file_list, columns_to_keep)
    if bool_file_list not in call_dictionary.keys():
        print("Invalid file combination")
        return False
    
    # using the call dictionary, call the correct controller function with the associated parameters
    function, parameters, default_multiIndex = call_dictionary[bool_file_list]
    if type(parameters) == list:
        user_dataframe = function(*parameters)
    else:
        user_dataframe = function(parameters)

    # columns to keep
    user_dataframe = select_columns_to_keep(user_dataframe=user_dataframe, columns_to_keep=columns_to_keep)

    # proteins to keep
    if proteins_to_keep != None:
        user_dataframe = user_dataframe.loc[user_dataframe['Protein Accession'].isin(proteins_to_keep)]

    # peptides to keep
    if peptides_to_keep != None:
        user_dataframe = user_dataframe.loc[user_dataframe['Peptide'].isin(peptides_to_keep)]
    
    # multiIndexing
    user_dataframe = select_multiIndex(user_dataframe=user_dataframe, multiIndex=multiIndex, default_multiIndex=default_multiIndex)
    
    save_df(joined_dataframe=user_dataframe, file_path=output_file_path)

    return user_dataframe