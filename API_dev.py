# Get all the relevant mzML data from a file
import pandas as pd
import numpy as np
import os
import pyteomics.mzml
import spectrum_utils.spectrum as sus
from pathlib import Path

# hmm, maybe when we load the data we should create a dictionary containing the 
# spectrum_dict for each scan. That way we won't repeat that part of the work.
# it also wouldn't be a bad idea to group all of that into a function that 
# loads the data (it can reference other loading functions).

mzml_file_paths = ["C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\Ex_Auto_J3_30umTB_2ngQC_60m_1.mzML"]
peptide_file_paths = ["C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\Ex_Auto_J3_30umTB_2ngQC_60m_1-calib_Peptides.psmtsv"]
protein_file_paths = ["C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\Ex_Auto_J3_30umTB_2ngQC_60m_1-calib_ProteinGroups.tsv"]

def load_tsv_data(file_paths):
    dataframe_dict = {}
    for file_path in file_paths:
        dataframe_dict[file_path] = pd.read_table(file_path, delimiter='\t')
    return dataframe_dict

# create the peptide and protein dataframe dictionaires
peptide_dataframe_dict = load_tsv_data(peptide_file_paths)
protein_dataframe_dict = load_tsv_data(protein_file_paths)

def load_mzml_data(mzml_file_paths, peptide_file_paths):
    dict_of_spectrum_dict = {}
    # here we will be assuming that the peptide file paths list is aligned with
    # the mzml path list
    for index, mzml_file_path in enumerate(mzml_file_paths):
        peptide_dataframe = peptide_dataframe_dict[peptide_file_paths[index]]
        for index, row in peptide_dataframe.iterrows():
            scan = row['Scan Number']
            dict_of_spectrum_dict[scan] = get_mzml_ms2_scan_data(scan=scan, mzML_file_path=mzml_file_path) 
    return dict_of_spectrum_dict

def get_mzml_ms2_scan_data(scan, mzML_file_path):
    # create the mzml object
    mzml = pyteomics.mzml.MzML(mzML_file_path)
    my_id = 'controllerType=0 controllerNumber=1 scan='+ str(scan)
    spectrum_dict = mzml.get_by_id(my_id)
    return spectrum_dict

dict_of_scan_data = load_mzml_data(mzml_file_paths=mzml_file_paths, peptide_file_paths=peptide_file_paths)

def printDictionary(myDict):
    for key in myDict.keys():
        print(f"{key}: {myDict[key]}")

printDictionary(dict_of_scan_data)
#def graph_ms2_spectrum(scan, mzML_file_path):




def mzml_helper(scan, mzml):

    my_id = 'controllerType=0 controllerNumber=1 scan='+ str(scan)
    spectrum_dict = mzml.get_by_id(my_id)

    if spectrum_dict['ms level'] != 2:
        print('This is a MS1 scan. The information will not be processed.')
        return

    spectrum_id = spectrum_dict['id']
    retention_time = (spectrum_dict['scanList']['scan'][0].get('scan start time', -1))
    mz_array = list(spectrum_dict['m/z array'])
    intensity_array = list(spectrum_dict['intensity array'])

    if 'precursorList' in spectrum_dict.keys():
        precursor = spectrum_dict['precursorList']['precursor'][0]
        precursor_ion = precursor['selectedIonList']['selectedIon'][0]
        precursor_mz = precursor_ion['selected ion m/z']
        if 'peak intensity' in precursor_ion:
            precursor_intensity =  precursor_ion['peak intensity']
        else:
            precursor_intensity = None
        if 'charge state' in precursor_ion:
            precursor_charge = int(precursor_ion['charge state'])
        elif 'possible charge state' in precursor_ion:
            precursor_charge = int(precursor_ion['possible charge state'])
        else:
            precursor_charge = 'NAN'
    else:
        precursor_intensity = None
        
    all_info = [mz_array,intensity_array,precursor_intensity]

    return all_info



