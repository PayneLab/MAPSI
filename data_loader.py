import pandas as pd
import numpy as np
import os
import pyteomics.mzml
import spectrum_utils.spectrum as sus
from pathlib import Path




def download_parsed_psm():

    list_of_file_names = [ 'bulk_rep1',"bulk_rep2","bulk_rep3",
                            '2ng_rep1','2ng_rep2','2ng_rep3','2ng_rep4','2ng_rep5','2ng_rep6',
                            '0.2ng_rep1','0.2ng_rep2','0.2ng_rep3','0.2ng_rep4','0.2ng_rep5','0.2ng_rep6',
                            'sc_rep1','sc_rep2','sc_rep3','sc_rep4','sc_rep5']

    # list_of_file_names = ['2ng_rep1']
    mm_files = {}
    #bulk
    mm_files["bulk_rep1"] = "data/MetaMorpheus_output/Project_PXD011163/OR11_20160122_PG_HeLa_CVB3_CT_A-calib_PSMs.psmtsv.gz"
    mm_files["bulk_rep2"] = "data/MetaMorpheus_output/Project_PXD011163/OR11_20160122_PG_HeLa_CVB3_CT_B-calib_PSMs.psmtsv.gz"
    mm_files["bulk_rep3"] = "data/MetaMorpheus_output/Project_PXD011163/OR11_20160122_PG_HeLa_CVB3_CT_C-calib_PSMs.psmtsv.gz"

    #2ng
    mm_files["2ng_rep1"] = "data/MetaMorpheus_output/2ng/Ex_Auto_J3_30umTB_2ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep2"] = "data/MetaMorpheus_output/2ng/Ex_Auto_J3_30umTB_2ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep3"] = "data/MetaMorpheus_output/2ng/Ex_Auto_K13_30umTA_2ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep4"] = "data/MetaMorpheus_output/2ng/Ex_Auto_K13_30umTA_2ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep5"] = "data/MetaMorpheus_output/2ng/Ex_Auto_W17_30umTB_2ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["2ng_rep6"] = "data/MetaMorpheus_output/2ng/Ex_Auto_W17_30umTB_2ngQC_60m_2-calib_PSMs.psmtsv.gz"

    #.2ng
    mm_files["0.2ng_rep1"] = "data/MetaMorpheus_output/02ng/Ex_Auto_J3_30umTB_02ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep2"] = "data/MetaMorpheus_output/02ng/Ex_Auto_J3_30umTB_02ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep3"] = "data/MetaMorpheus_output/02ng/Ex_Auto_K13_30umTA_02ngQC_60m_1-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep4"] = "data/MetaMorpheus_output/02ng/Ex_Auto_K13_30umTA_02ngQC_60m_2-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep5"] = "data/MetaMorpheus_output/02ng/Ex_Auto_W17_30umTA_02ngQC_60m_3-calib_PSMs.psmtsv.gz"
    mm_files["0.2ng_rep6"] = "data/MetaMorpheus_output/02ng/Ex_Auto_W17_30umTA_02ngQC_60m_4-calib_PSMs.psmtsv.gz"
    #
    #True SC
    mm_files["sc_rep1"] = "data/MetaMorpheus_output/true_SC/D19_15um30cm_SC1-calib_PSMs.psmtsv.gz"
    mm_files["sc_rep2"] = "data/MetaMorpheus_output/true_SC/D19_15um30cm_SC2-calib_PSMs.psmtsv.gz"
    mm_files["sc_rep3"] = "data/MetaMorpheus_output/true_SC/D19_15um30cm_SC3-calib_PSMs.psmtsv.gz"
    mm_files["sc_rep4"] = "data/MetaMorpheus_output/true_SC/D19_15um30cm_SC4-calib_PSMs.psmtsv.gz"
    mm_files["sc_rep5"] = "data/MetaMorpheus_output/true_SC/D19_15um30cm_SC5-calib_PSMs.psmtsv.gz"

    mzml_files = {}
    mzml_files["bulk_rep1"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_A.mzML"
    mzml_files["bulk_rep2"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_B.mzML"
    mzml_files["bulk_rep3"] = "data/mzMLs/OR11_20160122_PG_HeLa_CVB3_CT_C.mzML"

    mzml_files["2ng_rep1"] = "data/mzMLs/Ex_Auto_J3_30umTB_2ngQC_60m_1.mzML"
    mzml_files["2ng_rep2"] = "data/mzMLs/Ex_Auto_J3_30umTB_2ngQC_60m_2.mzML"
    mzml_files["2ng_rep3"] = "data/mzMLs/Ex_Auto_K13_30umTA_2ngQC_60m_1.mzML"
    mzml_files["2ng_rep4"] = "data/mzMLs/Ex_Auto_K13_30umTA_2ngQC_60m_2.mzML"
    mzml_files["2ng_rep5"] = "data/mzMLs/Ex_Auto_W17_30umTB_2ngQC_60m_1.mzML"
    mzml_files["2ng_rep6"] = "data/mzMLs/Ex_Auto_W17_30umTB_2ngQC_60m_2.mzML"
    #
    mzml_files["0.2ng_rep1"] = "data/mzMLs/Ex_Auto_J3_30umTB_02ngQC_60m_1.mzML"
    mzml_files["0.2ng_rep2"] = "data/mzMLs/Ex_Auto_J3_30umTB_02ngQC_60m_2.mzML"
    mzml_files["0.2ng_rep3"] = "data/mzMLs/Ex_Auto_K13_30umTA_02ngQC_60m_1.mzML"
    mzml_files["0.2ng_rep4"] = "data/mzMLs/Ex_Auto_K13_30umTA_02ngQC_60m_2.mzML"
    mzml_files["0.2ng_rep5"] = "data/mzMLs/Ex_Auto_W17_30umTA_02ngQC_60m_3.mzML"
    mzml_files["0.2ng_rep6"] = "data/mzMLs/Ex_Auto_W17_30umTA_02ngQC_60m_4.mzML"

    mzml_files["sc_rep1"] = "data/mzMLs/D19_15um30cm_SC1.mzML"
    mzml_files["sc_rep2"] = "data/mzMLs/D19_15um30cm_SC2.mzML"
    mzml_files["sc_rep3"] = "data/mzMLs/D19_15um30cm_SC3.mzML"
    mzml_files["sc_rep4"] = "data/mzMLs/D19_15um30cm_SC4.mzML"
    mzml_files["sc_rep5"] = "data/mzMLs/D19_15um30cm_SC5.mzML"




    Path('data/parsed_psm/').mkdir(parents=True, exist_ok=True)#make the folder that we'll store all parsed psms in

    for file_name in list_of_file_names:
        path_to_data_loader = os.path.abspath(os.path.dirname(__file__))
        complete_path_to_psm = os.path.join(path_to_data_loader, mm_files.get(file_name)) # We then append the relative path to the data file
        complete_path_to_mzml = os.path.join(path_to_data_loader, mzml_files.get(file_name))

        psm = pd.read_csv(complete_path_to_psm, sep = '\t', dtype={'Base Sequence': str, 'Missed Cleavages': str, 'Peptide Monoisotopic Mass':str, 'Mass Diff (Da)':str,'Mass Diff (ppm)':str,'	Protein Accession':str,'Peptide Description':str, 'Notch':str, 'Num Variable Mods':str, 'Decoy': str})

        #drop columns that we don't use
        unused = ['Score', 'Delta Score', 'Notch','Base Sequence','Essential Sequence',
        'Mods','Mods Chemical Formulas', 'Mods Combined Chemical Formula',
         'Num Variable Mods', 'Missed Cleavages','Organism Name','Identified Sequence Variations',
         'Splice Sites', 'Contaminant','Peptide Description', 'Start and End Residues In Protein',
          'Previous Amino Acid', 'Next Amino Acid', 'Theoreticals Searched', 'Decoy/Contaminant/Target',
          'Normalized Spectral Angle', 'Localized Scores',
       'Improvement Possible', 'Cumulative Target', 'Cumulative Decoy',
       'Cumulative Target Notch', 'Cumulative Decoy Notch',
       'QValue Notch', 'PEP', 'PEP_QValue']
        psm = psm.drop(columns = unused)

        psm = psm.rename({"Scan Number": "scan", "Full Sequence": "peptide"}, axis=1)

        #remove duplicate scans
        psm = psm.sort_values("QValue")
        psm = psm.drop_duplicates(subset=["scan"], keep="first")

        psm[["scan"]] = psm[["scan"]].apply(pd.to_numeric)
        psm["Scan Retention Time"] = psm["Scan Retention Time"].astype(str)

        #split the time column so you get one that's just minutes
        psm["temp_minute"] = psm["Scan Retention Time"].str.split("\.")
        psm.loc[:, 'minute'] = psm['temp_minute'].map(lambda x: x[0])

        psm[["minute"]] = psm[["minute"]].apply(pd.to_numeric)
        psm = psm.drop(columns='temp_minute')

        #append info from mzML
        # import pdb;pdb.set_trace()
        mzml = pyteomics.mzml.MzML(complete_path_to_mzml)
        psm['mzml_info'] = psm.apply(lambda row: _mzml_helper(row, mzml),axis=1)
        psm[['mz_array', 'intensity_array', 'precursor_intenisty']] = pd.DataFrame(psm['mzml_info'].tolist(), index= psm.index)
        psm=psm.drop(columns='mzml_info')
        write_file_path = "data/parsed_psm/" + file_name + ".csv"
        psm.to_csv(write_file_path, sep="\t")

def load_joined_psm_mzml(file_name):
    #return the joined table (we need the start_scan time)
    mm_files = {}
    mm_files["bulk_rep1"] = "data/parsed_psm/bulk_rep1.csv"
    mm_files["bulk_rep2"] = "data/parsed_psm/bulk_rep2.csv"
    mm_files["bulk_rep3"] = "data/parsed_psm/bulk_rep3.csv"

    mm_files["2ng_rep1"] = "data/parsed_psm/2ng_rep1.csv"
    mm_files["2ng_rep2"] = "data/parsed_psm/2ng_rep2.csv"
    mm_files["2ng_rep3"] = "data/parsed_psm/2ng_rep3.csv"
    mm_files["2ng_rep4"] = "data/parsed_psm/2ng_rep4.csv"
    mm_files["2ng_rep5"] = "data/parsed_psm/2ng_rep5.csv"
    mm_files["2ng_rep6"] = "data/parsed_psm/2ng_rep6.csv"

    mm_files["0.2ng_rep1"] = "data/parsed_psm/0.2ng_rep1.csv"
    mm_files["0.2ng_rep2"] = "data/parsed_psm/0.2ng_rep2.csv"
    mm_files["0.2ng_rep3"] = "data/parsed_psm/0.2ng_rep3.csv"
    mm_files["0.2ng_rep4"] = "data/parsed_psm/0.2ng_rep4.csv"
    mm_files["0.2ng_rep5"] = "data/parsed_psm/0.2ng_rep5.csv"
    mm_files["0.2ng_rep6"] = "data/parsed_psm/0.2ng_rep6.csv"

    mm_files["sc_rep1"] = "data/parsed_psm/sc_rep1.csv"
    mm_files["sc_rep2"] = "data/parsed_psm/sc_rep2.csv"
    mm_files["sc_rep3"] = "data/parsed_psm/sc_rep3.csv"
    mm_files["sc_rep4"] = "data/parsed_psm/sc_rep4.csv"
    mm_files["sc_rep5"] = "data/parsed_psm/sc_rep5.csv"

    path_to_data_loader = os.path.abspath(os.path.dirname(__file__)) # This gets the absolute path to the location of the data_loader.py file
    complete_path_to_data = os.path.join(path_to_data_loader, mm_files.get(file_name))


    df = pd.read_csv(complete_path_to_data, sep='\t', index_col=0, dtype={'Base Sequence': str, 'Missed Cleavages': str, 'Peptide Monoisotopic Mass':str, 'Mass Diff (Da)':str,'Mass Diff (ppm)':str,'	Protein Accession':str,'Peptide Description':str, 'Notch':str, 'Num Variable Mods':str, 'Peptide Description':str, 'Decoy': str})


    return(df)

def _mzml_helper(row, mzml):
    scan = str(row['scan'])
    my_id = 'controllerType=0 controllerNumber=1 scan='+ str(scan)
    spectrum_dict = mzml.get_by_id(my_id)

    spectrum_id = spectrum_dict['id']
    retention_time = (spectrum_dict['scanList']['scan'][0].get('scan start time', -1))
    mz_array = list(spectrum_dict['m/z array'])
    intensity_array = list(spectrum_dict['intensity array'])

    #precursor information
    precursor = spectrum_dict['precursorList']['precursor'][0]
    precursor_ion = precursor['selectedIonList']['selectedIon'][0]
    precursor_mz = precursor_ion['selected ion m/z']
    if 'peak intensity' in precursor_ion:
        precursor_intenisty =  precursor_ion['peak intensity']
    else:
        precursor_intenisty = None
    if 'charge state' in precursor_ion:
        precursor_charge = int(precursor_ion['charge state'])
    elif 'possible charge state' in precursor_ion:
        precursor_charge = int(precursor_ion['possible charge state'])
    else:
        precursor_charge = 'NAN'

    all_info = [mz_array,intensity_array,precursor_intenisty]

    return all_info