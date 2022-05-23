from email import header
import unittest
from MSFragger_Parser import *

# mzML file
mzml_file_path = "C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\MetaM\\2ng\\Ex_Auto_J3_30umTB_2ngQC_60m_1.mzML"

# MSFragger files
msfragger_psm_file_path = "C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\msfragger\\psm1.tsv"
msfragger_peptide_file_path = "C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\msfragger\\combined_peptide.tsv"
msfragger_protein_file_path = "C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\msfragger\\combined_protein.tsv"

outfile = "C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\MetaM\\2ng\\MMtester.tsv"

class TestParser(unittest.TestCase):

    # test loaders

    

    def test_load_psm_df_msfragger(self):
        # POSITIVE TEST
        # note that the file types are checked before this function is called so only a positive test is needed
        psm_df = load_psm_df_msfragger(msfragger_psm_file_path)
        with open(msfragger_psm_file_path, 'r') as psm_file:
            headers = psm_file.readline()    
        headers = headers.replace("Protein ID", "Protein Accession").replace('Spectrum File', 'File Name').replace('Spectrum', 'Scan Number')   
        expected_columns = headers.strip('\n').split('\t')
        expected_columns.remove('Mapped Genes')
        expected_columns.remove('Mapped Proteins')
        obs_columns = psm_df.columns.values.tolist()
        self.assertEqual(sorted(obs_columns), sorted(expected_columns))
    
    def test_load_msfragger_peptideQ(self):
        # POSITIVE TEST
        # note that the file types are checked before this function is called so only a positive test is needed
        peptideQ_df = load_msfragger_peptideQ(msfragger_peptide_file_path)
        with open(msfragger_peptide_file_path, 'r') as peptideQ_file:
            headers = peptideQ_file.readline()
        headers = headers.replace("Protein ID", "Protein Accession").replace('Peptide Sequence', 'Peptide')
        expected_columns = headers.strip('\n').split('\t')
        obs_columns = peptideQ_df.columns.values.tolist()
        self.assertEqual(obs_columns, expected_columns)

    def test_load_msfragger_protein(self):
        # POSITIVE TEST
        # note that the file types are checked before this function is called so only a positive test is needed
        protein_df = load_msfragger_protein(msfragger_protein_file_path)
        with open(msfragger_protein_file_path, 'r') as protein_file:
            headers = protein_file.readline()
        headers = headers.replace("Protein ID", "Protein Accession").replace('Peptide Sequence', 'Peptide')
        expected_columns = headers.strip('\n').split('\t')
        obs_columns = protein_df.columns.values.tolist()
        self.assertEqual(obs_columns, expected_columns)

    # test parser helper functions

    def test_check_user_inputs(self):
        # NEGATIVE TESTS
        # test that it will user inputs that are not the right dtype
        with self.assertRaises(Exception):
            check_user_inputs(input_files=45, output_file_path=outfile, columns_to_keep=None, multiIndex=None, proteins_to_keep=None, peptides_to_keep=None, scans_to_keep=None)
        with self.assertRaises(Exception):
            check_user_inputs(input_files=msfragger_peptide_file_path, output_file_path=45, columns_to_keep=None, multiIndex=None, proteins_to_keep=None, peptides_to_keep=None, scans_to_keep=None)
        with self.assertRaises(Exception):
            check_user_inputs(input_files=msfragger_peptide_file_path, output_file_path=outfile, columns_to_keep=45, multiIndex=None, proteins_to_keep=None, peptides_to_keep=None, scans_to_keep=None)
        
        # POSITIVE TESTS
        # test that it will not raise an exception with correct file types
        try:
            check_user_inputs(input_files=msfragger_peptide_file_path, output_file_path=outfile, columns_to_keep="Peptide", multiIndex="Peptide", proteins_to_keep="Q25364", peptides_to_keep="AAAGGGCC", scans_to_keep="54263")
        except Exception:
            self.fail(f"{check_user_inputs} raised an Exception unexpectedly!")
        # test that it will convert a str into a list if it was inputted as such
        input_files, output_file_path, columns_to_keep, multiIndex, proteins_to_keep, peptides_to_keep, scans_to_keep = check_user_inputs(input_files=msfragger_peptide_file_path, output_file_path=outfile, columns_to_keep=None, multiIndex=None, proteins_to_keep=None, peptides_to_keep=None, scans_to_keep=None)
        self.assertTrue(type(input_files) == list, "check_input_files did not convert the input_files from a str into a list")

    def test_assign_file_types(self):
        # POSITIVE TESTS

        # test mzml inputs
        ex_file_list = [mzml_file_path]
        expected_output_file_list = [mzml_file_path, None, None, None]
        obs_output_file_list = assign_file_types(ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)
        # test psm inputs
        ex_file_list = [msfragger_psm_file_path]
        expected_output_file_list = [None, msfragger_psm_file_path, None, None]
        obs_output_file_list = assign_file_types(ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)
        # test peptideQ inputs
        ex_file_list = [msfragger_peptide_file_path]
        expected_output_file_list = [None, None, msfragger_peptide_file_path, None]
        obs_output_file_list = assign_file_types(ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)
        # test protein inputs
        ex_file_list = [msfragger_protein_file_path]
        expected_output_file_list = [None, None, None, msfragger_protein_file_path]
        obs_output_file_list = assign_file_types(ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)
        # test a combination of inputs
        ex_file_list = [msfragger_protein_file_path, mzml_file_path, msfragger_psm_file_path, msfragger_peptide_file_path]
        expected_output_file_list = [mzml_file_path, msfragger_psm_file_path, msfragger_peptide_file_path, msfragger_protein_file_path]
        obs_output_file_list = assign_file_types(ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)

        # NEGATIVE TESTS

        # test no valid inputs 
        ex_file_list = ['asfkjl']
        with self.assertRaises(Exception):
            assign_file_types(ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)

    def test_generate_bool_file_list(self):
        # POSITIVE TEST
        # Note: because the input file types are checked in the 'assign_file_types' function
        #       I did not include checks for file types
        ex_file_list = [None, msfragger_psm_file_path, msfragger_peptide_file_path, None]
        expected_output = '[0,1,1,0]'
        obs_output = generate_bool_file_list(ex_file_list)
        self.assertEqual(obs_output, expected_output)

    def test_load_call_dictionary(self):
        # POSITIVE TEST
        # test that the dictionary is returned correctly
        interpreted_file_list = [None, msfragger_psm_file_path, msfragger_peptide_file_path, None]
        expected_value = [psm_and_peptideQ_controller, [interpreted_file_list[1], interpreted_file_list[2], None], ['Protein Accession','Peptide', 'Scan Number']]
        call_dict = load_call_dictionary(interpreted_file_list=interpreted_file_list, columns_to_keep=None)
        obs_value = call_dict['[0,1,1,0]']
        self.assertEqual(obs_value, expected_value)

    def test_select_columns_to_keep(self):
        # POSITIVE TESTS
        # test that the correct columns are selected
        psm_df = load_psm_df_msfragger(msfragger_psm_file_path)
        columns_to_keep = ['Peptide', 'Scan Number', 'Protein Accession']
        expected_columns = str(columns_to_keep)
        mod_psm_df = select_columns_to_keep(user_dataframe=psm_df, columns_to_keep=columns_to_keep)
        obs_columns = str(mod_psm_df.columns.values.tolist())
        self.assertEqual(obs_columns, expected_columns)
        # test that if columns_to_keep = None, the original dataframe is returned
        mod_psm_df = select_columns_to_keep(user_dataframe=psm_df, columns_to_keep=None)
        expected_columns = str(psm_df.columns.values.tolist())
        obs_columns = str(mod_psm_df.columns.values.tolist())
        self.assertEqual(obs_columns, expected_columns)

        # NEGATIVE TEST
        # test that function catches incorrect columns and returns the original dataframe
        mod_psm_df = select_columns_to_keep(user_dataframe=psm_df, columns_to_keep=['asdhjk'])
        expected_columns = str(psm_df.columns.values.tolist())
        obs_columns = str(mod_psm_df.columns.values.tolist())
        self.assertEqual(obs_columns, expected_columns)

    def test_select_rows_to_keep(self):
        # NEGATIVE TESTS
        # check that the function checks that the correct ___to_keep is selected for the dataframe
        protein_df = load_msfragger_protein(msfragger_protein_file_path)
        scans_to_keep = ['6000', '5000']
        mod_protein_df = select_rows_to_keep(user_dataframe=protein_df, proteins_to_keep=None, peptides_to_keep=None, scans_to_keep=scans_to_keep)
        self.assertEqual(len(protein_df), len(mod_protein_df))
        # check that the function will return the entire dataframe if an element of a 
        # rows_to_keep parameter was not found
        psm_df = load_psm_df_msfragger(msfragger_psm_file_path)
        proteins_to_keep = ['60000000', '500000000000000']
        peptides_to_keep = ['60000000', '500000000000000']
        mod_psm_df = select_rows_to_keep(user_dataframe=psm_df, proteins_to_keep=proteins_to_keep, peptides_to_keep=peptides_to_keep, scans_to_keep=None)
        self.assertEqual(len(psm_df), len(mod_psm_df))
        # POSITIVE TESTS
        # check that the function recognizes all of the instances where the value in ___to_keep is found
        # in the given column
        peptideQ_df = load_msfragger_peptideQ(msfragger_peptide_file_path)
        proteins_to_keep = ['Q86U42', 'P37108']
        mod_peptideQ_df = select_rows_to_keep(user_dataframe=peptideQ_df, proteins_to_keep=proteins_to_keep, peptides_to_keep=None, scans_to_keep=None)
        expected_num_rows = 0
        protein_accessions = peptideQ_df['Protein Accession'].tolist()
        for protein in proteins_to_keep:
            if protein in protein_accessions:
                expected_num_rows += protein_accessions.count(protein)
        self.assertEqual(len(mod_peptideQ_df), expected_num_rows)

    def test_select_multiIndex(self):
        # POSITIVE TESTS
        # test that the function works with a correct multiIndex
        psm_df = load_psm_df_msfragger(msfragger_psm_file_path)
        multiIndex = ['Protein Accession', 'Peptide']
        mod_psm_df = select_multiIndex(user_dataframe=psm_df, multiIndex=multiIndex, default_multiIndex=['Protein Accession','Peptide', 'Scan Number'])
        expected_index = multiIndex
        obs_index = mod_psm_df.index.names
        self.assertEqual(obs_index, expected_index)
        # test that the function works with multiIndex=None and returns the dataframe with the default multiIndex for the dataframe
        psm_df = load_psm_df_msfragger(msfragger_psm_file_path)
        mod_psm_df = select_multiIndex(user_dataframe=psm_df, multiIndex=None, default_multiIndex=['Protein Accession','Peptide', 'Scan Number'])
        expected_index = ['Protein Accession','Peptide', 'Scan Number']
        obs_index = mod_psm_df.index.names
        self.assertEqual(obs_index, expected_index)

        # NEGATIVE TEST
        # test that the function catches incorrect multiIndices
        psm_df = load_psm_df_msfragger(msfragger_psm_file_path)
        multiIndex = ['Protein Accession', 'Peptide', 'asdjk']
        mod_psm_df = select_multiIndex(user_dataframe=psm_df, multiIndex=multiIndex, default_multiIndex=['Protein Accession','Peptide', 'Scan Number'])
        expected_index = ['Protein Accession','Peptide']
        obs_index = mod_psm_df.index.names
        self.assertEqual(obs_index, expected_index)

    def test_join_psm_and_peptideQ_dataframes(self):
        # POSITIVE TESTS 
        psm_df = load_psm_df_msfragger(msfragger_psm_file_path)
        peptideQ_df = load_msfragger_peptideQ(msfragger_peptide_file_path)
        joined_columns = psm_df.columns.values.tolist()
        for column in peptideQ_df.columns.values.tolist():
            if column not in joined_columns:
                joined_columns.append(column)
        joined_dataframe = join_psm_and_peptideQ_dataframes(psm_df=psm_df, peptideQ_df=peptideQ_df)
        obs_columns = joined_dataframe.columns.values.tolist()
        self.assertEqual(sorted(joined_columns), sorted(obs_columns))
        # test that the number of rows looks right
        psm_peptides = psm_df['Peptide'].values.tolist()
        peptideQ_peptides = peptideQ_df['Peptide'].values.tolist()
        num_of_overlaps = 0
        for peptide in psm_peptides:
            if peptide in peptideQ_peptides:
                num_of_overlaps += 1
        self.assertEqual(len(joined_dataframe), num_of_overlaps)
    
    def test_join_peptideQ_and_protein_dataframes(self):
        # POSITIVE TESTS 
        # test that the right number of duplicate columns are there
        peptideQ_df = load_msfragger_peptideQ(msfragger_peptide_file_path)
        protein_df = load_msfragger_protein(msfragger_protein_file_path)
        joined_dataframe = join_peptideQ_and_protein_dataframes(protein_df=protein_df, peptideQ_df=peptideQ_df)
        obs_columns = joined_dataframe.columns.tolist()
        obs_duplicates = 0
        for column in obs_columns:
            if "_peptide" in column or "_protein" in column:
                obs_duplicates += 1
        expected_num_duplicates = 0
        for column in peptideQ_df.columns.tolist():
            if column in protein_df and column != "Protein Accession" and 'Spectral Count' in column:
                expected_num_duplicates += 2
        self.assertEqual(expected_num_duplicates, obs_duplicates)
        # test that there are the right number of rows
        peptideQ_protein_accessions = peptideQ_df['Protein Accession'].values.tolist()
        protein_protein_accessions = protein_df['Protein Accession'].values.tolist()
        num_of_overlaps = 0
        for protein_accession in peptideQ_protein_accessions:
            if protein_accession in protein_protein_accessions:
                num_of_overlaps += 1
        self.assertEqual(len(joined_dataframe), num_of_overlaps)

if __name__ == '__main__':
    unittest.main()