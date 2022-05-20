import unittest
from MetaMorpheus_Parser import *

# mzML file
mzml_file_path = "C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\MetaM\\2ng\\Ex_Auto_J3_30umTB_2ngQC_60m_1.mzML"

# MetaMorpheus files
mm_psm_file_path = "C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\MetaM\\2ng\\Ex_Auto_J3_30umTB_2ngQC_60m_1-calib_PSMs.psmtsv"
mm_peptideQ_file_path = "C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\MetaM\\2ng\\AllQuantifiedPeptides.tsv"
mm_protein_file_path = "C:\\Users\\Sarah Curtis\\OneDrive - BYU\\Documents\\Single Cell Team Documents\\API_dev\\MetaM\\2ng\\AllQuantifiedProteinGroups.tsv"

class TestParser(unittest.TestCase):

    # test loaders

    def test_load_mzml_df(self):
        # add this one last
        pass

    def test_load_psm(self):
        # POSITIVE TEST
        # note that the file types are checked before this function is called so only a positive test is needed
        psm_df = load_psm(mm_psm_file_path)
        with open(mm_psm_file_path, 'r') as psm_file:
            headers = psm_file.readline()    
            
        expected_columns = headers.replace("Full Sequence", "Peptide").strip('\n').split('\t')
        obs_columns = psm_df.columns.values.tolist()
        self.assertEqual(obs_columns, expected_columns)
    
    def test_load_peptideQ(self):
        # POSITIVE TEST
        # note that the file types are checked before this function is called so only a positive test is needed
        peptideQ_df = load_peptideQ(mm_peptideQ_file_path)
        with open(mm_peptideQ_file_path, 'r') as peptideQ_file:
            headers = peptideQ_file.readline()
        expected_columns = headers.replace("Protein Groups", "Protein Accession").strip('\n').split('\t')
        expected_columns[expected_columns.index('Sequence')] = 'Peptide'
        obs_columns = peptideQ_df.columns.values.tolist()
        self.assertEqual(obs_columns, expected_columns)

    def test_load_protein(self):
        # POSITIVE TEST
        # note that the file types are checked before this function is called so only a positive test is needed
        protein_df = load_protein(mm_protein_file_path)
        with open(mm_protein_file_path, 'r') as protein_file:
            headers = protein_file.readline()
        expected_columns = headers.strip('\n').split('\t')
        obs_columns = protein_df.columns.values.tolist()
        self.assertEqual(obs_columns, expected_columns)

    # test parser helper functions

    def test_assign_file_types(self):
        # POSITIVE TESTS

        # test mzml inputs
        ex_file_list = [mzml_file_path]
        expected_output_file_list = [mzml_file_path, None, None, None]
        obs_output_file_list = assign_file_types(ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)
        # test psm inputs
        ex_file_list = [mm_psm_file_path]
        expected_output_file_list = [None, mm_psm_file_path, None, None]
        obs_output_file_list = assign_file_types(ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)
        # test peptideQ inputs
        ex_file_list = [mm_peptideQ_file_path]
        expected_output_file_list = [None, None, mm_peptideQ_file_path, None]
        obs_output_file_list = assign_file_types(ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)
        # test protein inputs
        ex_file_list = [mm_protein_file_path]
        expected_output_file_list = [None, None, None, mm_protein_file_path]
        obs_output_file_list = assign_file_types(ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)
        # test a combination of inputs
        ex_file_list = [mm_protein_file_path, mzml_file_path, mm_psm_file_path, mm_peptideQ_file_path]
        expected_output_file_list = [mzml_file_path, mm_psm_file_path, mm_peptideQ_file_path, mm_protein_file_path]
        obs_output_file_list = assign_file_types(ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)

        # NEGATIVE TESTS

        # test no valid inputs 
        ex_file_list = ['asfkjl']
        expected_output_file_list = ['asfkjl', None, None, None]
        obs_output_file_list = assign_file_types(ex_file_list)
        print(obs_output_file_list)
        self.assertWarns(Warning, assign_file_types, ex_file_list)
        self.assertEqual(obs_output_file_list, expected_output_file_list)

    def test_generate_bool_file_list(self):
        # POSITIVE TEST
        # Note: because the input file types are checked in the 'assign_file_types' function
        #       I did not include checks for file types
        ex_file_list = [None, mm_psm_file_path, mm_peptideQ_file_path, None]
        expected_output = '[0,1,1,0]'
        obs_output = generate_bool_file_list(ex_file_list)
        self.assertEqual(obs_output, expected_output)

    def test_load_call_dictionary(self):
        # POSITIVE TEST
        # test that the dictionary is returned correctly
        interpreted_file_list = [None, mm_psm_file_path, mm_peptideQ_file_path, None]
        expected_value = [psm_and_peptideQ_controller, [interpreted_file_list[1], interpreted_file_list[2], None], ['Protein Accession','Peptide', 'Scan Number']]
        call_dict = load_call_dictionary(interpreted_file_list=interpreted_file_list, columns_to_keep=None)
        obs_value = call_dict['[0,1,1,0]']
        self.assertEqual(obs_value, expected_value)

    def test_select_columns_to_keep(self):
        # POSITIVE TESTS
        # test that the correct columns are selected
        psm_df = load_psm(mm_psm_file_path)
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

    def test_select_multiIndex(self):
        # POSITIVE TESTS
        # test that the function works with a correct multiIndex
        psm_df = load_psm(mm_psm_file_path)
        multiIndex = ['Protein Accession', 'Peptide']
        mod_psm_df = select_multiIndex(user_dataframe=psm_df, multiIndex=multiIndex, default_multiIndex=['Protein Accession','Peptide', 'Scan Number'])
        expected_index = multiIndex
        obs_index = mod_psm_df.index.names
        self.assertEqual(obs_index, expected_index)
        # test that the function works with multiIndex=None and returns the dataframe with the default multiIndex for the dataframe
        psm_df = load_psm(mm_psm_file_path)
        mod_psm_df = select_multiIndex(user_dataframe=psm_df, multiIndex=None, default_multiIndex=['Protein Accession','Peptide', 'Scan Number'])
        expected_index = ['Protein Accession','Peptide', 'Scan Number']
        obs_index = mod_psm_df.index.names
        self.assertEqual(obs_index, expected_index)

        # NEGATIVE TEST
        # test that the function catches incorrect multiIndices
        psm_df = load_psm(mm_psm_file_path)
        multiIndex = ['Protein Accession', 'Peptide', 'asdjk']
        mod_psm_df = select_multiIndex(user_dataframe=psm_df, multiIndex=multiIndex, default_multiIndex=['Protein Accession','Peptide', 'Scan Number'])
        expected_index = ['Protein Accession','Peptide']
        obs_index = mod_psm_df.index.names
        self.assertEqual(obs_index, expected_index)

    # test joining functions
    def test_join_psm_and_peptideQ_dataframes(self):
        # positive test 
        psm_df = load_psm(mm_psm_file_path)
        peptideQ_df = load_peptideQ(mm_peptideQ_file_path)
        joined_columns = psm_df.columns.values.tolist()
        for column in peptideQ_df.columns.values.tolist():
            if column not in joined_columns:
                joined_columns.append(column)
        joined_dataframe = join_psm_and_peptideQ_dataframes(psm_df=psm_df, peptideQ_df=peptideQ_df)
        obs_columns = joined_dataframe.columns.values.tolist()
        self.assertEqual(sorted(joined_columns), sorted(obs_columns))
    # should I test the other ones? 
    # I'm guessing that the mzml one would be kinda because I'm not sure how to create the positive control





if __name__ == '__main__':
    unittest.main()