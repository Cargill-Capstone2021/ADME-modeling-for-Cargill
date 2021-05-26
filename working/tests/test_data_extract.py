
import sqlite3
import pandas as pd
import unittest

import sys
sys.path.append('../')
import data_extract


class Testdataextract(unittest.TestCase):
    def test_dbfile_to_dataframe(self):
        result = list(data_extract.dbfile_to_dataframe('source.db').columns)
        self.assertEqual(result, ['organism', 'standard_type', 'standard_value', 'standard_units', 'chembl_id', 'canonical_smiles']), 'unexpected result'
        return

    def test_merged_data(self):
        data1 = pd.read_csv("Oryctolagus_cuniculus.csv", sep=';')
        data2 = pd.read_csv("Ovis_Aries.csv", sep=';')
        df = data_extract.merged_data(data1, data2)
        result = list(df['organism'].unique())
        self.assertEqual(result,['Oryctolagus cuniculus', 'Ovis aries']), 'unexpected result: the dataframes are not combined properly'
        return

