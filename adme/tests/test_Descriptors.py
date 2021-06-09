from adme.modules import Descriptors
import numpy as np
import pandas as pd


def test_generate():
    """
    Unit test for generating molecular descriptors.
    """
    df = pd.read_csv('tests/df_Cmax_test.csv', sep=';')
    input = df['canonical_smiles']
    columnNames = ["MolWt", "NumRotatableBonds"]
    baseData = np.arange(1, 1)
    output = Descriptors.generate(input)
    result = pd.DataFrame(data=baseData, columns=columnNames)
    expected = pd.DataFrame
    assert result.equals(expected), "The function is broken!"


def test_generate_1():
    df = pd.read_csv('tests/df_Cmax_test.csv', sep=';')
    input = df['canonical_smiles']
    columnNames = ["MolMR", "HeavyAtomCount", "HeavyAtomMolwt"]
    baseData = np.arange(1, 1)
    output = Descriptors.generate(input)
    result = pd.DataFrame(data=baseData, columns=columnNames)
    expected = pd.DataFrame
    assert result.equals(expected), "The function is broken!"
