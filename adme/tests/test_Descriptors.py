from adme import Descriptors
import pandas as pd
import unittest

def test_generate(self):
  """
  Unit test for generating molecular descriptors.
  """
  df = pd.read('Cmax_final.csv')
  input = df.canonical_smiles
  columnNames=["MolWt","NumRotatableBonds"] 
  baseData= np.arange(1,1)
  output = generate(input)
  result = pd.DataFrame(data=baseData,columns=columnNames)
  expected = pd.DataFrame
  self.assertTrue(result, excpected)

  return
   

def test_generate_1(self):

  df = pd.read('Cmax_final.csv')
  input = df.canonical_smiles
  columnNames=["MolMR","HeavyAtomCount","HeavyAtomMolwt"]   
  baseData= np.arange(1,1)
  output = generate(input)
  result = pd.DataFrame(data=baseData,columns=columnNames)
  expected = pd.DataFrame
  self.assertTrue(result, excpected)

  return
 
