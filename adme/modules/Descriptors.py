import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors


def generate(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_MolWt = Descriptors.MolWt(mol)
        desc_NumRotatableBonds = Descriptors.NumRotatableBonds(mol)
        row = np.array([
                        desc_MolWt,
                        desc_NumRotatableBonds])
        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1
    columnNames = ["MolWt", "NumRotatableBonds"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_1(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)
    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_NumHAcceptors = Descriptors.NumHAcceptors(mol)
        desc_NumHDonors = Descriptors.NumHDonors(mol)
        desc_NumHeteroatoms = Descriptors.NumHeteroatoms(mol)
        row = np.array([desc_NumHAcceptors,
                        desc_NumHDonors,
                        desc_NumHeteroatoms])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["NumHAcceptors", "NumHDonors", "NumHeteroatoms"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_2(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_MolMR = Descriptors.MolMR(mol)
        desc_HeavyAtomCount = Descriptors.HeavyAtomCount(mol)
        desc_HeavyAtomMolwt = Descriptors.HeavyAtomMolWt(mol)

        row = np.array([desc_MolMR,
                        desc_HeavyAtomCount,
                        desc_HeavyAtomMolwt])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["MolMR", "HeavyAtomCount", "HeavyAtomMolwt"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_3(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)
    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_NOCount = Descriptors.NOCount(mol)
        desc_NumValenceElectrons = Descriptors.NumValenceElectrons(mol)
        row = np.array([desc_NOCount,
                        desc_NumValenceElectrons])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["NOCount", "NumValenceElectrons"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_4(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_RingCount = Descriptors.RingCount(mol)
        desc_qed = Descriptors.qed(mol)
        desc_MaxAbsEStateIndex = Descriptors.MaxAbsEStateIndex(mol)
        row = np.array([desc_RingCount,
                        desc_qed,
                        desc_MaxAbsEStateIndex])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["RingCount", "qed", "MaxAbsEStateIndex"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_5(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_TPSA = Descriptors.TPSA(mol)
        desc_LabuteASA = Descriptors.LabuteASA(mol)
        desc_PEOE_VSA1 = Descriptors.PEOE_VSA1(mol)
        row = np.array([desc_TPSA,
                        desc_LabuteASA,
                        desc_PEOE_VSA1])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["TPSA", "LabuteASA", "PEOE_VSA1"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_7(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)
    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_FractionCSP3 = Descriptors.FractionCSP3(mol)
        desc_MaxPartialCharge = Descriptors.MaxPartialCharge(mol)
        desc_MinPartialCharge = Descriptors.MinPartialCharge(mol)
        row = np.array([desc_FractionCSP3,
                        desc_MaxPartialCharge,
                        desc_MinPartialCharge])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["FractionCSP3", "MaxPartialCharge", "MinPartialCharge"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_8(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)
    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_FpDensityMorgan1 = Descriptors.FpDensityMorgan1(mol)
        desc_FpDensityMorgan2 = Descriptors.FpDensityMorgan2(mol)
        desc_FpDensityMorgan3 = Descriptors.FpDensityMorgan3(mol)
        row = np.array([desc_FpDensityMorgan1,
                        desc_FpDensityMorgan2,
                        desc_FpDensityMorgan3])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["FpDensityMorgan1", "FpDensityMorgan2", "FpDensityMorgan3"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_9(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_Chi3n = Descriptors.Chi3n(mol)
        desc_Chi3v = Descriptors.Chi3v(mol)
        row = np.array([desc_Chi3n,
                        desc_Chi3v])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["Chi3n", "Chi3v"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_10(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_Kappa3 = Descriptors.Kappa3(mol)
        desc_PEOE_VSA10 = Descriptors.PEOE_VSA10(mol)
        desc_SMR_VSA4 = Descriptors.SMR_VSA4(mol)
        row = np.array([desc_Kappa3,
                        desc_PEOE_VSA10,
                        desc_SMR_VSA4])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["Kappa3", "PEOE_VSA10", "SMR_VSA4"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_11(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_PEOE_VSA12 = Descriptors.PEOE_VSA12(mol)
        desc_PEOE_VSA2 = Descriptors.PEOE_VSA2(mol)
        desc_PEOE_VSA5 = Descriptors.PEOE_VSA5(mol)
        row = np.array([desc_PEOE_VSA12,
                        desc_PEOE_VSA2,
                        desc_PEOE_VSA5])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1
    columnNames = ["PEOE_VSA12", "PEOE_VSA2", "PEOE_VSA5"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_11(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_SMR_VSA10 = Descriptors.SMR_VSA10(mol)
        desc_PEOE_VSA9 = Descriptors.PEOE_VSA9(mol)
        desc_PEOE_VSA6 = Descriptors.PEOE_VSA6(mol)
        row = np.array([desc_SMR_VSA10,
                        desc_PEOE_VSA9,
                        desc_PEOE_VSA6])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1
    columnNames = ["SMR_VSA10", "PEOE_VSA9", "PEOE_VSA6"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_14(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_SlogP_VSA4 = Descriptors.SlogP_VSA4(mol)
        desc_SlogP_VSA5 = Descriptors.SlogP_VSA5(mol)
        desc_SlogP_VSA2 = Descriptors.SlogP_VSA2(mol)
        row = np.array([desc_SlogP_VSA4,
                        desc_SlogP_VSA5,
                        desc_SlogP_VSA2])
        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["SlogP_VSA4", "SlogP_VSA5", "SlogP_VSA2"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_15(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_SlogP_VSA7 = Descriptors.SlogP_VSA7(mol)
        desc_Estate_VSA1 = Descriptors.EState_VSA1(mol)
        row = np.array([desc_SlogP_VSA7,
                        desc_Estate_VSA1])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["SlogP_VSA7", "Estate_VSA1"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_16(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_Estate_VSA10 = Descriptors.EState_VSA10(mol)
        desc_Estate_VSA2 = Descriptors.EState_VSA2(mol)
        desc_Estate_VSA3 = Descriptors.EState_VSA3(mol)
        row = np.array([desc_Estate_VSA10,
                        desc_Estate_VSA2,
                        desc_Estate_VSA3])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["Estate_VSA10", "Estate_VSA2", "Estate_VSA3"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_17(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_Estate_VSA4 = Descriptors.EState_VSA4(mol)
        desc_Estate_VSA7 = Descriptors.EState_VSA7(mol)
        desc_Estate_VSA8 = Descriptors.EState_VSA8(mol)
        row = np.array([desc_Estate_VSA4,
                        desc_Estate_VSA7,
                        desc_Estate_VSA8])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["Estate_VSA4", "Estate_VSA7", "Estate_VSA8"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def generate_18(smiles, verbose=False):
    moldata = []
    for elem in smiles:
        mol = Chem.MolFromSmiles(elem)
        moldata.append(mol)

    baseData = np.arange(1, 1)
    i = 0
    for mol in moldata:
        desc_EState_VSA9 = Descriptors.EState_VSA9(mol)
        desc_VSA_Estate1 = Descriptors.VSA_EState1(mol)
        row = np.array([desc_EState_VSA9,
                        desc_VSA_Estate1])

        if(i == 0):
            baseData = row
        else:
            baseData = np.vstack([baseData, row])
        i = i+1

    columnNames = ["EState_VSA9", "VSA_Estate1"]
    descriptors = pd.DataFrame(data=baseData, columns=columnNames)

    return descriptors


def AromaticAtoms(m):
    aromatic_atoms = [m.GetAtomWithIdx(i).GetIsAromatic()
                      for i in range(m.GetNumAtoms())]
    aa_count = []
    for i in aromatic_atoms:
        if i is True:
            aa_count.append(1)

    sum_aa_count = sum(aa_count)

    return sum_aa_count
