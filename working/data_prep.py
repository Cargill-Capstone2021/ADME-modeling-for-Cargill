import pandas as pd
import numpy as np


def pre_clean_data(df):
    """
    Pre-clean the DataFrame directly extracted from ChemBL data base by 
    dropping unnecessary columns and rows.

    Parameters
    ----------
    df: DataFrame
        DataFrame directly extracted from ChemBL data base with sql.
    
    Returns
    -------
    DataFrame
        The cleaned DataFrame.
    """
    columns_to_drop = [
        'alogp', 'hba', 'hbd', 'psa', 'rtb', 'ro3_pass', 
        'num_ro5_violations', 'cx_most_apka', 'cx_most_bpka', 'full_mwt', 
        'cx_logp', 'cx_logd', 'molecular_species', 'mw_freebase', 
        'aromatic_rings', 'heavy_atoms', 'qed_weighted', 'mw_monoisotopic', 
        'full_molformula', 'hba_lipinski', 'hbd_lipinski', 
        'num_lipinski_ro5_violations', 'assay_type', 'target_type', 'tax_id', 
        'molregno', 'molregno.1', 'structure_type', 'molecule_type', 'chirality', 
        'natural_product', 'relation', 'standard_inchi', 
        'standard_relation', 'type', 'value', 'units']
    
    df = df.drop(columns=columns_to_drop)
    
    types_to_drop = [
        'Tissue Severity Score', 'WEIGHT', 'BUN', 'ALT', 'ALP', 'CHOL', 
        'GLUC', 'CREAT', 'PHOS', 'ALB', 'AST', 'RBC', 'HCT', 'MCV', 
        'PROT', 'WBC', 'HGB', 'MCH', 'MCHC', 'SODIUM', 'POTASSIUM', 
        'PLAT', 'NEUTLE', 'EOSLE', 'BASOLE', 'MONOLE', 'LYMLE', 'BILI',
        'Stabilty','Stability', 'Cl', 'RETIRBC', 'TRIG', 'PHOSLPD', 
        'CALCIUM', 'ALBGLOB', 'GGT', 'TERMBW', 'BILDIR', 'PT', 'APTT', 'FIBRINO', 
        'CHLORIDE', 'CK', 'LIPASE', 'URATE', 'Fu', 'CO2', 'BASO', 'EOS' 'LYM', 
        'MONO', 'NEUTSG', 'RBCNUC', 'TD50', 'ID/g', 'Vd', 'CLH', 'CC50', 'TIME',
        'MTD', 'UI', 'Dose/organ', 'Concentration', 'Survival', 
        'Radiolabel recovery', 'GI50', 'Dose/g', 'Radioactivity', 'Activity', 
        'Cp', 'LDH', 'Tmax', 'EOS', 'LYM', 'Biodistribution', 'Injected dose/g',
        'Cp(f)', 'Ratio TD50/ED50', 'ERH', 'Distribution of radioactivity', 
        'Plasma level']

    for i in types_to_drop:
        index_names = df[df['standard_type'] == i].index
        df.drop(index_names, inplace = True)

    df = df.reset_index(drop=True)
    
    return df


def rearrange_and_clean_data(df):
    """
    Rearrange the DataFrame returned by the pre_clean_data(df) function. 
    Add new empty columns with 'standard_types/standard_units' as the column 
    headers and then fill the table cells with standard values. 
    Drop unnecessary rows and columns.
    
    Parameters
    ----------
    df: DataFrame
        DataFrame returned by the pre_clean_data(df) function.
    
    Returns
    -------
    DataFrame
        The rearranged and cleaned DataFrame.
    """
    # The features we want to predict.
    header_list = [
        'T1/2 /hr', 'F /%', 'AUC /ng.hr.mL-1', 'Drug uptake /% ID/g', 'IC50 /nM', 
        'Cmax /ug.mL-1', 'Vdss /L.kg-1', 'CL /mL.min-1.kg-1', 'LD50 /mg.kg-1', 
        'PPB /%', 'Drug metabolism /nan', 'Distribution /%', 'ED50 /mg.kg-1', 
        'Permeability /nan', 'Ratio /nan', 'Papp /ucm/s']
    
    # Add empty columns with the features as the column headers.
    for i in header_list:
        df[i] = np.nan

    rows_to_drop = []
    
    # Fill the values in the table cells and drop unnecessary rows.
    for index, row in df.iterrows():
        types = str(row['standard_type']) + ' /' + str(row['standard_units'])
        if types in header_list:
            df.loc[index, types] = row['standard_value']
        else:
            rows_to_drop.append(index)

    df.drop(rows_to_drop, inplace = True)
   
    df = df.drop(columns=['standard_type', 'standard_value', 'standard_units'])
    df = df.reset_index(drop=True)    
    
    return df


def combine_rows(df):
    """
    Combine the rows that represent the same molecule (same 'chembl_id').
    
    Parameters
    ----------
    df: DataFrame
        DataFrame returned by rearrange_and_clean_data(df) function.

    Returns
    -------
    DataFrame
        The final DataFrame ready to be trained.
    """
    # Rearrange the rows based on their chembl_ids.
    df = df.sort_values(by=['chembl_id'], ascending=False).reset_index(drop=True)
   
    df_final = pd.DataFrame(columns=df.columns)

    # Combine the rows that represent the same molecule and put the combined row 
    # into a new DataFrame df_final.
    for i, row in df.iterrows():
        if i==0:
            continue
        if row['chembl_id'] == df.loc[i-1, 'chembl_id']:
            df.iloc[i] = df.iloc[i-1].combine_first(df.iloc[i])
        else:
            df_final = df_final.append(df.iloc[i-1], ignore_index=True)

    df_final = df_final.append(df.iloc[len(df)-1], ignore_index=True)
    
    return df_final


def data_prep(df):
    """
    Clean and reorganize the data extracted directly from ChemBL data base for 
    the purpose of training the machine learning model.

    Parameters
    ----------
    df: DataFrame
        DataFrame extracted directly from ChemBL data base.
    Returns
    -------
    DataFrame
        The final DataFrame ready to be trained.
    """
    df1 = pre_clean_data(df)
    df2 = rearrange_and_clean_data(df1)
    df_final = combine_rows(df2)
    
    return df_final
