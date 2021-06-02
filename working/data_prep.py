import pandas as pd
import numpy as np
import math


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
        'T1/2 /hr', 'F /%', 'AUC /ng.hr.mL-1', 'IC50 /nM', 'Cmax /ug.mL-1', 
        'Vdss /L.kg-1', 'CL /mL.min-1.kg-1', 'LD50 /mg.kg-1', 'PPB /%', 
        'Papp /ucm/s']
    
    new_columns = ['chembl_id', 'organism', 'canonical_smiles', 'T1/2 /hr', 
        'F /%', 'AUC /ng.hr.mL-1', 'IC50 /nM', 'Cmax /ug.mL-1', 
        'Vdss /L.kg-1', 'CL /mL.min-1.kg-1', 'LD50 /mg.kg-1', 'PPB /%', 
        'Papp /ucm/s']

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
    df = df[new_columns]
    
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



def remove_NaN(df_final):
    """
    Remove all the NaN values for individual feature.
    
    Parameters
    ----------
    df: DataFrame
        DataFrame obtained from combine_rows.
        
    Returns
    -------
    DataFrame
        The final DataFrames ready to be trained.
    """
    columns_T_half = ['chembl_id', 'organism', 'canonical_smiles', 'T1/2 /hr']
    df_T_half_pre = df_final[columns_T_half]

    columns_F = ['chembl_id', 'organism', 'canonical_smiles', 'F /%']
    df_F_pre = df_final[columns_F] 

    columns_AUC = ['chembl_id', 'organism', 'canonical_smiles', 'AUC /ng.hr.mL-1']
    df_AUC_pre = df_final[columns_AUC] 

    columns_Cmax = ['chembl_id', 'organism', 'canonical_smiles', 'Cmax /ug.mL-1']
    df_Cmax_pre = df_final[columns_Cmax] 

    columns_Vdss = ['chembl_id', 'organism', 'canonical_smiles', 'Vdss /L.kg-1']
    df_Vdss_pre = df_final[columns_Vdss] 

    columns_CL = ['chembl_id', 'organism', 'canonical_smiles', 'CL /mL.min-1.kg-1']
    df_CL_pre = df_final[columns_CL]  

    columns_LD50 = ['chembl_id', 'organism', 'canonical_smiles', 'LD50 /mg.kg-1']
    df_LD50_pre = df_final[columns_LD50]

    columns_IC50 = ['chembl_id', 'organism', 'canonical_smiles', 'IC50 /nM']
    df_IC50_pre = df_final[columns_IC50]

    columns_PPB = ['chembl_id', 'organism', 'canonical_smiles', 'PPB /%']
    df_PPB_pre = df_final[columns_PPB]
    
    rows_to_drop_T_half = []
    for index, row in df_T_half_pre.iterrows():
        if math.isnan(row['T1/2 /hr']) == True:
            rows_to_drop_T_half.append(index)

    df_T_half_final = df_T_half_pre.drop(rows_to_drop_T_half)
    df_T_half_final = df_T_half_final.reset_index(drop=True)


    rows_to_drop_F = []
    for index, row in df_F_pre.iterrows():
        if math.isnan(row['F /%']) == True:
            rows_to_drop_F.append(index)

    df_F_final = df_F_pre.drop(rows_to_drop_F)
    df_F_final = df_F_final.reset_index(drop=True)


    rows_to_drop_AUC = []
    for index, row in df_AUC_pre.iterrows():
        if math.isnan(row['AUC /ng.hr.mL-1']) == True:
            rows_to_drop_AUC.append(index)

    df_AUC_final = df_AUC_pre.drop(rows_to_drop_AUC)
    df_AUC_final = df_AUC_final.reset_index(drop=True)


    rows_to_drop_Cmax = []
    for index, row in df_Cmax_pre.iterrows():
        if math.isnan(row['Cmax /ug.mL-1']) == True:
            rows_to_drop_Cmax.append(index)

    df_Cmax_final = df_Cmax_pre.drop(rows_to_drop_Cmax)
    df_Cmax_final = df_Cmax_final.reset_index(drop=True)


    rows_to_drop_Vdss = []
    for index, row in df_Vdss_pre.iterrows():
        if math.isnan(row['Vdss /L.kg-1']) == True:
            rows_to_drop_Vdss.append(index)

    df_Vdss_final = df_Vdss_pre.drop(rows_to_drop_Vdss)
    df_Vdss_final = df_Vdss_final.reset_index(drop=True)


    rows_to_drop_CL = []
    for index, row in df_CL_pre.iterrows():
        if math.isnan(row['CL /mL.min-1.kg-1']) == True:
            rows_to_drop_CL.append(index)

    df_CL_final = df_CL_pre.drop(rows_to_drop_CL)
    df_CL_final = df_CL_final.reset_index(drop=True)


    rows_to_drop_LD50 = []
    for index, row in df_LD50_pre.iterrows():
        if math.isnan(row['LD50 /mg.kg-1']) == True:
            rows_to_drop_LD50.append(index)

    df_LD50_final = df_LD50_pre.drop(rows_to_drop_LD50)
    df_LD50_final = df_LD50_final.reset_index(drop=True)


    rows_to_drop_IC50 = []
    for index, row in df_IC50_pre.iterrows():
        if math.isnan(row['IC50 /nM']) == True:
            rows_to_drop_IC50.append(index)

    df_IC50_final = df_IC50_pre.drop(rows_to_drop_IC50)
    df_IC50_final = df_IC50_final.reset_index(drop=True)


    rows_to_drop_PPB = []
    for index, row in df_PPB_pre.iterrows():
        if math.isnan(row['PPB /%']) == True:
            rows_to_drop_PPB.append(index)

    df_PPB_final = df_PPB_pre.drop(rows_to_drop_PPB)
    df_PPB_final = df_PPB_final.reset_index(drop=True)
    
    return df_T_half_final, df_F_final, df_AUC_final, df_Cmax_final,\
           df_Vdss_final, df_CL_final, df_LD50_final, df_IC50_final,\
           df_PPB_final
