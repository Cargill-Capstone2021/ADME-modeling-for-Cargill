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
        'molregno', 'molregno.1', 'structure_type', 'molecule_type',
        'chirality', 'natural_product', 'relation', 'standard_inchi',
        'standard_relation', 'type', 'value', 'units']

    df = df.drop(columns=columns_to_drop)

    types_to_drop = [
        'Tissue Severity Score', 'WEIGHT', 'BUN', 'ALT', 'ALP', 'CHOL',
        'GLUC', 'CREAT', 'PHOS', 'ALB', 'AST', 'RBC', 'HCT', 'MCV',
        'PROT', 'WBC', 'HGB', 'MCH', 'MCHC', 'SODIUM', 'POTASSIUM',
        'PLAT', 'NEUTLE', 'EOSLE', 'BASOLE', 'MONOLE', 'LYMLE', 'BILI',
        'Stabilty', 'Stability', 'Cl', 'RETIRBC', 'TRIG', 'PHOSLPD',
        'CALCIUM', 'ALBGLOB', 'GGT', 'TERMBW', 'BILDIR', 'PT', 'APTT',
        'FIBRINO', 'CHLORIDE', 'CK', 'LIPASE', 'URATE', 'Fu', 'CO2',
        'BASO', 'EOS' 'LYM', 'MONO', 'NEUTSG', 'RBCNUC', 'TD50', 'ID/g',
        'Vd', 'CLH', 'CC50', 'TIME', 'MTD', 'UI', 'Dose/organ',
        'Concentration', 'Survival', 'Radiolabel recovery', 'GI50',
        'Dose/g', 'Radioactivity', 'Activity', 'Cp', 'LDH', 'Tmax',
        'EOS', 'LYM', 'Biodistribution', 'Injected dose/g', 'Cp(f)',
        'Ratio TD50/ED50', 'ERH', 'Distribution of radioactivity',
        'Plasma level']

    for i in types_to_drop:
        index_names = df[df['standard_type'] == i].index
        df.drop(index_names, inplace=True)

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
                   'Vdss /L.kg-1', 'CL /mL.min-1.kg-1', 'LD50 /mg.kg-1',
                   'PPB /%', 'Papp /ucm/s']

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

    df.drop(rows_to_drop, inplace=True)

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
    df = df.sort_values(by=['chembl_id'], ascending=False).\
        reset_index(drop=True)

    df_final = pd.DataFrame(columns=df.columns)

    # Combine the rows that represent the same molecule and put the
    # combined row into a new DataFrame df_final.
    for i, row in df.iterrows():
        if i == 0:
            continue
        if row['chembl_id'] == df.loc[i-1, 'chembl_id']:
            df.iloc[i] = df.iloc[i-1].combine_first(df.iloc[i])
        else:
            df_final = df_final.append(df.iloc[i-1], ignore_index=True)

    df_final = df_final.append(df.iloc[len(df)-1], ignore_index=True)

    return df_final


def remove_nan(df_final):
    """
    Remove all the NaN values for individual feature.

    Parameters
    ----------
    df_final: DataFrame
        DataFrame obtained from combine_rows.

    Returns
    -------
    DataFrame
        The final DataFrames ready to be trained.
    """
    columns_t_half = ['chembl_id', 'organism', 'canonical_smiles', 'T1/2 /hr']
    df_t_half_pre = df_final[columns_t_half]

    columns_f = ['chembl_id', 'organism', 'canonical_smiles', 'F /%']
    df_f_pre = df_final[columns_f]

    columns_auc = ['chembl_id', 'organism', 'canonical_smiles',
                   'AUC /ng.hr.mL-1']
    df_auc_pre = df_final[columns_auc]

    columns_cmax = ['chembl_id', 'organism', 'canonical_smiles',
                    'Cmax /ug.mL-1']
    df_cmax_pre = df_final[columns_cmax]

    columns_vdss = ['chembl_id', 'organism', 'canonical_smiles',
                    'Vdss /L.kg-1']
    df_vdss_pre = df_final[columns_vdss]

    columns_cl = ['chembl_id', 'organism', 'canonical_smiles',
                  'CL /mL.min-1.kg-1']
    df_cl_pre = df_final[columns_cl]

    columns_ld50 = ['chembl_id', 'organism', 'canonical_smiles',
                    'LD50 /mg.kg-1']
    df_ld50_pre = df_final[columns_ld50]

    columns_ic50 = ['chembl_id', 'organism', 'canonical_smiles',
                    'IC50 /nM']
    df_ic50_pre = df_final[columns_ic50]

    columns_ppb = ['chembl_id', 'organism', 'canonical_smiles',
                   'PPB /%']
    df_ppb_pre = df_final[columns_ppb]

    rows_to_drop_t_half = []
    for index, row in df_t_half_pre.iterrows():
        if math.isnan(row['T1/2 /hr']):
            rows_to_drop_t_half.append(index)

    df_t_half_final = df_t_half_pre.drop(rows_to_drop_t_half)
    df_t_half_final = df_t_half_final.reset_index(drop=True)

    rows_to_drop_f = []
    for index, row in df_f_pre.iterrows():
        if math.isnan(row['F /%']):
            rows_to_drop_f.append(index)

    df_f_final = df_f_pre.drop(rows_to_drop_f)
    df_f_final = df_f_final.reset_index(drop=True)

    rows_to_drop_auc = []
    for index, row in df_auc_pre.iterrows():
        if math.isnan(row['AUC /ng.hr.mL-1']):
            rows_to_drop_auc.append(index)

    df_auc_final = df_auc_pre.drop(rows_to_drop_auc)
    df_auc_final = df_auc_final.reset_index(drop=True)

    rows_to_drop_cmax = []
    for index, row in df_cmax_pre.iterrows():
        if math.isnan(row['Cmax /ug.mL-1']):
            rows_to_drop_cmax.append(index)

    df_cmax_final = df_cmax_pre.drop(rows_to_drop_cmax)
    df_cmax_final = df_cmax_final.reset_index(drop=True)

    rows_to_drop_vdss = []
    for index, row in df_vdss_pre.iterrows():
        if math.isnan(row['Vdss /L.kg-1']):
            rows_to_drop_vdss.append(index)

    df_vdss_final = df_vdss_pre.drop(rows_to_drop_vdss)
    df_vdss_final = df_vdss_final.reset_index(drop=True)

    rows_to_drop_cl = []
    for index, row in df_cl_pre.iterrows():
        if math.isnan(row['CL /mL.min-1.kg-1']):
            rows_to_drop_cl.append(index)

    df_cl_final = df_cl_pre.drop(rows_to_drop_cl)
    df_cl_final = df_cl_final.reset_index(drop=True)

    rows_to_drop_ld50 = []
    for index, row in df_ld50_pre.iterrows():
        if math.isnan(row['LD50 /mg.kg-1']):
            rows_to_drop_ld50.append(index)

    df_ld50_final = df_ld50_pre.drop(rows_to_drop_ld50)
    df_ld50_final = df_ld50_final.reset_index(drop=True)

    rows_to_drop_ic50 = []
    for index, row in df_ic50_pre.iterrows():
        if math.isnan(row['IC50 /nM']):
            rows_to_drop_ic50.append(index)

    df_ic50_final = df_ic50_pre.drop(rows_to_drop_ic50)
    df_ic50_final = df_ic50_final.reset_index(drop=True)

    rows_to_drop_ppb = []
    for index, row in df_ppb_pre.iterrows():
        if math.isnan(row['PPB /%']):
            rows_to_drop_ppb.append(index)

    df_ppb_final = df_ppb_pre.drop(rows_to_drop_ppb)
    df_ppb_final = df_ppb_final.reset_index(drop=True)

    return df_t_half_final, df_f_final, df_auc_final, df_cmax_final,\
        df_vdss_final, df_cl_final, df_ld50_final, df_ic50_final,\
        df_ppb_final
