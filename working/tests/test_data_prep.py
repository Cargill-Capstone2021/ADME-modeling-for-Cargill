import sys
sys.path.append('../')
from data_prep import pre_clean_data
from data_prep import rearrange_and_clean_data
from data_prep import combine_rows
from data_prep import remove_nan
import pandas as pd
import numpy as np
import math


def test_pre_clean_data():
    """
    Unit test for pre_clean_data(df) function to see if the function
    produces correct output.
    """
    input_df = pd.read_csv('input_df_1.csv', sep=';')
    output_df = pd.read_csv('output_df_1.csv', sep=';')
    df = pre_clean_data(input_df)
    assert output_df.equals(df),\
        "The function pre_clean_data is broken!"


def test_rearrange_and_clean_data():
    """
    Unit test for rearrange_and_clean_data(df) function to see if the function
    produces correct output.
    """
    input_df = pd.read_csv('output_df_1.csv', sep=';')
    output_df = pd.read_csv('output_df_2.csv', sep=';')
    df = rearrange_and_clean_data(input_df)
    assert output_df.equals(df),\
        "The function rearrange_and_clean_data is broken!"


def test_combine_rows():
    """
    Unit test for combine_rows(df) function to see if the function
    produces correct output.
    """
    input_df = pd.read_csv('output_df_2.csv', sep=';')
    output_df = pd.read_csv('output_df_3.csv', sep=';')
    df = combine_rows(input_df)
    assert output_df.equals(df),\
        "The function combine_rows is broken!"


def test_remove_nan():
    """
    Unit test for remove_nan(df_final) function to see if the function
    produces correct output.
    """
    input_df = pd.read_csv('output_df_3.csv', sep=';')
    df_t_half_output = pd.read_csv('df_T_half_test.csv', sep=';')
    df_f_output = pd.read_csv('df_F_test.csv', sep=';')
    df_auc_output = pd.read_csv('df_AUC_test.csv', sep=';')
    df_cmax_output = pd.read_csv('df_Cmax_test.csv', sep=';')
    df_vdss_output = pd.read_csv('df_Vdss_test.csv', sep=';')
    df_cl_output = pd.read_csv('df_CL_test.csv', sep=';')
    df_ld50_output = pd.read_csv('df_LD50_test.csv', sep=';')
    df_ic50_output = pd.read_csv('df_IC50_test.csv', sep=';')
    df_ppb_output = pd.read_csv('df_PPB_test.csv', sep=';')

    df_t_half, df_f, df_auc, df_cmax, df_vdss, df_cl, df_ld50,\
        df_ic50, df_ppb = remove_nan(input_df)
    assert df_t_half.equals(df_t_half_output) and df_f.equals(df_f_output)\
        and df_auc.equals(df_auc_output) and df_cmax.equals(df_cmax_output)\
        and df_vdss.equals(df_vdss_output) and df_cl.equals(df_cl_output)\
        and df_ld50.equals(df_ld50_output) and df_ic50.equals(df_ic50_output)\
        and df_ppb.equals(df_ppb_output), "The function remove_NaN is broken!"
