import pandas as pd
import numpy as np
import math

import sys
sys.path.append('../')
from data_prep import pre_clean_data
from data_prep import rearrange_and_clean_data
from data_prep import combine_rows
from data_prep import remove_NaN


def test_pre_clean_data():
    """
    Unit test for pre_clean_data(df) function to see if the function 
    produces correct output.
    """
    input_df = pd.read_csv('input_df_1.csv',sep=';') 
    output_df = pd.read_csv('output_df_1.csv',sep=';') 
    df = pre_clean_data(input_df)
    assert output_df.equals(df),\
        "The function pre_clean_data is broken!"  


def test_rearrange_and_clean_data():
    """
    Unit test for rearrange_and_clean_data(df) function to see if the function 
    produces correct output.
    """
    input_df = pd.read_csv('output_df_1.csv',sep=';') 
    output_df = pd.read_csv('output_df_2.csv',sep=';') 
    df = rearrange_and_clean_data(input_df)
    assert output_df.equals(df),\
        "The function rearrange_and_clean_data is broken!"  


def test_combine_rows():
    """
    Unit test for combine_rows(df) function to see if the function 
    produces correct output.
    """
    input_df = pd.read_csv('output_df_2.csv',sep=';') 
    output_df = pd.read_csv('output_df_3.csv',sep=';') 
    df = combine_rows(input_df)
    assert output_df.equals(df),\
        "The function combine_rows is broken!" 


def test_remove_NaN():
    """
    Unit test for remove_NaN(df_final) function to see if the function 
    produces correct output.
    """
    input_df = pd.read_csv('output_df_3.csv',sep=';') 
    df_T_half_output = pd.read_csv('df_T_half_test.csv', sep=';')
    df_F_output = pd.read_csv('df_F_test.csv', sep=';')
    df_AUC_output = pd.read_csv('df_AUC_test.csv', sep=';')
    df_Cmax_output = pd.read_csv('df_Cmax_test.csv', sep=';')
    df_Vdss_output = pd.read_csv('df_Vdss_test.csv', sep=';')
    df_CL_output = pd.read_csv('df_CL_test.csv', sep=';')
    df_LD50_output = pd.read_csv('df_LD50_test.csv', sep=';')
    df_IC50_output = pd.read_csv('df_IC50_test.csv', sep=';')
    df_PPB_output = pd.read_csv('df_PPB_test.csv', sep=';')
 
    df_T_half, df_F, df_AUC, df_Cmax, df_Vdss, df_CL, df_LD50,\
    df_IC50, df_PPB = remove_NaN(input_df)
    assert df_T_half.equals(df_T_half_output) and df_F.equals(df_F_output)\
           and df_AUC.equals(df_AUC_output) and df_Cmax.equals(df_Cmax_output)\
           and df_Vdss.equals(df_Vdss_output) and df_CL.equals(df_CL_output)\
           and df_LD50.equals(df_LD50_output) and df_IC50.equals(df_IC50_output)\
           and df_PPB.equals(df_PPB_output), "The function remove_NaN is broken!"   
