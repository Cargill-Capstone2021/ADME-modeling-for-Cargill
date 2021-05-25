import pandas as pd
import numpy as np

import sys
sys.path.append('../')
from data_prep import pre_clean_data
from data_prep import rearrange_and_clean_data
from data_prep import combine_rows
from data_prep import data_prep


def test_pre_clean_data():
    """
    Unit test for pre_clean_data(df) function to see if the function 
    produces correct output.
    """
    input_df = pd.read_csv('input_df_1.csv',sep=';') 
    output_df = pd.read_csv('output_df_1.csv',sep=';') 
    df = pre_clean_data(input_df)
    assert output_df.equals(df),\
        "The function encode_column is broken!"  


def test_rearrange_and_clean_data():
    """
    Unit test for rearrange_and_clean_data(df) function to see if the function 
    produces correct output.
    """
    input_df = pd.read_csv('output_df_1.csv',sep=';') 
    output_df = pd.read_csv('output_df_2.csv',sep=';') 
    df = rearrange_and_clean_data(input_df)
    assert output_df.equals(df),\
        "The function encode_column is broken!"  


def test_combine_rows():
    """
    Unit test for combine_rows(df) function to see if the function 
    produces correct output.
    """
    input_df = pd.read_csv('output_df_2.csv',sep=';') 
    output_df = pd.read_csv('output_df_3.csv',sep=';') 
    df = combine_rows(input_df)
    assert output_df.equals(df),\
        "The function encode_column is broken!" 


def test_data_prep():
    """
    Unit test for data_prep(df) function to see if the function 
    produces correct output.
    """
    input_df = pd.read_csv('input_df_1.csv',sep=';')
    output_df = pd.read_csv('final_output_df.csv',sep=';') 
    df = data_prep(input_df)
    assert output_df.equals(df),\
        "The function encode_column is broken!"   
