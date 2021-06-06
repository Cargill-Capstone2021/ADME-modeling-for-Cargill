#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sqlite3
import pandas as pd


# In[2]:


def database_file_to_data_frame(db):
    """
    Creates a SQL connection to the SQLite. This function reads SQL database files,
    executes SQL queries in python platform and extracts data for specific table and
    returns a pandas DataFrame
    
    Parameters
    ----------
    - db: database file/ full path of database file

    
    Returns
    -------
    - modified DataFrame based on the user selected tables and columns
    """

    con = sqlite3.connect(db)
    # The desired sql query could be written in """ """ to join the desired table and columns
    df = pd.read_sql_query("""SELECT * FROM source""", con)

    return df


# In[3]:


def merged_data(df1, df2):
    """This function combines two DataFrames
    
    Parameters
    ----------
    - df1,df2: DataFrames to be combined
    
    Returns
    -------
    - merged DataFrame
    """
    df_merged = pd.concat([df1, df2])

    return df_merged
