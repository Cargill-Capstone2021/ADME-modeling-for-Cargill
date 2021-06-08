import Descriptors.py

import numpy as np
import pandas as pd

from rdkit import Chem

import sklearn
from sklearn.model_selection import train_test_split
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn import neighbors, svm
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix
from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVC, SVR
from sklearn.ensemble import BaggingClassifier,\
    RandomForestRegresson,\
    GradientBoostingClassifier
from sklearn.model_selection import GridSearchCV

import matplotlib
# matplotlib.use('GTKAgg')
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap
%matplotlib inline
matplotlib.rcParams.update({'font.size': 18})


df = data_IC_50
# half-maximal inhibitory concentration (how much drug is
# needed to inhibit a biological process by half, thus providing
# a measure of potency of an antagonist drug.
df = df[df['IC50 /nM'] < 500]  # removing outliers


df1 = generate(df.canonical_smiles)
df2 = generate_1(df.canonical_smiles)
df3 = generate_2(df.canonical_smiles)
df4 = generate_3(df.canonical_smiles)
df5 = generate_4(df.canonical_smiles)
df6 = generate_5(df.canonical_smiles)
# df7 = generate_6(df.canonical_smiles) # descriptors not relevant
df8 = generate_7(df.canonical_smiles)
df9 = generate_8(df.canonical_smiles)
df10 = generate_9(df.canonical_smiles)
df11 = generate_10(df.canonical_smiles)
df12 = generate_11(df.canonical_smiles)
# df13 = generate_12(df.canonical_smiles) # descriptors not relevant
# df14 = generate_13(df.canonical_smiles) * descriptors not relevant
df15 = generate_14(df.canonical_smiles)
df16 = generate_15(df.canonical_smiles)
df17 = generate_16(df.canonical_smiles)
df18 = generate_17(df.canonical_smiles)
df19 = generate_18(df.canonical_smiles)


mol_list = []
for element in df.canonical_smiles:
    mol = Chem.MolFromSmiles(element)
    mol_list.append(mol)

desc_AromaticProportion = [AromaticAtoms(element) /
                           Descriptors.HeavyAtomCount(element)
                           for element in mol_list]
df_desc_AromaticProportion = pd.DataFrame(desc_AromaticProportion,
                                          columns=['AromaticProportion'])

X = pd.concat([df1, df2, df3, df4, df5, df6,
               df8, df9, df10, df11, df12,
               df15, df16, df17, df18, df19], axis=1)

y = df.iloc[:, 3]

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2)

X_train = StandardScaler().fit(X_train).transform(X_train)
X_test = StandardScaler().fit(X_test).transform(X_test)

regressor = linear_model.LinearRegression()
regressor.fit(X_train, y_train)

y_pred_test = regressor.predict(X_test)

print('Coefficients:', model.coef_)
print('Intercept:', model.intercept_)
print('Mean squared error (MSE): %.2f'
      % mean_squared_error(y_test, y_pred_test))
print('Coefficient of determination (R^2): %.2f'
      % r2_score(y_test, y_pred_test))
