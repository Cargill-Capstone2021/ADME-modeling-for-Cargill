import Descriptors
import pandas as pd

from rdkit import Chem

from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

import matplotlib


matplotlib.rcParams.update({'font.size': 18})


data_IC_50 = pd.read_csv('adme/data/IC_50_final.csv')
df = data_IC_50
# half-maximal inhibitory concentration (how much drug is
# needed to inhibit a biological process by half, thus providing
# a measure of potency of an antagonist drug.
df = df[df['IC50 /nM'] < 500]  # removing outliers


df1 = Descriptors.generate(df.canonical_smiles)
df2 = Descriptors.generate_1(df.canonical_smiles)
df3 = Descriptors.generate_2(df.canonical_smiles)
df4 = Descriptors.generate_3(df.canonical_smiles)
df5 = Descriptors.generate_4(df.canonical_smiles)
df6 = Descriptors.generate_5(df.canonical_smiles)
# df7 = Descriptors.generate_6(df.canonical_smiles) descriptors not relevant
df8 = Descriptors.generate_7(df.canonical_smiles)
df9 = Descriptors.generate_8(df.canonical_smiles)
df10 = Descriptors.generate_9(df.canonical_smiles)
df11 = Descriptors.generate_10(df.canonical_smiles)
df12 = Descriptors.generate_11(df.canonical_smiles)
# df13 = Descriptors.generate_12(df.canonical_smiles) descriptors not relevant
# df14 = Descriptors.generate_13(df.canonical_smiles) descriptors not relevant
df15 = Descriptors.generate_14(df.canonical_smiles)
df16 = Descriptors.generate_15(df.canonical_smiles)
df17 = Descriptors.generate_16(df.canonical_smiles)
df18 = Descriptors.generate_17(df.canonical_smiles)
df19 = Descriptors.generate_18(df.canonical_smiles)


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

print('Coefficients:', regressor.coef_)
print('Intercept:', regressor.intercept_)
print('Mean squared error (MSE): %.2f'
      % mean_squared_error(y_test, y_pred_test))
print('Coefficient of determination (R^2): %.2f'
      % r2_score(y_test, y_pred_test))
