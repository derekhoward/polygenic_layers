import pandas as pd

### 1st sheet: compares CP vs SVZ
df1 = pd.read_excel('./Genes_216738_284_17896_v1.xls', sheetname=0)
df1 = df1[df1.pval < 0.05]

# compare cp with svz
# get genes expressed in CP
cp_genes = df1[df1.loc[:, ['b1-cp1', 'b1-cp2', 'b2-cp1', 'b2-cp2']].mean(axis=1) > df1.loc[:, ['b1-svz1', 'b1-svz2', 'b2-svz1', 'b2-svz2']].mean(axis=1)]
svz_genes = df1[df1.loc[:, ['b1-cp1', 'b1-cp2', 'b2-cp1', 'b2-cp2']].mean(axis=1) < df1.loc[:, ['b1-svz1', 'b1-svz2', 'b2-svz1', 'b2-svz2']].mean(axis=1)]

print(cp_genes.shape)
print(svz_genes.shape)

cp_genes.to_csv('./mouse_CP_genes_vs_svz.csv', index=None)
svz_genes.to_csv('./mouse_SVZ_genes_vs_cp.csv', index=None)

### 2nd sheet compares CP vs VZ
df2 = pd.read_excel('./Genes_216738_284_17896_v1.xls', sheetname=1)
df2 = df2[df2.pval < 0.05]

cp_genes = df2[df2.loc[:, ['b1-cp1', 'b1-cp2', 'b2-cp1', 'b2-cp2']].mean(axis=1) > df2.loc[:, ['b1-vz1', 'b1-vz2', 'b2-vz1', 'b2-vz2']].mean(axis=1)]
vz_genes = df2[df2.loc[:, ['b1-cp1', 'b1-cp2', 'b2-cp1', 'b2-cp2']].mean(axis=1) < df2.loc[:, ['b1-vz1', 'b1-vz2', 'b2-vz1', 'b2-vz2']].mean(axis=1)]

print(cp_genes.shape)
print(vz_genes.shape)

cp_genes.to_csv('./mouse_CP_genes_vs_vz.csv', index=None)
svz_genes.to_csv('./mouse_VZ_genes_vs_cp.csv', index=None)

### 3rd sheet: compares VZ vs SVZ
df3 = pd.read_excel('./Genes_216738_284_17896_v1.xls', sheetname=2)
df3 = df3[df3.pval < 0.05]

vz_genes = df3[df3.loc[:, ['b1-vz1', 'b1-vz2', 'b2-vz1', 'b2-vz2']].mean(axis=1) > df3.loc[:, ['b1-svz1', 'b1-svz2', 'b2-svz1', 'b2-svz2']].mean(axis=1)]
svz_genes = df3[df3.loc[:, ['b1-vz1', 'b1-vz2', 'b2-vz1', 'b2-vz2']].mean(axis=1) < df3.loc[:, ['b1-svz1', 'b1-svz2', 'b2-svz1', 'b2-svz2']].mean(axis=1)]

cp_genes.to_csv('./mouse_VZ_genes_vs_svz.csv', index=None)
svz_genes.to_csv('./mouse_SVZ_genes_vs_vz.csv', index=None)
