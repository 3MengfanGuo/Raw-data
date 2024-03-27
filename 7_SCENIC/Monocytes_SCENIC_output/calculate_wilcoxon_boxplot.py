import numpy as np
import pandas as pd
from scipy.stats import ranksums

df1 = pd.read_csv("./meta.txt", header=0, index_col=0, sep="\t")
df2 = pd.read_csv("./AUCmatrix_regulon.txt", header=0, index_col=0, sep="\t")
df = pd.concat([df1, df2], axis=1, join="inner")
total_arr = []
for i in range(1, df.shape[1]):
	sub_df = df.iloc[:, [0, i]]
	sub_df["TF"] = sub_df.columns[1]
	sub_df.columns = ["Tissue", "Score", "TF"]
	sub_df = sub_df.loc[sub_df["Score"]>0, :]
	x = sub_df.loc[sub_df["Tissue"]=="primary tumor", "Score"]
	y = sub_df.loc[sub_df["Tissue"]=="metastasis", "Score"]
	t, p = ranksums(x, y)
	if p < 0.05:
		total_arr.append(sub_df)
total_df = pd.concat(total_arr, axis=0)
total_df.to_csv("./Monocytes_diff_TF_boxplot.txt", header=True, index=True, sep="\t")





