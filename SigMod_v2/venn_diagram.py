from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn2
from pandas import DataFrame

df = DataFrame.from_csv("/home/vcabeli/stage_save/SigMod_v2/tab3b.tsv", sep="\t")
liste_bp=df.index.tolist()
ds = DataFrame.from_csv("/home/vcabeli/stage_save/SigMod_v2/tab3b_schizophrenia.tsv", sep="\t")
liste_sc=ds.index.tolist()

set1 = set(liste_bp)
set2 = set(liste_sc)

venn2([set1, set2], ('bipolarity', 'Schizophrenia'))
plt.show()