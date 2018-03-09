from matplotlib import pyplot as plt
import numpy as np
from matplotlib_venn import venn2
from pandas import DataFrame

df = DataFrame.from_csv("file:///home/vcabeli/stage/GRIV/HIV_interaction_intersect.tab", sep="\t")
liste_pr=df.index.tolist()

ds = DataFrame.from_csv("file:///home/vcabeli/stage/GRIV/np/HIV_interaction_intersect.tab", sep="\t")
liste_np=ds.index.tolist()

set1 = set(liste_pr)
set2 = set(liste_np)

v=venn2([set1, set2], ('PR', 'NP'))
v.get_label_by_id('11').set_text(list(set1.intersection(set2)))
plt.show()