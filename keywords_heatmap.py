from NCBI_scrape_key import NCBI_scrape_key
from collections import Counter
from fuzzywuzzy import fuzz
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sn

MIA = ["ajmalicine","ajmaline","akuammine","alstonine","arbophyllidine","burnamine","cadambine",\
    "camptothecin","conolidine","corynanthine","ibogaine","irinotecan","lochnerinine","melotenuine D","mitragynine",\
    "mitraphylline","pericine","psychollatine","umbellatine","quinidine","quinine","rauwolscine","reserpine",\
    "rhazyaminine","serpentine","tabersonine","topotecan","vallesiachotamine","vinblastine","vincamine","vincristine",\
    "vindesine","vinflunine","vinorelbine","yohimbine"]

map_data = {}
for i in range(len(MIA)):
    key = NCBI_scrape_key(MIA[i],"[title]")
    if not key:
        continue
    else:
        map_data.setdefault(MIA[i],key)
    print("Downloaded keywords of " + MIA[i])

tot_key = []
for i in map_data:
    for j in range(len(map_data[i])):
        if fuzz.token_set_ratio("roseus",map_data[i][j]) == 100:
            map_data[i][j] = "catharanthus roseus"
        if fuzz.token_set_ratio("alkaloid",map_data[i][j]) >= 90:
            map_data[i][j] = "alkaloids"
        tot_key.append(map_data[i][j])

counter = Counter(tot_key).most_common()

key_count = []
for c in counter:
    key_count.append(c[0])

print("Constructing counting matrix")
data = np.zeros((len(map_data),len(key_count)))
for i in range(len(key_count)):
    for j in map_data:
        for k in range(len(map_data[j])):
            if key_count[i] == map_data[j][k]:
                data[list(map_data).index(j),i] = 1
    if len(key_count[i]) > 25:
        for l in range(24,len(key_count[i])):
            if key_count[i][l] == " ":
                key_count[i] = key_count[i][0:l]
                break

ax = sn.clustermap(data[:,0:49],yticklabels=[*map_data],xticklabels=key_count[0:49],figsize=(14,18),cmap="winter")
plt.setp(ax.ax_heatmap.get_yticklabels(),rotation=0)
plt.setp(ax.ax_heatmap.get_xticklabels(),rotation=45,ha="right")
plt.savefig("heatmap_complete.png",bbox_inches="tight",dpi=300)