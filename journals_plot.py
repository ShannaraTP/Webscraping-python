from NCBI_scrape_jour import NCBI_scrape_jour
from fuzzywuzzy import fuzz
import numpy as np
import matplotlib.pyplot as plt

MIA = ["ajmalicine","ajmaline","akuammine","alstonine","arbophyllidine","burnamine","cadambine",\
    "camptothecin","conolidine","corynanthine","ibogaine","irinotecan","lochnerinine","melotenuine D","mitragynine",\
    "mitraphylline","pericine","psychollatine","umbellatine","quinidine","quinine","rauwolscine","reserpine",\
    "rhazyaminine","serpentine","tabersonine","topotecan","vallesiachotamine","vinblastine","vincamine","vincristine",\
    "vindesine","vinflunine","vinorelbine","yohimbine"]

MIA_jour = {}
MIA_found, count = [], []
for i in range(len(MIA)):
    jour, cnt = NCBI_scrape_jour(MIA[i],"[title]")
    if not jour:
        continue
    else:
        MIA_jour.setdefault(MIA[i],jour)
        MIA_found.append(MIA[i])
        count.append(cnt)
    print("Downloaded journals of " + MIA[i])

j_key = {}
cat = ["general","medical","plant","microbial","chemical"]
for i in range(len(cat)):
    jour = open(r"Journal_keywords\keywords_"+cat[i]+"_journals.txt","r").readlines()
    temp = []
    for k in range(len(jour)):
        j_key.setdefault(cat[i],[]).append(jour[k].replace("\n",""))

data_plot = np.zeros((len(cat)+1,len(MIA_found)))
for i in MIA_jour:
    for j in range(len(MIA_jour[i])):
        out = True
        for k in j_key:
            if k == "general":
                sim = 100
            else:
                sim = 75
            for l in range(len(j_key[k])):
                s = fuzz.token_set_ratio(MIA_jour[i][j],j_key[k][l])
                if s >= sim:
                    data_plot[cat.index(k),MIA_found.index(i)] += 1
                    out = False
                    break
        if out:
            data_plot[5,MIA_found.index(i)] += 1
print("Counting complete")

data_plot_per = np.zeros((len(cat)+1,len(MIA_found)))
for i in range(len(cat)+1):
    for j in range(len(MIA_found)):
        data_plot_per[i,j] = data_plot[i,j]/np.sum(data_plot[:,j])*100

cat.append("other")
r = range(len(MIA_found))
plt.figure(figsize=(18,8))
plt.bar(r,data_plot_per[0,:],width=0.9,label=cat[0])
for i in range(1,len(cat)):
    plt.bar(r,data_plot_per[i,:],bottom=np.sum(data_plot_per[0:i],axis=0),width=0.9,label=cat[i])
plt.xticks(r,MIA_found,rotation=45,ha="right",fontsize=12)
plt.ylabel("Percentage (%)")
plt.legend(bbox_to_anchor=(0,1.02,1,0.2),loc="lower center",ncol=len(cat))
for i in range(len(count)):
    plt.text(i,50,str(count[i]),fontsize=9,ha="center")
plt.savefig("barplot_complete.png",bbox_inches="tight",dpi=300)