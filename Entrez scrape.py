from Bio import Entrez
from Bio import Medline
from collections import Counter
from fuzzywuzzy import fuzz
import matplotlib.pyplot as plt

############## CHANGE THIS ################
key_word = "vindoline"
place = "[title]"
mc = 10
###########################################

Entrez.email = "shtapa@biosustain.dtu.dk"
search = Entrez.read(Entrez.esearch(db="pubmed", term=key_word+place, usehistory="y"))
count = int(search["Count"])
print("%i results found" %(count))

keywords = []
journal = []
batch_size = 100
for s in range(0, count, batch_size):
    end = min(count, s+batch_size)
    print("Downloading record %i to %i" %(s+1, end))
    fetch = Entrez.efetch(db="pubmed", rettype="medline", retmode="text", retstart=s, retmax=batch_size, webenv=search["WebEnv"], query_key=search["QueryKey"])
    records = Medline.parse(fetch)
    records = list(records)

    for r in records:
        l = r.get("OT", ["none"])
        for i in range(len(l)):
            if l[i] == "none":
                continue
            elif fuzz.token_set_ratio(l[i], key_word) >= 90:
                continue
            else:
                word = l[i].replace("*", "")
                keywords.append(word.lower())
        j = r.get("JT", "none")
        if j == "none":
            continue
        else:
            journal.append(j.lower())
    fetch.close

key_nonred = []
for i in range(len(keywords)):
    key_nonred.append(keywords[i])

for i in range(len(keywords)):
    for j in range(len(key_nonred)):
        if keywords[i] == key_nonred[j]:
            continue
        else:
            sim = fuzz.token_set_ratio(keywords[i],keywords[j])
            if sim >= 85:
                if len(keywords[i]) <= len(keywords[j]):
                    key_nonred[j] = keywords[i]

counter = Counter(key_nonred).most_common(mc)

key_abun = []
for c in counter:
    key_abun.append(c[0])
print(key_abun)

j_key = {}
cat = ["general","medical","plant","microbial","chemical"]
for i in range(len(cat)):
    jour = open(r"Journal_keywords\keywords_"+cat[i]+"_journals.txt","r").readlines()
    for k in range(len(jour)):
        j_key.setdefault(cat[i],[]).append(jour[k].replace("\n",""))

pie_plot = {}
for i in range(len(journal)):
    out = True
    for j in j_key:
        if j == "general":
            sim = 100
        else:
            sim = 75
        for k in range(len(j_key[j])):
            s = fuzz.token_set_ratio(journal[i],j_key[j][k])
            if s >= sim:
                pie_plot.setdefault(j,[]).append(journal[i])
                out = False
                break
    if out:
        pie_plot.setdefault("other",[]).append(journal[i])

num,labels = [],[]
for i in pie_plot:
    num.append(len(pie_plot[i]))
    labels.append(i)

for i in range(len(num)):
    if num[i] == 0:
        labels[i] = ""
plt.pie(num,labels=labels,autopct="%1.0f%%",textprops={"fontsize": 12})
# plt.savefig("xxx_pieplot.png",bbox_inches="tight",dpi=300)
plt.show()

#for i in range(0,mc-1):
#    search_key = Entrez.read(Entrez.esearch(db="pubmed", term=key_abun[i]+place, usehistory="y"))