# Biopython has to be installed separately (http://biopython.org/DIST/docs/tutorial/Tutorial.html)
# fuzzywuzzy has to be installed separately (https://github.com/seatgeek/fuzzywuzzy)
from Bio import Entrez
from Bio import Medline
from collections import Counter
from fuzzywuzzy import fuzz
import matplotlib.pyplot as plt

############## CHANGE THIS ################
key_word = "vindoline"
place = "[title]"
mc = 10 # returned number of most common keywords
###########################################

# please enter a valid email address (is really necessary)
Entrez.email = "xxxx@biosustain.dtu.dk"
# searches PubMed database for publications with the keyword in defined place
# always turn history on to only query once
search = Entrez.read(Entrez.esearch(db="pubmed", term=key_word+place, usehistory="y"))
count = int(search["Count"])
print("%i results found" %(count))
# batch size of 100 is advised when scraping for publications
# retrieving the metadata for each PubMed ID found
# extract the keywords and the journal title of each publication
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
# remove redundancy in list of keywords with string matching
# removes search keyword from list and changes similar spelling to the same word
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
# sorts the keyword list from high to low abundance
counter = Counter(key_nonred).most_common(mc)
# prints x most common keywords
key_abun = []
for c in counter:
    key_abun.append(c[0])
print(key_abun)
# use manually created lists of keywords to categorise the journal titles and thus the publications
# lists are available on the GitHub page in the Journal_kewords folder
j_key = {}
cat = ["general","medical","plant","microbial","chemical"]
for i in range(len(cat)):
    jour = open(r"Journal_keywords/keywords_"+cat[i]+"_journals.txt","r").readlines()
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
# the categorisation is plotted in a pie chart
num,labels = [],[]
for i in pie_plot:
    num.append(len(pie_plot[i]))
    labels.append(i)

for i in range(len(num)):
    if num[i] == 0:
        labels[i] = ""
plt.pie(num,labels=labels,autopct="%1.0f%%",textprops={"fontsize": 12})
# plt.savefig("xxx_pieplot.png",bbox_inches="tight",dpi=300) # to save figure to file
plt.show()

# optional extension of the script to search PubMed again with the most common keywords
#for i in range(0,mc-1):
#    search_key = Entrez.read(Entrez.esearch(db="pubmed", term=key_abun[i]+place, usehistory="y"))
