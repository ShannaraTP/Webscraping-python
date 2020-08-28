# input is the keyword and the place of search ([title] or [abstract]) as a string

def NCBI_scrape_key(key_word,place):
    from Bio import Entrez
    from Bio import Medline
    from collections import Counter
    from fuzzywuzzy import fuzz

    Entrez.email = "shtapa@biosustain.dtu.dk"
    search = Entrez.read(Entrez.esearch(db="pubmed", term=key_word+place, usehistory="y"))
    count = int(search["Count"])

    keywords = []
    batch_size = 500
    for s in range(0, count, batch_size):
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
        fetch.close

    key_nonred = []
    for i in range(len(keywords)):
        key_nonred.append(keywords[i])

    for i in range(len(keywords)):
        for j in range(len(key_nonred)):
            if keywords[i] == key_nonred[j]:
                continue
            else:
                if fuzz.token_set_ratio(keywords[i], keywords[j]) >= 85:
                    if len(keywords[i]) <= len(keywords[j]):
                        key_nonred[j] = keywords[i]

    counter = Counter(key_nonred).most_common()

    key_count = []
    for c in counter:
        key_count.append(c[0])

    return key_count