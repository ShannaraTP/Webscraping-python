# input is the keyword and the place of search ([title] or [abstract]) as a string

def NCBI_scrape_jour(key_word,place):
    from Bio import Entrez
    from Bio import Medline

    Entrez.email = "shtapa@biosustain.dtu.dk"
    search = Entrez.read(Entrez.esearch(db="pubmed", term=key_word+place, usehistory="y"))
    count = int(search["Count"])

    journal = []
    batch_size = 500
    for s in range(0, count, batch_size):
        fetch = Entrez.efetch(db="pubmed", rettype="medline", retmode="text", retstart=s, retmax=batch_size, webenv=search["WebEnv"], query_key=search["QueryKey"])
        records = Medline.parse(fetch)
        records = list(records)

        for r in records:
            j = r.get("JT", "none")
            if j == "none":
                continue
            else:
                journal.append(j.lower())
        fetch.close

    return journal, count