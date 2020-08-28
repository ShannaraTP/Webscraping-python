# input is the keyword and the place of search ([title] or [abstract]) as string
# Biopython has to be installed separately (http://biopython.org/DIST/docs/tutorial/Tutorial.html)

def NCBI_scrape_jour(key_word,place):
    from Bio import Entrez
    from Bio import Medline
    # please enter a valid email address (is really necessary)
    Entrez.email = "xxxx@biosustain.dtu.dk"
    # searches PubMed database for publications with the keyword in defined place
    # always turn history on to only query once
    search = Entrez.read(Entrez.esearch(db="pubmed", term=key_word+place, usehistory="y"))
    count = int(search["Count"])
    # retrieving the metadata for each PubMed ID found
    # extract the journal title of each publication
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
    # returns the list of journal titles and the total count of publications found
    return journal, count