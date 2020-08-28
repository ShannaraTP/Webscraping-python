# input is the keyword and the place of search ([title] or [abstract]) as string
# and number of output publications required as value
# Biopython has to be installed separately (http://biopython.org/DIST/docs/tutorial/Tutorial.html)

def citation_scrape(key_word,place,num):
    from Bio import Entrez
    from Bio import Medline
    from collections import Counter

    # please enter a valid email address (is really necessary)
    Entrez.email = "xxxx@biosustain.dtu.dk"
    # searches PubMed database for publications with the keyword in defined place
    # always turn history on to only query once
    # returns a list of PubMed ID's (not complete for some reason, could not find why)
    search = Entrez.read(Entrez.esearch(db="pubmed", term=key_word+place, usehistory="y"))
    count = int(search["Count"])

    pmid = []
    # batch size of 100 is advised when scraping for publications
    batch_size = 100
    for s in range(0, count, batch_size):
        # retrieving the metadata for each PubMed ID found
        fetch = Entrez.efetch(db="pubmed",rettype="medline",retmode="text",retstart=s,retmax=batch_size,webenv=search["WebEnv"],query_key=search["QueryKey"])
        records = Medline.parse(fetch)
        records = list(records)
        # extract the PubMed ID's (complete this time (matches the count))
        for r in records:
            p = r.get("PMID", "none")
            if p == "none":
                continue
            else:
                pmid.append(p)
        fetch.close

    # retrieve PMC ID's of publication which cite a publication found from the search
    # convert PMC ID's back to PubMed ID's afterwards
    pmid_cit = {}
    for i in pmid:
        r1 = Entrez.read(Entrez.elink(dbfrom="pubmed",db="pmc",LinkName="pubmed_pmc_refs",id=i))
        if not r1[0]["LinkSetDb"]:
            continue
        else:
            pmc = [l["Id"] for l in r1[0]["LinkSetDb"][0]["Link"]]
            r2 = Entrez.read(Entrez.elink(dbfrom="pmc",db="pubmed",LinkName="pmc_pubmed",id=",".join(pmc)))
            if not r2[0]["LinkSetDb"]:
                continue
            else:
                p_cit = [l["Id"] for l in r2[0]["LinkSetDb"][0]["Link"]]
                pmid_cit.setdefault(i,p_cit)

    # find the x number of publications which are cited the most
    most_cit,cit = [],[]
    for i in range(num):
        m = max(pmid_cit, key = lambda x: len(set(pmid_cit[x])))
        most_cit.append(m)
        cit.append(len(pmid_cit[m]))
        del(pmid_cit[m])
    # retrieve metadata from the x most cited and the x newest publications
    # the ID list from the search is in chronological order, so the top is always the newest
    new_cit = list(pmid_cit)
    title_new,title_most,date_new,date_most = [],[],[],[]
    for i in range(num):
        f1 = Entrez.efetch(db="pubmed",rettype="medline",retmode="text",id=most_cit[i])
        record1 = Medline.parse(f1)
        record1 = list(record1)
        f2 = Entrez.efetch(db="pubmed",rettype="medline",retmode="text",id=new_cit[i])
        record2 = Medline.parse(f2)
        record2 = list(record2)

        for r in record1:
            t = r.get("TI","none")
            title_most.append(t)
            d = r.get("DP","none")
            date_most.append(d)
        for r in record2:
            t = r.get("TI","none")
            title_new.append(t)
            d = r.get("DP","none")
            date_new.append(d)
    # returns title and publication date for x most cited and newest publications
    # and number of times cited of the most cited publications
    return date_new, title_new, date_most, title_most, cit