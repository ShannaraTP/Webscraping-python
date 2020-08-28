# input is the keyword as string and number of months between posted and today as value
# biorxiv_retriever has to be installed separately (https://github.com/TalSchuster/BiorxivRetriever)

def biorxiv_scrape(key_word,months):
    from biorxiv_retriever import BiorxivRetriever
    import datetime as dt

    # to retrieve the date x months before today in the correct format
    today = dt.date.today()
    diff = dt.timedelta(weeks=months*4)
    date_thres = (today - diff).strftime("%Y-%m-%d")

    # to scrape the biorxiv database (only metadate retrieved)
    br = BiorxivRetriever()
    papers = br.query(key_word, full_text=False)

    # papers is a list of dictionaries, so essential information extracted from each dictionary
    # of a publication posted within the last x months
    pub_title,pub_date,pub_link = [],[],[]
    for p in papers:
        date = dt.datetime.strptime(p["posted"],"%B %d, %Y.").strftime("%Y-%m-%d")
        if date_thres < date:
            pub_title.append(p["title"])
            pub_date.append(p["posted"])
            pub_link.append(p["biorxiv_url"])
        else:
            break

    # output of the function is the title, date posted and URL of each publication found
    # within the last x months
    return pub_title, pub_date, pub_link