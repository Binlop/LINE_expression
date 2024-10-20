import requests
from bs4 import BeautifulSoup
import schedule
import time
from urllib.parse import urlencode


class BlastProtein:

    NCBI_BLAST_URL = "https://blast.ncbi.nlm.nih.gov/Blast.cgi"

    def __init__(self, seq: str):
        self.seq = seq
        self.blast_data = []

    def blast(self):
        params = {
            "DATABASE": "nr",
            "ENTREZ_QUERY": "(none)",
            "EXPECT": "10.0",
            "HITLIST_SIZE": "50",
            "PROGRAM": "blastx",
            "QUERY": self.seq,
            "CMD": "Put"
        }
            
        rid = None
        i = [0] 
        message = urlencode(params).encode()

        headers = {"User-Agent": "BiopythonClient"}
        response = requests.put(self.NCBI_BLAST_URL, data=message, headers=headers)

        response.raise_for_status()

        soup = BeautifulSoup(response.content, "html.parser")
        rid = soup.find("input", {"name": "RID"})['value']

        finished = [False]

        schedule.every(10).seconds.do(lambda: self.check_status(rid, i, finished))

        while not finished[0]:
            schedule.run_pending()
            time.sleep(5)

    def check_status(self, rid, i, finished):
        i[0] += 1 
        params = {
            'RID': rid,
            'CMD': 'Get',
            # 'FORMAT_TYPE': 'XML'
        }
        response = requests.get(self.NCBI_BLAST_URL, params=params)
        soup = BeautifulSoup(response.content, "html.parser")
        table_stat = soup.find("table", {"id": "statInfo"})
        table_data = soup.find("table", {"id": "dscTable"})

        if table_data:
            self.extract_data(table_data) 
            finished[0] = True

        elif table_stat:
            rows = table_stat.find_all("tr")

            for row in rows:
                cells = row.find_all("td")
                if len(cells) > 1 and cells[0].get_text(strip=True) == "Status":
                    status_value = cells[1].get_text(strip=True)
                    break

            print(f"Status: {status_value}, запрос {i[0]}")
            if status_value == "Searching":
                schedule.every(10).seconds.do(lambda: self.check_status(rid, i, finished))


    def extract_data(self, table_data):
        rows = table_data.find_all("tr")
        for row in rows:
            id_cell = row.find("td", {"class": "l c0"})
            function_cell = row.find("td", {"class": "ellipsis c2"})
            coincidence_cell = row.find("td", {"class": "c6"})
            if id_cell and function_cell and coincidence_cell is not None:
                id_value = id_cell.find("input", {"name":"getSeqGi"})["value"]
                function_value = function_cell.find("a", class_="deflnDesc").get_text(strip=True)
                coincidence_value = coincidence_cell.get_text(strip=True)
                data = {'id': id_value, 'fucntion': function_value, 'coincidence': coincidence_value}
                self.blast_data.append(data)

        return self.blast_data


if __name__ == "__main__":
    sequence = "CTTTGCTGGGAAGAAAAGCAGAGGGCCTTTCCTGTCGCTTGCATCATCTTTCCTGCTCTGGTTTTATCTTGGAAACCTTGGGATCGAAGAAATTTTCCCAGTTCGTTTTCTTCTTGAGACAAGAA"
    blast = BlastProtein(seq=sequence)
    blast.blast()
    print(blast.blast_data[0])