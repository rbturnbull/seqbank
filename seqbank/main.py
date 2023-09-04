


def add_file(self, path:Path, format="fasta", all_accessions=None):
    with open_gzip(path) as f:
        for record in SeqIO.parse(f, format):
            self.add(record, record.id, all_accessions=all_accessions
                     

def add_url(self, url:str, format="fasta", all_accessions=None):
    with open_url(url) as f:
        for record in SeqIO.parse(f, format):
            self.add(record, record.id, all_accessions=all_accessions)