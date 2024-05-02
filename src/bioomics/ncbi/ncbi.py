"""
download data from NCIB FTP
"""
import os

from ..connector.conn_ftp import ConnFTP
from ..connector.conn_ftplib import ConnFTPlib

# lock the scope of groups in local version
ANATOMY_GROUPS = ['archaea', 'bacteria', 'fungi', 'invertebrate', 'plant',
    'protozoa', 'vertebrate_mammalian', 'vertebrate_other', 'viral', ]


class NCBI(ConnFTP):
    url = 'ftp.ncbi.nlm.nih.gov'

    def __init__(self, local_path:str, overwrite:bool=None):
        super().__init__(url=self.url, overwrite=overwrite)
        self.local_path = os.path.join(local_path, "NCBI")

    def download_assembly_summary(self, groups:list=None):
        '''
        download assembly_summary.txt. That is genome metadata
        '''
        if not groups:
            groups = ANATOMY_GROUPS
        groups = [i for i in groups if i in ANATOMY_GROUPS]
        
        res = {}
        for antonomy in groups:
            outdir = os.path.join(self.local_path, 'assembly_summary', antonomy)
            local_file = self.download_file(
                endpoint = f'genomes/refseq/{antonomy}/',
                file_name = 'assembly_summary.txt',
                local_path = outdir
            )
            res[antonomy] = local_file
        return self.local_path, res

    def download_genome(self, ftp_path:str, specie:str, version:str=None):
        '''
        download genome including subdirectories and files
        '''
        local_path = os.path.join(self.local_path, 'genome', specie)
        if version:
            local_path = os.path.join(local_path, version)

        # download sequences and annotations
        local_files = self.download_files(
            endpoint = ftp_path.replace('https://ftp.ncbi.nlm.nih.gov/', ''),
            match = '.gz',
            local_path = local_path,
        )
        return local_path, local_files
    

    def download_id_mapping(self):
        '''
        '''
        local_files = []
        file_names = ['gene2accession.gz', 'gene2ensembl.gz', 'gene2go.gz', \
            'gene2pubmed.gz', 'gene2refseq.gz', 'gene_group.gz', 'gene_history.gz', \
            'gene_info.gz', 'gene_neighbors.gz', 'gene_orthologs.gz', \
            'gene_refseq_uniprotkb_collab.gz',]
        for file_name in file_names:
            local_file = self.download_file(
                endpoint='gene/DATA',
                file_name=file_name,
                local_path=os.path.join(self.local_path, 'gene', 'DATA'),
            )
            if local_file:
                local_files.append(local_file)
        return local_files


    def download_refseq_uniprotkb(self):
        '''
        download gene_refseq_uniprotkb_collab.gz from 
        https://ftp.ncbi.nlm.nih.gov/refseq/uniprotkb/
        map refeseq ~ uniprotkb
        '''
        local_file = self.download_file(
            endpoint='refseq/uniprotkb/',
            file_name='gene_refseq_uniprotkb_collab.gz',
            local_path=os.path.join(self.local_path, 'refseq'),
        )
        return local_file
    
    def download_refseq_gpff(self):
        local_files = []
        species = ['H_sapiens', 'D_rerio', 'B_taurus', 'M_musculus',\
            'R_norvegicus', 'S_scrofa', 'X_tropicalis']
        for sub in species:
            local_files += self.download_files(
                local_path = os.path.join(self.local_path, 'refseq', 'mRNA_Prot'),
                endpoint = f'refseq/{sub}/mRNA_Prot',
                match = '.gpff.gz$'
            )
        return local_files

    def download_refseq_complete_gpff(self):
        local_files = []
        species = ['vertebrate_mammalian',]
        for sub in species:
            local_files += self.download_files(
            local_path = os.path.join(self.local_path, 'refseq', 'release', 'gpff'),
            endpoint = f'refseq/release/{sub}/',
            match = '.gpff.gz$'
        )
        return local_files

    def download_protein_fasta(self):
        local_files = []
        for i in range(1, 260):
            local_file = self.download_file(
                endpoint='ncbi-asn1/protein_fasta/',
                file_name=f"gbbct{i}.fsa_aa.gz",
                local_path=os.path.join(self.local_path, 'protein_fasta'),
            )
            if local_file:
                local_files.append(local_file)
        return local_files


    def download_gene_data(self):
        '''
        download /gene/DATA including subdirectories and files
        '''
        local_files = self.download_tree(
            local_path = os.path.join(self.local_path, 'gene', 'DATA'),
            endpoint = 'gene/DATA',
            match = '.gz$'
        )
        return local_files

    def download_pubmed(self):
        '''
        download /PubMed including subdirectories and files
        '''
        res = ConnFTPlib.download_tree(
            ftp_endpoint = self.url,
            ftp_path = '/pubmed',
            match = '.gz',
            local_path = self.local_path
        )
        return res