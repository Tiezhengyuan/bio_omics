'''
GENE-DB

fasta files:
* The file IMGTGENE-DB-ReferenceSequences.fasta-nt-WithGaps-F+ORF+inframeP 
    includes IMGT/GENE-DB nucleotide reference sequences for functional, 
    open reading frame and in-frame pseudogene genes and alleles, 
    with IMGT gaps for V and C genes and alleles.
* The file IMGTGENE-DB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+inframeP
    includes IMGT/GENE-DB nucleotide reference sequences for functional, 
    open reading frame and in-frame pseudogene genes and alleles, without IMGT gaps.
* The file IMGTGENE-DB-ReferenceSequences.fasta-nt-WithoutGaps-F+ORF+allP
    includes IMGT/GENE-DB nucleotide reference sequences for functional, 
    open reading frame and all pseudogene genes and alleles, without IMGT gaps.
* The file IMGTGENE-DB-ReferenceSequences.fasta-AA-WithGaps-F+ORF+inframeP
    includes IMGT/GENE-DB amino acid reference sequences for functional, 
    open reading frame and in-frame pseudogene genes and alleles, 
    with IMGT gaps for V and C genes and alleles.
* The file IMGTGENE-DB-ReferenceSequences.fasta-AA-WithoutGaps-F+ORF+inframeP
    includes IMGT/GENE-DB amino acid reference sequences for functional, 
    open reading frame and in-frame pseudogene genes and alleles, without IMGT gaps.
'''
from Bio import SeqIO
import re
import os
import pandas as pd

from .imgt import Imgt
from ..protein.process_igblast import ProcessIgblast
from ..common.dir import Dir

class ImgtGenedb(Imgt):

    def __init__(self, data_dir:str, verbose:bool=False):
        super().__init__(data_dir, verbose)
        self.data_dir = os.path.join(self.data_dir, 'GENE-DB')
        # keys: specie, gene_name, region_name
        self.data = {}
    
    def build_gene_records(self):
        # read fasta and integrate sequences to self.data
        self.data = {}
        self.parse_fasta()

        germline_id, records, seq_records = 0, [], []
        for specie, sub1 in self.data.items():
            for gene_name, sub2 in sub1.items():
                for region_name, sub3 in sub2.items():
                    for allele_name, sub4 in sub3.items():
                        germline_id += 1
                        rec = {
                            'germline_id': germline_id,
                            'specie': specie,
                            'isotype': gene_name[:2].upper(),
                            'chain_type': gene_name[2].upper(),
                            'gene_name': gene_name,
                            'allele_name': allele_name,
                            'region_name': region_name,
                        }
                        records.append(rec)
                        for k,v in sub4['seq'].items():
                            seq_rec = {
                                'germline_id': germline_id,
                                'seq_type': k,
                                'seq': v.upper(),
                            }
                            seq_records.append(seq_rec)
        return records, seq_records

    def build_specie_genelist_json(self):
        self.data = {}
        # parse gene list to self.data
        df = self.parse_genelist(True)
        g = df.groupby(['Species', 'IMGT/GENE-DB'])
        for (specie, gene_name), sub in g:
            if specie not in self.data:
                self.data[specie] = {}
            if gene_name not in self.data[specie]:
                self.data[specie][gene_name] = {}
            sub = sub.to_dict(orient='records')
            self.data[specie][gene_name]['genelist'] = sub
        # read fasta and integrate sequences to self.data
        self.parse_fasta()

        # save data
        outdir = os.path.join(self.data_dir, 'genelist')
        Dir(outdir).init_dir()
        for specie, data1 in self.data.items():
            json_file = os.path.join(outdir, f"{specie}.json")
            self.save_json(data1, json_file)
        if self.verbose:
            print('Export genes to ', outdir)
    
    def parse_genelist(self, to_df=False):
        '''
        # retrieve data from GeneList
        IMGTGENEDB-GeneList.txt includes the list of genes 
            with the following fields separated by ";":
        - the species
        - the IMGT/GENE-DB gene name
        - the IMGT gene functionality 
        - the IMGT and HGNC gene definition
        - the number of alleles 
        - the chromosome
        - the chromosomal localization
        - the IMGT/LIGM-DB reference sequence for the allele *01
        - the NCBI Gene ID
        '''
        data = []
        infile = os.path.join(self.data_dir, 'IMGTGENEDB-GeneList.txt')
        with open(infile, 'r') as f:
            header = next(f).rstrip()
            header = header.split(';')
            for line in f:
                line = line.rstrip()
                items = line.split(';')
                if len(items) > len(header):
                    items = items[:7]+ [';'.join(items[7:-2]),] + items[-2:]
                # format specie
                items[0] = items[0].replace(' ', '_')
                data.append(items)
        #
        if to_df:
            df = pd.DataFrame(data, columns=header)
            df = df.drop(['',], axis=1)
            df = self.label_chain(df)
            return df
        return data

    def label_chain(self, df):
        def func(x):
            col = 'IMGT and HGNC gene definition'
            # IG chain:heavy, domain: V C and J
            s = str(x[col]).lower()
            if 'immunoglobulin heavy variable' in s:
                return pd.Series(['heavy', 'v'])
            if 'immunoglobulin heavy constant' in s:
                return pd.Series(['heavy', 'c'])
            if 'immunoglobulin heavy joining' in s:
                return pd.Series(['heavy', 'j'])
            if 'immunoglobulin heavy diversity' in s:
                return pd.Series(['heavy', 'd'])
            
            # light chain
            if 'immunoglobulin lambda constant' in s:
                return pd.Series(['lambda', 'c'])
            col = 'IMGT/GENE-DB'
            if 'IGLV' in str(x[col]):
                return pd.Series(['lambda', 'v'])
            if 'IGLJ' in str(x[col]):
                return pd.Series(['lambda', 'j'])
            # IG kappa, V and J
            if 'immunoglobulin kappa constant' in s:
                return pd.Series(['kappa', 'c'])
            if 'IGKV' in str(x[col]):
                return pd.Series(['kappa', 'v'])
            if 'IGKJ' in str(x[col]):
                return pd.Series(['kappa', 'j'])
            # t cell receptor
            if 't cell receptor' in s:
                return pd.Series(['tcr', ''])
            return pd.Series(['other', ''])
        
        df[['chain', 'domain']] = df.apply(func, axis=1)
        return df
    
    def parse_fasta(self):
        '''
        update self.data
        '''
        for postfix, parser in self.parse_fasta_file():
            for record in parser:
                items = record.description.split('|')
                allele_name, region_name = items[1], items[4]
                gene_name = allele_name.split('*')[0]
                organism = re.split(r'[ |_]', items[2])
                organism = f'{organism[0].title()}_{organism[1].lower()}'
                if organism not in self.data:
                    self.data[organism] = {}
                if gene_name not in self.data[organism]:
                    self.data[organism][gene_name] = {}
                if region_name not in self.data[organism][gene_name]:
                    self.data[organism][gene_name][region_name] = {}
                if allele_name not in self.data[organism][gene_name][region_name]:
                    self.data[organism][gene_name][region_name][allele_name] = {
                        'specie': items[2],
                        'gene_name': gene_name,
                        'region': region_name,
                        'allele_name': allele_name,
                        'seq': {},
                    }
                self.data[organism][gene_name][region_name][allele_name]['seq'][postfix] = str(record.seq)

    def parse_fasta_file(self):
        prefix = 'IMGTGENEDB-ReferenceSequences.fasta'
        names = [
            'nt-WithGaps-F+ORF+inframeP',
            'nt-WithoutGaps-F+ORF+inframeP',
            'nt-WithoutGaps-F+ORF+allP',
            'AA-WithGaps-F+ORF+inframeP',
            'AA-WithoutGaps-F+ORF+inframeP',
        ]
        for postfix in names:
            data1 = {}
            file_name = f"{prefix}-{postfix}"
            infile = os.path.join(self.data_dir, file_name)
            with open(infile, 'r') as f:
                parser = SeqIO.parse(f, 'fasta')
                yield postfix, parser

    def region_fasta_file(self, ref_type, region_name):
        outdir = os.path.join(self.data_dir, ref_type, 'region')
        Dir(outdir).init_dir()
        return os.path.join(outdir, f"{region_name}.fasta")

    def build_region_fasta(self):
        for postfix, parser in self.parse_fasta_file():
            data = {}
            for record in parser:
                items = record.description.split('|')
                allele_name, region_name = items[1], items[4]
                species = re.split(r'[ |_]', items[2])
                specie = f'{species[0].title()}_{species[1].lower()}'
                if region_name not in data:
                    data[region_name] = []
                record.id = allele_name + '|' + specie
                record.description = ''
                data[region_name].append(record)
            # export
            for region_name, records in data.items():
                outfile = self.region_fasta_file(postfix, region_name)
                with open(outfile, 'w') as f:
                    SeqIO.write(records, f, 'fasta')
                if self.verbose:
                    print(outfile)

    def build_igblastp_region(self, bin_dir=None):
        '''
        build igblastp database by V/D/J/C regions
        '''
        # define outdir
        ref_type = 'AA-WithoutGaps-F+ORF+inframeP'
        outdir = os.path.join(self.data_dir, ref_type, 'region_igblastdb')
        p = ProcessIgblast(outdir, self.verbose)

        # build db by region
        regions = ['V-REGION', 'D-REGION', 'J-REGION', 'C-REGION']
        meta = {r:[] for r in regions}
        for region in regions:
            fa_file = self.region_fasta_file(ref_type, region)
            # build database of igblastp
            result = p.build_prot_db(fa_file, bin_dir)
            meta[region].append(result)
        return meta

    def organism_region_fasta_file(self, ref_type, organism, region_name):
        outdir = os.path.join(self.data_dir, ref_type, 'organism', organism)
        Dir(outdir).init_dir()
        return os.path.join(outdir, f"{region_name}.fasta")

    def build_organism_fasta(self):
        '''
        The IMGT/GENE-DB FASTA header contains 15 fields separated by '|':
        1. IMGT/LIGM-DB accession number(s)
        2. IMGT gene and allele name
        3. species (may be followed by an "_" and the name of the strain, breed or isolate, if defined)
        4. IMGT gene and allele functionality
        5. exon(s), region name(s), or extracted label(s)
        6. start and end positions in the IMGT/LIGM-DB accession number(s)
        7. number of nucleotides in the IMGT/LIGM-DB accession number(s)
        8. codon start, or 'NR' (not relevant) for non coding labels
        9. +n: number of nucleotides (nt) added in 5' compared to the corresponding label extracted from IMGT/LIGM-DB
        10. +n or -n: number of nucleotides (nt) added or removed in 3' compared to the corresponding label extracted from IMGT/LIGM-DB
        11. +n, -n, and/or nS: number of added, deleted, and/or substituted nucleotides to correct sequencing errors, or 'not corrected' if non corrected sequencing errors
        12. number of amino acids (AA): this field indicates that the sequence is in amino acids
        13. number of characters in the sequence: nt (or AA)+IMGT gaps=total
        14. partial (if it is)
        15. reverse complementary (if it is)
        '''
        for postfix, parser in self.parse_fasta_file():
            data1 = {}
            for record in parser:
                items = record.description.split('|')
                allele_name, region_name = items[1], items[4]
                species = re.split(r'[ |_]', items[2])
                specie = f'{species[0].title()}_{species[1].lower()}'
                if specie not in data1:
                    data1[specie] = {}
                if region_name not in data1[specie]:
                    data1[specie][region_name] = []
                record.id = allele_name
                record.description = ''
                data1[specie][region_name].append(record)
            # export
            for specie, data2 in data1.items():
                for region, records in data2.items():
                    outfile = self.organism_region_fasta_file(postfix, specie, region)
                    with open(outfile, 'w') as f:
                        SeqIO.write(records, f, 'fasta')

    def build_igblastp_organism(self, bin_dir=None):
        ref_type = 'AA-WithoutGaps-F+ORF+inframeP'
        species = ('Mus_musculus', 'Macaca_mulatta', 'Homo_sapiens',
            'Rattus_norvegicus', 'Oryctolagus_cuniculus',)
        regions = ['V-REGION', 'D-REGION', 'J-REGION', 'C-REGION']
        meta = {r:[] for r in regions}
        for specie in species:
            outdir = os.path.join(self.data_dir, ref_type, 'organism_igblastdb', specie)
            p = ProcessIgblast(outdir, self.verbose)
            for region in regions:
                fa_file = self.organism_region_fasta_file(ref_type, specie, region)
                result = p.build_prot_db(fa_file, bin_dir)
                meta[region].append(result)
        return meta

    def genedb_region(self, ref_type, specie, region):
        v = []
        infile = os.path.join(self.data_dir, ref_type, specie, f'{region}.fasta')
        with open(infile, 'r') as f:
            parser = SeqIO.parse(f, 'fasta')
            for record in parser:
                _id = record.id
                seq = str(record.seq)
                v.append({
                    'isotype': re.findall(r'^(\w*)V', _id)[0],
                    'imgt_db': re.findall(r'(.*)\*', _id)[0],
                    'gene_name': _id,
                    'region': re.findall(r'V\d*', _id)[0],
                    'allele': re.findall(r'\*(\d*)', _id)[0],
                    'seq': seq,
                    'len': len(seq),
                })
        df = pd.DataFrame(v)
        return df
