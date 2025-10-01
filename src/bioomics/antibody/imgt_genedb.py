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
from copy import deepcopy
import re
import os
import numpy as np
import pandas as pd
import requests
import subprocess
import sys

from .imgt import Imgt
from .process_data import ProcessData

class ImgtGenedb(Imgt):

    def __init__(self, data_dir:str):
        super().__init__(data_dir)
        self.data_dir = os.path.join(self.data_dir, 'GENE-DB')
        # keys: specie, gene_name, region_name
        self.data = {}
    
    def __call__(self):
        # gene list
        df = self.parse_genelist(True)
        g = df.groupby(['Species', 'IMGT/GENE-DB'])
        for (specie, gene_name), sub in g:
            if specie not in self.data:
                self.data[specie] = {}
            if gene_name not in self.data[specie]:
                self.data[specie][gene_name] = {}
            sub = sub.to_dict(orient='records')
            self.data[specie][gene_name]['genelist'] = sub
        
        # read fasta and integrate sequences
        self.parse_fasta()

        # save data
        outdir = os.path.join(self.data_dir, 'gene')
        self.init_dir(outdir)
        for specie, data1 in self.data.items():
            json_file = os.path.join(outdir, f"{specie}.json")
            self.save_json(data1, json_file)
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
            # print(len(header), header)
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
                return pd.Series(['tcr', np.nan])
            return pd.Series(['other', np.nan])
        
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

    def build_specie_fasta(self):
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
                outdir = os.path.join(self.data_dir, postfix, specie)
                self.init_dir(outdir)
                for region, records in data2.items():
                    outfile = os.path.join(outdir, f"{region}.fasta")
                    with open(outfile, 'w') as f:
                        SeqIO.write(records, f, 'fasta')
        
    def build_genedb_igblastdb(self):
        bin_dir = '/home/yuan/data/ncbi-igblast-1.22.0/bin'
        if bin_dir not in sys.path:
            sys.path.append(bin_dir)

        
        ref_types = ('AA-WithoutGaps-F+ORF+inframeP',)
        species = ('Mus_musculus', 'Macaca_mulatta', 'Homo_sapiens',
            'Rattus_norvegicus', 'Oryctolagus_cuniculus',)
        regions = ['V-REGION', 'D-REGION', 'J-REGION', 'C-REGION']
        meta = {r:[] for r in regions}
        for ref_type in ref_types:
            for specie in species:
                for region in regions:
                    outdir = os.path.join(self.data_dir, ref_type, specie)
                    fasta_file = os.path.join(outdir, f'{region}.fasta')
                    try:
                        outprefix = os.path.join(outdir, region)
                        cmd = [
                            f'{bin_dir}/makeblastdb',
                            '-parse_seqids -dbtype prot',
                            f'-in {fasta_file}',
                            f'-out {outprefix}',
                        ]
                        cmd = ' '.join(cmd)
                        result = subprocess.run(cmd, capture_output=True, \
                            shell=True, text=True, check=True)
                        # print(result.stdout)
                        meta[region].append(outprefix)
                    except Exception as e:
                        pass
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
