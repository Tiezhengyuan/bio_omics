# bio_omics
Download, retrieve and process omics data, or biological informatics data from public database

Comprehensive Databases
- NCBI: genome database, 
- UniProt: protein database, https://www.expasy.org/resources/uniprotkb-swiss-prot

Sepecific Databases
- miRBase: mircoRNA database, https://www.mirbase.org/
- RNACentral: non-coding RNA squence database, https://rnacentral.org/
- IEDB: immune epitope database, https://www.iedb.org/


See the help documents of example coding at https://www.fbridges.com/pipeline/bio_omics.

https://www.iedb.org/downloader.php?file_name=doc/epitope_full_v3.zip


## data model
ETL data processing is composed of some steps including downloads, retrieval, organization, combination, integration, enrichment, formation. This packages focus on downloads, retrieval, and combination of omics data.
It is suggested that the data model would be consistent. Data are organized by entity namely protein, or antigen. An example of data is showed as the below. Here the pair 'key' defines unique identifier of this entity. "ID" is automatically created. Retrieved data are pushed as one key-value. 
Note:
- Abundant data are possible and to be allowed.
- The key-value is defined by this corresponding database source.
- Used for Integration rather than enrichment. Therefore, data combination or aggregation is not recommended.
- Data from various source could be different or invalid. Those would be validated in the afterwards step rather than this step.
```
{
    "ID": "79541",
    "key": "H0YED9",
    "Uni-Prot": {
        ....
    },
    "NCBI": {
        ....
    },
    "PDB": {
        ....
    },
}
```