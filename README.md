# SAFPredDB

Bacterial synteny database

SAFPredDB is a bacterial synteny database built for the gene function prediction tool SAFPred, Synteny Aware Function Predictor. The database is a collection of conserved synteny and operons found across the bacterial kingdom. First, we formulated a synteny model based on experimentally known operons and the genomic features common in bacteria. We designed a bottoms-up, purely computational approach to build our database based on the proposed synteny model using complete bacterial genome assemblies from the Genome Taxonomy Database (GTDB).

Although we initially built SAFPred for our prediction tool only, it can be used for other purposes where such catalog is needed. As a standalone database, it can be queried to mine information about conserved genomic patterns in bacteria. In this repository, we provide the Python scripts we used to build SAFPredDB, so it can be reproduced or updated as newer assemblies are added to GTDB. In addition, you can run the scripts on any given genomic data to build your own database, tailored to a specific task downstream.

## Installation

You can create a new conda environment and install all dependencies with the environment file in this repository.

```bash
conda env create -f environment.yml
```

## Dependencies

- biopython

- pandas

- numpy

## Building SAFPredDB - Synteny database

As input, you will need the following 4 mapping dictionaries (stored as pickle files) in a single directory:

- `gene_dict`: Mapping genes to their clusters (if the genes were clustered to reduce redundancy before building the database), positions on the genome and the contig they are located on
  
  `{gene_id: {'cluster_id': cluster_id, 'contig': contig, 'pos': [start_pos, end_pos, strand]}}`

- `contig_dict`: Mapping contigs to the contig length and the genes located on the contig
  
  `{contig_id: {'contig_len': contig_len, 'genes': genes}}`

- `cluster_dict`: Mapping cluster IDs to the genes within the cluster and the genomes these genes are located on
  
  `{cluster_id: {'genes': genes, 'genomes': genomes}}`

- `genome_dict`: Mapping genome to the contigs it contains
  
  `{genome: contigs}`

You can run the python script `build_safpreddb.py` as follows

```bash
python safpreddb.py -i data -o out --max_dist 5000 --max_intergenic_dist 300
```

This command will build a synteny database and write the database contents as a pandas dataframe into the directory `out` in a compressed pickle file titled `safpreddb.pkl.gz`

```bash
python safpreddb.py -h 
usage: safpreddb.py [-h] --input_dir INPUT_DIR --output_dir OUTPUT_DIR [--max_dist MAX_DIST]
                    [--max_intergenic_dist MAX_INTERGENIC_DIST]

build SAFPredDB

optional arguments:
  -h, --help            show this help message and exit
  --input_dir INPUT_DIR, -i INPUT_DIR
                        Directory containing the input mapping files: gene_dict, contig_dict,
                        cluster_dict, genome_dict
  --output_dir OUTPUT_DIR, -o OUTPUT_DIR
                        Output directory where the database will be stored
  --max_dist MAX_DIST, -r MAX_DIST
                        Maximum distance allowed within a syntenic region. Default is 5000
  --max_intergenic_dist MAX_INTERGENIC_DIST, -g MAX_INTERGENIC_DIST
                        Maximum intergenic distance allowed within a syntenic region. Default is
                        300
```

## Edit SAFPredDB

Additionally, you can edit an existing database if you want to remove any clusters and the syntenic regions they are associated with. We used this script to generate the additional databases for our benchmark study where we reduced the maximum sequence similarity between the training and the test set. To have a database also consistent with the sequence similarity threshold, we removed clusters from the existing database.

```python
from database_utils import edit_database

outdb_path = 'out/safpreddb_edited.pkl.gz'
indb_path = 'out/safpreddb.pkl.gz'
# If you want to add back the singleton regions (regions with only one gene) to the database
add_singletons = True

with open('data/keep_clusters.txt', 'r') as f:
    keep_clusters = [int(line.atrip()) for line in f]

out_db_df = edit_database(indb_path, outdb_path, keep_clusters, max_intergenic_dist=300,
                          add_singletons=add_singletons)
```

## Citation

If you find our method or any of the original script in this repository useful, please cite the manuscript:

```python
@article {urhan2024safpred,
    author = {Aysun Urhan and Bianca-Maria Cosma and Ashlee M. Earl and Abigail L. Manson and Thomas Abeel},
    title = {SAFPred: Synteny-aware gene function prediction for bacteria using protein embeddings},
    year = {2024},
    volume = {40},
    number = {6},
    pages = {btae328},
    month = {05},
    doi = {10.1093/bioinformatics/btae328},
    journal = {Bioinformatics}
}
```

If you use the full SAFPredDB database from the [4TU.ResearchData](https://doi.org/10.4121/ac84802e-853f-46f1-9786-b9d29c0f7557.v1) please cite

```python
@misc{urhan2024safpreddb,
  doi = {10.4121/AC84802E-853F-46F1-9786-B9D29C0F7557},
  url = {https://data.4tu.nl/datasets/ac84802e-853f-46f1-9786-b9d29c0f7557},
  author = {Urhan, Aysun and Cosma, Bianca-Maria and Earl, Ashlee M. and Manson, Abigail L. and Abeel, Thomas},
  keywords = {Microbiology, FOS: Biological sciences, Genetics, Biological Sciences, bionformatics, microbial genomics, genomics, protein language model, bacterial genomics, comparative genomics, protein embeddings, sequence analysis, bacterial synteny},
  language = {en},
  title = {SAFPredDB: Bacterial synteny database},
  publisher = {4TU.ResearchData},
  year = {2024},
  copyright = {Creative Commons Attribution Non Commercial 4.0 International}
}
```
