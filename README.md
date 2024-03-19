# SAFPredDB
Bacterial synteny database

SAFPredDB is a bacterial synteny database built for the gene function prediction tool SAFPred, Synteny Aware Function Predictor. The database is a collection of conserved synteny and operons found across the bacterial kingdom. First, we formulated a synteny model based on experimentally known operons and the genomic features common in bacteria. We designed a bottoms-up, purely computational approach to build our database based on the proposed synteny model using complete bacterial genome assemblies from the Genome Taxonomy Database (GTDB).

Although we initially built SAFPred for our prediction tool only, it can be used for other purposes where such catalog is needed. As a standalone database, it can be queried to mine information about conserved genomic patterns in bacteria. In this repository, we provide the Python scripts we used to build SAFPredDB, so it can be reproduced or updated as newer assemblies are added to GTDB. In addition, you can run the scripts on any given genomic data to build your own database, tailored to a specific task downstream.
