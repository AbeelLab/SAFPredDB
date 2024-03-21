import pickle
import pandas as pd
import input_utils

from tqdm import tqdm

def find_gene_neighbors_df(gene_df, curpos, contig_len, gene, max_dist=5000):
    """ 
    Helper function to quickly calculate the neighbors of a gene in a pandas dataframe
    Parameters
    ----------
    gene_df : pandas dataframe
        A pandas dataframe that contains gene ID and gene position vector (start, end, strand)
    curpos : list
        Current position vector of the seed gene (start, end , strand)
    contig_len : int
        Length of the contig the gene is located on. Needed only if the contig is circular
    gene : str
        Locus tag or ID of the seed gene. Should match the IDs in the gene_df
    max_dist : int
        Maximum distance between the neighbor genes. Default is 5000
    Returns
    -------
    neighbors : list
        All neighbors of the seed gene
    """
    res = gene_df.apply(lambda x: input_utils.calc_intergenic_dist(curpos, x.pos, 
                                                                   contig_len = contig_len), axis=1)
    neighbors = list(res[res < max_dist].index.difference([gene]))
    return neighbors

def find_neighborhood(gene_dict, contig_dict, cluster_dict, genome_dict, max_dist=5000, out_path=None):
    """
    Find the neighbors in a given set of genes after running clustering algorithm 
    Parameters
    ----------
    gene_dict : dict 
        A mapping of gene IDs to its contig, position on the genome and its cluster ID 
        {gene_id: {'cluster_id': cluster_id, 'contig': contig, 'pos': [start_pos, end_pos, strand]}}
    contig_dict : dict
        A mapping of contigs to the contig length and the genes located on the contig
        {contig_id: {'contig_len': contig_len, 'genes': genes}}
    cluster_dict : dict
        A mapping of cluster IDs to the genes within and the genomes these genes are located on
        {cluster_id: {'genes': genes, 'genomes': genomes}}
    genome_dict : dict
        A mapping of genome to the contigs it contains
        {genome: contigs}
    max_dist : int
        Maximum distance (bp) between two genes in a neighborhood. Default: 5000
    out_path : str
        Path to write the neighbors out    
    Returns
    -------
    res_dict : dict
        A dictionary with final results mapping genes to their neighbors, also reports
        all contigs and genomes processed
    """
    tqdm.pandas()
    
    neighbor_dict = dict()
    res_dict = dict()
    genomes_done = []
    contigs_done = []
    
    gene_list = gene_dict.keys()
    for cluster_id, cluster_vals in cluster_dict.items():
        genomes = cluster_vals['genomes']
        num_genomes = len(genomes) 
        for i, genome in enumerate(genomes):
            saveflag = False
            print("{} genomes left in cluster {}".format(num_genomes - i, cluster_id))
            if genome in genomes_done: # This genome was processed before: skip it
                continue
            contigs = genome_dict[genome]
            for contig in contigs:
                try: # To check if there are actually any genes on the contig
                    genes = contig_dict[contig]['genes']
                    genes = set(genes).intersection(gene_list)
                except KeyError:
                    #cluster_dict()print("There are no genes on contig {}!!".format(contig))
                    continue
                contig_len = contig_dict[contig]['contig_len']
                gene_df = pd.DataFrame({'contig_id': {gene: contig for gene in genes}, 
                           'pos': {gene: gene_dict[gene]['pos'] for gene in genes}, 
                           'gene': {gene: gene for gene in genes}})
                neighbors = gene_df.progress_apply(lambda x: find_gene_neighbors_df(gene_df, x.pos, 
                                                                                    contig_len, 
                                                                                    x.name, max_dist), axis=1)
                saveflag = True
                for gene, neighbors in neighbors.iteritems(): 
                    neighbors_passed = [n for n in neighbors if n in gene_list]
                    neighbor_dict[gene] = {'neighbors': neighbors, 'neighbors_passed': neighbors_passed}
                contigs_done.append(contig)
            genomes_done.append(genome)
            if not saveflag: continue # Didn't find any neighbors for any new genomes and/or contigs
            # Save the neighbors into the results dict  
            res_dict['neighbor_dict'] = neighbor_dict
            res_dict['genomes_done'] = genomes_done
            res_dict['contigs_done'] = contigs_done
    if out_path:
        with open(out_path, 'wb') as f: 
            pickle.dump(res_dict, f)

    return res_dict

def find_synteny(neighbor_dict, gene_dict, contig_dict):
    """
    Find the initial syntenic regions, based on the gene neighborhoods collected
    ----------
    neighbor_dict : dict 
        A mapping of gene IDs to its neighbor genes 
        {gene_id: {'neighbors': neighbors, 'neighbors_passed': neighbors_passed}}
    gene_dict : dict 
        A mapping of gene IDs to its contig, position on the genome and its cluster ID 
        {gene_id: {'cluster_id': cluster_id, 'contig': contig, 'pos': [start_pos, end_pos, strand]}}        
    contig_dict : dict
        A mapping of contigs to the contig length and the genes located on the contig
        {contig_id: {'contig_len': contig_len, 'genes': genes}}  
    Returns
    -------
    synteny_dict : dict
        A dictionary with final results listing initial synteny vectors mapped to their ID and
        intergenic distance
    synteny_vec : list
        A list of all synteny vectors generated from the neighborhoods
    intergenic_dist_vec : list
        A list of all intergenic distances within a synteny vector
    """
    num_seed_genes = len(neighbor_dict)
    gene_list = gene_dict.keys()
    synteny_vec = []
    intergenic_dist_vec = []

    for i, (seed_gene, seed_gene_vals) in enumerate(neighbor_dict.items()):
        synteny_genes = []
        add_synteny = []
        add_intergenic_dist = []
        if seed_gene not in gene_list:
            continue
        if i % 1e3 == 0:
            print("Finished {:.2f}% of the seed genes".format(i/num_seed_genes*100))
        neighbors = seed_gene_vals['neighbors_passed']
        # We have neighbors: find the synteny vector
        if len(neighbors) > 0:
            synteny_genes = deepcopy(neighbors)
            synteny_genes.append(seed_gene)
            synteny_genes = sorted(synteny_genes)
            if gene_dict[seed_gene]['pos'][-1] == '-': # Seed gene is on the opposite strand: reverse
                synteny_genes = synteny_genes[::-1]
            add_synteny = [gene_dict[g]['cluster_id'] for g in synteny_genes]
            synteny_vec.append(add_synteny)

            contig_len = contig_dict[gene_dict[seed_gene]['contig']]['contig_len']
            if len(synteny_genes) == 2:
                gene1, gene2 = synteny_genes
                pos1, pos2 = gene_dict[gene1]['pos'], gene_dict[gene2]['pos']
                intergenic_dist_vec.append([input_utils.calc_intergenic_dist(pos1, pos2, contig_len)])
            else: # Calculate all vs all intergenic distance
                add_intergenic_dist = []
                for gene1, gene2 in zip(synteny_genes, synteny_genes[1:]):
                    pos1, pos2 = gene_dict[gene1]['pos'], gene_dict[gene2]['pos']
                    add_intergenic_dist.append(input_utils.calc_intergenic_dist(pos1, pos2, contig_len))
                intergenic_dist_vec.append(add_intergenic_dist)

    # Remove redundant synteny vectors and create a dict to store the results
    synteny_collected = {}
    synteny_count = 0
    for synteny, intergenic_dist in zip(synteny_vec, intergenic_dist_vec):
        synteny_tuple = tuple(synteny)
        if synteny_tuple in synteny_collected.keys():
            synteny_collected[synteny_tuple]['intergenic_dist'].append(intergenic_dist)
        else:
            synteny_collected[synteny_tuple] = {'synteny_id': synteny_count, 'intergenic_dist': [intergenic_dist]}
            synteny_count = synteny_count + 1

    synteny_dict = {v['synteny_id']: {'region': k, 'intergenic_dist': v['intergenic_dist']}
                                      for k, v in synteny_collected.items()}
    
    return synteny_dict, synteny_vec, intergenic_dist_vec

def break_synteny(synteny_dict, max_intergenic_dist=300):
    """
    Break syntenic regions with genes too far apart
    ----------
    synteny_dict : dict
        A dictionary with final results listing initial synteny vectors mapped to their ID and
        intergenic distance
    max_intergenic_dist : int
        Maximum intergenic distance allowed within a syntenic region. Default is 300
    Returns
    -------
    db_df : pandas dataframe
        A dataframe that lists all the final syntenic regions, mapping region ID to gene clusters,
        intergenic distances, and the region length (number of gene clusters in a region)
    """
    synteny_df = pd.DataFrame(synteny_dict.values(), index=synteny_dict.keys())
    synteny_df.loc[:,'region_len'] = synteny_df.region.apply(len)
    synteny_df.loc[:,'min_intergenic_dist'] = synteny_df.intergenic_dist.apply(lambda x: np.array(x).min(axis=0))
    
    num_synteny = len(synteny_df)
    synteny_count = 0
    synteny_dict = {}
    for i, (synteny_id, row) in enumerate(synteny_df.iterrows()):
        if i % 1e3 == 0:
            print("Finished {:.2f}% of the syntenic regions".format(i/num_synteny*100))
        add_synteny = []
        add_intergenic_dist = []
        for j, (c1, c2) in enumerate(zip(row.region, row.region[1:])):
            intergenic_dist = row.min_intergenic_dist[j]
            if intergenic_dist <= max_intergenic_dist:
                if len(add_synteny) == 0:
                    add_synteny.append(c1)
                add_synteny.append(c2)
                add_intergenic_dist.append(intergenic_dist)
            else:
                if len(add_synteny) > 0:
                    synteny_dict[synteny_count] = {'region': add_synteny, 'intergenic_dist': add_intergenic_dist}
                    add_synteny = []
                    add_intergenic_dist = []
                    synteny_count = synteny_count + 1
        if len(add_synteny) > 0:
            synteny_dict[synteny_count] = {'region': add_synteny, 'intergenic_dist': add_intergenic_dist}
            synteny_count = synteny_count + 1

    # Write the final synteny regions into a dataframe
    db_df = pd.DataFrame(synteny_dict.values(), index=synteny_dict.keys())
    db_df.loc[:,'region_len'] = db_df.region.apply(lambda x: len(x))

    return db_df

    

