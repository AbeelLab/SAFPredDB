import pickle
import pandas as pd
import input_utils


def find_gene_neighbors_df(gene_df, curpos, contig_len, gene, max_dist=5000):
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
        A mapping of gene IDs to its position on the genome and its cluster ID 
        {gene_id: {'cluster_id': cluster_id, 'pos': [start_pos, end_pos, strand]}}
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
                except KeyError:
                    print("There are no genes on contig {}!!".format(contig))
                    continue
                contig_len = contig_dict[contig]['contig_len']
                gene_df = pd.DataFrame({'contig_id': {gene: contig for gene in genes}, 
                           'pos': {gene: gene_dict[gene] for gene in genes}, 
                           'gene': {gene: gene for gene in genes}})
                neighbors = gene_df.progress_apply(lambda x: find_gene_neighbors_df(genes, x.pos, 
                                                                                    contig_len, 
                                                                                    x.name, max_dist), axis=1)
                saveflag = True
                for gene, neighbors in neighbors.iteritems(): 
                    neighbors_passed = [n for n in neighbors if n in gene_list]
                    neighbor_dict[idx] = {'neighbors': neighbors, 'neighbors_passed': neighbors_passed}
                contigs_done.append(contig)
            genomes_done.append(genome)
            if not saveflag: continue # Didn't find any neighbors for any new genomes and/or contigs
            # Save the neighbors into the results dict  
            res_dict['neighbors_dict'] = neighbors_dict
            res_dict['genomes_done'] = genomes_done
            res_dict['contigs_done'] = contigs_done
  
      if out_path:
          with open(out_path, 'wb') as f: 
              pickle.dump(res_dict, f)
