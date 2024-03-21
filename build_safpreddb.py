import argparse
import os
import pickle
from build_database import find_neighborhood, find_synteny, break_synteny

def main():
    parser = argparse.ArgumentParser(description = 'build SAFPredDB')
    parser.add_argument('--input_dir', '-i', type=str, required=True, 
                        help = 'Directory containing the input mapping files: gene_dict, contig_dict, cluster_dict, genome_dict')
    parser.add_argument('--output_dir', '-o', type=str, required=True, 
                        help='Output directory where the database will be stored')
    parser.add_argument('--max_dist', '-r', type=int, default = 300, 
                        help='Maximum distance allowed within a syntenic region. Default is 5000')
    parser.add_argument('--max_intergenic_dist', '-g', type=int, default=300, 
                        help='Maximum intergenic distance allowed within a syntenic region. Default is 300')

    args = parser.parse_args()
    args.out = args.out[0]
    print(args.out)
    input_dir = args.input_dir
    output_dir = args.output_dir
    max_dist = args.max_dist
    max_intergenic_dist = args.max_intergenic_dist

    print("Loading input data")
    with open(os.path.join(input_dir, 'gene_dict.pkl'), 'rb') as f:
        gene_dict = pickle.load(f)
    with open(os.path.join(input_dir, 'contig_dict.pkl'), 'rb') as f:
        contig_dict = pickle.load(f)    
    with open(os.path.join(input_dir, 'cluster_dict.pkl'), 'rb') as f:
        cluster_dict = pickle.load(f)
    with open(os.path.join(input_dir, 'genome_dict.pkl'), 'rb') as f:
        genome_dict = pickle.load(f) 

    print("Finding gene neighboods")   
    neighbor_res_dict = find_neighborhood(gene_dict, contig_dict, cluster_dict, genome_dict, max_dist=max_dist)
    neighbor_dict = neighbor_res_dict['neighbor_dict']
    print("Finding itial syntenic regions")
    synteny_dict, _, _ = find_synteny(neighbor_dict, gene_dict, contig_dict)
    print("Finalizing the syntenic regions -- breaking regions that are too far apart")
    db_df = break_synteny(synteny_dict, max_intergenic_dist=max_intergenic_dist)

    output_file = os.path.join(output_dir, 'safpreddb.pkl.gz')
    db_df.to_pickle(output_file, compression='gzip')

if __name__ == '__main__':
    main()
