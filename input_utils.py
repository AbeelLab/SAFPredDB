#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Utility functions to generate and preprocess input to build synteny database

import os
import numpy as np

def calc_ham_dist(seq1, seq2, to_ref=True):
    """
    Calculate Hamming distance between 2 sequences of varying lengths, trim if larger
    Input: 2 sequences, flag to use the first sequence as reference
    Returns: Hamming distance
    """
    d = sum(s1 != s2 for (s1, s2) in zip(seq1, seq2))
    if to_ref: d = d/len(seq1)
    return d
    
def calc_intergenic_dist(pos1, pos2, contig_len=None):
    """
    Calculate intergenic distance between two genes on the same contig, also 
    takes strandedness into account
    Input: positions of 2 genes (as intervals), circular flag (default: True)
    and contig length (necessary only for circular contigs)
    Returns the intergenic distance - or None if the genes are on different strands
    ref: Mihelčić, M., Šmuc, T. & Supek, F. Patterns of diverse gene functions 
    in genomic neighborhoods predict gene function and phenotype. 
    Sci Rep 9, 19537 (2019). https://doi.org/10.1038/s41598-019-55984-0
    """
    if pos1[-1] != pos2[-1]:
        return None
    if pos1[0] > pos2[1]:
        corpos1 = pos2
        corpos2 = pos1
    elif pos2[0] > pos1[1]:
        corpos1 = pos1
        corpos2 = pos2
    else:
        return 0
    if contig_len:
        d = min(abs(corpos2[0] - corpos1[1]), abs(corpos2[1] - contig_len - corpos1[0]))
    else:
        d = abs(corpos2[0] - corpos1[1])

    return d

def extract_cdhit_clusters(cluster_file_path, rep_file_path, out_path=None):
    """
    Helper function to preprocess and clean CD-HIT output files
    Input: cluster file output from CD-HIT run with extention .clstr, cluster
    representatives file output from CD-HIT run (should have no extensions),
    output file to write the resulting clusters in a human readable form (optional)
    Returns a DataFrame object to describe the clusters
    """
    import pandas as pd
    from Bio import SeqIO
  
    rep_list = []
    prot2cluster = {}
    for i, rec in enumerate(SeqIO.parse(rep_file_path, 'fasta')):
        rep_id = rec.id
        rep_seq = str(rec.seq)
        rep_len = len(rec.seq)
        rep_list.append({'rep_id': rep_id, 'rep_seq': rep_seq, 'rep_len': rep_len})
        prot2cluster[rep_id] = i
    
    cluster_list = []
    prot_list = []
    len_list = []   
    with open(cluster_file_path, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('>Cluster'):
                c = int(line.strip().split()[-1])
                cluster_list.append(c)
                if c != 0: 
                    prot_list.append(cur_prot)
                    len_list.append(cur_len)
                cur_len = []
                cur_prot = []
            else:
               f1, f2 = line.strip().split('\t')[-1].split(', ')
               prot_len = int(f1.replace('aa',''))
               prot_id = f2.split('...')[0].replace('>','')
               cur_len.append(prot_len)
               cur_prot.append(prot_id)
    prot_list.append(cur_prot)
    len_list.append(cur_len)
    
    cluster_data = [{'members': m, 'seq_len': s} for m, s in zip(prot_list, len_list)]
    cluster_df = pd.DataFrame(cluster_data, index = cluster_list)
    cluster_dict = {c: {} for c in cluster_df.index}
    for i, rep in enumerate(rep_list):
        if i % 1e5 == 0: print("Porocessing representative {}".format(i))
        c = prot2cluster[rep['rep_id']]
        cluster_dict[c] = {'rep_id': rep['rep_id'], 'rep_seq': rep['rep_seq'], 'rep_seqlen': rep['rep_len']}
    rep_df = pd.DataFrame(cluster_dict.values(), index=cluster_dict.keys())    
    cluster_df = pd.concat([cluster_df, rep_df], axis=1)
    
    if out_path: 
      cluster_df.to_pickle(out_path, compression='gzip')
    return cluster_df
