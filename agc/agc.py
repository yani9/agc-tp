#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter, defaultdict
#import numpy as np
import time
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw 

__author__ = "REN Yani & MIZINIAK Hippolyte"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["REN Yani & MIZINIAK Hippolyte"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "REN Yani & MIZINIAK Hippolyte"
__email__ = "yani.ren@etu.u-paris.fr & hippolyte.miziniak@gmail.com"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication  (default 100)")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()

def read_fasta(amplicon_file, minseqlen):
    """Read the amplicon fasta format file and return a list of sequences where respect
    the minimum sequence length for dereplication condition.
    :Parameters: 
        amplicon_file : str, the path of the amplicon fasta format file.
        minseqlen : int, the minimum sequence length for dereplication
    :Return:
        seq : list, contains sequences in amplicon file which have minimum length imposed by the 'minseqlen' value.
    """
    seq = []
    with gzip.open(amplicon_file, "rt") as  monfich:
        lines=monfich.readlines()
        for i in range(1,len(lines)): 
            if lines[i].startswith(">") or int(i) == len(lines)-1 :
                if len("".join(seq))>=minseqlen:
                    yield "".join(seq) 
                seq = []
            elif not lines[i].startswith(">"):
                seq.append(lines[i].strip())
    

def dereplication_fulllength(amplicon_file, minseqlen, mincount):
    """Read amplicon file and return the most commun seq.
    :Parameters: 
        amplicon_file : str, the path of the amplicon file.
        minseqlen : int, the value of the minimum sequence length.
        mincount : int, the value of the minimum sequence count. 
    :Return: 
        seq_occ : generator of tuple, the sequence and its count.
    """
    start = time.time()

    sequences = list(read_fasta(amplicon_file, minseqlen))
    end = time.time() - start 
    print("time read fasta  ",end)
    print("nombre de seq sortie par read_fast() : ",len(sequences))
 
    seq_set= list(sequences)
    count = Counter()
    for seq in sequences: 
            count[seq]+=1 
    
    count_sorted = count.most_common()
    for seq_occ in count_sorted: 
        #print("SEQ and OCCURENCE ", seq_occ)
        if seq_occ[1]>=mincount:
            yield seq_occ


def get_unique(ids):
    return {}.fromkeys(ids).keys()

def get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size):
    """Create (or update) and return a dictionnary of kmers for every non chimera sequence
    :Parameters: 
        kmer_dict : dict, contains kmers of non chimera sequence or empty dict. 
        sequence : str, non chimeric sequence to cut into kmers and add to kmer_dict.
        id_seq : str, the index of the sequence in the original order in the amplicon file. 
        kmer_size : int, the size of kmers.
    :Return: 
        kmer_dict : dict, updated dict that contains kmers of non chimera sequence.
    """
    kmer_list = list(cut_kmer(sequence, kmer_size))
    for kmer in kmer_list: 
        if not kmer in kmer_dict: 
            kmer_dict[kmer]=[id_seq]  
        else: 
            kmer_dict[kmer].append(id_seq)
    return kmer_dict 

def search_mates(kmer_dict, sequence, kmer_size):
    """Search and return the id of parent sequences for the input sequence.
    :Parameters:
        kmer_dict : dict, contains kmers of non chimera sequence. 
        sequence : str, the sequence for which we search the id of the parents sequences. 
        kmer_size : int, the size of kmers.
    :Return: 
        id_parents : list, contains the id of parent sequences for the input sequence.
    """
    id_seq_list = list()
    kmer_candidat = list(cut_kmer(sequence, kmer_size))
    for kmer in kmer_candidat:
        if kmer in list(kmer_dict.keys()): 
            id_seq_list.extend(kmer_dict[kmer])

    two_most_common = Counter(id_seq_list).most_common(2)
    id_parents = [id[0] for id in two_most_common]
    #print("id parents", id_parents)
    return id_parents

            

def detect_chimera(perc_identity_matrix):
    """Analyze the percentage identity matrix : compute the standard deviation value of each pair of percentage identity of alignment
    between each segment of the sequence and the segment of theirs two parents. 
    If the mean of all standard deviation value of percentage identity is superior to 5 and if there are difference of similarity of sequence with the two parents, 
    we consider the sequence as non chimeric. 
    Return a boolean to indidates if the sequence is chimeric or not. 
    :Parameter: 
        perc_identity_matrix : array of array, each segment contains two percentage identity (with the two parents). 
    """
    all_stdev = [statistics.stdev(seg_array_id) for seg_array_id in perc_identity_matrix]
    simil_0 = 0
    simil_1 = 0
    # if the mean of standard deviation value of identity percentage is superior to 5 and ..
    if statistics.mean(all_stdev)>5:
        for seg_array_id in perc_identity_matrix: 
            if seg_array_id[0]>seg_array_id[1]: 
                simil_0+=1
            elif seg_array_id[0]<seg_array_id[1]:
                simil_1+=1    
        # if we detect at leat one segment with higher similarity to sequence 1 and one segment having a higher similarity to sequence 2   
        if (simil_0>=1 and simil_1>=1): 
            return True # the seq is chimeric
        else: 
            return False  
    else: 
        return False    
    

def common(lst1, lst2): 
    """Return common elements between two lists."""
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """Divides the sequence into segment and return a list of not overlapping segment (chunk).
    :Parameters:
        sequence : str, sequence to divide into segment.
        chunk_size : int, the size of the segment.
    """
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size] 
              for i in range(0, len_seq, chunk_size) 
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers. Return all kmers of input sequence."""
    for i in range(0, len(sequence) - kmer_size + 1):
        yield sequence[i:i+kmer_size]

def get_identity(alignment_list):
    """Prend en une liste de séquences alignées au format ["SE-QUENCE1", "SE-QUENCE2"]
    Retourne le pourcentage d'identite entre les deux."""
    id_nu = 0
    for i in range(len(alignment_list[0])):
        if alignment_list[0][i] == alignment_list[1][i]:
            id_nu += 1
    return round(100.0 * id_nu / len(alignment_list[0]), 2)


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size): 
    """Read amplicon file and filter dereplication and chimeric sequence and return non chimeric sequences.
    :Parameters: 
        amplicon_file : str, the path of the amplicon file.
        minseqlen : int, the value of the minimum sequence length.
        mincount : int, the value of the minimum sequence count. 
        chunk_size : int, the size of segments.
        kmer_size : int, the size of kmers.
    :Return: 
        seq : generator of list, non chimeric sequences and its count.
    """
    start = time.time()
    seq_list = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    print(seq_list)
    print("nombre seq sortie par dereplication() ", len(seq_list))
    end = time.time() - start
    print("time dereplication ", end)
    seq_nonchim = seq_list[:2]
    kmer_dict = defaultdict(list)
    for id_seq, sequence in enumerate(seq_nonchim): 
        kmer_dict.update(get_unique_kmer(kmer_dict, sequence, id_seq, kmer_size))
    for i in range(2, len(seq_list)):
        id_parents = search_mates(kmer_dict, seq_list[i][0], kmer_size)
        seq_parents = [seq_list[id][0] for id in id_parents]
        if len(seq_parents)<2: 
            # The candidat sequence is considered non chimeric if it has less than 2 parent sequences
            seq_nonchim.append(seq_list[i][0])
            for chunk in get_chunks(seq_list[i][0], chunk_size):
                kmer_dict.update(get_unique_kmer(kmer_dict, chunk, i, kmer_size))
                
        else : 
            chunk_candidat = list(get_chunks(seq_list[i][0], chunk_size))
            chunk_parents = []
            for seq in seq_parents: 
                chunk_parents.append(get_chunks(seq, chunk_size))

            perc_identity_matrix = [[] for c in range(4)]         

            # Compute the identity percentage alignment between every 4 segments of the 2 parent sequences and candidat sequence.
            for k in range(2): 
                for j in range(4):
                    alignment_list = nw.global_align(chunk_candidat[j],chunk_parents[k][j], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))

                    identity = get_identity(alignment_list)   
                    perc_identity_matrix[j].append(identity)

            #print("matrice ", perc_identity_matrix)
            ischimera = detect_chimera(perc_identity_matrix)
            #print(ischimera)
            if not ischimera: 
                seq_nonchim.append(seq_list[i])
                kmer_dict.update(get_unique_kmer(kmer_dict, seq_list[i][0], id_seq, kmer_size))
    for seq in seq_nonchim:  
        yield seq 

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """Check for each non chimeric sequence if it is an OTU.
    :Parameters: 
        amplicon_file : str, the path of the amplicon file.
        minseqlen : int, the value of the minimum sequence length.
        mincount : int, the value of the minimum sequence count. 
        chunk_size : int, the size of segments.
        kmer_size : int, the size of kmers.
    :Return: 
        otu_list : list, contains all sequences detected as an OTU.
    """
   
    start = time.time()
    seq_nonchim = list(chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size))
    end = time.time() - start 
    print("time chimera ", end)
    print("CHIM",seq_nonchim)
    otu_list = []
    otu_list.append(seq_nonchim[0])

    for i in range(1,len(seq_nonchim)): 
        flag = True
        for j in range(len(otu_list)): 
            alignment_list = nw.global_align(seq_nonchim[i][0], otu_list[j][0], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            identity = get_identity(alignment_list)
            #print("identity value ", identity) 
            if identity > 97:
                flag = False
                break
        if flag:
           otu_list.append(seq_nonchim[i])
    
    return otu_list

  

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    """Write an output file with all detected OTU sequences, in fasta format file.
    :Parameters: 
        OTU_list : list, contains OTU sequences.
        output_file : str, the path of the output file.
    """
    with open(output_file, "wt") as filout: 
        for i, otu in enumerate(OTU_list): 
            filout.write(">OTU_{} occurrence:{}\n".format(i+1, otu[1])) 
            seq = otu[0]+"\n" 
            filout.write(str(fill(seq))) 


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici
    OTU_list = abundance_greedy_clustering(args.amplicon_file, args.minseqlen, args.mincount, args.chunk_size, args.kmer_size)
    print("Number of OTU : {}, the list OTU sequences are write in {} file ", len(OTU_list), args.output_file)
    write_OTU(OTU_list, args.output_file)

if __name__ == '__main__':
    main()


