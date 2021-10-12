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
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
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
"""
def read_fasta(amplicon_file, minseqlen):
    with gzip.open(amplicon_file, "rt") as  monfich:
        lines=monfich.readlines()
        for i in range(2,len(lines),3):
            if len(lines[i])>=minseqlen:
                yield lines[i].strip()
""" 
def read_fasta(amplicon_file, minseqlen):
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
    """Read amplicon file and return the most commun seq."""
    sequences = list(read_fasta(amplicon_file, minseqlen))
    count = Counter()
    print("ICI ",len(sequences))
    for seq in sequences: 
            count[seq]+=1 
    
    count_sorted = count.most_common()
    for seq_occ in count_sorted: 
        print("OCCURENCE ", seq_occ, seq_occ[1])
        if seq_occ[1]>=mincount:
            yield seq_occ


def get_unique(ids):
    return {}.fromkeys(ids).keys()


def common(lst1, lst2): 
    return list(set(lst1) & set(lst2))


def get_chunks(sequence, chunk_size):
    """"""
    len_seq = len(sequence)
    if len_seq < chunk_size * 4:
        raise ValueError("Sequence length ({}) is too short to be splitted in 4"
                         " chunk of size {}".format(len_seq, chunk_size))
    return [sequence[i:i+chunk_size] 
              for i in range(0, len_seq, chunk_size) 
                if i+chunk_size <= len_seq - 1]


def cut_kmer(sequence, kmer_size):
    """Cut sequence into kmers"""
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
    pass

def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    seq_list = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    otu_list = []
    otu_list.append(seq_list[0])
    align = nw 
    for i in range(len(seq_list)): 
        for j in range(len(otu_list)): 
            alignment_list = nw.global_align(seq_list[i][0], otu_list[j][0], gap_open=-1, gap_extend=-1, matrix=os.path.abspath(os.path.join(os.path.dirname(__file__),"MATCH")))
            identity = get_identity(alignment_list)
            if identity <=97:
                otu_list.append([seq_list[i][0], seq_list[i][1]])
            elif identity > 97 : 
                if otu_list[j][1]<seq_list[i][1]:
                    otu_list[j]
                
    return(otu_list)

  

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))

def write_OTU(OTU_list, output_file):
    pass

#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    #args = get_arguments()
    # Votre programme ici

    amplicon_file = "./tests/test_sequences.fasta.gz"
    minseqlen = 0
    mincount = 0
    chunk_size = 100
    kmer_size = 8
    seq = read_fasta(amplicon_file, minseqlen)
    print(next(seq))
    print(next(seq))
    print(next(seq))

    
    #seq_info = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    #print("SEQ INFO ",seq_info)
    seq_list = dereplication_fulllength(amplicon_file, minseqlen, mincount)
    #print(next(seq_list))

    otu_list = abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size)
    print("OTU ", otu_list)
    #print(next(seq_info))
    #print(next(seq_info))
    
if __name__ == '__main__':
    main()
