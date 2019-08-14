from __future__ import print_function
from dnapilib.apred import adapter_prediction

import argparse

def auto_detect(fasta_file):
    adapter = adapter_prediction(fasta_file, 1.4, 9, 50000)
    if adapter:
        if adapter[0][1] == "TGGAATTCTCGGGTGCCAAGG":
            return "illumina"
        elif adapter[0][1] == "AACTGTAGGCACCATCAAT":
            return "qiaseq"
        elif adapter[0][1] == "AAAAAAAAA":
            return "cats"
    return protocol


parser = argparse.ArgumentParser()
parser.add_argument("--fa",
                    help="File with sequences.", required=True)
parser.add_argument("-o", "--out", default="spikeins.fa",
                    help="Name used for output files.")
args = parser.parse_args()

protocol = auto_detect(args.fa)

if protocol == "illumina":
    clip_R1 = 0
    three_prime_clip_R1 = 0
    three_prime_adapter = "TGGAATTCTCGGGTGCCAAGG"
elif protocol == "nextflex":
    clip_R1 = 4
    three_prime_clip_R1 = 4
    three_prime_adapter = "TGGAATTCTCGGGTGCCAAGG"
elif protocol == "qiaseq":
    clip_R1 = 0
    three_prime_clip_R1 = 0
    three_prime_adapter = "AACTGTAGGCACCATCAAT"
elif protocol == "cats":
    clip_R1 = 3
    three_prime_clip_R1 = 0
    # three_prime_adapter = "GATCGGAAGAGCACACGTCTG"
    three_prime_adapter = "AAAAAAAA"