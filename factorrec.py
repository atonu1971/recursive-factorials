#! /usr/bin/env python
#import required packages & modules

from Bio.Seq import Seq
from Bio.Seq import translate
import re
import sys

#defining recursive factorial 


        
class Sort:
    def __init__(self, ID, Date, Phenotype, Sequence):
        self.ID = ID
        self.date = Date
        self.pheno = Phenotype
        self.seq = Sequence


def factorrec(n):
    if n == 1:
        return 1
    else:
        return n * factorrec(n-1)
        
     


def results(fasta_file, p, output, ln, n, blue, orange):
    """
This script will produce result.txt file which will contain the frequency of orange phenotype in specific population by using recursive function.    
    Input of the argument will contain the FASTA ID (fasta_file), input for orange population frequency (p),list for blue and orange phenotype 
   output contains a txt file determining the frequency of orange phenotype in a population.
    n=total number of individuals in sample set
   k=list of all individuals with the "orange" trait
   [f(k;n,p)]= the probability of finding your number of observed orange traits given the sample size and the frequency of orange within the population 
    Output Samples will be:
     Results
    
        p (the frequency of "orange" in the population) = 0.3
        n (the number of sampled individuals) = 100
        k (the number of "orange" individuals in the sample set) = 26
        Probability of collecting 32 individuals with 5 being "orange" (given a population frequency of 0.3) = 0.06126913528290234

    """

    print('Opening' + fasta_file)
    with open(fasta_file, 'r') as in_stream:
        print('opening' + output)
        with open(output, 'w') as out_stream:


            for line in in_stream:
                ln = ln + 1

                if  re.match('>',line):
                    n = n + 1
                    ls = re.split('_|\ |>', line)
                    ID = ls[1]
                    Date = ls[2]


                else:
                    Sequence = line
                    ps = translate(line)
                    if ps[3] == 'R':
                        Phenotype = "orange"
                        orange.append(ID)
                    if ps[3] == 'S':
                        Phenotype = "blue"
                        blue.append(ID)
                    globals()[f"{ID}_Aubie"] = Sort(ID, Date, Phenotype, Sequence)

            q = 1-p
            k = len(orange)
            bern = (factorrec(n)/(factorrec(n-k)*factorrec(k)))*((p**k)*(q**(n-k)))

        #closing the opened files from earlier
            out_stream.write("Results\n\np (the frequency of \"orange\" in the population) = " + str(p))
            out_stream.write("\nn (the number of sampled individuals) = " + str(n))
            out_stream.write("\nk (the number of \"orange\" individuals in the sample set) = " + str(k))
            out_stream.write("\n\nProbability of collecting 32 individuals with 5 being \"orange\" (given a population frequency of 0.3) = " + str(bern) + "\n")

    print("Done!")
    print(sys.argv[1] + ' is closed?', in_stream.closed)
    print(sys.argv[3] + ' is closed?', out_stream.closed)

if __name__ == '__main__':
    fasta_file = (sys.argv[1]) 
    p = float(sys.argv[2])
    output = (sys.argv[3])
    ln = 0
    n = 0
    blue = []
    orange = []
    results(fasta_file, p, output, ln, n, blue, orange)
