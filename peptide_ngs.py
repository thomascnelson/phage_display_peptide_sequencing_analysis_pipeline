"""
Peptide NGS backend data processing
Created Dec 2020
@author: Tom Nelson
"""
import sys
import pandas as pd
import re
import regex
import gzip
from collections import Counter
import pickle
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
from Bio.SeqUtils.ProtParam import ProteinAnalysis

#### Change the location of these files and directories when installing ####
MIN_LEN = 3 # minimum allowable peptide length, can be overridden
CUSTOM_CODON_TABLE_FILE = 'my_custom_codon_tables.pkl' # can be overridden
DATA_DIR = '/path/to/sequencing/data/'
OUT_DIR = '/path/to/output/'

# Functions

def load_custom_codon_table(name, table_file = CUSTOM_CODON_TABLE_FILE):
    '''Function to load a custom codon table from pickle file database.
    Takes the codon table name and the pickle file name as arguments.
    Pickle file name defaults to CUSTOM_CODON_TABLE_FILE but can be overridden.
    Returns a Biopython codon table object.'''
    tables = []
    with (open(table_file, "rb")) as openfile:
        while True:
            try:
                tables.append(pickle.load(openfile))
            except EOFError:
                break
    for codon_table in tables:
        if str(codon_table.names[0]) == name:
            my_table = codon_table
            break
        else:
            my_table = None
    if my_table:
        return my_table
    else:
        raise ValueError('The codon table name was not found in the file')  

def get_insert(sequence, primer1, primer2):
    '''Function to extract the sequence between 2 specified primer or adapter sequences.
    Uses new and improved regex package to allow mismatches and improve number of matches returned.
    The current parameters are a reasonable balance between reducing false negatives
    while minimizing false positives (allows 2 substitutions, 1 deletion, 1 insertion in total).
    Takes the read sequence and 2 flanking primer sequences as peamaters.  Returns the segment
    of sequence in between the flanking sequences'''
    pattern = "(" + primer1 + "(.+?)" + primer2 + ")" + "{s<=2,d<=1,i<=1}" 
    match = regex.search(pattern, sequence, regex.BESTMATCH)
    if match:
        return match.group(2)
    else:
        return 0

def process_fastq_file(file_name, left_primer, right_primer, codon_table, sample, min_len = MIN_LEN):
    '''Open and parse through a gzipped fastq file.
    Extract the insert sequences, translate into peptides, create list of peptides and count up summary. 
    Takes a file name, left and right flanking sequences and a codon table object as arguments.
    Optionally, you can pass a minimum allowable peptide length (default is 3).
    Returns a table of peptides, a peptides counts table and summary statistics dictionary.'''
    print("Processing: ", file_name)
    count_reads = 0
    count_inserts = 0
    count_no_match = 0
    invalid_length = 0
    right_primer = str(Seq(right_primer).reverse_complement())
    my_sample_name = []
    my_names = []
    my_reads = []
    my_qual_score = []
    my_inserts = []
    my_peptides = []
    with gzip.open(file_name,'rt') as file_handle:
        for name, sequence, qual_score in FastqGeneralIterator(file_handle):
            count_reads += 1
            insert = get_insert(sequence, primer1 = left_primer, primer2 = right_primer)
            if insert:
                if len(insert)%3: #a valid insert length is a multiple of 3
                    invalid_length += 1
                else:
                    count_inserts += 1
                    mySeq = Seq(insert).translate(table=codon_table, to_stop=True)
                    if len(mySeq) > min_len:
                        my_sample_name.append(sample)
                        my_names.append(name)
                        my_reads.append(sequence)
                        my_qual_score.append(qual_score)
                        my_inserts.append(insert)
                        my_peptides.append(str(mySeq))
            else:
                count_no_match += 1
    peptide_counts = pd.DataFrame(Counter(my_peptides).most_common(), 
                                  columns=(['peptide', 'count']))
    peptide_counts.insert(0, 'sample', sample)
    results_table = pd.DataFrame(list(zip(my_sample_name, my_names, my_reads, my_qual_score, my_inserts, my_peptides)), 
                                 columns =['sample', 'name', 'read', 'qual_score', 'insert', 'peptide'])
    sample_stats_dict = {'sample':sample,
                  'reads':count_reads,
                  'inserts':count_inserts,
                  'no_match':count_no_match,
                  'invalid_length':invalid_length,
                  'peptides':len(my_peptides),
                  'uniques':len(peptide_counts),
                  'left':left_primer,
                  'right':right_primer,
                  'codon_table': codon_table}
    rows = []
    rows.append(sample_stats_dict)
    sample_stats = pd.DataFrame(rows)
    return results_table, peptide_counts, sample_stats

def get_cyclization_positions(sequence, residues):
    """Function to find the positions of cyclization resudues.
    These positions are used to determine the loop structure, and hence the
    HELM notation.
    
    Pass the function the peptide sequence in our formal notation, and a list
    of cyclization residues (usually ["C"] or ["[Caf]", "C"]
    """
    # split residues into a list (parse on single letters and residues in brackets)
    my_sequence = list(filter(None, re.split("(?=[^\]]*(?:\[|$))", sequence)))
    # return the positons in the sequence where the cyclization should occur
    return [pos + 1 for pos, char in enumerate(my_sequence) if char in residues]

# this should be robustified with a connection to the monomer database
def peptide2helm(peptide, is_cyclized=False, residues=["[Caf]", "C"]):
    '''Function to convert peptide sequence to helm.  Takes a peptide sequence 
    as input.  Returns a string with peptide in HELM format.'''
    if is_cyclized:
        positions = get_cyclization_positions(peptide, residues)
        if len(positions) == 2:
            myPeptide = ".".join(list(filter(None, re.split("(?=[^\]]*(?:\[|$))", peptide))))
            myPeptide = "PEPTIDE1" + "{" + myPeptide + "}$PEPTIDE1,PEPTIDE1," + str(positions[1]) + ":R3-" + str(positions[0]) + ":R3$$$"
            return myPeptide
        else:
            myPeptide = ".".join(list(filter(None, re.split("(?=[^\]]*(?:\[|$))", peptide))))
            myPeptide = "PEPTIDE1" + "{" + myPeptide + "}"
            return myPeptide
    else:
        myPeptide = ".".join(list(filter(None, re.split("(?=[^\]]*(?:\[|$))", peptide))))
        myPeptide = "PEPTIDE1" + "{" + myPeptide + "}"
        return myPeptide

def extract_features(sequence):
    '''Function to compute peptide features.  Takes a peptide sequence and returns
    the features.  Right now only works for natural amino acids.
    '''
    protein_properties = []
    protein_properties.append(ProteinAnalysis(sequence).aromaticity())
    protein_properties.append(ProteinAnalysis(sequence).gravy())
    protein_properties.append(ProteinAnalysis(sequence).instability_index())
    protein_properties.append(ProteinAnalysis(sequence).isoelectric_point())
    protein_properties.extend([value for value in ProteinAnalysis(sequence).molar_extinction_coefficient()]) # reduced and not
    protein_properties.append(ProteinAnalysis(sequence).molecular_weight())
    protein_properties.extend([value for value in ProteinAnalysis(sequence).secondary_structure_fraction()]) # helix, turn, sheet
    return protein_properties

def peptide_parameters(unique_peptides):
    '''Computes the peptide features for all of the peptides in a unique peptides
    list from a sample.  Takes the unique_peptides dataframe as input.  Output
    is a results dataframe.
    '''
    results_list = []
    for row in range(len(unique_peptides)):
        results_list.append(extract_features(unique_peptides.iloc[row].peptide))
    columns = ['aromaticity','gravy','instability_index','isoelectric_point','ext_coeff_reduced','ext_coeff_not_reduced',
           'molecular_weight','helix_fraction','turn_fraction','sheet_fraction']
    results = pd.DataFrame(results_list, columns=columns)
    return results

def process_samples(sample_sheet):
    '''Function to process samples.  It takes a sample sheet as input.  The sample sheet contains
    all of the relavant information for each sample.  It is a csv file format.'''
    samples = pd.read_csv(sample_sheet)
    for index, sample in samples.iterrows():
        my_codon_table = load_custom_codon_table(sample.codon_table)
        results_table, unique_peptides, sample_stats = process_fastq_file(DATA_DIR + sample.file_name, 
                                                      sample.left_primer, 
                                                      sample.right_primer, 
                                                      my_codon_table,
                                                      sample['sample'])
        # add helm notation column
        unique_peptides['helm'] = unique_peptides['peptide'].apply(lambda x : peptide2helm(x, is_cyclized=sample.cyclized))

        # add protein parameters to unique_peptides
        # right now only works for standard codon and natural amino acids
        if sample.codon_table == "Standard":
            parameters = peptide_parameters(unique_peptides)
            unique_peptides = pd.concat([unique_peptides, parameters], axis=1)
        else:
            parameters = pd.DataFrame(columns=['aromaticity','gravy','instability_index','isoelectric_point','ext_coeff_reduced','ext_coeff_not_reduced',
           'molecular_weight','helix_fraction','turn_fraction','sheet_fraction'])
            unique_peptides = pd.concat([unique_peptides, parameters], axis=1)

        # add cyclization info to sample_stats
        sample_stats['cyclized'] = sample.cyclized
        
        ### code to write output to a database will go here... ###
        # for now just writing to files
        results_table.to_csv(OUT_DIR + sample['sample'] + '_data.csv', index=False)
        unique_peptides.to_csv(OUT_DIR + sample['sample'] + '_unique_peptides.csv', index=False)
        sample_stats.to_csv(OUT_DIR + sample['sample'] + '_samples.csv', index=False)



if __name__ == '__main__':
    if len(sys.argv) != 2:
        print('Usage: python peptide_ngs.py sample_sheet.csv')
    else:
        print('Running pipeline on:', sys.argv[1])
        process_samples(sys.argv[1])


#process_samples('sample_sheet.csv')






