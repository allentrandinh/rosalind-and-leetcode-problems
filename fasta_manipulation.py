#functions to process input from files downloaded from rosalind

def single_seq_noheading(file_path):
    '''input is path to plain txt, no fasta header
     :return: dna sequence with no new line'''
    file_open=open(file_path,"r")
    dna_seq = file_open.read()
    dna_seq = dna_seq.replace("\n", "")
    return(dna_seq)

def multiple_fasta_seq(file_path):
    '''
    input txt file contain multiple fasta sequences
    :return: dna sequence in a list without header
    '''
    file_open = open(file_path, "r")
    dna_seq = file_open.read()
    splitted = dna_seq.split(">")[1:]
    fasta_dna = {}
    for each_fasta in splitted:
        content = each_fasta.split("\n")
        header = content[0]
        dna_content = ''.join(content[1:])
        fasta_dna[header]=dna_content
    return(fasta_dna)

def multiple_seq_noheading(file_path):
    '''
    input text file contain multiple sequences with no headers separated by \n
    :param file_path
    :return: list contains sequences
    '''
    file_open = open(file_path, "r")
    dna_seq = file_open.read()
    splitted = dna_seq.split("\n")
    return(splitted)


