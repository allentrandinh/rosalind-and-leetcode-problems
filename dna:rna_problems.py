import fasta_manipulation as fm
import time

def count_each_nucleotide(DNA):
    '''
    :return: list with count of A,C,G,T respectively
    '''
    count_num = []
    for each_nu in ["A","C","G","T"]:
        count_num.append(DNA.count(each_nu))
    return(count_num)

def dna_T2U(DNA):
    '''
    :param DNA: 5'-3' DNA
    :return: 5'-3' RNA
    '''
    return(DNA.replace("T","U"))


#transcription (3'-5' DNA to RNA)

def dna_2_(DNA):
    '''
    :param DNA: 3'-5' DNA
    :return: 5'-3' RNA
    '''
    def base_DNA_to_RNA(one_base):
        DNA_base = ["A", "T", "G", "C"]
        RNA_base = ["U", "A", "C", "G"]
        return (RNA_base[DNA_base.index(one_base)])
    RNA = list(map(base_DNA_to_RNA,DNA))
    return(''.join(RNA))

#complementary strand
def complementary_strand(DNA):
    ''':return: complementary strand in same orientation'''
    def com_base(one_base):
        base_strand_1 = ["A", "T", "G", "C"]
        base_strand_2 = ["T", "A", "C", "G"]
        return(base_strand_2[base_strand_1.index(one_base)])
    com_dna = list(map(com_base,DNA))
    return(''.join(com_dna)[::-1])

#gc content
def gc_content(DNA):
    '''
    :return: GC content of DNA
    '''
    return((DNA.count("G")+DNA.count("C"))/len(DNA))

#input fasta seqs, output fasta of seq with higest GC content
def highest_gc_content(dna):
    '''
    :param dna: dictionary with key = fasta header, value = dna seq
    dictionary can be created with multiple_fasta_seq in fasta_manipulation.py
    :return: fasta header whose sequence have highest gc_content and gc_content of it
    '''
    fasta = list(dna.keys())
    dna_seqs = list(dna.values())
    calculated_gc = list(map(gc_content,dna_seqs))
    index = calculated_gc.index(max(calculated_gc))
    return(fasta[index],max(calculated_gc)*100)

#input 2 sequences of equal length, finds number of point mutations
def num_point_mutation(dna_list):
    '''
    :param dna_list: list containing 2 dna seqs of equals length
    :return: number of point mutations
    '''
    point_mutation_count=0
    for each_base in range(len(dna_list[0])):
        if dna_list[0][each_base] != dna_list[1][each_base]:
            point_mutation_count += 1
    return(point_mutation_count)

def motif_find(DNA,motif):
    '''
    :param DNA: DNA sequence to look for motif
    :param motif: motif sequence
    :return: position of motif
    '''
    position=[]
    for i in range(len(DNA)-len(motif)):
        if DNA[i:i+len(motif)]==motif:
            position.append(i+1)
    return(position)

def consensus_matrix(DNA_list):
    '''
    :param DNA_list: list of dna sequence of equal length
    :return: A consensus string and profile matrix for the collection, order are A,T,G,C
    '''
    ATGC_profile=["","","",""]
    dna_base = ["A","C","G","T"]
    consensus = ""
    for i in range(len(DNA_list[0])):
        position = [DNA_list[k][i] for k in range(len(DNA_list))]
        for m in range(4):
            ATGC_profile[m] += str(position.count(dna_base[m]))
    for n in range(len(DNA_list[0])):
        each_base_count = [ATGC_profile[k][n] for k in range(4)]
        consensus += dna_base[each_base_count.index(max(each_base_count))]
    profile = f"A: {ATGC_profile[0].replace('',' ')}\n" \
              f"C: {ATGC_profile[1].replace('',' ')}\n" \
              f"G: {ATGC_profile[2].replace('',' ')}\n" \
              f"T: {ATGC_profile[3].replace('',' ')}\n"
    return(consensus, profile)

gencode = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
    'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W'}
def candidate_protein(DNA):
    '''
    :param DNA
    :return: candidate protein on all 6 possible frames, candidate without stop codon does not count
    '''
    sequence = [DNA,complementary_strand(DNA)]
    protein_list=[]
    for seq in sequence:
        possible_start_site=motif_find(seq,"ATG")
        #-1 to switch back to 0-based system
        for each_start_site in possible_start_site:
            protein=""
            while each_start_site < len(DNA)-2:
                codon = seq[each_start_site-1:each_start_site+2]
                if codon in ["TAA","TAG","TGA"]:
                    protein_list.append(protein)
                    break
                else:
                    protein += gencode[codon]
                    each_start_site += 3
    return(set(protein_list))

def is_palindromic(DNA):
    '''
    :param DNA:
    :return: TRUE if the DNA seq is palindromic, FALSE if no
    '''
    return(DNA == complementary_strand(DNA))

def palindromic_site(DNA,min_motif_length,max_motif_length):
    '''
    :param DNA:
    :return: position of site, length of site
    '''
    palindromic_sites = {}
    for i in range(min_motif_length,max_motif_length+1):
        for k in range(len(DNA)-i+1):
            subset= DNA[k:(k+i)]
            if is_palindromic(subset):
                palindromic_sites[k+1]=i
    return(palindromic_sites)

def rna_splicing(DNA_list):
    '''
    :param DNA_list: dna seqs in a lits, first sequence: full length dna start with ATG, other sequences: introns
    :return: protein resulting from splicing
    '''
    full_length = DNA_list[0]
    intron_list = DNA_list[1:]
    for intron in intron_list:
        full_length = full_length.replace(intron,"")
    protein=""
    for i in range(0,len(full_length)-2,3):
        codon = full_length[i:i+3]
        if codon in ["TAA", "TAG", "TGA"]:
            break
        else:
            protein += gencode[codon]
    return(protein)

codon_code = {
    'AUA':'I', 'AUC':'I', 'AUU':'I', 'AUG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACU':'T',
    'AAC':'N', 'AAU':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGU':'S', 'AGA':'R', 'AGG':'R',
    'CUA':'L', 'CUC':'L', 'CUG':'L', 'CUU':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCU':'P',
    'CAC':'H', 'CAU':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGU':'R',
    'GUA':'V', 'GUC':'V', 'GUG':'V', 'GUU':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCU':'A',
    'GAC':'D', 'GAU':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGU':'G',
    'UCA':'S', 'UCC':'S', 'UCG':'S', 'UCU':'S',
    'UUC':'F', 'UUU':'F', 'UUA':'L', 'UUG':'L',
    'UAC':'Y', 'UAU':'Y', 'UAA':'_', 'UAG':'_',
    'UGC':'C', 'UGU':'C', 'UGA':'_', 'UGG':'W'}

def rna_2_protein(RNA):
    '''
    :param RNA: CDS (length must be divisible by 3)
    :return: protein seq
    '''
    def codon_to_aa(codon):
        return(codon_code[codon])
    codon = [RNA[3*i:3*i+3] for i in range(int(len(RNA)/3))]
    protein = ''.join(list(map(codon_to_aa,codon)))
    return(protein)

dna = fm.multiple_fasta_seq("/Users/apd20500/Desktop/rosalind_lcsm.txt")
dna_seq = list(dna.values())

def shared_motif(DNA_list):
    '''
    :param DNA_list: DNA seqs in a list
    :return: longest motif found (min_lenght = 2)
    '''
    #min motif = 2
    #extract possible motif from first DNA, check if it appears in other motif, if no, breaks immediately
    for i in range(len(DNA_list[0]),1,-1):
        for k in range(0,int(len(DNA_list[0])-i+1)):
            candidate_motif = DNA_list[0][k:k+i]
            count = 0
            for each_other_sequence in DNA_list[1:]:
                if candidate_motif in each_other_sequence:
                    count = count +1
                else:
                    break
            if count == len(DNA_list[1:]):
                return(candidate_motif)



#dna = fm.single_seq_noheading("/Users/apd20500/Desktop/rosalind_revc.txt")