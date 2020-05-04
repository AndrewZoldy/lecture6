import sys
import argparse

class Base:
    minimal_orf_lenght = 100
    translation_dict = {'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'UUA': 'L', 'UUG': 'L',
                        'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L', 'CGU': 'R', 'CGC': 'R',
                        'CGA': 'R', 'CGG': 'R', 'AGA': 'R', 'AGG': 'R', 'AAA': 'K', 'AAG': 'K',
                        'AAU': 'N', 'AAC': 'N', 'AUG': 'M', 'GAU': 'D', 'GAC': 'D', 'UUU': 'F',
                        'UUC': 'F', 'UGU': 'C', 'UGC': 'C', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P',
                        'CCG': 'P', 'CAA': 'Q', 'CAG': 'Q', 'UCU': 'S', 'UCC': 'S', 'UCA': 'S',
                        'UCG': 'S', 'AGU': 'S', 'AGC': 'S', 'GAA': 'E', 'GAG': 'E', 'ACU': 'T',
                        'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G',
                        'GGG': 'G', 'UGG': 'W', 'CAU': 'H', 'CAC': 'H', 'UAU': 'Y', 'UAC': 'Y',
                        'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V',
                        'GUG': 'V', 'UAG': '*', 'UGA': '*', 'UAA': '*'}
    weights = {}
    water_weight = 18

    def __init__(self, s):
        self.s = s.upper()

    def __str__(self):
        return self.s

    def as_dna(self):
        raise NotImplementedError()

    def as_rna(self):
        raise NotImplementedError()

    def as_protein(self):
        raise NotImplementedError()

    def weight(self):
        return sum([self.weights[m] - self.water_weight for m in str(self)]) + self.water_weight

    def show_type(self):
        raise NotImplementedError()

    @staticmethod
    def check_sequence(s):
        raise NotImplementedError()

    @staticmethod
    def find_longest_orf(s):
        raise NotImplementedError()


class DNA(Base):
    weights = {'A': 331.2218, 'T': 322.2085, 'G': 347.2212, 'C': 307.1971}

    def as_dna(self):
        return self

    def show_type(self):
        return 'dna'


class RNA(Base):
    back_transcription_dict = {'U': 'A', 'A': 'T', 'G': 'C', 'C': 'G'}
    start_codon = 'AUG'
    stop_codons = ['UAG', 'UGA', 'UAA']
    weights = {'A': 347.2212, 'G': 363.2206, 'C': 323.1965, 'U': 324.1813}

    def as_dna(self):
        return CodingDNA(''.join([self.back_transcription_dict[letter] for letter in self.s]))

    def as_rna(self):
        return self

    def as_protein(self):
        s = str(self.as_rna())
        protein = str()
        for n in range(0, len(s), 3):
            protein += self.translation_dict[s[n:n + 3]]
        return Protein(protein)

    def show_type(self):
        return 'rna'

    @staticmethod
    def check_sequence(s):
        if 'U' in s:
            return True
        else:
            return False

    @staticmethod
    def find_longest_orf(s):
        indexes = None
        max_len = RNA.minimal_orf_lenght
        for i in range(0, len(s) - 3):
            if s[i:i + 3] == RNA.start_codon:
                for j in range(i + 3, len(s), 3):
                    if s[j:j + 3] in RNA.stop_codons:
                        if len(s[i:j + 3]) > max_len:
                            max_len = len(s[i:j + 3])
                            indexes = [i, j + 3]
                        break
        if indexes:
            return s[indexes[0]:indexes[1]]
        else:
            return None

    @staticmethod
    def find_orfs(s):
        orf_list = []
        max_len = RNA.minimal_orf_lenght
        for i in range(0, len(s) - 3):
            if s[i:i + 3] == RNA.start_codon:
                for j in range(i + 3, len(s), 3):
                    if s[j:j + 3] in RNA.stop_codons:
                        if len(s[i:j + 3]) > max_len:
                            orf_list.append(s[i:j + 3])
                        break
        if len(orf_list) > 0:
            return orf_list
        else:
            return None


class CodingDNA(DNA):
    transcription_dict = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
    start_codon = 'TAC'
    stop_codons = ['ATC', 'ACT', 'ATT']

    def as_rna(self):
        return RNA(''.join([self.transcription_dict[letter] for letter in self.s]))

    def as_protein(self):
        s = str(self.as_rna())
        protein = ''
        for n in range(0, len(s), 3):
            protein += self.translation_dict[s[n:n + 3]]
        return Protein(protein)

    @staticmethod
    def check_sequence(s):
        if 'U' not in s:
            return True
        else:
            return False

    @staticmethod
    def find_longest_orf(s):
        indexes = None
        max_len = CodingDNA.minimal_orf_lenght
        for i in range(0, len(s) - 3):
            if s[i:i + 3] == CodingDNA.start_codon:
                for j in range(i + 3, len(s), 3):
                    if s[j:j + 3] in CodingDNA.stop_codons:
                        if len(s[i:j + 3]) > max_len:
                            max_len = len(s[i:j + 3])
                            indexes = [i, j + 3]
                        break
        if indexes:
            return s[indexes[0]:indexes[1]]
        else:
            return None

    @staticmethod
    def find_orfs(s):
        orf_list = []
        max_len = CodingDNA.minimal_orf_lenght
        for i in range(0, len(s) - 3):
            if s[i:i + 3] == CodingDNA.start_codon:
                for j in range(i + 3, len(s), 3):
                    if s[j:j + 3] in CodingDNA.stop_codons:
                        if len(s[i:j + 3]) > max_len:
                            orf_list.append(s[i:j + 3])
                        break
        if len(orf_list) > 0:
            return orf_list
        else:
            return None

class NoncodingDNA(DNA):
    transcription_dict = {'A': 'U', 'T': 'A', 'G': 'C', 'C': 'G'}
    start_codon = 'ATG'
    stop_codons = ['TAG', 'TGA', 'TAA']

    def as_rna(self):
        return RNA(self.s.replace('T', 'U'))

    def as_protein(self):
        s = str(self.as_rna())
        protein = ''
        for n in range(0, len(s), 3):
            protein += self.translation_dict[s[n:n + 3]]
        return Protein(protein)

    @staticmethod
    def check_sequence(s):
        if 'U' not in s:
            return True
        else:
            return False

    @staticmethod
    def find_longest_orf(s):
        indexes = None
        max_len = NoncodingDNA.minimal_orf_lenght
        for i in range(0, len(s) - 3):
            if s[i:i + 3] == NoncodingDNA.start_codon:
                for j in range(i + 3, len(s), 3):
                    if s[j:j + 3] in NoncodingDNA.stop_codons:
                        if len(s[i:j + 3]) > max_len:
                            max_len = len(s[i:j + 3])
                            indexes = [i, j + 3]
                        break
        if indexes:
            return s[indexes[0]:indexes[1]]
        else:
            return None

    @staticmethod
    def find_orfs(s):
        orf_list = []
        max_len = NoncodingDNA.minimal_orf_lenght
        for i in range(0, len(s) - 3):
            if s[i:i + 3] == NoncodingDNA.start_codon:
                for j in range(i + 3, len(s), 3):
                    if s[j:j + 3] in NoncodingDNA.stop_codons:
                        if len(s[i:j + 3]) > max_len:
                            orf_list.append(s[i:j+3])
                        break

        if len(orf_list) > 0:
            return orf_list
        else:
            return None

class Protein(Base):
    weights = {'A': 89.0932, 'L': 131.1729, 'R': 174.201, 'K': 146.1876, 'N': 132.1179, 'M': 149.2113,
               'D': 133.1027, 'F': 165.1891, 'C': 121.1582, 'P': 115.1305, 'Q': 146.1445, 'S': 105.0926,
               'E': 147.1293, 'T': 119.1192, 'G': 75.0666, 'W': 204.2252, 'H': 155.1546, 'Y': 181.1885,
               'I': 131.1729, 'V': 117.1463, '*': 0}

    def as_protein(self):
        return self

    def check_sequence(s):
        for letter in ['D' , 'E' , 'I' , 'L' , 'F' , 'V' , 'R' ,
                       '*' , 'K' , 'P' , 'W' , 'N' , 'Q' , 'H' ,
                       'M' , 'S' , 'Y']:
            if letter in s:
                return True
        return False

    def show_type(self):
        return 'protein'

def create_longest_orf_obj(sequence):
    possible = (NoncodingDNA, CodingDNA, RNA)
    for t in possible:
        ORF = t.find_longest_orf(sequence)
        if t.check_sequence(sequence) and ORF:
            return t(ORF)

def create_orf_objs(sequence):
    possible = (NoncodingDNA, CodingDNA, RNA)
    for t in possible:
        status=True
        ORFs = t.find_orfs(sequence)
        if t.check_sequence(sequence) and ORFs:
            return [t(ORF) for ORF in ORFs]

def create_obj(sequence):
    possible = (Protein, NoncodingDNA, CodingDNA, RNA)
    for t in possible:
        if t.check_sequence(sequence):
            return t(sequence)


#data = sys.stdin.read()
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='tool for working around biological sequence data')
    parser.add_argument('-m', '--mode', type=str, default='all_orfs', help='possible work mods: view,'
                                                                           'longest_orf, all_orfs')
    args = parser.parse_args()
    work_modes= ['view', 'longest_orf', 'all_orfs']
    current_mode = args.mode
    with open('lec6_hw.fsa', 'r') as inp:
        data = inp.read()

    if current_mode == 'longest_orf':
        longest_ORFs = [create_longest_orf_obj(''.join(seq.split('\n')[1:])) for seq in data.split('>')[1:]]
        for orf in longest_ORFs:
            dna = orf.as_dna()
            rna = orf.as_rna()
            protein = orf.as_protein()
            print('DNA Sequence: {:>5}  Weight: {:>5}\n'.format(str(dna), dna.weight()))
            print('RNA Sequence: {:>5}  Weight: {:>5}\n'.format(str(rna), rna.weight()))
            print('Protein Sequence: {:>5}  Weight: {:>5}\n'.format(str(protein), protein.weight()))

    if current_mode == 'all_orfs':
        ORF_dict = dict()
        for seq in data.split('>')[1:]:
            ORF_dict.update({seq.split('\n')[0]: create_orf_objs(''.join(seq.split('\n')[1:]))})
        counter = 0
        for seq in ORF_dict:
            for orf in ORF_dict[seq]:
                counter+=1
                dna = orf.as_dna()
                rna = orf.as_rna()
                protein = orf.as_protein()
                print('Gene: {:>5}  DNA Sequence: {:>5}  Weight: {:>5}\n'.format(seq, str(dna), dna.weight()))
                print('Gene: {:>5}  RNA Sequence: {:>5}  Weight: {:>5}\n'.format(seq, str(rna), rna.weight()))
                print('Gene: {:>5}  Protein Sequence: {:>5}  Weight: {:>5}\n'.format(seq, str(protein), protein.weight()))
    print(counter)

