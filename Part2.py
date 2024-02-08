import numpy as np
from abc import ABC, abstractmethod
import pandas as pd

import matplotlib.pyplot as plt
import seaborn as sns

from io import BytesIO
import base64


class GeneticEntity(ABC):

    @abstractmethod
    def __init__(self):
        pass


class OrganicElement(GeneticEntity):

    pass


class Nucleotide(OrganicElement):
    __nt_dict = {'A': False, 'T': False, 'U': False, 'G': False, 'C': False, 'N': False}
    __pairs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

    def __init__(self, base):
        self._base = base

    @property
    def base(self):
        return self._base

    @classmethod
    def find_instance(cls, nt):
        if not cls.__nt_dict[nt]:
            cls.__nt_dict[nt] = Nucleotide(nt)
        return cls.__nt_dict[nt]

    @classmethod
    def get_complementary_nt(cls, nt):
        return cls.find_instance(cls.__pairs[nt])

    @classmethod
    def get_nts(cls):
        return list(cls.__nt_dict.keys())

    @classmethod
    def new_nt(cls, nt, pair=False):
        cls.__nt_dict[nt] = False

        if pair:
            cls.__pairs[nt] = pair
        else:
            cls.__pairs[nt] = nt

    def __repr__(self):
        return self._base


class AminoAcid(OrganicElement):
    __AA_dict = {'Y': False, 'N': False, 'A': False, 'F': False, 'L': False,
                 'K': False, 'H': False, 'T': False, 'S': False, 'V': False,
                 'C': False, 'R': False, 'G': False, 'E': False, 'P': False,
                 'M': False, 'I': False, 'Q': False, 'D': False, 'W': False}

    def __init__(self, AA):
        self._AA = AA

    @property
    def AA(self):
        return self._AA

    @classmethod
    def find_instance(cls, AA):
        if not cls.__AA_dict.get(AA):
            cls.__AA_dict[AA] = AminoAcid(AA)
        return cls.__AA_dict[AA]

    @classmethod
    def get_AAs(cls):
        return list(cls.__AA_dict.keys())

    @classmethod
    def new_AA(cls, aa):
        cls.__AA_dict[aa] = False

    def __repr__(self):
        return self._AA


class Sequence(GeneticEntity):
    @abstractmethod
    def __init__(self, name, seq, uom, kind):
        self._name = name
        self._type = kind
        self._stats = Statistics(self._name, seq, uom)

    @property
    def name(self):
        return self._name

    @name.setter
    def name(self, new_name):
        self._name = new_name

    @property
    def type(self):
        return self._type

    @property
    def stats(self):
        return self._stats

    def transcribe(self):
        return self


class AAChain(Sequence):
    __dataset = pd.DataFrame(columns=["Name", "Length", "mRNA of origin", "Seq"])

    def __init__(self, seq, name, origin_mRNA, kind):
        self.__seq = np.array([AminoAcid.find_instance(AA) for AA in seq])
        super().__init__(name, self.__seq, "aa", kind)
        self.origin_mRNA = origin_mRNA
        AAChain.__dataset.loc[len(AAChain.__dataset)] = [self._name, len(self.__seq), self.origin_mRNA, self]

    @property
    def seq(self):
        return self.__seq.copy()

    @property
    def str_seq(self):
        return ''.join([str(AA) for AA in self.__seq])

    @classmethod
    def create_chain(cls, seq, name=False, origin_mRNA='unknown'):

        if len(seq) > 20:
            f = Protein
        else:
            f = Oligopeptide
        return f(seq, name, origin_mRNA)

    @classmethod
    def find_dataset(cls, filter_by_mRNA=False, filter_by_name=False):
        if filter_by_name:
            return AAChain.__dataset.loc[AAChain.__dataset["Name"] == filter_by_name]["Seq"].values[0]

        if filter_by_mRNA:
            return AAChain.__dataset.loc[AAChain.__dataset["mRNA of origin"] == filter_by_mRNA]
        return AAChain.__dataset

    @classmethod
    def order_dataset(cls, column_name, how):
        if column_name == "Name":
            return AAChain.__dataset.sort_values("Seq", ascending=how)
        else:
            return AAChain.__dataset.sort_values(column_name, ascending=how)

    def __gt__(self, other):
        if self.name.split("_")[0] > other.name.split("_")[0]:
            return True
        elif self.name.split("_")[0] < other.name.split("_")[0]:
            return False
        else:
            return int(self.name.split("_")[1]) > int(other.name.split("_")[1])


    def __str__(self):
        s = "N- "
        w = self.str_seq
        for i in range(0, self.stats.length, 40):
            s += w[i:i+40] + "\n\n   "
        s = s[:-5]+" -C"
        return s


class Oligopeptide(AAChain):
    __counter = 0

    def __init__(self, seq, name, origin_mRNA):
        if not name:
            name = f"Oligo_{Oligopeptide.__counter}"
            Oligopeptide.__counter += 1

        super().__init__(seq, name, origin_mRNA, "Oligopeptide")


class Protein(AAChain):
    __counter = 0

    def __init__(self, seq, name, origin_mRNA):
        if not name:
            name = f"Protein_{Protein.__counter}"
            Protein.__counter += 1

        super().__init__(seq, name, origin_mRNA, "Protein")


class DNA(Sequence):

    def __init__(self, seq, name):
        self.__fwd_strand = np.array([Nucleotide.find_instance(nt) for nt in seq])
        self.__rev_strand = np.array([Nucleotide.get_complementary_nt(nt) for nt in seq])
        super().__init__(name, self.__fwd_strand, "bp", "DNA")

    @property
    def fwd_strand(self):
        return self.__fwd_strand.copy()

    @property
    def rev_strand(self):
        return self.__rev_strand.copy()

    @staticmethod
    def __str_array(array):
        return ''.join([str(nt) for nt in array])

    @property
    def str_seq(self):
        return self.__str_array(self.__fwd_strand)

    def transcribe(self, reverse=False):

        if reverse:
            mRNA_seq = self.rev_strand[::-1]
            mRNA_seq[mRNA_seq == Nucleotide.find_instance("T")] = Nucleotide.find_instance("U")
            strand = "-"
        else:
            mRNA_seq = self.fwd_strand
            mRNA_seq[mRNA_seq == Nucleotide.find_instance("T")] = Nucleotide.find_instance("U")
            strand = "+"

        return mRNA(mRNA_seq, origin_DNA=self._name, strand=strand)

    def __str__(self):
        s = ""

        if self.stats.length <= 40:
            s += (f"5'- {DNA.__str_array(self.__fwd_strand[:40])} -3'\n"
                  f"3'- {DNA.__str_array(self.__rev_strand[:40])} -5'")
            return s

        for i in range(0, self.stats.length, 40):
            if i == 0:
                s += (f"5'- {DNA.__str_array(self.__fwd_strand[i:i+40])}\n"
                      f"3'- {DNA.__str_array(self.__rev_strand[i:i+40])}\n\n")
            elif i+40 >= self.stats.length:
                s += (f"    {DNA.__str_array(self.__fwd_strand[i:i + 40])} -3'\n"
                      f"    {DNA.__str_array(self.__rev_strand[i:i + 40])} -5'")
            else:
                s += (f"    {DNA.__str_array(self.__fwd_strand[i:i + 40])}\n"
                      f"    {DNA.__str_array(self.__rev_strand[i:i+40])}\n\n")
        return s


class mRNA(Sequence):
    __translation = {'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L',
                     'CUA': 'L', 'CUG': 'L', 'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M',
                     'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V', 'UCU': 'S', 'UCC': 'S',
                     'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
                     'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A',
                     'GCA': 'A', 'GCG': 'A', 'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*',
                     'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 'AAU': 'N', 'AAC': 'N',
                     'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
                     'UGU': 'C', 'UGC': 'C', 'UGA': '*', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R',
                     'CGA': 'R', 'CGG': 'R', 'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
                     'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

    def __init__(self, seq, name=False, origin_DNA=False, strand=False):
        if isinstance(seq, np.ndarray):
            self.__transcript = seq
        else:
            self.__transcript = np.array([Nucleotide.find_instance(nt) for nt in seq])

        if not name:
            name = "mRNA transcript of " + origin_DNA + f' [{strand}]'

        super().__init__(name,self.__transcript,"nt", "mRNA")
        self.__string = str(self)[4:-4].replace("\n\n    ", "")

    @property
    def transcript(self):
        return self.__transcript.copy()

    @property
    def str_seq(self):
        return self.__string

    def translate(self):
        for ORF in range(0, 3):
            codon_seq = np.array([mRNA.__translation.get(self.__string[i:i + 3], 'invalid')
                                  for i in range(ORF, self.stats.length, 3)])
            stops = np.where(codon_seq == '*')[0]
            starts = np.where(codon_seq == 'M')[0]

            for start in starts:
                for stop in stops:
                    if start < stop:
                        if len(np.where(codon_seq[start:stop] == 'invalid')[0]) == 0:
                            AAChain.create_chain(codon_seq[start:stop], origin_mRNA=self.name)
                            break
                    else:
                        stops = stops[1:]

    def __str__(self):
        s = "5'- "
        w = ''.join([str(nt) for nt in self.__transcript])
        for i in range(0, self.stats.length, 40):
            s += w[i:i+40] + "\n\n    "
        s = s[:-6]+" -3'"
        return s


class Statistics:
    __kv_dict = {"aa": [AminoAcid.get_AAs(), AminoAcid.find_instance],
                 "bp": [Nucleotide.get_nts()[:2] + Nucleotide.get_nts()[3:], Nucleotide.find_instance],
                 "nt": [Nucleotide.get_nts()[:1] + Nucleotide.get_nts()[2:], Nucleotide.find_instance]}
    def __init__(self, name, seq, unit):
        self.name = name
        self.seq = seq
        self.length = len(self.seq)
        self.unit = unit

    @staticmethod
    def percent(val, length):
        return round(100*val/length, 2)

    def frequencies(self, rel=True):
        freq = {}
        keys = Statistics.__kv_dict[self.unit][0]
        values_fn = Statistics.__kv_dict[self.unit][1]

        for key in keys:
            freq[key] = np.count_nonzero(self.seq == values_fn(key))

        if rel:
            for key, value in freq.items():
                freq[key] = [value,
                             Statistics.percent(value, self.length)]
        return freq

    def GC_content(self):
        values_fn = Statistics.__kv_dict[self.unit][1]
        gc_count = (np.count_nonzero(self.seq == values_fn("G")) +
                    np.count_nonzero(self.seq == values_fn("C")))
        return Statistics.percent(gc_count, self.length)

    def make_histogram(self):

        if self.unit == "aa":
            cols = ["#0AFFFF", "#FF1DCE", "#C7EA46", "#50C878", "#4CBB17",
                    "#800080", "#005A9C", "#FF5F1F", "#F04A00", "#00A86B",
                    "#FFF700", "#000080", "#FF7518", "#DC143C", "#FFD800",
                    "#0B6623", "#4CBB17", "#E0115F", "#FF0800", "#007FFF"]
            title = "Amino Acid"
        else:
            cols = ["#E68994", "#337167", "#96031A", "#FFC759", "#28231D"]
            title = "Nucleotide"

        freqs = self.frequencies(rel=False)

        sns.set_style(rc={'axes.facecolor': (0, 0, 0, 0), 'figure.facecolor': (0, 0, 0, 0)})
        fig = sns.barplot(freqs, palette=cols)

        sns.despine(right=True, top=True, left=True)

        plt.yticks([])
        plt.xlabel(f'{title}', fontsize=12, fontweight='semibold')
        plt.title(f'{title} Absolute Frequencies', fontsize=14, fontweight='semibold')

        for i in range(len(fig.containers)):
            fig.bar_label(fig.containers[i], fmt='%d', label_type='edge')

        image_stream = BytesIO()
        plt.savefig(image_stream, format='png')
        image_stream.seek(0)
        plt.clf()

        encode_image = base64.b64encode(image_stream.getvalue()).decode('utf-8')

        return encode_image

    def create_file(self):
        file_content = f"Name: {self.name}"
        file_content += f"\r\n\r\nLength: {self.length}\r\n"
        for key, value in self.frequencies(rel=False).items():
            file_content += f"\r\n{key}: {value}"
        if self.unit != "aa":
            file_content += f"\r\n\r\nGC_content: {self.GC_content()}"

        in_memory_file = BytesIO()
        in_memory_file.write(file_content.encode('utf-8'))
        in_memory_file.seek(0)

        return in_memory_file
