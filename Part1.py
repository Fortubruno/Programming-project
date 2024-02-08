import pandas as pd
from io import BytesIO


class DATASET:
    def __init__(self):
        self.__Dataset = pd.DataFrame(columns=["Type", "Name", "Seq"])

    def add_FASTA(self, file):
        l = Fasta.parse_FASTA(file)

        if not self.__is_in(l[1]):
            self.__Dataset.loc[len(self.__Dataset)] = l
        return l

    def add_obj(self, obj):
        if not self.__is_in(obj.name):
            self.__Dataset.loc[len(self.__Dataset)] = [obj.type, obj.name, obj]

    def __find_indexes(self, name):
        indexes = []
        for row in self.__Dataset.iterrows():
            if name in row[1]["Name"] and name != row[1]["Name"]:
                indexes.append(row[0])
        return indexes

    def find(self, name="optional", exactly=True, each=False):
        if each:
            return self.__Dataset.groupby("Type")

        if exactly:
            x = self.__Dataset.loc[self.__Dataset["Name"] == name, "Seq"]
            if len(x) == 0:
                return False
            return x.values[0]

        return self.__Dataset.loc[self.__find_indexes(name)].groupby("Type")

    def replace_FASTA(self, obj):
        self.__Dataset.loc[self.__Dataset['Name'] == obj.name, ['Seq']] = obj

    def remove(self, *names):
        for name in names:
            self.__Dataset.drop(self.__Dataset.index[self.__Dataset["Name"] == name], inplace=True)

    def __is_in(self, name):
        if len(self.__Dataset.loc[self.__Dataset['Name'] == name]["Seq"]) > 0:
            return True
        return False

    def is_empty(self):
        if len(self.__Dataset) == 0:
            return True
        return False


class Fasta:

    @staticmethod
    def create_FASTA(name, seq):
        file_content = f'>{name}'
        for i in range(0, len(seq), 80):
            file_content += f"\r\n{seq[i:i+80]}"

        in_memory_file = BytesIO()
        in_memory_file.write(file_content.encode('utf-8'))
        in_memory_file.seek(0)

        return in_memory_file

    @staticmethod
    def parse_FASTA(fasta_file):
        name, *seq = fasta_file.split('\r\n')
        name = name[1:].replace("/", "-")
        seq = ''.join(seq).replace('\r\n', '')

        if 'U' in seq:
            kind = 'mRNA'
        elif 'M' in seq:
            kind = 'AA_chain'
        else:
            kind = 'DNA'

        return [kind, name, seq]
