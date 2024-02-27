import itertools
import os

from constants.constants import (
    DATA_PATH,
    MOTIFF_BASE_LENGTH,
    NEGATIVE,
    POSITIVE_HIGH,
    THRESHOLD,
    VALID_AMINOACIDS,
)


class Parser:
    # Class to parse all the used files
    def __init__(
        self,
        IEDB_file=None,
        csv_file=None,
        aligment_file=None,
        valid_letters=VALID_AMINOACIDS,
        length=MOTIFF_BASE_LENGTH,
    ):
        """
        :param IEDB_file: IEDB database in csv format
        :param csv_file: Additional file with user defined data format(peptide;hla;qualitative_value)
        :param aligment_file: file with the HLA sequences in selex format
        :param valid_letters: List of possible aminoacids
        :param length: length of the peptides to parse
        """
        # Files
        self.IEDB_file = IEDB_file
        self.csv_file = csv_file
        self.aligment_file = aligment_file
        self.positions_file = os.path.join(DATA_PATH, "key_positions_%s.txt" % length)
        self.valid_letters = set(valid_letters)
        self.length = length
        self.minimum_data_signal = 50
        self.minimum_data_background = 10

        # Containers
        self.hla_set = set()
        self.hla_sequences = {}

    def set_minimum_data_signal(self, minimum_data):
        self.minimum_data_signal = minimum_data

    def set_minimum_data_background(self, minimum_data):
        self.minimum_data_background = minimum_data

    @staticmethod
    def process_csv_data_generator(file):
        peptide = None
        hla = None
        qual = None
        try:
            with open(file, "r") as inn:
                for line in inn:
                    line = line.rstrip()
                    line = line.split(";")
                    peptide = line[0]
                    hla = line[1]
                    qual = line[2]
                    yield peptide, hla, qual
        except TypeError or IOError:
            yield peptide, hla, qual

    @staticmethod
    def process_IEDB_data_generator(file):
        hla = None
        qualitativeValue = None
        peptide = None
        with open(file, "r") as f:
            next(f)
            next(f)
            for line in f:
                h_line = line.split('","')
                try:
                    reference = h_line[3]
                    peptide = h_line[11]
                    method = h_line[79]
                    assay = h_line[80]
                    units = h_line[81]
                    qualitativeValue = h_line[83]
                    quantitativeValue = h_line[85]
                    hla = h_line[95]  # .rstrip()
                except IndexError:
                    continue
                yield peptide, hla, qualitativeValue

    @staticmethod
    def load_csv_matrix(matrix_file, invert=False):
        matrix_dict = {}
        index_position = {}
        with open(matrix_file, "r") as inn:
            first_line = inn.readline().rstrip()
            splitted_line = first_line.split(",")
            splitted_line.pop(0)
            index_position = {i: x for i, x in enumerate(splitted_line)}
            matrix_dict = {x: {} for x in splitted_line}
            for line in inn:
                line = line.rstrip()
                line = line.split(",")
                letter = line.pop(0)
                for i, entry in enumerate(line):
                    entry = float(entry)
                    if invert:
                        entry *= -1
                    matrix_dict[index_position[i]].setdefault(letter, entry)
        return matrix_dict

    def parse_data(self):

        # Initialize variables
        final_dict = {}
        test_data_dict = {}
        hla_with_signal = set()
        hla_with_background = set()
        used_peptides = {}

        # Extract and verify data
        print("     Reading data files")
        for element in itertools.chain(
            self.process_IEDB_data_generator(self.IEDB_file),
            self.process_csv_data_generator(self.csv_file),
        ):
            peptide, hla, qual = element
            if not peptide or not hla or not qual:
                continue
            peptide_set = set(peptide)
            if (
                hla in self.hla_set
                and not peptide_set.difference(self.valid_letters)
                and len(peptide) == self.length
            ):
                # sets remove identical repeated entry's without removing all of them
                used_peptides.setdefault(hla, {}).setdefault(peptide, set()).add(qual)
                final_dict.setdefault(qual, {}).setdefault(hla, set()).add(peptide)
                test_data_dict.setdefault(hla, {}).setdefault(qual, set()).add(peptide)

        # Filter repeated entry's
        # (those peptides that appear more than once on a given hla and have different experimental characterization)
        print("     Curating data")
        for hla in used_peptides:
            for peptide in used_peptides[hla]:
                different_qualitatives = used_peptides[hla][peptide]
                if len(different_qualitatives) > 1:
                    for qual in different_qualitatives:
                        final_dict[qual][hla].remove(peptide)
                        test_data_dict[hla][qual].remove(peptide)

        # Remove HLAs that do not reach minimum data requirements
        for hla, peptides_set in final_dict[POSITIVE_HIGH].items():
            if len(peptides_set) >= self.minimum_data_signal:
                hla_with_signal.add(hla)
        for hla, peptides_set in final_dict[NEGATIVE].items():
            if len(peptides_set) >= self.minimum_data_background:
                hla_with_background.add(hla)
        correct_hlas = list(
            self.hla_set.intersection(hla_with_signal, hla_with_background)
        )

        return final_dict, test_data_dict, correct_hlas

    def parse_key_position_file(self):
        key_pos = {}
        sim_weight = {}
        try:
            threshold = THRESHOLD[self.length]
        except KeyError:
            threshold = "default"
        with open(self.positions_file, "r") as inn:
            for line in inn:
                line = line.split()
                pos = int(line[0])
                residue = int(line[1])
                crystal_count = int(line[2])
                if crystal_count >= threshold:
                    key_pos.setdefault(pos, []).append(residue)
                    sim_weight.setdefault(pos, {}).setdefault(residue, crystal_count)
                # env = [int(x) for x in line[1].split(',')]
                # key_pos[pos] = env
        return key_pos, sim_weight

    def parse_aligment_file(self, alignment_file):
        with open(alignment_file, "r") as inn:
            for line in inn:
                line = line.rstrip()
                line = line.split()
                hla = line[0]
                sequence = line[1]
                self.hla_sequences.setdefault(hla, sequence)
                self.hla_set.add(hla)
        return self.hla_sequences

    def parse_all(self):

        # Process data (the order of execution is important)
        print("Reading Aligment file %s" % self.aligment_file)
        hla_aligment = self.parse_aligment_file(self.aligment_file)
        print("Extracting environment")
        key_positions, sim_weight = self.parse_key_position_file()
        print("Parsing Data")
        data, test_data, hla_list = self.parse_data()

        # Purge incorrect data
        hla_in_aligment = list(hla_aligment.keys())
        for hla in hla_in_aligment:
            if hla not in hla_list:
                hla_aligment.pop(hla, None)

        print("Parsing Finished")
        return data, test_data, hla_list, hla_aligment, key_positions, sim_weight


if __name__ == "__main__":
    print("PARSING")
    parser = Parser(
        os.path.join(__file__, "../data/IEDB_data.csv"),
        None,
        aligment_file=os.path.join(__file__, "../data/hla_aligment.txt"),
    )
    data, test_data, hla_list, hla_aligment, key_positions, sim_weight = (
        parser.parse_all()
    )
    print(list(test_data.keys()))
    with open(os.path.join(__file__, "../test/data_to_predict.txt"), "w") as inn:
        for hla in hla_list:
            for qual in test_data[hla]:
                for peptide in test_data[hla][qual]:
                    inn.write("%s\t%s\t%s\n" % (peptide, hla, qual))

    print(hla_list)
