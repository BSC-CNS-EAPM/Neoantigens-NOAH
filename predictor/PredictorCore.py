from NOAH.constants.constants import *
import numpy as np


class PredictorCore:
    def __init__(self, hla_list, motif_length, key_positions, sim_weight, valid_letters=VALID_AMINOACIDS):
        """
        Abstract class with all the base properties of the Motif
        :param hla_list: list of hlas to use
        :param motif_length: length of the motif
        :param key_positions: dictionary of the aa positions that define the environment
        :param sim_weight: dictionary of the weight that each environment position has based on crystallographic data
        :param valid_letters: list of valid amino acids
        """

        # Base Parameters
        self.hla_list = hla_list
        self.valid_letters = valid_letters
        self.motif_length = motif_length

        # Constants
        self.similarity_matrix = None
        self.similarity_weight = sim_weight
        self.key_positions = key_positions

        # Mapping
        self.hla_to_env = None  # {pos:{hla:env}}  None for the abstract class
        self.env_to_hla = None  # {pos:{hla:[hlas_to_use]}}
        self.letters_to_nums = {letter: i for i, letter in enumerate(valid_letters)}  # {letter: num}
        self.hla_to_num = {hla: i for i, hla in enumerate(hla_list)}  # {hla: num}

    def set_similarity_matrix(self, matrix):
        self.similarity_matrix = matrix

    def compare_hla_envs(self, hla, function_to_use):
        """
        Function that compares the binding environment of a given HLA with the environments
        of all the other reported HLA's.
        :param hla: HLA to analyse
        :param function_to_use: criteria to use to determine the similarity of the environments
        :return: sorted list of the environments by similarity
        """
        hla_similarities = {}
        for position in range(self.motif_length):
            hla_similarities.setdefault(position, [])
            for hla_2 in self.hla_list:
                if hla_2 != hla:
                    similarity = function_to_use(hla, hla_2, position)
                    hla_similarities.setdefault(position, []).append((similarity, hla_2))
            hla_similarities[position].sort(reverse=True)
        return hla_similarities

    def compare_two_global_environemnts(self,  hla_1, hla_2, position):
        """
        Compares two HLAs taking into account all the environments used to build the model (the hla sequence and the ones fused to it)
        Note: This is one of the posible functions used by compare hla envs
        :param hla_1: hla that you want to compare (str)
        :param hla_2: second hla that you want to compare (str)
        :param position: position of the motif that you want to compare
        :return: similarity of the given position
        """
        env_1 = self.hla_to_env[position][hla_1]
        similarity_list = []
        for hla in self.env_to_hla[position][hla_2]:
            env_2 = self.hla_to_env[position][hla]
            similarity = np.array(list(map(self.compare_one_letter, env_1, env_2)), dtype=np.float64)
            position_weight = np.array([self.similarity_weight[position][x] for x in self.key_positions[position]], dtype=np.float64)
            position_weight /= np.sum(position_weight)
            similarity *= position_weight
            similarity = np.sum(similarity)
            similarity_list.append(similarity)
        similarity = sum(similarity_list)/len(similarity_list)
        return similarity

    def compare_two_HLAs(self, hla_1, hla_2, position):
        """
        Compares two HLAs using only their sequences
        Note: This is one of the posible functions used by compare hla envs
        :param hla_1: hla that you want to compare (str)
        :param hla_2: second hla that you want to compare (str)
        :param position: position of the motif that you want to compare
        :return: similarity of the given position
        """
        env_1 = self.hla_to_env[position][hla_1]
        env_2 = self.hla_to_env[position][hla_2]
        similarity = np.array(list(map(self.compare_one_letter, env_1, env_2)), dtype=np.float64)
        position_weight = np.array([self.similarity_weight[position][x] for x in self.key_positions[position]], dtype=np.float64)
        position_weight /= np.sum(position_weight)
        similarity *= position_weight
        similarity = np.sum(similarity)
        return similarity

    def compare_one_letter(self, letter_1, letter_2):
        # Compares two amino acids using a similarity matrix (higher numbers mean more similarity)
        try:
            value = self.similarity_matrix[letter_1][letter_2]
        except KeyError:
            print("WARNING UNKNOW LETTERS")
            value = 0
        return value

    def extract_binding_environment(self, aligment_dict):
        """
        Method that given a dictionary with the hla and their sequence extracts the aa that are relevant for the binding
        based on the key positions file
        :param aligment_dict: dictionary with the HLA's and their sequences
        :return:dict {position: {hla: aa_environment}}
        """
        hla_to_env = {}
        for hla, sequence in aligment_dict.items():
            data = self.find_environment(sequence)
            for i in range(self.motif_length):
                hla_to_env.setdefault(i, {}).setdefault(hla, data[i])
        return hla_to_env

    def find_environment(self, sequence):
        # function that maps the positions of the key positions file to the alignment
        envs = {}
        for i in range(self.motif_length):
            env_string = ""
            for key_position in self.key_positions[i]:
                env_string += sequence[key_position]
            envs[i] = env_string
        return envs

    def prepare_env_to_hla(self):
        # method that initiazes the dictionary required for the mapping
        env_to_hla = {}
        for i in range(self.motif_length):
            for hla in self.hla_list:
                env_to_hla.setdefault(i, {}).setdefault(hla, [hla])
        return env_to_hla


