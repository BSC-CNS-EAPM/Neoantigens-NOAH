import numpy as np
import multiprocessing as mp
import copy
from NOAH.constants.constants import *
from NOAH.predictor import Scorer
from NOAH.predictor.PredictorCore import PredictorCore

# ignore warning of 0 division errors (infinit positions are expected if no pseudocounts are used)
np.seterr(divide='ignore')


class MotifMaker(PredictorCore):
    """
    Class That constructs and refines the model
    """

    def __init__(self, data, test_data, hla_list, motif_length, key_positions, sim_weight, hla_aligment,
                 valid_letters=VALID_AMINOACIDS, pseudocounts=1):
        """

        :param data: data to use to build the model
        :param test_data: data to use to test the model
        :param hla_list: list of valid hlas
        :param motif_length: length of the motif
        :param key_positions: dictionary of the aa positions that define the environment
        :param sim_weight: dictionary of the weight that each environment position has based on crystallographic data
        :param hla_aligment: Alignment of the different provided HLAs
        :param valid_letters: list of valid amino acids
        :param pseudocounts: value to use when using pseudocounts
        """

        # Initialize Model
        PredictorCore.__init__(self, hla_list, motif_length, key_positions, sim_weight, valid_letters)

        # Base Parameters
        self.data = data  # {qualitative_value : {hla:[peptides]}}
        self.test_data = test_data  # {hla: {qualitative_value : [peptides]}}
        self.pseudocounts = pseudocounts
        self.hla_aligment = hla_aligment
        self.total_hla = len(self.hla_list)
        self.motif_length = motif_length

        # Background Model Parameters
        self.random_model_type = "negative"
        self.valid_random = ["random", "all", "negative", "unique", "fused"] # already implemented random models
        self.random_model = None

        # Mappings
        self.hla_to_env = self.extract_binding_environment(self.hla_aligment)  # {pos:{hla:env}}
        self.env_to_hla = self.prepare_env_to_hla()   # {pos:{hla:[hlas_to_use]}}

        # Model Parameters
        self.hla_data = {}
        self.count_matrix = None
        self.freq_matrix = None
        self.likelihood_matrix = None
        self.environments_mapping = {}

    def set_random_model_type(self, random_type):

        if random_type in self.valid_random:
            self.random_model_type = random_type
        else:
            print("Random model type %s not found. Valid random types: %s" % (random_type, " ".join(self.valid_random)))
            print("Using the default one")

    def set_random_model(self, random_model):
        self.random_model = random_model

    def set_data(self, data):
        self.data = data

    def _compute_random_model(self):
        # Private method that computes the random model to use
        final_random_model = None
        self.random_model += 1.0  # Add pseudocounts to avoid 0 division
        if self.random_model_type in ["all", "negative"]:
            # Background model using all the data available
            final_random_model = np.sum(self.random_model, axis=0)  # Add all positions
            final_random_model = np.sum(final_random_model, axis=0)  # Add all envs
            final_random_model /= np.sum(final_random_model)  # Compute frequencies using Total

        elif self.random_model_type == "unique":
            # Background model with only the frequencies for its HLA
            final_random_model = np.sum(self.random_model, axis=0)
            for hla_matrix in final_random_model:
                hla_matrix /= np.sum(hla_matrix)

        elif self.random_model_type == "random":
            # Background model using same frequencies for all the amino acids
            final_random_model = np.array([(1.0)/len(self.valid_letters) for x in range(len(self.valid_letters))])

        elif self.random_model_type == "fused":
            # Background model using the frequencies of all the used HLA to build the model
            final_random_model = np.sum(self.random_model, axis=0)
            final_random_model_copy = copy.deepcopy(final_random_model)
            for hla in self.hla_list:
                index_1 = self.hla_to_num[hla]
                hla_in_env = set()
                for position in range(self.motif_length):
                    for hla_2 in self.env_to_hla[position][hla]:
                        if hla_2 != hla:
                            hla_in_env.add(hla_2)

                for hla_2 in hla_in_env:
                    index_2 = self.hla_to_num[hla_2]
                    final_random_model[index_1] += final_random_model_copy[index_2]
                final_random_model[index_1] /= np.sum(final_random_model[index_1])

        return final_random_model

    def build(self):
        # Method that builds the model
        self.freq_matrix = self._create_matrix_skeleton()
        self._compute_frequencies()
        final_random_model = self._compute_random_model()
        self.likelihood_matrix = self._compute_likelihood(self.freq_matrix, final_random_model)
        return self._instantiate_model()

    def initialize(self):
        # Method that initializes all the required numpy matrices
        self.count_matrix = self._create_matrix_skeleton()
        self.random_model = self._create_matrix_skeleton()

        if self.random_model_type in ["negative", "unique", "fused"]:
            self._count(self.data[NEGATIVE], self.random_model)

        elif self.random_model_type in ["all"]:
            all_data = []
            for qualitative_value, data in self.data.items():
                all_data += data
            self._count(all_data, self.random_model)

        self._count(self.data[POSITIVE_HIGH], self.count_matrix, verbose=True)
        # Empty the data to decrease RAM usage
        self.data = None

    def _create_matrix_skeleton(self):
        # Private method that computes the size of the multidimensional array that will hold the data
        r_1 = self.motif_length
        r_2 = self.total_hla
        r_3 = len(self.valid_letters)
        matrix = np.zeros(shape=(r_1, r_2, r_3), dtype=float)
        return matrix

    def _count(self, data, matrix, verbose=False):
        """
        Private method that counts the number of entries
        :param data: data to count
        :param matrix: matrix to load the data into
        :param verbose: Boolean, whether to print missing data warnings or not
        :return:
        """
        hla_warnings = set()
        for hla in data:
            for peptide in data[hla]:
                for i, letter in enumerate(peptide):
                    try:
                        matrix[i][self.hla_to_num[hla]][self.letters_to_nums[letter]] += 1
                    except KeyError:
                        # Hla without enough data
                        hla_warnings.add(hla)
                        continue
        if hla_warnings and verbose:
            print("Warning: Hlas without enough data: ", list(hla_warnings))

    def _compute_frequencies(self):
        # Computes the frequencies
        matrix = self.count_matrix + self.pseudocounts
        for i in range(self.motif_length):
            for hla in self.hla_list:
                for uhla in self.env_to_hla[i][hla]:
                    uhla_num = self.hla_to_num[uhla]
                    self.freq_matrix[i][self.hla_to_num[hla]] += matrix[i][uhla_num]
                self.freq_matrix[i][self.hla_to_num[hla]] /= np.sum(self.freq_matrix[i][self.hla_to_num[hla]])

    @staticmethod
    def _compute_likelihood(matrix_1, matrix_2):
        final_matrix = np.log2(matrix_1/matrix_2)
        final_matrix *= -1  # Because pep likes negatives values as the good ones
        return final_matrix

    def _instantiate_model(self):
        return Scorer.Scorer(self.hla_to_env, self.env_to_hla, self.hla_list, self.motif_length, self.key_positions,
                             self.similarity_weight, self.valid_letters, self.likelihood_matrix)

    def compare_all_envs(self):
        # compares all the environments with each other
        hla_dict = {}
        for hla in self.hla_list:
            hla_dict[hla] = self.compare_hla_envs(hla, self.compare_two_HLAs)
        return hla_dict

    def _compare_models(self, motif, hla, hla_dict, threshold=-1):
        # Builds and compares two models with different environments to select the best one
        good_fusions = {x: {} for x in range(self.motif_length)}
        refined_motif = copy.deepcopy(motif)
        original_model = motif.build()
        original_model.hla_list = [hla]
        MCC_1 = original_model.score_MCC(threshold, self.test_data)
        for i in range(self.motif_length):
            similarity = self.compare_two_HLAs(hla, hla, i)
            #self.env_to_hla = copy.deepcopy(old_env)
            for similarity_tuple in hla_dict[hla][i]:
                new_similarity = similarity_tuple[0]
                refined_motif.env_to_hla = copy.deepcopy(self.env_to_hla)
                refined_motif.env_to_hla[i][hla].append(similarity_tuple[1])
                refined_model = refined_motif.build()
                refined_model.hla_list = [hla]
                MCC_2 = refined_model.score_MCC(threshold, self.test_data)
                if (similarity == new_similarity and MCC_1 - MCC_2 < 0.1) or MCC_1 - MCC_2 < -0.1:
                    good_fusions.setdefault(i, {}).setdefault(hla, [hla]).append(similarity_tuple[1])
                else:
                    continue
        return good_fusions

    def refine_model(self, processors=1):
        # Compares the HLA envs to determine which ones to fuse, based on a loss function
        hla_dict = self.compare_all_envs()
        pool = mp.Pool(processors)
        workers = []
        print("Refining Model")
        for hla in self.hla_list:
            workers.append(pool.apply_async(self._compare_models, (copy.deepcopy(self), hla, hla_dict)))
        for worker in workers:
            good_fusions = worker.get()
            for i in range(self.motif_length):
                self.env_to_hla[i].update(good_fusions[i])
        print("building final model")
        pool.terminate()
        model = self.build()
        model.set_similarity_matrix(self.similarity_matrix)
        return model


