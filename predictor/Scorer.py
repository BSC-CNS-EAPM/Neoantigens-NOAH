from NOAH.predictor.PredictorCore import PredictorCore
from NOAH.constants.constants import *
import numpy as np
import pickle
import sys


def rounder(function):
    # Decorator that rounds with 3 decimals the returned value
    def wrapper(*args, **kwargs):
        result = np.round(function(*args, **kwargs), 3)
        return result
    return wrapper


def score_MCC(TP, FP, TN, FN):
    # Computes the MCC given a confusion matrix
    try:
        MCC = ((TP * TN) - (FP * FN)) / ((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN)) ** 0.5
    except ZeroDivisionError:
        MCC = ((TP * TN) - (FP * FN))
    return MCC


def confusion_matrix_from_predictet_data(data_scored, data_qual, threshold):
    # Builds a confusion matrix from the predicted data
    count = 0
    TP, FP, TN, FN = 0, 0, 0, 0
    for peptide, score in data_scored.items():
        if score <= threshold:
            if data_qual[peptide] in [POSITIVE_HIGH, POSITIVE_INTERMEDIATE]:
                TP += 1
                count += 1
            elif data_qual[peptide] in [NEGATIVE]:
                FP += 1
                count += 1
        elif score >= threshold:
            if data_qual[peptide] in [POSITIVE_HIGH, POSITIVE_INTERMEDIATE]:
                FN += 1
                count += 1
            elif data_qual[peptide] in [NEGATIVE]:
                TN += 1
                count += 1
    return TP, FP, TN, FN, count


class Scorer(PredictorCore):
    # Trained model

    def __init__(self, hla_to_env, env_to_hla, hla_list, motif_length, key_positions, sim_weight, valid_letters, likelyhood_matrix):
        """
        :param hla_to_env: dict, Mapping of the HLA to their environments
        :param env_to_hla: dict, Mapping of the env to all the HLA associated to it
        :param hla_list: list of valid HLAs
        :param motif_length: length of the motif
        :param key_positions: dictionary of the aa positions that define the environment
        :param sim_weight: dictionary of the weight that each environment position has based on crystallographic data
        :param valid_letters: list of valid amino acids
        :param likelyhood_matrix: numpy ndarray with all the likelihoods
        """
        PredictorCore.__init__(self, hla_list, motif_length, key_positions, sim_weight, valid_letters)

        self.hla_to_env = hla_to_env
        self.likelihood_matrix = likelyhood_matrix
        self.env_to_hla = env_to_hla
        self.unknown_hlas = []
        self.unknown_hla_map = {}  # {hla: {position: likelyhood matrix}}

    def _score(self, peptide, hla, verb=True):
        """
        Private method to score a peptide for a given HLA
        :param peptide: Peptide to score
        :param hla: Hla for which the peptide has to be scored
        :return: loglikelihood score
        """
        score = 0
        # hla_num = self.hla_to_num[hla]
        for i, letter in enumerate(peptide):
            hla_num = self.hla_to_num[self.env_to_hla[i][hla][0]]
            if letter in self.valid_letters:
                score += self.likelihood_matrix[i][hla_num][self.letters_to_nums[letter]]
            elif verb:
                print("WARNING: %s is not a valid character, skipping position %s of peptide %s" % (letter, i, peptide))
        return score

    def _score_unknown(self, peptide, hla, verb=True):
        """
        Private method to score a peptide for a HLA that does not form part of the model
        :param peptide: Peptide to score
        :param hla: Hla for which the peptide has to be scored
        :return: loglikelihood score
        """
        score = 0
        for i, letter in enumerate(peptide):
            if letter in self.valid_letters:
                score += self.unknown_hla_map[hla][i][self.letters_to_nums[letter]]
            elif verb:
                print("WARNING: %s is not a valid character, skipping position %s of peptide %s" % (letter, i, peptide))
        return score

    @rounder
    def _diffsize_score(self, peptide, hla):
        """
        Score peptides for sizes different than the one used to create the motif.
        This method returns the best score base don an sliding window
        :param peptide: Peptide to score
        :param hla: Hla for which the peptide has to be scored
        :return: loglikelihood score
        """
        pept_template_list = []
        extra_length = len(peptide) - self.motif_length
        # create template peptides based on sliding window
        lwindow = max(len(peptide), self.motif_length)
        for position in range(lwindow):
            new_peptide = list(peptide)
            warp_count = 0
            alter_positions = []
            for i in range(abs(extra_length)):
                pos = position + i
                if pos >= lwindow:
                    pos = 0 + warp_count
                    warp_count += 1
                alter_positions.append(pos)
            alter_positions.sort()
            if extra_length > 0:
                new_peptide = list(peptide)
                for j, alter in enumerate(alter_positions):
                    new_peptide.pop(alter - j)
            else:
                for j, alter in enumerate(alter_positions):
                    new_peptide.insert(alter, "X")
            pept_template_list.append("".join(new_peptide))
        # now score each of the peptides
        scores = []
        for n_peptide in pept_template_list:
            scores.append(self._score_rounded(n_peptide, hla, False))
        return min(scores)

    @rounder
    def _score_rounded(self, peptide, hla, verb=True):
        if hla in self.unknown_hla_map.keys():
            return self._score_unknown(peptide, hla, verb)
        else:
            return self._score(peptide, hla, verb)

    @rounder
    def score_MCC(self, threshold, test_data):
        TP, FP, TN, FN = self._build_confusion_matrix(threshold, test_data)
        MCC = score_MCC(TP, FP, TN, FN)
        return MCC

    def _build_confusion_matrix(self, threshold, test_data):
        # Builds a confusion matrix from the predicted data
        # Should be migrated and to utilities module and abstracted
        TP, FP, TN, FN = 0.0, 0.0, 0.0, 0.0
        for hla in test_data:
            if hla in self.hla_list:
                for qualitative_value in test_data[hla]:
                    for peptide in test_data[hla][qualitative_value]:
                        if qualitative_value == POSITIVE:
                            continue
                        score = self._score(peptide, hla)
                        if score <= threshold:
                            if qualitative_value in [POSITIVE_HIGH, POSITIVE_INTERMEDIATE]:
                                TP += 1
                            elif qualitative_value == NEGATIVE:
                                FP += 1
                        elif score > threshold:
                            if qualitative_value in [POSITIVE_HIGH, POSITIVE_INTERMEDIATE]:
                                FN += 1
                            elif qualitative_value == NEGATIVE:
                                TN += 1
        return TP, FP, TN, FN

    def save_pickle(self, name):
        print(name)
        with open(name, 'wb') as inn:
            pickle.dump(self, inn)

    def score_peptide(self, sequence, *args):
        """
        Main method to score peptides.
        It can predict peptides individually and for multiple HLA at the same time.
        If an exact hla name is provided it computes the score for that hla, if a substring is provided,
        all the matching HLA will be used.
        :param sequence: sequence to predict
        :param args: HLA's to use
        :return: dictionary with {hla:peptide:score}
        """
        hla_to_evaluate = []
        for hla in args:
            if type(hla) is list:
                hla_to_evaluate = hla_to_evaluate + hla
                continue
            else:
                hla_to_evaluate.append(hla)
        for hla in hla_to_evaluate:
            unknow_hla = True
            for valid_hla in self.hla_list:
                if hla in valid_hla:
                    unknow_hla = False
            if unknow_hla:
                print("HLA %s not known\nMaking deNovo prediction" % hla)
                self._prepare_for_denovo(hla)
        if len(hla_to_evaluate) == 1:
            result = {}
            if len(sequence) == self.motif_length:
                score = self._score_rounded(sequence, hla_to_evaluate[0])
            else:
                score = self._diffsize_score(sequence, hla_to_evaluate[0])
            result.setdefault(hla_to_evaluate[0], {}).setdefault(sequence, score)
            return result
        elif len(hla_to_evaluate) > 1:
            results = {}
            for hla in hla_to_evaluate:
                if len(sequence) == self.motif_length:
                    score = self._score_rounded(sequence, hla)
                else:
                    score = self._diffsize_score(sequence, hla)
                results.setdefault(hla, {}).setdefault(sequence, score)
            return results
        else:
            print("No HLAs Provided")
            return {}

    def _prepare_for_denovo(self, hla):
        # Loads and prepares the model to be able to do deNovo predictions
        if hla in self.unknown_hlas[0]:
            self.hla_list.append(hla)
            similarities = self.compare_hla_envs(hla, self.compare_two_global_environemnts)
            self.unknown_hla_map.setdefault(hla, self._create_unknownhla_matrix_skeleton())
            for position in range(self.motif_length):
                bestscore = similarities[position][0][0]
                score = similarities[position][0][0]
                count = 0
                while bestscore == score:
                    hla_num = self.hla_to_num[self.env_to_hla[position][similarities[position][0][1]][0]]
                    self.unknown_hla_map[hla][position] += self.likelihood_matrix[position][hla_num]
                    count += 1
                    try:
                        score = similarities[position][count][0]
                    except:
                        break
                self.unknown_hla_map[hla][position] /= count
        else:
            sys.stderr.write("Error no sequence provided for HLA %s\n" % hla)
            sys.stderr.write("Please load a file with the sequences in Selex format\n")
            exit(1)
        return

    def load_sequences(self, aligment_dict):
        # loads new HLAs for the deNovo prediction
        self.unknown_hlas += self.extract_binding_environment(aligment_dict)

    def _create_unknownhla_matrix_skeleton(self):
        # Private method that computes the size of the multidimensional array that will hold the data
        r_1 = self.motif_length
        r_2 = len(self.valid_letters)
        matrix = np.zeros(shape=(r_1, r_2), dtype=float)
        return matrix


