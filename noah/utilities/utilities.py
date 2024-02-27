import multiprocessing as mp
import pickle
import sys

import numpy as np


def load_model(file_path):
    try:
        with open(file_path, "rb") as inn:
            motif = pickle.load(inn)
        return motif
        # file = open(file_path, "rb")
        # motif = pickle.load(file)
        # file.close()
        # return motif
    except IOError:
        raise Exception(
            "ERROR UNABLE TO LOAD MOTIF, try using another python version\n"
        )


def load_data(file_path):
    data = []
    data_qual = {}
    try:
        with open(file_path, "r") as inn:
            for line in inn:
                if line.startswith("peptide,HLA"):
                    continue
                line = line.rstrip()
                line = line.split(",")
                peptide = line[0]
                hla = line[1]
                data.append((peptide, hla))
                if len(line) == 3:
                    qualitativeValue = line[2]
                    data_qual.setdefault(hla, {}).setdefault(peptide, qualitativeValue)
    except IOError:
        raise Exception("Error: input file not found\n")
    if data_qual:
        return data, data_qual
    else:
        return data


def process_peptides(motif, data):
    results = {}
    for element in data:
        peptide = element[0]
        hla = element[1]
        try:
            result = motif.score_peptide(peptide, hla)
            for processed_hla in result:
                results.setdefault(processed_hla, {}).setdefault(
                    peptide, result[processed_hla][peptide]
                )
        except:
            raise Exception(
                "Error ocurred while processing %s for HLA %s\n" % (peptide, hla)
            )
    return results


def score_peptides_paralleled(processors, loaded_data, motif):
    pool = mp.Pool(processors)
    workers = []
    splited_data = np.array_split(loaded_data, processors)
    output_data = {}
    for i in range(processors):
        workers.append(pool.apply_async(process_peptides, (motif, splited_data[i])))
    for worker in workers:
        result = worker.get()
        for hla in result:
            for peptide in result[hla]:
                output_data.setdefault(hla, {}).setdefault(
                    peptide, result[hla][peptide]
                )
    pool.terminate()
    return output_data


def print_motif(motif, hla):
    motif_length = motif.motif_length
    motif_dict = {}
    for i in range(motif_length):
        hla_num = motif.hla_to_num[motif.env_to_hla[i][hla][0]]
        for letter in motif.valid_letters:
            value = motif.likelihood_matrix[i][hla_num][motif.letters_to_nums[letter]]
            motif_dict.setdefault(i, {}).setdefault(letter, value)
    return motif_dict
