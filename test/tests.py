from predictor.Model_builder import MotifMaker
from utilities import utilities as utilities
import copy
from hlaizer.parser import Parser


def IEDB_prediction_test(data, test_data, hla_list, hla_aligment, key_positions, sim_weight, data_to_score, processors,
                         hla_seq, similarity_matrix):
    """
    Test performance of the model
    :param data: data to use to build the model
    :param test_data: test data to use to build the model
    :param hla_list: list of hlas for which the model will be build
    :param hla_aligment: hla aligment
    :param key_positions: dict with the sequences positions that define the binding environment
    :param sim_weight: weight that each position of the binding environment has
    :param data_to_score: data to score to test performance
    :param processors: processors to use
    :param hla_seq: file with the sequences of the HLAs
    :param similarity_matrix: similarity matrix to use to compare binding environments
    :return: predictions
    """
    parser = Parser()
    print("Making Model with fused background")
    motif = MotifMaker(data=data, test_data=test_data, hla_list=hla_list, motif_length=9,
                       key_positions=key_positions, sim_weight=sim_weight, hla_aligment=hla_aligment)
    motif.set_similarity_matrix(similarity_matrix)
    motif.set_random_model_type("fused")
    motif.initialize()
    motif.build()
    model_fused = motif.refine_model(int(processors))
    print("Making Model with Fused background")
    motif.set_random_model_type("unique")
    model_unique = motif.build()
    print("Making Model with Unique background")
    hla_align = parser.parse_aligment_file(hla_seq)
    model_fused.load_sequences(hla_align)
    model_unique.load_sequences(hla_align)
    prediction_negative = utilities.score_peptides_paralleled(processors, data_to_score, model_fused)
    prediction_unique = utilities.score_peptides_paralleled(processors, data_to_score, model_unique)
    return prediction_negative, prediction_unique


def deNovo_prediction_test(data, test_data, hla_list, hla_aligment, key_positions, sim_weight, data_to_score, processors,
                           hla_seq, similarity_matrix):
    """
    Test deNovo performance of the model
    :param data: data to use to build the model
    :param test_data: test data to use to build the model
    :param hla_list: list of hlas for which the model will be build
    :param hla_aligment: hla aligment
    :param key_positions: dict with the sequences positions that define the binding environment
    :param sim_weight: weight that each position of the binding environment has
    :param data_to_score: data to score to test performance
    :param processors: processors to use
    :param hla_seq: file with the sequences of the HLAs
    :param similarity_matrix: similarity matrix to use to compare binding environments
    :return: predictions
    """
    deNovo_fused = {}
    deNovo_unique = {}
    data_to_score_formated = {}
    for element in data_to_score:
        peptide = element[0]
        hla = element[1]
        data_to_score_formated.setdefault(hla, []).append(element)
    print("Making deNovo Motifs")
    for hla in hla_list:
        if hla not in data_to_score_formated:
            print("HLA  %s not present in the scoring list, skipping allele" % hla)
            continue
        print("Making deNovo Motif for HLA: %s" % hla)
        tmp_data = copy.deepcopy(data)
        tmp_test_data = copy.deepcopy(test_data)
        tmp_hla_list = copy.deepcopy(hla_list)
        for qualvalue in tmp_data:
            tmp_data[qualvalue].pop(hla, None)
        tmp_test_data[hla].pop(hla, None)
        tmp_hla_list.remove(hla)
        motif = MotifMaker(data=tmp_data, test_data=tmp_test_data, hla_list=tmp_hla_list, motif_length=9,
                           key_positions=key_positions, sim_weight=sim_weight, hla_aligment=hla_aligment)
        motif.set_similarity_matrix(similarity_matrix)
        motif.set_random_model_type("fused")
        motif.initialize()
        motif.build()
        model_fused = motif.refine_model(int(processors))
        motif = MotifMaker(data=tmp_data, test_data=tmp_test_data, hla_list=tmp_hla_list, motif_length=9,
                           key_positions=key_positions, sim_weight=sim_weight, hla_aligment=hla_aligment)
        motif.set_similarity_matrix(similarity_matrix)
        motif.set_random_model_type("unique")
        motif.initialize()
        motif.build()
        model_unique = motif.refine_model(int(processors))
        print("Making deNovo prediction for HLA: %s" % hla)
        parser = Parser()
        hla_align = parser.parse_aligment_file(hla_seq)
        model_fused.load_sequences(hla_align)
        model_unique.load_sequences(hla_align)
        prediction_fused = utilities.score_peptides_paralleled(processors, data_to_score_formated[hla], model_fused)
        prediction_unique = utilities.score_peptides_paralleled(processors, data_to_score_formated[hla], model_unique)
        deNovo_fused.setdefault(hla, prediction_fused[hla])
        deNovo_unique.setdefault(hla, prediction_unique[hla])

    return deNovo_fused, deNovo_unique


def Individual_prediction_test(data, test_data, hla_list, hla_aligment, key_positions, sim_weight, data_to_score,
                               processors, hla_seq, similarity_matrix):
    """
    Test performance of the model without merging environments
    :param data: data to use to build the model
    :param test_data: test data to use to build the model
    :param hla_list: list of hlas for which the model will be build
    :param hla_aligment: hla aligment
    :param key_positions: dict with the sequences positions that define the binding environment
    :param sim_weight: weight that each position of the binding environment has
    :param data_to_score: data to score to test performance
    :param processors: processors to use
    :param hla_seq: file with the sequences of the HLAs
    :param similarity_matrix: similarity matrix to use to compare binding environments
    :return: predictions
    """
    print("Making Individual Model with Fused background")
    motif = MotifMaker(data=data, test_data=test_data, hla_list=hla_list, motif_length=9,
                       key_positions=key_positions, sim_weight=sim_weight, hla_aligment=hla_aligment)
    motif.set_similarity_matrix(similarity_matrix)
    motif.set_random_model_type("fused")
    motif.initialize()
    model_fused = motif.build()
    print("Making Individual Model with Unique background")
    motif = MotifMaker(data=data, test_data=test_data, hla_list=hla_list, motif_length=9,
                       key_positions=key_positions, sim_weight=sim_weight, hla_aligment=hla_aligment)
    motif.set_similarity_matrix(similarity_matrix)
    motif.set_random_model_type("unique")
    motif.initialize()
    model_unique = motif.build()
    print("Making predictions")
    parser = Parser()
    hla_align = parser.parse_aligment_file(hla_seq)
    model_fused.load_sequences(hla_align)
    model_unique.load_sequences(hla_align)
    prediction_fused = utilities.score_peptides_paralleled(processors, data_to_score, model_fused)
    prediction_unique = utilities.score_peptides_paralleled(processors, data_to_score, model_unique)
    return prediction_fused, prediction_unique
