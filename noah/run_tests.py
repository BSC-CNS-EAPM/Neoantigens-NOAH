from test import tests

import predictor.Scorer as scorer
from constants.constants import GRANTHAMS, HLA_ALIGMENT_FILE, IEDB_DATA_FILE, TEST_DATA
from hlaizer.parser import Parser
from utilities import utilities


def main(runIEDB, rundeNovo, runIndividual, processors, similarity_tuple, output):
    """
    Script that creates and runs models to test their performance.
    :param runIEDB: test models using the IEDB data
    :param rundeNovo: test deNovo performance
    :param runIndividual: test models without using merged HLAs
    :param processors: number of processors to run
    :param similarity_tuple: similarity tuple to use
    :param output: file to save the results
    :return:
    """

    print("Loading Data")
    parser = Parser(
        IEDB_file=IEDB_DATA_FILE,
        csv_file=None,
        aligment_file=HLA_ALIGMENT_FILE,
        length=9,
    )
    similarity_matrix = parser.load_csv_matrix(similarity_tuple[0], similarity_tuple[1])
    parser.set_minimum_data_background(10)
    parser.set_minimum_data_signal(50)
    data, test_data, hla_list, hla_aligment, key_positions, sim_weight = (
        parser.parse_all()
    )
    data_to_score, data_qual = utilities.load_data(TEST_DATA)
    predictions = []

    if runIEDB:
        pred_fused, pred_unique = tests.IEDB_prediction_test(
            data,
            test_data,
            hla_list,
            hla_aligment,
            key_positions,
            sim_weight,
            data_to_score,
            processors,
            HLA_ALIGMENT_FILE,
            similarity_matrix,
        )
        predictions.append(("IEDB Fused Background", pred_fused))
        predictions.append(("IEDB Unique Background", pred_unique))

    if rundeNovo:
        pred_fused, pred_unique = tests.deNovo_prediction_test(
            data,
            test_data,
            hla_list,
            hla_aligment,
            key_positions,
            sim_weight,
            data_to_score,
            processors,
            HLA_ALIGMENT_FILE,
            similarity_matrix,
        )
        predictions.append(("deNovo Fused Background", pred_fused))
        predictions.append(("deNovo Unique Background", pred_unique))

    if runIndividual:
        pred_fused, pred_unique = tests.Individual_prediction_test(
            data,
            test_data,
            hla_list,
            hla_aligment,
            key_positions,
            sim_weight,
            data_to_score,
            processors,
            HLA_ALIGMENT_FILE,
            similarity_matrix,
        )
        predictions.append(("Individual Fused Background", pred_fused))
        predictions.append(("Individual Unique Background", pred_unique))

    MCC_dict = {}
    hla_predicted = [hla for hla in sorted(predictions[0][1])]
    set_hla_with_errors = set()
    print("Computing MCCs")
    for prediction in predictions:
        model = prediction[0]
        results = prediction[1]
        for hla in hla_predicted:
            try:
                TP, FP, TN, FN, count = scorer.confusion_matrix_from_predictet_data(
                    results[hla], data_qual[hla], -1
                )
                MCC = scorer.score_MCC(TP, FP, TN, FN)
                MCC_dict.setdefault(model, {}).setdefault(hla, MCC)
            except KeyError:
                print("Model %s does not have data for HLA %s" % (model, hla))
                set_hla_with_errors.add(hla)
                MCC_dict.setdefault(model, {}).setdefault(hla, 0)
    if set_hla_with_errors:
        print("HLA without data:", set_hla_with_errors)

    print("Saving results")
    with open(output, "w") as inn:
        hla_string = ";".join(hla_predicted)
        inn.write("MODEL;%s\n" % hla_string)
        for model in MCC_dict:
            MCC_string = ";".join(
                [str(round(MCC_dict[model][hla], 3)) for hla in hla_predicted]
            )
            inn.write("%s;%s\n" % (model, MCC_string))

    print("Test finished")


if __name__ == "__main__":
    # Valid similarity_matrices = BLOSUM62, PAM250, GRANTHAMS, SNEATH
    main(
        runIEDB=True,
        rundeNovo=True,
        runIndividual=True,
        processors=2,
        similarity_tuple=GRANTHAMS,
        output="csv_results_prova.csv",
    )
