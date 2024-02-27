import argparse
import multiprocessing as mp
import sys

import numpy as np
import utilities
from hlaizer.parser import Parser


def parse_args():
    """
    Parse command line arguments
    :returns:
    """
    desc = """Main Script that loads a given model and a file with the peptides to predict. 
    Check the Documentation "README.md" for more information about the usage"""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-i",
        required=True,
        help="File with the peptides to Predict. The file must have the following"
        "structure; peptide  HLA",
    )
    parser.add_argument(
        "-seq",
        default=None,
        help="File with the proteic sequences for the unknown HLAs (Selex format)",
    )
    parser.add_argument("-o", required=True, help="Output file")
    parser.add_argument(
        "-model", required=True, help="Path to where the models are stored"
    )
    parser.add_argument("-processors", default=1, help="Number of processors to use")
    args = parser.parse_args()
    return args.i, args.seq, args.o, args.model, args.processors


def process_peptides(model, data):
    """
    Scores each peptide in the data array
    :param model: model to use to score the peptides
    :param data: list of tuples containing the data to score (peptide, HLA)
    :return: scored data
    """
    results = {}
    for element in data:
        print(data)
        peptide = element[0]
        hla = element[1]
        try:
            result = model.score_peptide(peptide, hla)
            for processed_hla in result:
                results.setdefault(processed_hla, {}).setdefault(
                    peptide, result[processed_hla][peptide]
                )
        except:
            sys.stderr.write(
                "Error ocurred while processing %s for HLA %s\n" % (peptide, hla)
            )
    return results


def main(input_file, hla_seq, output, model, processors):
    print("Starting NOAH")
    try:
        scorer = utilities.load_model(model)
    except:
        raise Exception("Error: Unable to load model %s\n" % model)

    if hla_seq:
        parser = Parser()
        hla_align = parser.parse_aligment_file(hla_seq)
        scorer.load_sequences(hla_align)

    pool = mp.Pool(processors)
    workers = []
    loaded_data = utilities.load_data(input_file)
    splited_data = np.array_split(loaded_data, processors)
    output_data = {}

    for i in range(processors):
        workers.append(pool.apply_async(process_peptides, (scorer, splited_data[i])))
    for worker in workers:
        result = worker.get()
        for hla in result:
            for peptide in result[hla]:
                output_data.setdefault(hla, {}).setdefault(
                    peptide, result[hla][peptide]
                )

    print("Saving results")
    try:
        file = open(output, "w")
    except IOError:
        sys.stderr.write(
            "WARNING: Can't open outputfile %s, using _tmp_result.txt instead\n"
            % output
        )
        file = open("_tmp_result.txt", "w")

    for hla in output_data:
        for peptide in output_data[hla]:
            file.write("%s\t%s\t%s\n" % (hla, peptide, output_data[hla][peptide]))
    file.close()
    print("Prediction finished")


if __name__ == "__main__":
    # main("data_proba.txt", os.path.join(DATA_PATH, "HLA-A.pfam"), "resu_random.txt",
    # os.path.join(DATA_PATH, "NOAH_9.pkl"), 1)
    input_file, hla_seq, output, models, processors = parse_args()
    main(input_file, hla_seq, output, models, processors)
