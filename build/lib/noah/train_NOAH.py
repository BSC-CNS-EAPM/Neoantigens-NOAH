import argparse
import os
import sys

from constants.constants import DATA_PATH, SIMILARITY_DICT
from hlaizer.parser import Parser
from predictor.Model_builder import *


def parse_args():
    """
    Parse command line arguments
    :returns:
    """
    desc = """Script that automatically retrains the model with the given data"""
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument(
        "-o", required=True, help="Name for the output model (without file extension)"
    )
    parser.add_argument(
        "--length", required=True, type=int, help="Length of the peptides to use.\n"
    )
    parser.add_argument(
        "--iedb",
        default=os.path.join(DATA_PATH, "IEDB_data.csv"),
        help="Path to the IEDB data file.\n",
    )
    parser.add_argument(
        "--processors", default=1, help="Number of processors to use.\n"
    )
    parser.add_argument(
        "--data",
        default=None,
        help="Path to an additional data file. The file must be a csv with the "
        "following format: peptide;hla;qualitative_Value.\nSee the "
        "documentation for more information about the format.\n",
    )
    parser.add_argument(
        "--alignment",
        default=os.path.join(DATA_PATH, "HLA.pfam"),
        help="Path to the file with the HLAs" " aligment (Selex format)\n",
    )
    parser.add_argument(
        "-b",
        "--background",
        default="fused",
        help="Background model to use. Options are :"
        "(random, all, negative, unique, fused).\nfused "
        "is the recommended one.\n",
    )
    parser.add_argument(
        "--signal",
        default=50,
        type=int,
        help="Minimum number of binding peptides that an HLA must have to be modeled.\n",
    )
    parser.add_argument(
        "--noise",
        default=10,
        type=int,
        help="Minimum number of non-binding peptides that an HLA must have to be modeled.\n",
    )
    parser.add_argument(
        "--simMatrix",
        default="sneath",
        type=str,
        help="Similarity matrix to use. Options are: [blosum62, pam250, granthams, sneath]",
    )
    args = parser.parse_args()
    return (
        args.o,
        args.length,
        args.iedb,
        args.processors,
        args.data,
        args.alignment,
        args.background,
        args.signal,
        args.noise,
        args.simMatrix,
    )


def main(
    output,
    motif_length,
    iedb_data,
    processors,
    data,
    aligment_file,
    background_model,
    signal,
    noise,
    similarity_tuple,
):
    parser = Parser(
        IEDB_file=iedb_data,
        csv_file=data,
        aligment_file=aligment_file,
        length=motif_length,
    )
    parser.set_minimum_data_background(noise)
    parser.set_minimum_data_signal(signal)
    similarity_matrix = parser.load_csv_matrix(similarity_tuple[0], similarity_tuple[1])
    data, test_data, hla_list, hla_aligment, key_positions, sim_weight = (
        parser.parse_all()
    )
    motif = MotifMaker(
        data=data,
        test_data=test_data,
        hla_list=hla_list,
        motif_length=motif_length,
        key_positions=key_positions,
        sim_weight=sim_weight,
        hla_aligment=hla_aligment,
    )
    motif.set_random_model_type(background_model)
    motif.set_similarity_matrix(similarity_matrix)
    motif.initialize()
    motif.build()
    model = motif.refine_model(int(processors))
    model.save_pickle("%s.pkl" % output)
    return 0


if __name__ == "__main__":
    (
        output,
        motif_length,
        iedb_data,
        processors,
        data,
        aligment_file,
        background_model,
        signal,
        noise,
        simMatrix,
    ) = parse_args()
    try:
        sim_tuple = SIMILARITY_DICT[simMatrix]
    except KeyError:
        sys.stderr.write("Error: Invalid similarity matrix type selected\n")
        sys.stderr.write("Valid similarity matrices are:\n")
        sys.stderr.write("%s\n" % " | ".join(SIMILARITY_DICT.keys()))
        exit(1)
    main(
        output,
        motif_length,
        iedb_data,
        processors,
        data,
        aligment_file,
        background_model,
        signal,
        noise,
        sim_tuple,
    )
