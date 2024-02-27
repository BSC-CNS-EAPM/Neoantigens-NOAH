import os
import test

import data

DATA_PATH = os.path.dirname(data.__file__)
TEST_PATH = os.path.dirname(test.__file__)
CONSTANTS_PATH = os.path.dirname(__file__)

HLA_ALIGMENT_FILE = os.path.join(DATA_PATH, "HLA-A.pfam")

# BASE_DATA_FILE = os.path.join(DATA_PATH, "data.csv")

IEDB_DATA_FILE = os.path.join(DATA_PATH, "IEDB_data.csv")

POSITIVE_HIGH = "Positive-High"

POSITIVE_INTERMEDIATE = "Positive-Intermediate"

POSITIVE = "Positive"

NEGATIVE = "Negative"

POSITIVE_LOW = "Positive-Low"

MOTIFF_BASE_LENGTH = 9

TEST_DATA = os.path.join(TEST_PATH, "data_to_predict.txt")

VALID_AMINOACIDS = [
    "A",
    "C",
    "E",
    "D",
    "G",
    "F",
    "I",
    "H",
    "K",
    "M",
    "L",
    "N",
    "Q",
    "P",
    "S",
    "R",
    "T",
    "W",
    "V",
    "Y",
]

# minimum structural crystal evidence
THRESHOLD = {8: 0, 9: 0, 10: 0, "default": 0}

# similarity matrices
# (path, bool->"whether they have to be inverted or not. if true all the value will be multiplied by -1"
BLOSUM62 = (os.path.join(CONSTANTS_PATH, "blosum62.csv"), False)
PAM250 = (os.path.join(CONSTANTS_PATH, "pam250.csv"), False)
GRANTHAMS = (os.path.join(CONSTANTS_PATH, "granthams_distance.csv"), True)
SNEATH = (os.path.join(CONSTANTS_PATH, "sneath_index.csv"), True)
SIMILARITY_DICT = {
    "blosum62": BLOSUM62,
    "pam250": PAM250,
    "granthams": GRANTHAMS,
    "sneath": SNEATH,
}
