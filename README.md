# NOAH

To make the code work, the IEDB database should be putted into data/ with the name IEDB_data.csv

To use the code install it with pip or pip3:
    
    pip install .
    

## Algorithm

The module core is divided in the following parts:
1. Constants -> Variables that are constant and should not change, such as the valid amino acids or the similarity matrices

2. hlaizer/parser -> class that does the parsing of the text files that are required

3. predictor/PredictorCore -> class that contains the methods that both the model builder and the scorer require to work properly

4. predictor/scorer -> class that contains the function that scores the data using an existing motif.

5. predictor/Model_builder -> class that contains the methods that build the motif

## How Noah works

All the methods are documented inside the code. This is an overview of the most important parts of it.

The model builder uses a numpy matrix to hold all the likelihoods and to allow faster processing speeds of the model. 
This matrix has 3 dimensions: positions -> different hlas -> aminoacids.

You can access to any likelihood value by using the correct indexes.

To know which index corresponds to each hla, and to know which binding environment corresponds to each hla, the code
contains variables with the mappings (correspondences)

This variables are:

   hla_to_env -> dictionary that maps the hla names to their environments {position:{hla:environment}}
    
   env_to_hla -> dictionary that contains which hlas are similar {position:{hla:[list of similar hlas]}}
   it is important to take into account that this dictionary is not symmetric, so for example: if HLA-A-01:01 is similar to HLA-B-02:02, it dos not mean that HLA-B-02;02 is similar to HLA-A-01:01.
   
   hla_to_num -> maps the hla (string) to its index in the numpy matrix
   
   letters_to_nums -> maps the aminoacids to their corresponding index in the numpy matrix
   
So, basically, if you want to access to a given likelihood value you have to acces the numpy matrix in the following way:
     
     self.likelihood_matrix[position][self.hla_to_num[hla]][self.letters_to_nums[aminoacid]]
     
To fuse the model, NOAH does the following:
First it builds a model of all the hla individually. Then for each hla it creates a list of all the hlas sorted by their similarity, and finally it starts building a new model joining the data following a predefined criteria.

Right now the criteria is the following: if the sequence for that position is equal, NOAH fuses both models if there is a gain in MCC, if they are not different, NOAH only fuses  them if the gain is bigger than 0.1

### *Denovo* prediction
To do *denovo* prediction what noah does is just to compare the new sequence with the already known ones. Noah doesn't build or computes anything new (because there is no data). Noah just uses an average of the most similar known hlas for each position.

### HOW TO RUN NOAH

#### 1. Running tests:

In the main folder there is an script called run_tests.py. At the bottom of the script (under if __name__ == "__main__":)
you can change some parameters as well as specify which tests do you want to run. Then you just have to call the script.

######

#### 2. Creating a Model:

To automatically build a model you can use the train_NOAH.py providing the needed arguments. you can always check the available arguments with --help
the arguments that you can use are the following:
    
    # Required arguments
    -o : Name for the output model (without file extension)
    --length : Length of the peptides to use to build the model

    # Optional arguments
    --iedb : Path to the IEDB data file. (if you followed the instructions of how to install NOAH it is not required)
    --processors : Number of processors to use. (default is no multiprocessing)
    --data : Path to an additional data file. The file must be a csv with the following format: peptide;hla;qualitative_Value.
    --alignment : Path to the file with the HLAs (default is the file named HLA.pfam inn the data folder) (Selex format)
    --background : Background model to use. Options are : (random, all, negative, unique, fused).fused is the recommended one and also the default
    --signal : Minimum number of peptides with good binding that an HLA must have to be modeled. default 50.
    --noise : Minimum number of non-binding peptides that an HLA must have to be modelled. default 10
    --simMatrix : Similarity matrix to use. Options are: [blosum62, pam250, granthams, sneath]. default is sneath. I do not remember why, the latests test were performed with GRANTHAMS. However the matrix used does not have a big impact on the results.

>[!TIP]
>The command should look similar to:
> 
>     python noah/train_NOAH.py -o model_name --length 9

At the end you will have a pickel with the model.



#####

#### 3. Scoring a list of peptides

You have two options, you can use the script that scores the peptides automatically or you could import noah and do it manually.
The first option is automatic, (how the code is meant to be used). The second option is more flexible but requires you to know what you are doing.

##### Using the script:
you can use the script main_NOAH.py in command line with the required arguments. I have not tested  this script extensively and there has been a lot of modifications of the code so it would be good to use it more with different cases to see if it still works.
The arguments are the following: (you can always check them with --help)
    
    # Required arguments
    -i : File with the peptides to Predict. The file must have the following structure; peptide  HLA
    -o : Output file. (with file extension) the output is always a csv
    -model : Path to the model to use
    # Optional arguments
    -seq : File with the proteic sequences for the unknown HLAs (Selex format) (right now you must give a selex file if there is any HLA not modelled in your list, pending to be changed)
    -processors : Number of processors to use, default 1

>[!TIP]
>The command should look similar to:
> 
>     python noah/main_NOAH.py -i path_input_csv -o name_output.csv -model path_to_the_model

##### using NOAH directly:
first you must load the model using the utilities module of NOAH:

    import noah.utilities
    my_model = utilities.load_model("PATH/TO/MODEL/model.pkl")

next you have to load the sequences that are not part of the model but should be recognized (you can give just give all of them it doe snot matter if they are already known)

##### to load the file

    from hlaizer.parser import Parser 
    sequence_parser = Parser()
    hla_alignment = sequence_parser.parse_aligment_file("PATH/TO/ALIGMENT_FILE_WITH_SELEX_FORMAT.pfam")
##### to add the sequences to the model
    my_model.load_sequences(hla_alignment)

Finally you can just call the score function which will handle the rest (including peptides with different length).

If you don't want to deal with formatting you can use the function:
process_peptides from utilities which will handle the output formatting for you.

NOAH.utilities.process_peptides(my_model, data) # data is a list of tuples with (peptide, hla)

The output of this function is a dictionary with the following levels: {hla: {peptide: score}}

If you want you can also call the scorer function and score each peptide one by one.
    
    my_model.score_peptide(peptide, hla)

NOTE: Here what hla means is a bit tricky. If HLA is an string the peptide will be evaluated for this hla only. However if the string is part of an HLA identifier it will process all the hlas that match that string.
For example, if you write HLA-A. the peptide will be evaluated for all the known HLA-A's. (but not for the unknown hlas)
as a comment we may want to remove this feature or if we want to keep it we should also support the unknown ones ( the ones that do not have a model).
Finally, if hla is a list, the peptide will be evaluated for all the hla in that list. And again if one of the hla in that list is not a full hla identifier, all the hla that match that part of the identifier will be used which could lead to crazy outputs depending on how much repetitions you have.
This is maybe more reason to remove this feature. (this also applies to utilities.process_peptides)

The output of the function is a dictionary with the following levels: {hla: {peptide: score}}
Since this function can not process more than one peptide at once you should consider creating you own output handling method or just used the function in utilities.

Alternatively, if you do not want to deal with this none sense of how crazy the hla system works you can use the internal "private" scoring functions.
However, this requires more extra steps:

###### Processing a known peptide:

    my_model._score_rounded(peptide, hla)

Here hla is an string with one identifier, the output is just the score

processing an unknown peptide

first you have to prepare the *deNOVO* functionality with:

    my_model._prepare_for_denovo(hla)

Here hla is also an string matching only one identifier. this function returns nothing and should be called only once for each new hla.
Afterwards, you can score any number of peptides with that hla without calling _prepare_for_denovo again.

> [!WARNING]
> Do not call _prepare_for_denovo if the hla is part of the hla used to create the model or it will crash!

You can find which HLAs have a model assigned with: my_model.hla_list  which is a list of hla identifiers
You can find which unknown HLAS have been already prepared with: my_model.unknown_hla_map.keys() which is also a list of hla identifiers

Finally if you want to process peptides with different length than the model, you should use

    my_model._diffsize_score(peptide, hla)

Again hla is only one identifier, if the hla is not part of the hlas with a model assigned you should call _prepare_for_denovo before as previously explained.
The output is also an score. (right now is the average of all the possible binding modes that the peptide can have)



--------
Thanks to Roc Farriol and Albert Cañellas for testing and bugs solving
