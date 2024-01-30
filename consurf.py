# Add software executables here. If they are not in your PATH, then add the full path of the executable
soft_exec_paths = {
    "phmmer":"phmmer",
    "cd-hit":"cd-hit",
    "mafft":"mafft-linsi",
    "rate4site":"rate4site"
}

Alphabet = set(['M', 'I', 'L', 'V', 'A', 'C', 'D', 'E', 'N', 'Q',
                'F', 'W', 'Y', 'G', 'P', 'S', 'T', 'H', 'K', 'R',
                'B', 'Z', 'J',
                'O', 'U', 'X'])

# Imports
import os
import subprocess
import argparse
from collections import defaultdict
from linecache import getline
from concurrent.futures import ProcessPoolExecutor

# Set (global) directory
from pathlib import Path
working_dir = Path(os.getcwd())

def PHMMER_stdout_parser(phmmer_stdout_file, E_val, query_seq, min_seq_identity, min_seq_len):
    """
    Extracts homologous sequences from a PHMMER output file generated with the --notextw parameter.
    Only homologues with an E-value <= E_val, sequence identity >= min_seq_identity and sequence length >= min_seq_len are included.
    Returns a 'phmmer.fasta' file in the current directory containing the sequences in FASTA format.
    
    Partly inspired by Biopython's phmmer parser: https://biopython.org/DIST/docs/api/Bio.SearchIO.HmmerIO.hmmer3_text-pysrc.html#Hmmer3TextParser.__init__
    """    
    min_len = min_seq_len * len(query_seq) # ex: 0.6 * 60 = 36 => homologues must be at least 36 residues long
    output_file = 'phmmer.fasta'
    
    with open(output_directory_path / phmmer_stdout_file, 'rt') as file:
        phmmer_stdout_readlines = file.readlines()
    
    # Here all indexes correspond to line numbers in the phmmer_stdout_file
    # All lines corresponding to a new sequence/entry begin with '>>'
    sequence_index_list = [index for index, line in enumerate(phmmer_stdout_readlines) if line.strip().startswith('>>')]
    
    # All lines corresponding to a homologous domain of a given sequence/entry begins with '==' 
    domain_index_list = [index for index, line in enumerate(phmmer_stdout_readlines) if line.strip().startswith('==')]
    
    # Dictionary with sequences' index and their associated domain indexes. Simpler and faster to parse the index lists in reverse
    sequence_domain_index = defaultdict(list)
    
    for sequence_index in reversed(sequence_index_list):
        for domain_index in reversed(domain_index_list):

            if domain_index > sequence_index:
                sequence_domain_index[sequence_index].append(domain_index)
                del domain_index_list[-1]
            else:
                break
                
    del sequence_index_list, domain_index_list # No longer needed
    
    # Write the output FASTA file
    with open(output_directory_path / output_file, 'wt') as file:
        
        for sequence_index, domain_index_list in sequence_domain_index.items():
            # Only parse the sequence if there is a domain (in some rare cases there is none)
            if domain_index_list:
                
                assert(phmmer_stdout_readlines[sequence_index + 1].strip().startswith('#'))
                assert(phmmer_stdout_readlines[sequence_index + 2].strip().startswith('---'))
                
                # Dictionary with i-Evalue of each domain
                domain_e_val = {}
                for i in range(len(domain_index_list)):
                    line = phmmer_stdout_readlines[sequence_index + 3 + i].split()
                    domain_number, e_val = line[0], line[5]

                    domain_e_val[domain_number] = str(e_val)
                
                # Select domains
                for domain_index in domain_index_list:
                    domain_number = phmmer_stdout_readlines[domain_index].split()[2]
                    
                    alignment_len = len(phmmer_stdout_readlines[domain_index + 1].split()[-2])
                    
                    seq_identity = len([char for char in phmmer_stdout_readlines[domain_index + 2].strip() if char not in ('+', ' ')]) / alignment_len # + and space are the characters used by HMMER to signify lack of alignment
                    
                    homologous_id, begin_id, homologous_seq, end_id = phmmer_stdout_readlines[domain_index + 3].split()
                    homologous_seq = homologous_seq.replace('-', '').upper() # Remove indels in the target sequence to obtain original domain sequence

                    # Check validity of the domain
                    if float(domain_e_val[domain_number]) <= E_val and len(homologous_seq) >= min_len and seq_identity >= min_seq_identity:

                        header = f'> {homologous_id}_{begin_id}_{end_id} | E_val: {domain_e_val[domain_number]}\n'
                        # Transform sequence to FASTA format
                        homologous_seq = ''.join([homologous_seq[i : i + 60] + '\n' for i in range(0, len(homologous_seq), 60)])

                        file.write(header)
                        file.write(homologous_seq)

    return

def run_PHMMER(query_fasta_file_path, database_fasta_file_path,
               min_seq_identity = 0.35, min_seq_len = 0.60, E_val = 0.0001,
               N_CPU = 2):
    """
    Run phmmer with the query against the given database (in FASTA format both) with the given parameter settings.
    Returns two files:
     - Output file generated by phmmer (phmmer_stdout.txt)
     - A FASTA file with all the homologous sequences that fulfill the required parameters (phmmer.fasta).
    """
    # Check input FASTA file and extract sequence
    query_seq = str()
    with open(query_fasta_file_path, 'rt') as f:
        for index, line in enumerate(f):
            
            if index == 0: 
                assert line.startswith('>'), "Incorrect header format, first line must begin with '>'"
            else: 
                assert all(char in Alphabet for char in line.strip()), "Unknown character found in the query's FASTA sequence"
                query_seq += line.strip()

    # Run. Usage: phmmer [-options] <seqfile> <seqdb>
    output_file = 'phmmer_stdout.txt'

    cmd = [soft_exec_paths['phmmer'], '-o', str(output_directory_path / output_file), '--notextw', '--cpu', str(N_CPU), '--domE', str(E_val), 
            '--incE', str(E_val), '-E', str(E_val), str(query_fasta_file_path), str(database_fasta_file_path)]
    subprocess.run(cmd, check=True)

    # Extract homologues at given E-value, minimal sequence identity and minimum sequence length
    PHMMER_stdout_parser(phmmer_stdout_file = output_file,
                         E_val = E_val,
                         query_seq = query_seq,
                         min_seq_identity = min_seq_identity, min_seq_len = min_seq_len)

    # Make sure file is not empty (at least 1 homologous sequence)
    if os.stat(output_directory_path / output_file).st_size == 0:
        raise ValueError("No homologous sequences where found in the database with the given parameters.")

    return

def run_CDHIT(query_fasta_file_path, clustering_seq_identity = 0.95, max_homol = 150):
    """
    """
    input_file = 'phmmer.fasta'
    output_file = 'selected_homologues.fasta'
    
    # Run. Usage: cd-hit [Options]
    cmd = [soft_exec_paths['cd-hit'], '-i', str(output_directory_path / input_file), '-o', str(output_directory_path / output_file),
           '-c', str(clustering_seq_identity)]
    with open(str(output_directory_path / 'cdhit.log'), 'wt') as stdout_file_path:
        subprocess.run(cmd, check=True, stdout=stdout_file_path)

    # Check the number of homologoues sequences left after clustering by reading the cd-hit log file:
    # - If there are less than 5 -> throw error
    # - If there are between 5-10 -> print warning
    with open(output_directory_path / 'cdhit.log', 'rt') as file:
        readlines = file.readlines()

    # The line with the n° of clusters is written in the format 'n° ... finished ... n° ... clusters'
    for index, line in enumerate(readlines):
        line = line.split()
        if all(keyword in line for keyword in ('finished', 'clusters')):
            n_clusters = int(line[-2])
            break
    else:
        raise ValueError("Could not determine the n° of clusters generated by cd-hit")

    if n_clusters < 5:
        raise ValueError(f'Only {n_clusters} homologues were found which is insufficient (need at least 5)')

    if 5 <= n_clusters <= 10: 
        print('Only {n_clusters} homologues were found, calculations continue nonetheless')
       
    # If there are more than 150 homologous sequences after clustering:
    # - Sort homologues by increasing E-value
    # - Sample 150 homologoues sequences
    # - Overwrite output file with the selected sequences
    if n_clusters > max_homol:
        
        with open(output_directory_path / output_file, 'rt') as file:
            readlines = file.readlines()

        indexes = [index for index, line in enumerate(readlines) if line.startswith('>')]
        
        # Sort. Simpler and faster to use reversed index
        sorted_homologues = []
        for index in reversed(indexes):
            
            E_val = float(readlines[index].split()[-1])
            homologue_data = ''.join(readlines[index:]) # Includes header and sequence
            
            sorted_homologues.append((E_val, homologue_data))
            
            del readlines[index:]
        
        sorted_homologues.sort() # Sorts by increasing value by default

        # Sample
        selected_indexes = [int(i * (n_clusters/max_homol)) for i in range(max_homol)] # Equivalent of (int(i) for i in range(0, n_clusters, n_clusters/150)), but range doesn't allow float steps so use this trick instead. int(i) always rounds down (ex: int(1.9) = 1)
        selected_homologues = [sorted_homologues[index][1] for index in selected_indexes]

        # Overwrite
        with open(output_directory_path / output_file, 'wt') as file:
            for homologue in selected_homologues:
                file.write(homologue) # Header and sequence
                
    
    # Rewrite the output file so that the query is the first sequence in the file (MAFFT considers the first sequence as the reference to build the alignment)
    with open(query_fasta_file_path, 'rt') as file: 
        readlines = file.readlines()
    with open(output_directory_path / output_file, 'rt') as file:
        readlines += file.readlines()
        
    with open(output_directory_path / output_file, 'wt') as file:
        for line in readlines: 
            file.write(line)
    
    # Remove .clstr file generated by cd-hit
    os.remove(output_directory_path / (output_file + '.clstr'))

    return

def run_MAFFT(N_CPU = 2):
    """
    """
    input_file = 'selected_homologues.fasta'
    output_file = 'query_msa.aln'
    
    # Run. Usage: mafft [Options] in > out
    cmd = [soft_exec_paths['mafft'], '--thread', str(N_CPU), '--quiet', str(output_directory_path / input_file)]
    with open(str(output_directory_path / output_file), 'wt') as stdout_file_path:
        subprocess.run(cmd, check=True, stdout=stdout_file_path)

    return

def run_evolutionary_model(model, input_file):
    """
    """
    input_file_path             = str(output_directory_path / input_file)
    normalized_grades_file_path = str(output_directory_path / f'{model}_normalized_grades.txt')
    original_grades_file_path   = str(output_directory_path / f'{model}_original_grades.txt')
    tree_file_path              = str(output_directory_path / f'{model}_tree.txt')
    
    # Run. Usage: rate4site [-options]
    cmd = [soft_exec_paths['rate4site'], f'-{model}', '-s', input_file_path, '-o', normalized_grades_file_path, '-y', original_grades_file_path,
           '-x', tree_file_path]
    with open(os.devnull, 'wt') as stdout_file_path:
        subprocess.run(cmd, check=True, stdout=stdout_file_path) # Still does not remove some of rate4site's print messages to console

    # Get the negative log likelihood (-lnl) of the evolutionary model
    line = getline(filename = normalized_grades_file_path, lineno = 13)
    if not line.startswith('#LL='): # Just to make sure
        raise ValueError(f'Could not read the rate4site output file correctly, check {normalized_grades_file_path}') 
    
    likelihood = float(line.split('=')[-1])

    return likelihood

def run_rate4site(evolutionary_model = 'best', N_CPU = 2):
    """
    """
    input_file = 'query_msa.aln'
    
    evolutionary_models_dico = {'best':['Mj', 'Md', 'Mw', 'Ml'],
                                'JTT':['Mj'], 'Dayhoff':['Md'], 'WAG':['Mw'], 'LG':['Ml'],
                                'mitochondrial':['mtREV'],
                                'chloroplast':['cpREV']}
    
    selected_models = evolutionary_models_dico[evolutionary_model]
    
    # Run the models and select the one with maximum likelihood if there are multiple models to evaluate
    models_likelihood = {}
    with ProcessPoolExecutor(max_workers = N_CPU) as executor:
        for model, likelihood in zip(selected_models, executor.map(run_evolutionary_model, selected_models, [input_file]*len(selected_models))):
            models_likelihood[model] = likelihood

    best_model = max(models_likelihood, key=models_likelihood.get)
    
    evolutionary_models_translation = {'Mj':'JTT', 'Md':'Dayhoff', 'Mw':'WAG', 'Ml':'LG',
                                       'mtREV':'mitochondrial(mtREV)', 'cpREV':'chloroplast(cpREV)'}
    
    print(f'Selected model: {evolutionary_models_translation[best_model]}')

    
    # Remove files of unwanted models and rename the files of the selected model (remove prefix)
    for model in models_likelihood:
        if model != best_model:
            os.remove(output_directory_path / f'{model}_original_grades.txt')
            os.remove(output_directory_path / f'{model}_normalized_grades.txt')
            os.remove(output_directory_path / f'{model}_tree.txt')
    
    os.rename(src = output_directory_path / f'{best_model}_original_grades.txt',   dst = output_directory_path / f'original_grades.txt')
    os.rename(src = output_directory_path / f'{best_model}_normalized_grades.txt', dst = output_directory_path / f'normalized_grades.txt')
    os.rename(src = output_directory_path / f'{best_model}_tree.txt',              dst = output_directory_path / f'tree.txt')

    # A 'r4s.res' empty file is generated by rate4site (don't know why)
    if os.path.isfile(working_dir / 'r4s.res') and os.stat(working_dir / 'r4s.res').st_size == 0:
        os.remove(working_dir / 'r4s.res')
    
    return


def check_arguments(args):
    """
    Check validity of parameters/arguments given in the command-line. 
    """
    if not os.path.isfile(args.input_file):
        raise argparse.ArgumentTypeError(f"The input file path is invalid ({args.input_file})")

    if not os.path.isfile(args.database):
        raise argparse.ArgumentTypeError(f"The database file path is invalid ({args.database})")

    if not args.e_val > 0:
        raise argparse.ArgumentTypeError("E val must be > 0")

    if not 0 <= args.min_seq_identity <= 1:
        raise argparse.ArgumentTypeError("Minimal global sequence identity must be between 0-1")
        
    if not 0 <= args.min_seq_len <= 1:
        raise argparse.ArgumentTypeError("Minimal relative length must be between 0-1")
        
    if not 1 <= args.max_homol <= 199:
        raise argparse.ArgumentTypeError("The number of homologues must be between 1-199")
        
    return

def argument_parser():
    """
    """
    parser = argparse.ArgumentParser(description = "Description:\n  Python implementation of the ConSurf pipeline for the estimation of residue's evolutionary conservation.",
                                     usage = 'consurf <input_file> [parameters]',
                                     epilog = "Results are saved in a folder in the current directory\n ",
                                     formatter_class = argparse.RawTextHelpFormatter,
                                     add_help = False)

    parser._positionals.title = 'Input Data'
    parser._optionals.title = 'Parameters'

    # Input data
    parser.add_argument('input_file', type = str,
                        help = 'Full path of a FASTA file containing the query sequence (1 chain only)')

    parser.add_argument('-database', type = str, metavar = '',
                        help = "Full path of a FASTA formated database that is used to search for homologous sequences\n ")
    # Parameters
    parser.add_argument('-model', type = str, choices = ['best', 'JTT', 'Dayhoff', 'WAG', 'LG', 'mitochondrial', 'chloroplast'], default = 'best', metavar = '',
                        help = "Options: best(default), JTT, Dayhoff, WAG, LG, mitochondrial, chloroplast\nEvolutionary substitution model to be used by the rate4site program. By default it will\ndetermine which of JTT/Dayhoff/WAG/LG fits best with the analyzed sequences.\n ")

    parser.add_argument('-e_val', type = float, default = 0.0001, metavar = '',
                        help = "Default: 0.0001\nMinimal E-value of homologous sequences selected by phmmer\n ")

    parser.add_argument('-min_seq_identity', type = float, default = 0.35, metavar = '',
                        help = 'Default: 0.35\nMinimal global sequence identity of homologous sequences w.r.t the input query sequence\n ')

    parser.add_argument('-min_seq_len', type = float, default = 0.60, metavar = '',
                        help = 'Default: 0.60\nMinimal relative length of homologous sequence w.r.t the input query sequence\n ')

    parser.add_argument('-clust_seq_identity', type = float, default = 0.95, metavar = '',
                        help = 'Default: 0.95\nSequence identity cutoff used by cd-hit to cluster similar homologues and select a\nrepresentative sequence (avoids having redundant homologues)\n ')
    
    parser.add_argument('-max_homol', type = int, default = 150, metavar = '',
                        help = 'Default: 150\nMaximum number of homologues to use. If their number exceeds the limit then they are sorted\nby their E value and selected at a regular interval (ex: if there are 600 homologues in\ntotal, then homolgues 1,5,9,... will be selected)\n ')

    parser.add_argument('-n_cpu', type = int, choices = [1, 2], default = 1, metavar = '', ### Adapt choices for github release
                        help = "Options: 1(default) or 2\nNumber of cores\n ")

    parser.add_argument('-h', '--help', action = 'help', default= argparse.SUPPRESS,
                        help = "Show this help message and exit\n ")
    
    return parser

def main(query_fasta_file_path, database_fasta_file_path,
         E_val = 0.0001, min_seq_identity = 0.35, min_seq_len = 0.60,
         clustering_seq_identity = 0.95, max_homol = 150,
         evolutionary_model = 'best',
         N_CPU = 2):
    """
    """
    # Create output directory where all generated files will be saved. If a directory with the same name already exists
    # then create one with a '_n°' extension
    query_id = os.path.split(query_fasta_file_path)[-1] # Tail of path
    query_id = query_id.rsplit('.', maxsplit = 1)[0] if '.' in query_id else query_id
    
    database_id = os.path.split(database_fasta_file_path)[-1]
    database_id = database_id.rsplit('.', maxsplit = 1)[0] if '.' in database_id else database_id

    global output_directory_path # Avoids having to add the directory as as parameter in each function
    output_directory_path = working_dir / f'{query_id}_{database_id}'

    counter = 1
    while os.path.isdir(output_directory_path):
        output_directory_path = str(output_directory_path).replace(f'_{str(counter - 1)}', '') # Removes extension if there is one
        output_directory_path = Path(f'{output_directory_path}_{str(counter)}')

        counter += 1

    os.mkdir(output_directory_path)
    
    # Run pipeline
    print('Finding homologous sequences ...')
    run_PHMMER(query_fasta_file_path = query_fasta_file_path, database_fasta_file_path = database_fasta_file_path,
               E_val = E_val, min_seq_identity = min_seq_identity, min_seq_len = min_seq_len, N_CPU = N_CPU)
    
    print('Clustering sequences to remove redundancy ...')
    run_CDHIT(query_fasta_file_path = query_fasta_file_path, clustering_seq_identity = clustering_seq_identity,
              max_homol = max_homol)
    
    print('Building multiple sequence alignement ...')
    run_MAFFT(N_CPU = N_CPU)
    
    print('Calculating conservation scores ...')
    run_rate4site(evolutionary_model = evolutionary_model, N_CPU = N_CPU)

    print(f'Consurf run was successful. Results are in {output_directory_path}.')
    
    return

if __name__ == '__main__':
    # Parse and validate command-line arguments
    parser = argument_parser()
    args = parser.parse_args()

    check_arguments(args)

    # Run
    main(query_fasta_file_path = args.input_file, database_fasta_file_path = args.database,
         E_val = args.e_val, min_seq_identity = args.min_seq_identity, min_seq_len = args.min_seq_len,
         clustering_seq_identity = args.clust_seq_identity, max_homol = args.max_homol,
         evolutionary_model = args.model,
         N_CPU = args.n_cpu)

# Allow for pre-built MSA input ?
# Write log file in case of errors ?
# Perform mapping to PDB ATOM records ?