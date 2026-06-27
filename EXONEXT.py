from Bio import SearchIO
import pandas as pd
import os, shutil
import time
import re
from collections import defaultdict
    
print("\n", "Which directory has the out files to use? ")
direct = "Exonerate_output_all" #input()

print("\n", "Which method should be used? ")
print("  Method 2 is much faster, but not all cigar entries will be registered in the cigar_lines.xlsx output")
print("   (1) using raw out files")
print("   (2) using hit files (out files filtered by raw score)")
print("   (3) using hit files only for the large out files")

method = 2 #int(input())

print("\n", "Do you want to register the location and size of introns? (Y/N)")

intron_input = "Y" #input()
t0 = time.time()

#######################################
#%% Filtering out files (hit script)
#######################################
def extract_hits(file_path, threshold):
     """
     Retrieves all hits of an output file that have raw scores above a threshold
 
     Parameters
     ----------
     file_path : out file path
     value : raw score threshold
 
     Returns
     -------
     hits : List of the hits with raw scores above the defined threshold
 
     """
     hits = []
     hit_lines = []
     with open(file_path, 'r') as file:
         for line in file:
             if 'Raw score:' in line:
                 raw_score = float(line.split(':')[1].strip())
                 
                 if raw_score >= threshold:
                     hit_lines.append(line)
                     
             elif '# --- START OF GFF DUMP ---' in line:
                 if hit_lines:
                     hits.append(hit_lines)
                     hit_lines = []
             elif hit_lines:
                 hit_lines.append(line)
     
     return hits

# Get initial directory
cwd = os.getcwd()
out_d = os.path.join(cwd, direct)
    
if method not in [1, 2, 3]:
    raise TypeError("The selected method must be 1 or 2. %s was given" %method)
        
elif method == 2 or method == 3:
    # Define raw score threshold
    print("\n", "Define the raw score threshold: ")
    threshold = 200 #int(input())
    print("\n", "Selecting the hits with raw scores above %i" %threshold)
         
    # Change to the out directory
    os.chdir(out_d)
              
    # List the out files
    outfiles = [f for f in os.listdir() if f.endswith('.out')]
        
    # Process each output file individually
    for out_file in outfiles:            
        if method == 2:
            # Retrieve the hits corresponding to the desired raw score values  
            hits = extract_hits(out_file, threshold)
            
            # Write the hits to the output file
            #======================= modified at 12/11/2025 for the GRC gene extraction
            hit_file = f'{out_file.split("_")[0]}_{out_file.split("_")[-2]}_{out_file.split("_")[-1]}{str(threshold)}.hit'
            
            if "GRC" in out_file:
                hit_file = f'{out_file.split("_v3")[0]}.out{str(threshold)}.hit'
            #======================================================================================
                    
            with open(hit_file, 'w') as outfile:
                for hit in hits:
                    outfile.writelines(hit)
                    outfile.write('\n')
            
            print("%i hits extracted from %s" % (len(hits), outfile.name), "\n")
    
    os.chdir(cwd)

###################################################################
#%% List the cigar lines
###################################################################

results1 = pd.DataFrame()
hits = []
cwd = os.getcwd()
start_time = time.time()
total_time = 0
sum_time = 0
sum_files = 0
sum_size = 0
i = 0
used_entries = set()
    
def get_directory_size(directory):
    total_size = 0
    
    for dirpath, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            filepath = os.path.join(dirpath, filename)
            total_size += os.stat(filepath).st_size
    return total_size

# List cigar lines
dst = os.path.join(cwd, direct)
os.chdir(dst)

if method == 1:
    flist = [f for f in os.listdir() if f.endswith('.out')]
    
elif method == 2:
    flist = [f for f in os.listdir() if f.endswith(str(threshold) + '.hit')]
    
elif method == 3:
    flist = [f for f in os.listdir() if f.endswith('.out')] + [f for f in os.listdir() if f.endswith(str(threshold) +'.hit')]
       
# Create a dataframe with all the cigar lines in the out files present in the selected directory
print("Selecting the best exonerate hits:")

for i, fname in enumerate(flist):
    start_gene = time.time() # start the clock

    print(fname, ":", round(i/len(flist) * 100, 2), "%")
        
    # Open the file
    with open(fname, 'r') as f:
        for j, line in enumerate(f):
                       
            # Register the information in the cigar lines
            if line.startswith('cigar:'): 
                splits = line.split()
                    
                new_row = {
                    "taxon": direct,
                    "file": fname,
                    "cigar": line,
                    "species": fname.split("_")[1].split(".")[0],
                    "gene": splits[1].split("_")[-1],
                    "paralog": "GRC" if "GRC" in fname else "ACHR",
                    "query_start": splits[2],
                    "query_stop": splits[3],
                    "scaffold": splits[5],
                    "start": splits[6],
                    "stop": splits[7],
                    "frame": splits[8],
                    "score": int(splits[9])
                    }

                results1 = pd.concat(
                    [results1, pd.DataFrame([new_row])],
                    ignore_index=True
                )
        
# Final processing
results1.sort_values(["gene", "species", "score"],
                         ascending=[True, True, False], inplace=True)
         
#Write the DataFrame to an Excel file
output_path = os.path.join(cwd, f"{direct}_cigar_lines.xlsx")
results1.to_excel(output_path, index=False)
 
os.chdir(cwd) 
        
print("\n", direct, " cigar output created:")
print(results1.head())
print("\n")

###################################################################
#%% NEW: READING ANNOTATION FILES TO SELECT CIGAR LINES
def read_annotation(gtf_file):
    """
    Returns:
        gene_map = {
            gene_name: [
                (gene_id, full_line),
                (gene_id, full_line),
                ...
            ]
        }
    """
    gene_map = defaultdict(list)

    # opener for gz or plain text
    opener = open
    if gtf_file.endswith(".gz"):
        import gzip
        opener = gzip.open

    with opener(gtf_file, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue

            feature = cols[2].lower()
            
            if feature not in ["transcript"]:  
                continue
            attributes = cols[8]

            # Extract gene_id and gene_name from GTF attributes
            gene_id_match = re.search(r'gene_id "([^"]+)"', attributes)
            gene_name_match = re.search(r'gene_name "([^"]+)"', attributes)

            if not gene_id_match or not gene_name_match:
                continue

            gene_id = gene_id_match.group(1)
            gene_name = gene_name_match.group(1)

            # Store tuple: (gene_id, full GTF line)
            gene_map[gene_name].append((gene_id, line.rstrip("\n")))

    return gene_map

# Create a dictionary for the four GRC annotations
annotation_files = [f for f in os.listdir("../0-assemblies_annotations") if f.endswith('.gtf')]

all_annotations = {}

for gtf in annotation_files:
    gene_dict = read_annotation(f"../0-assemblies_annotations/{gtf}")
    
    # Use filename (without path) as the dictionary key
    key = gtf.split(".")[0]
    if "GRC" in gtf:
        key = gtf.split("cds_to")[-1].split("_")[0]
        
    all_annotations[key] = gene_dict
 
#####################################################################
#%% Sequence selection (select hits to extract and create cut files) MAX SCORE FOR ACHRs AND ANNOTATED SEQUENCE FOR GRCs
#####################################################################
def extract_alignment_info(filename, cigar_line):
    """
    Extract only the lines between the nearest preceding 'Query range:' 
    and the specified CIGAR line, inclusive.
    """
    cigar_line = cigar_line.strip()
    alignment_info = []
    query_index = None
    cigar_index = None

    # Read the file
    with open(filename, 'r') as f:
        lines = [line.rstrip() for line in f]

    # Find the CIGAR line
    for i, line in enumerate(lines):
        if line.strip() == cigar_line:
            cigar_index = i
            break

    if cigar_index is None:
        raise ValueError(f"CIGAR line not found in {filename}: {cigar_line}")

    # Find the nearest preceding "Query range:" before the CIGAR line
    for i in range(cigar_index, -1, -1):
        if "Query range:" in lines[i]:
            query_index = i
            break

    if query_index is None:
        raise ValueError(f"No preceding 'Query range:' found before CIGAR line in {filename}")

    # Extract only the block between them
    alignment_info = lines[query_index:cigar_index + 1]

    return alignment_info
  
#%% ACHRs  
# creating the cut files
hits = pd.read_excel(direct + '_cigar_lines.xlsx')
hits_ACHR = hits[hits["paralog"] == "ACHR"]
os.chdir(os.path.join(cwd, direct))
results = pd.DataFrame()

# Find the rows with the highest score for each gene/species combination
max_scores = hits_ACHR.groupby(['gene', 'species', "paralog"])['score'].idxmax()
df_max_scores_ACHR = hits_ACHR.loc[max_scores]
df_max_scores_ACHR["query"] = df_max_scores_ACHR["cigar"].str.split().str[1]
query_split = df_max_scores_ACHR["query"].str.split("_", n=2, expand=True) # Split the query
df_max_scores_ACHR["transcript_id"] = query_split[0] + "_" + query_split[1]

# Iterate over the rows with the highest score
for index, row in df_max_scores_ACHR.iterrows():
    # Access the required columns as needed
    file = row["file"]
    cigar = row["cigar"]
    gene = row['gene']
    specie = row['species']
    score = row['score']
    paralog = row['paralog']
    transcript_id = row['transcript_id']
    
    output = f"{gene}_{specie}_{paralog}_{transcript_id}_cut.txt"
           
    alignment_info = extract_alignment_info(file, cigar) # run the previous function with that info
    
    # Create cuts folder if it doesn't exist
    cuts_folder = os.path.join(cwd, "cuts")
    os.makedirs(cuts_folder, exist_ok=True)
    output_path = os.path.join(cuts_folder, output)

    # Save extracted region
    with open(output_path, 'w') as output_file:
        for line in alignment_info:
            output_file.write(line + '\n')

    print(f"Extracted information saved to {output_path}")
    print(cigar)
              
os.chdir(cwd)
    
#%% GRCs
hits_GRC = hits[hits["paralog"] == "GRC"]
df_max_scores_GRC = pd.DataFrame() 
genomes = hits["species"].unique() 

for genome in genomes: 
    df_species = hits_GRC[hits_GRC["species"] == genome] 
    genes = df_species["gene"].unique()    
    
    for gene in genes:
        df_gene = df_species[df_species["gene"] == gene]
        transcript_entries = all_annotations[genome][gene]  # this is a list of tuples
        df_gene["copy_number"] = len(transcript_entries)
        
        # add a query column (to help select the correct entry in the annotation)
        df_gene["query"] = df_gene.cigar.str.split().str[1]
        
        for transcript_id, full_line in transcript_entries:
            cols = full_line.strip().split("\t")
            scaffold = cols[0]
            frame = cols[6]
            
            if frame == "+":
                start = int(cols[3])
                end = int(cols[4])
                
            elif frame == "-":
                start = int(cols[4])
                end = int(cols[3])
                
            score = int(cols[5])
            ids = cols[8].strip().split('"')
            gene_id = ids[1]
            transcript_id = ids[3]
            gene_name = ids[5]
            
            #Filter the cigar line dataframe for entries in annotation
            query_id = f"{transcript_id.split('_')[1]}_{transcript_id.split('_')[2]}_{gene}"
            df_annotated = df_gene[df_gene["query"] == f"{query_id}"] #same query
            df_annotated0 = df_annotated[df_annotated["scaffold"] == scaffold.upper()] #same scaffold
            df_annotated1 = df_annotated0[df_annotated0["frame"] == frame] #same reading frame
            
            # filter the transcript line that has virtually the same starting position and score
            df_annotated2 = df_annotated1[ (df_annotated1["score"] - score).abs() <= 10 ] # similar scores
            df_annotated3 = df_annotated2[ (df_annotated2["start"] - start).abs() <= 10 ].iloc[:1] #very similar start
            df_annotated3["transcript_id"] = transcript_id
            
            # save entries that match the annotation
            df_max_scores_GRC = pd.concat(
                [df_max_scores_GRC, df_annotated3],
                ignore_index=True
            )

# Iterate over the rows that match the annotation
for index, row in df_max_scores_GRC.iterrows():
    file = row["file"]
    cigar = row["cigar"]
    gene = row['gene']
    specie = row['species']
    score = row['score']
    paralog = row['paralog']
    transcript_id = row['transcript_id']
    
    output = f"{gene}_{specie}_{paralog}_{transcript_id}_cut.txt"
           
    alignment_info = extract_alignment_info(f"{direct}/{file}", cigar) 
    
    # Create cuts folder if it doesn't exist
    cuts_folder = os.path.join(cwd, "cuts")
    os.makedirs(cuts_folder, exist_ok=True)
    output_path = os.path.join(cuts_folder, output)

    # Save extracted region
    with open(output_path, 'w') as output_file:
        for line in alignment_info:
            output_file.write(line + '\n')

    print(f"Extracted information saved to {output_path}")
    
    print(cigar)
              
os.chdir(cwd)

df_max_scores = pd.concat([df_max_scores_ACHR, df_max_scores_GRC])
df_max_scores.to_excel(direct + "_selected_sequences.xlsx")
 
######################################################
#%% Sequence extraction
######################################################

cuts_folder = os.path.join(cwd, "cuts")
os.chdir(cuts_folder)

for filename in os.listdir(cuts_folder): 
    if filename.endswith("cut.txt"):
        nome = filename
        f = open(nome, 'r', newline = '\n')
        i = 0
        j = 4
        fin = []
        
        lines = f.readlines()
    
        for line in lines[3:]:
            i += 1
            if i == j :
                fin.append(line)
                j += 5
        f.close()
        upper=[]
        for el in fin:
            for i in el:
                if i.isupper():
                    upper.append(i)
        
        result = ''.join(upper)
        
        num = nome.count('_')
        nome_final= nome.split('_')
        header =''
        title = ''
        for i in range(0, num):
            title += nome_final[i]
            if i < num-1:
                title+='_'
        
        fasta_folder = os.path.join(cwd, "fasta_sequences")
        
        if not os.path.exists(fasta_folder):
              os.mkdir(fasta_folder) 
              
        os.chdir(fasta_folder)
        title1 = title + '.fas'
        j = open(title1, 'w')
        j.write('>' + title + '\n')
        j.write(result)
        j.close()
        
        os.chdir(cuts_folder)
        
os.chdir(cwd) 

print("\n", "Fasta sequences created!", "\n")
        
##############################################################################
#%% See if the sequences are complete
##############################################################################
"""
    Checks if the fasta sequences are not partial
    Creates a list of incomplete sequences
    Returns the percentage of completeness
"""

from Bio import SeqIO

cwd = os.getcwd()
dst = os.path.join(cwd, direct)
fasta_folder = os.path.join(cwd, "fasta_sequences")
file_path = direct + '_selected_sequences.xlsx'
hits = pd.read_excel(file_path)
results = pd.DataFrame()
    
species_list = hits["species"].unique()
    
for species in species_list:
    species_df = hits[hits["species"] == species]
    genes = species_df["gene"].unique()
        
    for gene in genes:        
        gene_df = species_df[species_df["gene"] == gene].copy().reset_index()  # create a copy to avoid warning
        paralog = gene_df["paralog"].iloc[0]
        
        # Split transcript_id into parts
        split_transcripts = gene_df["transcript_id"].str.split("_", expand=True)
        
        # Combine parts 1 and 2 (second and third elements, 0-based)
        gene_df["transcript"] = gene_df["transcript_id"]
        
        if paralog == "GRC":
            gene_df["transcript"] = split_transcripts[1] + "_" + split_transcripts[2]
            
        transcripts = gene_df["transcript_id"].unique()
        
        for transcript_id in transcripts:
            df_transcript = gene_df[gene_df["transcript_id"] == transcript_id].copy()
            transcript = df_transcript["transcript"].iloc[0]
            query_id = df_transcript["query"].iloc[0]      
            records = SeqIO.to_dict(SeqIO.parse(f"Query_all/{gene}.fas", "fasta"))

            query_seq = str(records[query_id].seq)
            query_size = len(query_seq)
    
            info = f"{gene}_{species}_GRC_{transcript_id}.fas"
             
            df_transcript["fasta_file"] = info
            df_transcript["query_size"] = query_size
            df_transcript["sequence_size"] = df_transcript["query_stop"] - df_transcript["query_start"]
            df_transcript["completeness"] = df_transcript["sequence_size"]/query_size * 100
    
            results = pd.concat([results, df_transcript], ignore_index=True)  # add new row to the output
    
# write results to Excel after processing all species and genes
output_file_path = direct + "_selected_sequences_completeness.xlsx"
results.to_excel(output_file_path, index=False)
    
print("Info on completeness added to," , output_file_path)
print("\n")  

###################################
#%% Find pseudogenes (stop codons)
###################################

def detect_sequence_stops(fasta_file):
    """
    Detect potential pseudogenes in the fasta sequences within a specified region while maintaining reading frame.
    Creates a list of potential pseudogenes.
    
    """
    potential_stops = ""

    stop_codons = {"TAA": 0, "TAG": 0, "TGA": 0}

    with open(fasta_file, 'r') as file:
            sequences = file.read().split('>')[1:]  # Skip the empty first element

    for sequence in sequences:
            header, seq = sequence.split('\n', 1)
            seq = seq.replace('\n', '')
            pseudo_location = []

            for i in range(0, len(seq) - 3, 3):
                if seq[i] == "T" and seq[i + 1] == "A" and seq[i + 2] == "A":
                    stop_codons["TAA"] += 1
                    pseudo_location.append(i)
                    
                elif seq[i] == "T" and seq[i + 1] == "A" and seq[i + 2] == "G":
                    stop_codons["TAG"] += 1
                    pseudo_location.append(i)
                
                elif seq[i] == "T" and seq[i + 1] == "G" and seq[i + 2] == "A":
                    stop_codons["TGA"] += 1
                    pseudo_location.append(i)

            if (stop_codons["TAA"] + stop_codons["TAG"] + stop_codons["TGA"]) > 0:
                #stop_info = fasta_file + ", stop at " + str(pseudo_location) + " bps"
                potential_pseudogenes = pseudo_location
                potential_stops = [int(pos / 3) for pos in potential_pseudogenes]
                
                print("%s is pseudogenized (TAA: %i; TAG: %i; TGA: %i)" % (fasta_file,
                                                        stop_codons["TAA"],
                                                        stop_codons["TAG"],
                                                        stop_codons["TGA"]))
                
    return potential_stops  

def get_mutation_positions(cut_filepath):
    # Step 1: Read and extract specific lines
    f = open(f"{cut_filepath}", 'r', newline='\n')  # abrir o ficheiro

    lines = f.readlines()
    fin = lines[5:-1:5]  # start at 4th line, then every 5th line

    # Step 2: Clean lines
    new_fin = []

    for line in fin:
        line = re.sub(r"\s+", "", line)  # Remove all whitespace
        line = line.replace("<->", "ARROW")  # Temporarily protect "<->"
        line = re.sub(r"([A-Z])-([a-z])", r"\1\2", line)  # Remove "-" between uppercase and lowercase letters
        matches = re.findall(r"[A-Za-z\*\-#]|ARROW", line)  # Keep letters, *, -, # and special token
        cleaned_line = "".join(matches)
        cleaned_line = cleaned_line.replace("ARROW", "<->")  # Restore "<->"
        new_fin.append(cleaned_line)

    fin = new_fin

    # Step 3: Concatenate lines
    upper = "".join(fin).replace(" ", "")

    # Step 4: Parse codons / amino acids
    amino_acids = {
        "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu",
        "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
        "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"
    }

    sequence = []
    i = 0
    n = len(upper)

    while i < n:
        if upper[i] == "#":
            sequence.append(upper[i])
            
        if upper[i:i+3] in ["***", "---", "<->"]:
            sequence.append(upper[i:i+3])
            i += 3
            continue

        codon = upper[i:i+3]

        if re.search(r"[A-Za-z]", codon):
            letters = []
            j = i
            while len(letters) < 3 and j < n:
                if upper[j].isalpha():
                    letters.append(upper[j])
                j += 1
            candidate = "".join(letters)
            if candidate in amino_acids:
                sequence.append(candidate)
            i = j
        else:
            sequence.append(codon)
            i += 3

    # Step 5: Track stop codons
    stop_positions = []
    frameshift_positions = []
    deletions = []
    current_index = 0

    for i in sequence:
        current_index += 1
        if i == "***":
            stop_positions.append(current_index)
        
        if "#" in i:
            frameshift_positions.append(current_index)   
        
        if "---" in i:
            deletions.append(current_index)   

    return stop_positions, frameshift_positions, deletions

def get_insertion_positions(cut_filepath):
    # Step 1: Read and extract specific lines
    f = open(f"{cut_filepath}", 'r', newline='\n')  # abrir o ficheiro

    lines = f.readlines()
    fin = lines[3:-1:5]  # start at 4th line, then every 5th line

    # Step 2: Clean lines
    new_fin = []

    for line in fin:
        line = re.sub(r"\s+", "", line)  # Remove all whitespace
        line = line.replace("<->", "ARROW")  # Temporarily protect "<->"
        line = re.sub(r"([A-Z])-([a-z])", r"\1\2", line)  # Remove "-" between uppercase and lowercase letters
        matches = re.findall(r"[A-Za-z\*\-#]|ARROW", line)  # Keep letters, *, -, # and special token
        cleaned_line = "".join(matches)
        cleaned_line = cleaned_line.replace("ARROW", "<->")  # Restore "<->"
        new_fin.append(cleaned_line)

    fin = new_fin

    # Step 3: Concatenate lines
    upper = "".join(fin).replace(" ", "")

    # Step 4: Parse codons / amino acids
    amino_acids = {
        "Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu",
        "Gly", "His", "Ile", "Leu", "Lys", "Met", "Phe",
        "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"
    }

    sequence = []
    i = 0
    n = len(upper)

    while i < n:            
        if upper[i:i+3] in ["<->"]:
            sequence.append(upper[i:i+3])
            i += 3
            continue

        codon = upper[i:i+3]

        if re.search(r"[A-Za-z]", codon):
            letters = []
            j = i
            while len(letters) < 3 and j < n:
                if upper[j].isalpha():
                    letters.append(upper[j])
                j += 1
            candidate = "".join(letters)
            if candidate in amino_acids:
                sequence.append(candidate)
            i = j
        else:
            sequence.append(codon)
            i += 3

    # Step 5: Track stop codons
    insertions = []
    current_index = 0

    for i in sequence:
        current_index += 1        
        if "<->" in i:
            insertions.append(current_index) 

    return insertions

##############################################################
#%% Register the stop codons and frameshifts to the dataframe
##############################################################
results = []

file_path = direct + "_selected_sequences_completeness.xlsx"
hits = pd.read_excel(file_path)

for specie in hits["species"].unique():
    species_df = hits[hits["species"] == specie]

    for gene in species_df["gene"].unique():
        gene_df = species_df[species_df["gene"] == gene]
        
        for transcript in gene_df["transcript_id"].unique():
            df = gene_df[gene_df.transcript_id == transcript].copy()
            query = df["query"].iloc[0]
            query_start = df["query_start"].iloc[0]
            paralog = df["paralog"].iloc[0]

            info = f"{gene}_{specie}_{paralog}_{transcript}"
            cut_filepath = f"cuts/{info}_cut.txt"
            fasta_file = f"fasta_sequences/{info}.fas"
            query_file = f"Query_all/{gene}.fas"

            with open(cut_filepath, "r") as cut_file:
                cut = cut_file.read()

            stop_sites, frameshift_sites, deletions = get_mutation_positions(cut_filepath)
            insertions = get_insertion_positions(cut_filepath)
            
            # adjust positions to query start (so that they match the query, not the sequence)
            stop_sites = [x + query_start for x in stop_sites]
            frameshift_sites = [x + query_start for x in frameshift_sites]
            deletions = [x + query_start for x in deletions]
            insertions = [x + query_start for x in insertions]
            
            # register to dataframe
            df.loc[:, "STOP_sites"] = str(stop_sites)
            df.loc[:, "N_STOPs"] = len(stop_sites)
            
            # Also register the position of frameshift mutations
            df.loc[:, "frameshift_sites"] = str(frameshift_sites)
            df.loc[:, "N_frameshifs"] = len(frameshift_sites)
            
            #register deletions and insertions
            df.loc[:, "deletion_sites"] = str(deletions)
            df.loc[:, "N_deletions"] = len(deletions)
            
            df.loc[:, "insertion_sites"] = str(insertions)
            df.loc[:, "N_insertions"] = len(insertions)
            
            #register stop codons in the sequence
            sequence_stops = detect_sequence_stops(fasta_file)
            
            if len(sequence_stops) > 0:
                df.loc[:, "functional?"] = "pseudogene"
                df.loc[:, "pseudogene_type"] = (
                    "frameshift" if "#" in cut else "STOP"
                )
                
                df.loc[:, "sequence_STOP_sites"] = str(sequence_stops)
                df.loc[:, "N_sequence_STOPs"] = len(sequence_stops)

            else:
                df.loc[:, "functional?"] = "functional"

            results.append(df)

# Final output
results_df = pd.concat(results, ignore_index=True)
results_df.to_excel(direct + "_selected_sequences_final.xlsx", index=False)

print("Saved:", direct + "_selected_sequences_final.xlsx")

########################################
# Add a mutation site sheet
########################################
import ast

print("Registering pseudo sites to their own long dataframe ...")

# Prepare a new dataframe with one site per row
mutation_rows = []

for idx, row in results_df.iterrows():
    # Process STOP sites
    stop_sites = row.get("STOP_sites", [])
    
    if pd.isna(stop_sites):
        stop_sites_list = []
        
    else:
        stop_sites_list = ast.literal_eval(stop_sites)

    for site in stop_sites_list:
        mutation_rows.append({
            "species": row["species"],
            "gene": row["gene"],
            "transcript_id": row["transcript_id"],
            "site": site,
            "mutation_type": "STOP"
        })
    
    # Process frameshift codons
    frameshift_sites = row.get("frameshift_sites", [])
    if pd.isna(frameshift_sites):
        frameshift_sites_list = []
        
    else:
        frameshift_sites_list = ast.literal_eval(frameshift_sites)
        
    for site in frameshift_sites_list:
        mutation_rows.append({
            "species": row["species"],
            "gene": row["gene"],
            "transcript_id": row["transcript_id"],
            "site": site,
            "mutation_type": "frameshift"
        })
        
    # Process deletions
    deletions = row.get("deletion_sites", [])
    if pd.isna(deletions):
        deletions_list = []
        
    else:
        deletions_list = ast.literal_eval(deletions)
        
    for site in deletions_list:
        mutation_rows.append({
            "species": row["species"],
            "gene": row["gene"],
            "transcript_id": row["transcript_id"],
            "site": site,
            "mutation_type": "deletion"
        })    

# Create a new dataframe
mutations_df = pd.DataFrame(mutation_rows)

#########################################
#%% Register introns location and size
#########################################
if intron_input == "Y": 
    sel_seqs = pd.read_excel(output_file_path)
    os.chdir(out_d)
            
    os.chdir(cwd)
    results = pd.DataFrame()
    results1 = pd.DataFrame()
    previous_end = []
    previous_start = []
    
    # Changing directory
    cwd = os.getcwd()
    dst = os.path.join(cwd, direct)
    os.chdir(dst)
    
    print("Finding intron information: ", "\n")
    
    # Rebuild list of output files after cleanup
    flist = [f for f in os.listdir() if f.endswith(".out")]

    for fname in flist:
        all_qresult = list(SearchIO.parse(fname, 'exonerate-text'))
        print("\n")
        
        for query_result in all_qresult:
            for hit in query_result:
                for hsp in hit:    
                    previous_end = []
                    previous_start = []
                    hit_starts = []
                    intron_sizes = []
                    n = 0
                    
                    for hsp_fragment in hsp:
                        previous_end += [hsp_fragment.hit_end]
                        previous_start += [hsp_fragment.hit_start]
                        frame = hsp_fragment.hit_frame
                        n += 1
                        
                        if frame > 0:
                            frame_abs = "+"
                            if len(previous_end) > 1:
                                size_pos = hsp_fragment.hit_start - previous_end[-2]
                                
                                if size_pos > 10:
                                    hit_starts += [hsp_fragment.query_start]
                                    intron_sizes += [size_pos]
                                    
                                if n == 1 and len(intron_sizes) > 0:
                                    del intron_sizes[0]
                                    del hit_starts[0]
                                                            
                        elif frame < 0:
                            frame_abs = "-"
                            if len(previous_end) > 1:
                                size_neg = previous_start[-2] - hsp_fragment.hit_end
                                
                                if size_neg > 10:
                                    hit_starts += [hsp_fragment.query_start]
                                    intron_sizes += [size_neg]  
                                
                                if n == 1 and len(intron_sizes) > 0:
                                    del intron_sizes[0]
                                    del hit_starts[0]
                      
                    hit_min = min(previous_start)
                    hit_max = max(previous_end)
                    
                    if frame_abs == "-":
                        hit_max = min(previous_start)
                        hit_min = max(previous_end)
                        
                    for i in range(len(hit_starts)):
                        print("File %s, Gene %s, scaffold %s, Raw score %i, Intron %i: %i %i bp" % (fname, 
                                                                        query_result.id,
                                                                        hit.id,
                                                                        hsp.score,
                                                                        i + 1,
                                                                        hit_starts[i],
                                                                        intron_sizes[i]))
                        # save locations and sizes in long format
                        results1 = pd.concat([
                            results1,
                            pd.DataFrame({
                                "species": fname.split("_")[1],
                                "out_file": fname,
                                "query": query_result.id,
                                "gene": query_result.id.split("_")[2],
                                "score": hsp.score,
                                "scaffold": hit.id,
                                "start": hit_min,
                                "stop": hit_max,
                                "intron": i + 1,
                                "intron_location": hit_starts[i],
                                "intron_size": intron_sizes[i],
                                "N_introns": len(intron_sizes)
                            }, index=[0])
                        ], ignore_index=True)
                        
                    # save locations and sizes as a str(list)   
                    results = pd.concat([
                            results,
                            pd.DataFrame({
                                "species": fname.split("_")[1],
                                "out_file": fname,
                                "query": query_result.id,
                                "gene": query_result.id.split("_")[2],
                                "score": hsp.score,
                                "scaffold": hit.id,
                                "start": hit_min,
                                "stop": hit_max,
                                "intron_locations": str(hit_starts),
                                "intron_sizes": str(intron_sizes),
                                "N_introns": len(intron_sizes)
                            }, index=[0])
                        ], ignore_index=True)
                                           
    os.chdir(cwd)
    results.to_excel(direct + "_intron_information.xlsx")    
    results1.to_excel(direct + "_intron_information1.xlsx")    

    # Creation of the final output
    selected = pd.read_excel(direct + "_selected_sequences_final.xlsx")
    introns = pd.read_excel(direct + "_intron_information.xlsx")
    introns1 = pd.read_excel(direct + "_intron_information1.xlsx")

    final_introns = pd.merge(selected, introns,
                             on = ["species", "gene", "score", "scaffold", "query", "start", "stop"],
                             how = "inner")

    introns_df = pd.merge(selected, introns1,
                             on = ["species", "gene", "score", "scaffold", "query", "start", "stop"],
                             how = "inner")

    final_introns = final_introns.drop(["Unnamed: 0_x", "Unnamed: 0_y", "index"], axis=1).drop_duplicates()
    introns_df = introns_df.drop(["Unnamed: 0_x", "Unnamed: 0_y", "index"], axis=1).drop_duplicates()

    # create the ouput file and save it
    output = f"{direct}_final_table.xlsx"
    with pd.ExcelWriter(os.path.join(cwd, output), engine="xlsxwriter") as writer:
        final_introns.to_excel(writer, sheet_name="final_table")
        introns_df.to_excel(writer, sheet_name="introns")
        mutations_df.to_excel(writer, sheet_name = "pseudogenes")
        
    for f in [
        f"{direct}_cigar_lines.xlsx",
        f"{direct}_selected_sequences_completeness.xlsx",
        f"{direct}_selected_sequences.xlsx",
        f"{direct}_selected_sequences_final.xlsx",
        f"{direct}_intron_information1.xlsx",
        f"{direct}_intron_information.xlsx"
    ]:
        if os.path.exists(f):
            print(f)
            os.remove(f)

    print("\n", "Final output saved to %s" % output)
    print("Done!")
        
print("\n", "Final output saved to %s" % output)        
print("Done!")  
