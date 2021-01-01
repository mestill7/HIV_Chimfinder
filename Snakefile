import os
import sys
from os import listdir
from os.path import isfile, join, isdir
import re
from Bio import SeqIO
import gzip
from collections import Counter
from pathlib import Path

configfile: "config.yaml"
myfastqpath = "fastq/"

## HELPER FUNCTIONS
# Create function for creating rule sets
# Create read-pair inputs for sample processing
def create_endings(x):
    """
    Returns a list of likely fastq file endings
    Input Parameter: 
    x (int): Either 1 or 2, indicating the forward (1) or reverse (2) read.
    Returns: 
    list: A list of strings, representing the file endings a user might 
    use for denoting their fastq files.
    """
    return(["_R" + str(x) + "_001.fastq", "_R" + str(x) + "_001.fq",
            "_R" + str(x) + ".fastq", "_R" + str(x) + ".fq",
            "_" + str(x) + ".fastq", "_" + str(x) + ".fq",
            ".R" + str(x) + "_001.fastq", ".R" + str(x) + "_001.fq",
            ".R" + str(x) + ".fastq", ".R" + str(x) + ".fq",
            "." + str(x) + ".fastq", "." + str(x) + ".fq",
            "_r" + str(x) + "_001.fastq", "_r" + str(x) + "_001.fq",
            "_r" + str(x) + ".fastq", "_r" + str(x) + ".fq",
            ".r" + str(x) + "_001.fastq", ".r" + str(x) + "_001.fq",
            ".r" + str(x) + ".fastq", ".r" + str(x) + ".fq"])

# Function to list the fastq files present in the fastq folder


def getfilelist(myfastqpath):
    """
    Extracts fastq files from the files present in your fastq directory.
    Input Parameter: 
    myfastqpath (string): directory containing your fastq files.
    Returns:
    list: List containing two strings. 
    1st string is all non-metadata files in the fastq directory
    2nd string is all non-metadata files ending in '.gz'
    """
    onlyfiles = [f for f in listdir(myfastqpath) if
                 isfile(join(myfastqpath, f))]
    onlyfiles = [i for i in onlyfiles if
                 i.endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz"))]
    gzfiles = [i for i in onlyfiles if i.endswith((".gz"))]
    return([onlyfiles, gzfiles])


def rename_files(oldname, replacement, myfastqpath):
    [os.rename(join(myfastqpath, i), join(myfastqpath, y))
     for i, y in zip(oldname, replacement)]

# Unify fastq files to single file ending


def fix_input_files(file_suffix, input_fileset, myfastqpath):
    """
    Renames mixed input fastq files to the most common file ending and 
    returns the selected file ending. NOTE: This step permenantly 
    renames your fastq files from their original file ending. 
    Input Parameter: 
    file_suffix (string): ".gz" or ""; Gzipped fastq files are expected 
    to end with the suffix ".gz". If files are NOT gzipped, the input is "".
    input_fileset (list): List of fastq file names to be examined. 
    As written, gzipped files are listed within the variable 'gzfiles' 
    and non-gzipped files are listed within the variable 'onlyfiles'.
    Returns: 
    list: A list containing four strings, the selected Read1 (forward read) 
    file ending and the corresponding Read2 (reverse read) file ending, 
    a list of all fastq-like files, and a list of gzipped fastq-like files.
    """

    # file_suffix, input_fileset = [".gz", gzfiles]
    # Create the series of fastq file endings to search
    base_endings_r1, base_endings_r2 = [create_endings(i) for i in (1, 2)]
    # Define the R1 and R2 suffix pairs for reference
    ending_dictionary = dict(zip(base_endings_r1, base_endings_r2))
    mylist = list()  # Create empty list

    # Traverse the R1 base endings to find the common ending
    for x in base_endings_r1:
        matched_ends = [
            i for i in input_fileset if i.endswith(x + file_suffix)]
        if(len(matched_ends) > 0):
            mylist.extend([x]*len(matched_ends))

    # If all samples are single-end
    if len(mylist) == 0:
        print("Your dataset appears to be entirely single-end files.")
        odd_files = [i for i in input_fileset
                     if i.endswith(".fq" + file_suffix)]
        if len(odd_files) > 0:
            old_rep = [i.replace(".fq" + suffix, ".fastq" + suffix)
                       for i in odd_files]
            rename_files(odd_files, old_rep, myfastqpath)

        # Re-assess fastq directory content and return filenames
        return([".fastq", ".fastq", getfilelist(myfastqpath)[0], 
            getfilelist(myfastqpath)[1]])

    # If R1 endings are present, check values and correct file names
    else:
        # create dictionary of mixed file endings
        mylist_endings = list(Counter(mylist).keys())
        # Find most common file ending
        myR1_suffix = max(Counter(mylist).items(), key=lambda x: x[1])[0]
        # Match chosen R1 ending to correct R2 ending
        myR2_suffix = ending_dictionary[myR1_suffix]
        # remove main R1 suffix from dictionary
        mylist_endings.remove(myR1_suffix)

        # Process forward reads
        if len(mylist_endings) > 0:
            for x in mylist_endings:
                oldnames = [
                    i for i in input_fileset if i.endswith(x + file_suffix)]
                old_rep = [i.replace(x, myR1_suffix) for i in oldnames]
                rename_files(oldnames, old_rep, myfastqpath)

            mylist = list()  # Create empty list to hold R2 file endings
            # Traverse the R2 base endings to endings
            for x in base_endings_r2:
                matched_ends = [i for i in input_fileset if i.endswith(
                    x + file_suffix) and x != myR2_suffix]
                if(len(matched_ends) > 0):
                    mylist.append(x)  # Create list of R2 files to be renamed
            if len(mylist) > 0: # Rename R2 files that don't match desired ending
                for x in mylist:
                    oldnames = [
                        i for i in input_fileset if i.endswith(x + file_suffix)]
                    old_rep = [i.replace(x, myR2_suffix) for i in oldnames]
                    rename_files(oldnames, old_rep, myfastqpath)

        # Re-assess file names
        if file_suffix == ".gz":
            input_fileset = getfilelist(myfastqpath)[1]
        else:
            input_fileset = getfilelist(myfastqpath)[0]

        # Now process single end files
        # Identify files that do not match the current R1, R2 ending
        odd_files = [i for i in input_fileset if not
                     i.endswith(myR1_suffix + file_suffix) if not
                     i.endswith(myR2_suffix + file_suffix)]

        # Partition single end files according to ending
        fastq_odd_1 = [i for i in odd_files if i.endswith(
            ".fastq" + file_suffix)]
        fastq_odd_2 = [i for i in odd_files if i.endswith(".fq" + file_suffix)]

        # If any apparently single-end files exist, then rename them
        if len(odd_files) > 0:
            print("Now unifying " + str(len(odd_files)) +
                  " single-end files to \"" + myR1_suffix +
                  file_suffix + "\" ending")
            # rename 'fastq' single-end files to correct ending
            if len(fastq_odd_1) > 0:
                old_rep = [i.replace(".fastq" + file_suffix,
                            myR1_suffix + file_suffix) for i in fastq_odd_1]
                rename_files(fastq_odd_1, old_rep, myfastqpath)
            # rename 'fq' single-end files to correct ending
            if len(fastq_odd_2) > 0:
                old_rep = [i.replace(".fq" + file_suffix,
                            myR1_suffix + file_suffix) for i in fastq_odd_2]
                rename_files(fastq_odd_2, old_rep, myfastqpath)
        # Re-assess and return filenames and file endings
        return([myR1_suffix, myR2_suffix, getfilelist(myfastqpath)[0], 
            getfilelist(myfastqpath)[1]])

def create_inputs(config):
    """
    Creates the fastq file inputs needed for read trimming steps of 
    the snakemake pipeline
    Input Parameter: 
    config (dict): Dictionary derived from config.yaml and any 
    additional key:value pairs added during the file preperation steps. 
    Returns:
    list: List of two strings; 
    1st string denotes the forward read
    2nd string denotes the reverse read
    """
    return([os.path.join("fastq","{sample}"+expand("{ending}{suffix}", \
        ending=R1_file_ending, suffix=suffix)[0]+""),
        os.path.join("fastq","{sample}"+expand("{ending}{suffix}", \
            ending=R2_file_ending, suffix=suffix)[0]+"")])

## End helper functions

# Retrieve the list of fastq files
onlyfiles, gzfiles = getfilelist(myfastqpath)

# Raise exception if no fastq files present
if len(gzfiles) == 0 and len(onlyfiles) == 0:
    raise NameError(
        "You do not seem to have any fastq files present to process. Exiting.")

# Raise exception if fastq files are a mixture of gzipped and non-gzipped files
if len(gzfiles) > 0 and len(gzfiles) != len(onlyfiles):
    myinput = "You have a mixture of gzipped files and non-gzipped files\n \
                Only {} of total {} files are gzipped!"
    raise NameError(print(myinput.format(len(gzfiles), len(onlyfiles))))

# Unify fastq file endings and return the final ending to be used.
if len(gzfiles) > 0:
    R1_file_ending, R2_file_ending, onlyfiles, gzfiles = \
            fix_input_files(".gz", gzfiles, myfastqpath)
    suffix = ".gz"
else:
    R1_file_ending, R2_file_ending, onlyfiles, gzfiles = \
            fix_input_files("", onlyfiles, myfastqpath)
    suffix = ""

sample_string = os.path.join("fastq","{sample}"+R1_file_ending+suffix)
SAMPLES, = glob_wildcards(sample_string)

# Check the file pairing
# Raise exception for non-paired PE files
len_r1 = len([i for i in onlyfiles if i.endswith(R1_file_ending + suffix)])
if len_r1*2 != len(onlyfiles):
    myinput = "One or more samples do not have a read pair!\nIf using \
        paired-end samples, please ensure each sample has read 1 and \
        read 2 files\nAborting..."
    raise NameError(myinput)  # Raise exception to break workflow

if isfile(config["human_viral_fusion"]) and not isfile("gmap_index/gmap_fusion/gmap_fusion.ref061regiondb"):
    gmap_index_needed="TRUE"        
    print("GMAP index will be created in ./gmap_index/gmap_fusion")
elif isfile("gmap_index/gmap_fusion/gmap_fusion.ref061regiondb"):
    gmap_index_needed="FALSE"
    print("GMAP index available in ./gmap_index/gmap_fusion")
    Path("output/gmap_build_complete").touch()
else:
    print("Fasta file specified for creating GMAP index was invalid. Now exiting...")
    raise NameError(print(config["genome_fasta"]))


# myinput = []
# myinput.append(expand('output/spades_{sample}/transcripts.sam',sample=SAMPLES))
# myinput.append(expand('output/trim_fastq/{sample}_stage2_R1.fq.gz',sample=SAMPLES))
# print(myinput)

# Generate input rule for Snakemake
rule all:
    input:
        expand('output/spades_{sample}/transcripts.sam',sample=SAMPLES)

# Trim adaptors from paired-end fastq files
rule trim_fastq_fastqc:
    input:
        pair1 = create_inputs(config)[0],
        pair2 = create_inputs(config)[1]
    output:
        trimmed_pair1 = "output/trim_fastq/{sample}_R1_trimmed.fq.gz",
        trimmed_pair2 = "output/trim_fastq/{sample}_R2_trimmed.fq.gz",
        fastqc_zipfile1 = "output/fastqc/{sample}_R1_fastqc.zip",
        fastqc_zipfile2 = "output/fastqc/{sample}_R2_fastqc.zip"
    log:
        "output/logs/{sample}.trim_adapters.log"
    params:
        suffix
    run:
        shell("mkdir -p output/temp_dir")
        shell("cp {input.pair1} \
            output/temp_dir/{wildcards.sample}_R1.fq{params}")
        shell("cp {input.pair2} \
            output/temp_dir/{wildcards.sample}_R2.fq{params}")
        shell("trim_galore \
            output/temp_dir/{wildcards.sample}_R1.fq{params} \
            output/temp_dir/{wildcards.sample}_R2.fq{params} \
            --paired --gzip --basename {wildcards.sample} \
            -o output/trim_fastq")
        shell("fastqc output/temp_dir/{wildcards.sample}_R1.fq{params} \
            output/temp_dir/{wildcards.sample}_R2.fq{params} \
            -o ./output/fastqc")
        shell("mv output/trim_fastq/{wildcards.sample}_val_1.fq.gz \
            output/trim_fastq/{wildcards.sample}_R1_trimmed.fq.gz"),
        shell("mv output/trim_fastq/{wildcards.sample}_val_2.fq.gz \
            output/trim_fastq/{wildcards.sample}_R2_trimmed.fq.gz")

# Align paired-end trimmed reads to genome
rule fastq_to_bam_human:
    input:
        trimmed_pair1 = "output/trim_fastq/{sample}_R1_trimmed.fq.gz",
        trimmed_pair2 = "output/trim_fastq/{sample}_R2_trimmed.fq.gz"
    output:
        "output/bam/{sample}.sam"
    params:
        index = config["human_index"]
    threads: config["threads_for_alignment"]
    log:
        "output/logs/{sample}.alignment.log"
    run:
        shell("hisat2 -p {threads} -x {params.index} \
            -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} \
            -S output/bam/{wildcards.sample}.sam 2> {log}")

rule sort_human_bam:
    input:
        "output/bam/{sample}.sam"
    output:
        bam = "output/bam/{sample}_human.bam",
        bambai = "output/bam/{sample}_human.bam.bai"
    run:
        shell("samtools sort -O bam -o output/bam/{wildcards.sample}_human.bam \
            {input}")
        shell("samtools index {output.bam}")
        shell("rm output/bam/{wildcards.sample}.sam")        


rule fastq_to_bam_viral:
    input:
        trimmed_pair1 = "output/trim_fastq/{sample}_R1_trimmed.fq.gz",
        trimmed_pair2 = "output/trim_fastq/{sample}_R2_trimmed.fq.gz"
    params:
        index = config["viral_index"]
    output:
        bam = "output/bam/{sample}_viral.sam"
    threads: config["threads_for_alignment"]
    log:
        "output/logs/{sample}.alignment.log"
    run:
        shell("hisat2 -p {threads} -x {params.index} \
            -1 {input.trimmed_pair1} -2 {input.trimmed_pair2} \
            -S output/bam/{wildcards.sample}_viral.sam 2> {log}")

rule extract_human_reads:
    input:
        "output/bam/{sample}_human.bam"
    output:
        v1 = "output/temp_dir/{sample}_human_unmap_pair_reads.txt",
        v2 = "output/temp_dir/{sample}_human_mapped_singleton_reads.txt"
    run:
        shell("samtools view -F 4 {input} | awk '{{print $1}}' | sort | uniq > output/temp_dir/{wildcards.sample}_human_mapped_reads.txt") # identify reads that aligned concordant or discordantly to human genome
        shell("samtools view -F 12 {input} | awk '{{print $1}}' | sort | uniq > output/temp_dir/{wildcards.sample}_human_mapped_pair_reads.txt") # identify read pairs that aligned concordant or discordantly to human genome
        shell("samtools view -F 4 -f 8 {input} | awk '{{print $1}}' | sort | uniq > output/temp_dir/{wildcards.sample}_human_mapped_singleton_reads.txt") # identify singletons in human
        shell("samtools view -f 12 {input} | awk '{{print $1}}' | sort | uniq > output/temp_dir/{wildcards.sample}_human_unmap_pair_reads.txt") # identify unmapped reads with unmapped pairs in human

rule extract_viral_reads:
    input:
        "output/bam/{sample}_viral.sam"
    output:
        v1 = "output/temp_dir/{sample}_viral_unmap_pair_reads.txt",
        v2 = "output/temp_dir/{sample}_viral_mapped_singleton_reads.txt"
    run:
        shell("samtools view -F 4 {input} | awk '{{print $1}}' | sort | uniq > output/temp_dir/{wildcards.sample}_viral_mapped_reads.txt") # identify reads that aligned concordant or discordantly to human genome
        shell("samtools view -F 12 {input} | awk '{{print $1}}' | sort | uniq > output/temp_dir/{wildcards.sample}_viral_mapped_pair_reads.txt") # identify read pairs that aligned concordant or discordantly to human genome
        shell("samtools view -F 4 -f 8 {input} | awk '{{print $1}}' | sort | uniq > output/temp_dir/{wildcards.sample}_viral_mapped_singleton_reads.txt") # identify singletons in human
        shell("samtools view -f 12 {input} | awk '{{print $1}}' | sort | uniq > output/temp_dir/{wildcards.sample}_viral_unmap_pair_reads.txt") # identify unmapped reads with unmapped pairs in human

rule combine_unmapped_reads:
    input:
        human="output/temp_dir/{sample}_human_unmap_pair_reads.txt",
        viral="output/temp_dir/{sample}_viral_unmap_pair_reads.txt"
    output:
        "output/temp_dir/{sample}_both_unmap_pair_reads.txt"
    run:
        shell("comm -12 {input.human}  {input.viral} > output/temp_dir/{wildcards.sample}_both_unmap_pair_reads.txt")

rule create_whitelist:
    input:
        human="output/temp_dir/{sample}_human_mapped_singleton_reads.txt",
        viral="output/temp_dir/{sample}_viral_mapped_singleton_reads.txt",
        unmap="output/temp_dir/{sample}_both_unmap_pair_reads.txt"
    output:
        "output/temp_dir/{sample}_whitelist_sorted.txt"
    run:
        shell("cat {input.human} {input.viral} {input.unmap} > output/temp_dir/{wildcards.sample}_whitelist.txt")
        shell("cat output/temp_dir/{wildcards.sample}_whitelist.txt | sort | uniq > output/temp_dir/{wildcards.sample}_whitelist_sorted.txt")            

rule filter_whitelist:
    input:
        whitelist="output/temp_dir/{sample}_whitelist_sorted.txt",
        r1 = "output/trim_fastq/{sample}_R1_trimmed.fq.gz",
        r2 = "output/trim_fastq/{sample}_R2_trimmed.fq.gz"
    output:
        r1 = "output/trim_fastq/{sample}_stage2_R1.fq.gz",
        r2 = "output/trim_fastq/{sample}_stage2_R2.fq.gz"
    run:
        shell("filterbyname.sh in={input.r1} in2={input.r2} out={output.r1} out2={output.r2} names={input.whitelist} include=t")

rule spades:
    input:
        r1 = "output/trim_fastq/{sample}_stage2_R1.fq.gz",
        r2 = "output/trim_fastq/{sample}_stage2_R2.fq.gz"
    output:
        "output/spades_{sample}/transcripts.fasta"              
    run:
        shell("rnaspades.py -t 20 -m 140 -o output/spades_{wildcards.sample} \
            -1 output/trim_fastq/{wildcards.sample}_stage2_R1.fq.gz \
            -2 output/trim_fastq/{wildcards.sample}_stage2_R2.fq.gz")

if gmap_index_needed=="TRUE":
    rule create_gmap_index:
        input:
            fusion_fa=expand("{fa}", fa=config["human_viral_fusion"])[0]
        output:
            "output/gmap_build_complete"
        run:
            shell("gmap_build -d gmap_fusion -D gmap_index {input}")
            shell("touch output/gmap_build_complete")

rule map_denovo:
    input:
        fa="output/spades_{sample}/transcripts.fasta",
        gmap="output/gmap_build_complete"
    output:
        "output/spades_{sample}/transcripts.sam"
    params:
        config["threads_for_alignment"]
    run:
        shell("gmap -D gmap_index -d gmap_fusion -f samse -n 0 -t {params} --max-intronlength-ends 200000 -z sense_force {input.fa} > output/spades_{wildcards.sample}/transcripts.sam -x 20")

