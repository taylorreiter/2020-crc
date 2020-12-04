import pandas as pd

m = pd.read_csv("inputs/metadata.tsv", sep = "\t", header = 0)
ACCESSIONS = m.sort_values(by='read_count')['run_accession']
SAMPLES = m['sample_alias'].unique().tolist()
STUDY = m['study_accession'].unique().tolist()

rule all:
    input:
        expand("inputs/cat/{sample}_2.fastq.gz", sample = SAMPLES),
        expand("inputs/cat/{sample}_2.fastq.gz", sample = SAMPLES)

rule download_fastq_files_R1:
    output: 
        r1="inputs/raw/{accession}_1.fastq.gz",
    run:
        row = m.loc[m['run_accession'] == wildcards.accession]
        fastq_1 = row['fastq_ftp_1'].values
        fastq_1 = fastq_1[0]
        shell("wget -O {output.r1} {fastq_1}")


rule download_fastq_files_R2:
    output:
        r2="inputs/raw/{accession}_2.fastq.gz"
    run:
        row = m.loc[m['run_accession'] == wildcards.accession]
        fastq_2 = row['fastq_ftp_2'].values
        fastq_2 = fastq_2[0]
        shell("wget -O {output.r2} {fastq_2}")


rule cat_libraries_R1:
    input: expand("inputs/raw/{accession}_1.fastq.gz", accession = ACCESSIONS)
    output: expand("inputs/cat/{sample}_1.fastq.gz", sample = SAMPLES)
    run: 
        merge_df = m[['sample_alias','run_accession']]
        merge_df = copy.deepcopy(merge_df)
        merge_df['run_accession'] = merge_df['run_accession'].apply(lambda x: f"inputs/raw/{x}_1.fastq.gz")
        merge_dict = merge_df.groupby('sample_alias')['run_accession'].apply(lambda g: g.values.tolist()).to_dict()
        for sample in merge_dict.keys():
            # merge SRR files
            to_merge = merge_dict[sample]
            # Check if the merged file results from a single or multiple fastq files.
            # For n-to-1 merging, concatenate input files to produce the output file
            merge_nb = len(to_merge)
            if merge_nb > 1:
                cmd = "cat " + " ".join(to_merge) + " > " + "inputs/cat/" + sample + "_1.fastq.gz"
            else:
                cmd = "ln --relative --force -s " + " ".join(to_merge) + " inputs/cat/" + sample + "_1.fastq.gz"
            os.system(cmd)
    
rule cat_libraries_R2:
    input: expand("inputs/raw/{accession}_2.fastq.gz", accession = ACCESSIONS)
    output: expand("inputs/cat/{sample}_2.fastq.gz", sample = SAMPLES)
    run: 
        merge_df = m[['sample_alias','run_accession']]
        merge_df = copy.deepcopy(merge_df)
        merge_df['run_accession'] = merge_df['run_accession'].apply(lambda x: f"inputs/raw/{x}_2.fastq.gz")
        merge_dict = merge_df.groupby('sample_alias')['run_accession'].apply(lambda g: g.values.tolist()).to_dict()
        for sample in merge_dict.keys():
            # merge SRR files
            to_merge = merge_dict[sample]
            # Check if the merged file results from a single or multiple fastq files.
            # For n-to-1 merging, concatenate input files to produce the output file
            merge_nb = len(to_merge)
            if merge_nb > 1:
                cmd = "cat " + " ".join(to_merge) + " > " + "inputs/cat/" + sample + "_2.fastq.gz"
            else:
                cmd = "ln --relative --force -s " + " ".join(to_merge) + " inputs/cat/" + sample + "_2.fastq.gz"
            os.system(cmd)
