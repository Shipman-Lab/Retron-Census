#!/usr/bin/env python3
'''
short_insertion_caller.py

This script is designed to analyze sequencing data to identify and categorize genetic edits and errors from insertion editing data. It processes FASTQ files to detect genetic variants such as insertions, deletions, and substitutions, and categorizes them based on the analysis results.

Usage:
    python3 short_insertion_caller.py --prefix <prefix> --cores <cores> --params <params_file> --paths <paths_file> --minlen <minimum_length>

Arguments:
    --prefix    : Prefix for the output file names.
    --cores     : Number of CPU cores for parallel processing (default is 6).
    --params    : Path to the JSON file containing parameters for the analysis.
    --paths     : Path to the INI configuration file specifying various paths used in the script.
    --minlen    : Minimum read length cutoff for quality control (default is 100).

The script performs the following functions:
1. Loads and preprocesses FASTQ data, supporting gzip compressed files.
2. Employs fuzzy matching for variant calling to handle mismatches and indels.
3. Classifies reads into categories: wild type, edited, indel, and various undetermined outcomes.
4. Outputs detailed counts of outcomes, base counts, variant data, and categorized FASTA files.

Outputs:
- FASTA files of unique sequences categorized by read outcomes.
- Pickle files with detailed edit outcomes, base counts, and error information.
'''

import sys, os
import gzip
import multiprocessing
import time
import pickle
import argparse, configparser, json
from functools import partial
import uuid

import pandas as pd
import numpy as np

from collections import Counter
import fuzzysearch
import difflib
import math

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def trim_primer(qseq, primer, primer_pos):
    five_trim = 0
    three_trim = 0
    
    # Search and trim forward primer
    fwd_primer = primer['fwd']
    match_out = fuzzysearch.find_near_matches(fwd_primer, qseq, max_l_dist=math.floor(len(fwd_primer)/5))
    
    if len(match_out) == 0:
        pass
    elif len(match_out) == 1:
        match = match_out[0]
        five_trim = match.end
        qseq = qseq[five_trim:] # follows 0-base
    else:
        for match in match_out:
            if (match.end > primer_pos['fwd'][1]-5) and (match.end < primer_pos['fwd'][1]+15): # can change numbers later
                five_trim = match.end
                qseq = qseq[five_trim:]
     
    # If read length is longer than end of primer, then search and trim reverse primer
    if len(qseq) >= primer_pos['rev'][1]:
        pass
    else:
        # sometimes reads containing reverse primers are truncated so shorter than expected so turn this part off for now
#         return qseq, five_trim, three_trim
        pass
    
    # Search and trim reverse primer
    rev_primer = primer['rev']
    match_out = fuzzysearch.find_near_matches(rev_primer, qseq, max_l_dist=math.floor(len(rev_primer)/5)) 
    
    if len(match_out) == 0:
        pass
    elif len(match_out) == 1: 
        match = match_out[0]
        three_trim = match.end - match.start
        qseq = qseq[:match.start]
    else:
        for match in match_out:
            if (match.start > primer_pos['rev'][0]-15) and (match.start < primer_pos['rev'][0]+5):
                # these numbers are the most variable -- there are ones that start around 156 vs 189
                three_trim = match.end - match.start
                qseq = qseq[:match.start]
            else:
                pass
        
    return qseq, five_trim, three_trim


def find_fidelity(tseq, qseq):
    var_keys = ['type', 't_pos', 'q_pos', 't_idt', 'q_idt']
    var_dict = {k: [] for k in var_keys}

    s = difflib.SequenceMatcher(None, tseq, qseq)
    for tag, i1, i2, j1, j2 in s.get_opcodes(): # i is template position; j is query position
        if tag == 'equal':
            continue # add nothing to variants
        elif tag == 'delete':
            if j1 != j2:
                raise Exception(f'query start and stop of deletion must be equal')
            if j1 >= len(qseq):
                # deletions at the very end of qseq = read stops before end of tseq
                continue # stop counting
            elif j1 == 0:
                # deletions at the very beginning of qseq could be because 5' end of read got trimmed
                continue
            # add nothing to coverage
            for i in range(i1, i2):   
                var_dict['type'].append('del')
                var_dict['t_pos'].append(i)
                var_dict['q_pos'].append(j1) # j is the q_position AFTER the query-side deletion
                var_dict['t_idt'].append(tseq[i])
                var_dict['q_idt'].append('d') # placeholder string to represent deletion on q side
        elif tag == 'insert':
            if i1 != i2:
                raise Exception(f'template start and stop of insertion must be equal')
            if i1 >= len(tseq): # insertions at the very end of tseq = read length continues past tseq
                continue # stop counting
            for j in range(j1, j2):
                var_dict['type'].append('ins')
                var_dict['t_pos'].append(i1)
                var_dict['q_pos'].append(j)
                var_dict['t_idt'].append('i') # placeholder string to represent insertion on t side
                var_dict['q_idt'].append(qseq[j])
        elif tag == 'replace':
            for i,j in zip(range(i1, i2), range(j1, j2)):
                var_dict['type'].append('swp')
                var_dict['t_pos'].append(i)
                var_dict['q_pos'].append(j)
                var_dict['t_idt'].append(tseq[i])
                var_dict['q_idt'].append(qseq[j])
        else:
            raise Exception('no other tag options')
       
    variants = pd.DataFrame.from_dict(var_dict, orient='columns')
    return variants


def match_and_sort(qseq_in, min_length, params):
    # Trim primers
    qseq, five_trim, three_trim = trim_primer(qseq_in, params['primer'], params['primer_pos'])
    
    variants = None
    if len(qseq) < min_length: # min read length after trimming
        outcome = 'undetermined_bad_qual'
        return outcome, variants, qseq
    
    if five_trim == 0: # throw away reads without fwd primer
        outcome = 'undetermined_bad_qual'
        return outcome, variants, qseq
    
    # Match and sort to wt, indel, or ed
    left_flank = fuzzysearch.find_near_matches(params["L_search_seq"],
                                               qseq,
                                               max_l_dist=4)
    right_flank = fuzzysearch.find_near_matches(params["R_search_seq"],
                                                qseq,
                                                max_l_dist=4)
    
    if len(left_flank) != 1:
        outcome = 'undetermined_no_flanking_match'
        return outcome, variants, qseq
    if len(right_flank) != 1:
        outcome = 'undetermined_no_flanking_match'
        return outcome, variants, qseq
    
    left_match = left_flank[0]
    right_match = right_flank[0]
    region = qseq[left_match.end:right_match.start]

    match_out = [fuzzysearch.find_near_matches(params[var], 
                                               region,
                                               max_l_dist=math.floor(len(var)/5))
                                              for var in ['wt_seq', 'ed_seq']]
    if match_out[0] and match_out[1]: # fuzzysearch should not match both wt and ed sequences
        outcome = 'undetermined_no_site_match'
        return outcome, variants, qseq
    elif match_out[0]: # match wt
        if len(match_out[0]) > 1:
            print(f'fuzzysearch has more than 1 match', '\n', match_out[0], '\n', qseq, flush=True)
        outcome = 'wt'
        if len(region) > 25:
            outcome = 'wt_seq_artifact'
            return outcome, variants, qseq
    elif match_out[1]: # match ed
        if len(match_out[1]) > 1:
            print(f'fuzzysearch has more than 1 match', '\n', match_out[1], '\n', qseq, flush=True)
        outcome = 'ed'
        if len(region) > 35:
            outcome = 'ed_seq_artifact'
            return outcome, variants, qseq
    else: # if read doesn't match either, both empty
        outcome = 'undetermined_no_site_match'
        return outcome, variants, qseq
    
    variants = find_fidelity(params[f'{outcome}_seq'], region)
    variants['read_outcome'] = outcome
    if len(variants) == 0: # outcome cannot be indel if no variants
        variants = None
        return outcome, variants, qseq
    
    # Determine if wt matching reads are wt or indel
    if outcome == 'wt':
        if len(set(variants['t_pos'].to_list()) & set(params['cut_pos'])) > 0:
            outcome = 'indel'
            variants = None # don't keep variants for indel

    return outcome, variants, qseq
    

def globalize(func):
    def result(*args, **kwargs):
        return func(*args, **kwargs)
    result.__name__ = result.__qualname__ = uuid.uuid4().hex
    setattr(sys.modules[result.__module__], result.__name__, result)
    return result


def find_edits(min_length, params, outputs, inputs):
    read_counts, nbase_counts, seqs_to_fasta, list_variants = outputs
    reads, (seq, n) = inputs
    if len(reads) != n:
        raise Exception("Num of quality arrays doesn't match number of seq counts!")

    if len(seq) < min_length: # minimum read length
        outcome = 'undetermined_bad_qual'
    if 'A'*20 in seq: # A repeats
        outcome = 'undetermined_bad_qual'

    try:
        outcome
    except:
        outcome, variants, trimmed_seq = match_and_sort(seq, min_length, params)

    # Count read outcomes
    read_counts[outcome] = read_counts[outcome] + n
    # Count number of bases
    nbase = len(seq)
    nbase_counts[outcome] = nbase_counts[outcome] + nbase*n
    nbase_counts['total'] = nbase_counts['total'] + nbase*n
    
    seqs_to_fasta[outcome].append(SeqRecord(Seq(seq), id='', description=str(n)))

    if outcome.split('_')[0] == 'undetermined':
        return None
    if variants is None:
        return None
    
    variants['read_name'] = pd.Series([[read.id for read in reads]]*len(variants))
    variants['qual'] = variants["q_pos"].apply(lambda x: 
                                                [read.letter_annotations["phred_quality"][x-1] # for deletions, will take quality of preceding base
                                                for read in reads])
    try:
        variants = variants.explode(['read_name', 'qual'])
    except:
        print(len(reads), variants[['read_name', 'qual']], '\n', flush=True)
        raise Exception('variants not exploded')
    list_variants.append(variants)
    return None


def parallelize(reads, items, min_length, params, cores):
    # Define outputs
    edits_columns = [
        'total',
        'wt',
        'indel',
        'ed',
        'undetermined_bad_qual',
        'undetermined_no_site_match',
        'undetermined_no_flanking_match',
        'ed_seq_artifact',
        'wt_seq_artifact',
    ]
    read_counts = {column: 0 for column in edits_columns} # count number of reads
    nbase_counts = {column: 0 for column in edits_columns}
    seqs_to_fasta = {column: [] for column in edits_columns[1:]}
    list_variants = []

    # Paralellize analysis
    with multiprocessing.Manager() as manager:
        read_counts = manager.dict(read_counts)
        nbase_counts = manager.dict(nbase_counts)
        seqs_to_fasta = manager.dict({k: manager.list(v) for k, v in seqs_to_fasta.items()})
        list_variants = manager.list(list_variants)
        outputs = (read_counts, nbase_counts, seqs_to_fasta, list_variants)

        pool = multiprocessing.Pool(processes=cores)
        inputs = zip(reads, items)
        f = partial(find_edits, min_length, params, outputs)
        pool.map(f, inputs)
        pool.close()
        pool.join()

        read_counts = read_counts._getvalue()
        nbase_counts = nbase_counts._getvalue()
        seqs_to_fasta = {k: v._getvalue() for k, v in seqs_to_fasta._getvalue().items()}
        list_variants = list_variants._getvalue()
        
    # Concatenate all variant dataframes into a single dataframe
    all_variants = pd.concat(list_variants)

    return read_counts, nbase_counts, all_variants, seqs_to_fasta


def load_fastq(f_in):
    records = []
    if f_in.endswith('.fastq'):
        with open(f_in, 'r') as handle:
            for record in SeqIO.parse(handle, 'fastq'):
                records.append(record)
    else:
        with gzip.open(f_in, 'rt') as handle:
            for record in SeqIO.parse(handle, 'fastq'):
                records.append(record)

    total_n = len(records)
    read_seqs = np.array([str(i.seq) for i in records])
    read_counter = Counter(read_seqs)
    items = list(read_counter.items())
    reads_idxs = [np.where(read_seqs==seq)[0] for seq, n in items]
    reads = [[records[i] for i in idx] for idx in reads_idxs]
    return reads, items, total_n


def main(prefix, cores, min_length, paths, params):
    start = time.time()
    print(prefix, flush=True)

    # Load FASTQ
    reads, items, total_n = load_fastq(params['fastq_in'])
    # Run analysis
    edits_out, nbase_out, editerr_out, seqs_to_fasta = parallelize(reads, items, 
                                                                   min_length=min_length, params=params, cores=cores)
    edits_out['total'] = total_n
    print(edits_out, flush=True)

    # Make output drcts
    if not os.path.exists(f"{paths['notebook_drct']}/fasta_out/"):
        os.mkdir(f"{paths['notebook_drct']}/fasta_out/")
    if not os.path.exists(f"{paths['notebook_drct']}/results/"):
        os.mkdir(f"{paths['notebook_drct']}/results/")

    # Write unique sequences into FASTA by outcome
    for k,v in seqs_to_fasta.items():
        with open(f"{paths['notebook_drct']}/fasta_out/{prefix}_{k}.fasta", 'w') as handle:
            SeqIO.write(v, handle, 'fasta')

    # Write editing outcome and fidelity counts to pickle
    edits_out['total_mapped'] = edits_out['wt'] + edits_out['indel'] + edits_out['ed']
    editerr_out['sample'] = prefix
    with open(f"{paths['notebook_drct']}/results/{prefix}.edits_out.pickle", 'wb') as handle:
        pickle.dump(edits_out, handle, protocol=4) # pickle.HIGHEST_PROTOCOL 
    with open(f"{paths['notebook_drct']}/results/{prefix}.nbase_out.pickle", 'wb') as handle:
        pickle.dump(nbase_out, handle, protocol=4)
    with open(f"{paths['notebook_drct']}/results/{prefix}.editerr_out.pickle", 'wb') as handle:
        pickle.dump(editerr_out, handle, protocol=4)

    end = time.time()
    print(f"took {(end-start)/60} minutes!"+"\n"+
            "="*60 + "\n", flush=True)
    

if __name__ == '__main__':
    # Set up command-line argument parsing
    parser = argparse.ArgumentParser(description='Process arguments for short insertion caller')
    parser.add_argument('--prefix', type=str, help='Prefix for output files')
    parser.add_argument('--cores', type=int, default=6, help='Number of CPU cores to use')
    parser.add_argument('--params', type=str, default='params.json', help='Path to parameters JSON file')
    parser.add_argument('--paths', type=str, default='paths.ini',  help='Path to paths INI file')
    parser.add_argument('--minlen', type=int, default=100, help='Minimum read length for processing')
    args = parser.parse_args()

    PATHS = dict()
    config = configparser.ConfigParser(inline_comment_prefixes="#")
    config.read(args.paths)
    for key,value in config.items('paths'):
        PATHS[key] = value
    
    with open(args.params) as handle:
        PARAMS = json.load(handle)
        
    main(args.prefix, args.cores, args.minlen, PATHS, PARAMS)

