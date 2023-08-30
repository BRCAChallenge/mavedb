#!/usr/bin/env python

from cool_seq_tool.data_sources import AlignmentMapper, TranscriptMappings, \
    UTADatabase
from cool_seq_tool.schemas import Assembly, ResidueMode
from cool_seq_tool.data_sources import SeqRepoAccess
from cool_seq_tool.paths import TRANSCRIPT_MAPPINGS_PATH, LRG_REFSEQGENE_PATH, \
    MANE_SUMMARY_PATH
from cool_seq_tool.data_sources.uta_database import UTA_DB_URL

from biocommons.seqrepo import SeqRepo

import argparse
import asyncio
from glob import glob
import json
import os
import re
import requests

def alignment_mapper_setup():
    sr = SeqRepo(root_dir=os.environ.get('SEQREPO_ROOT_DIR'))
    seqrepo_access = SeqRepoAccess(sr)
    db_url: str = UTA_DB_URL
    db_pwd: str = "uta"
    uta_db = UTADatabase(db_url=db_url)
    transcript_mappings=TranscriptMappings(TRANSCRIPT_MAPPINGS_PATH)
    alignment_mapper = AlignmentMapper(seqrepo_access, transcript_mappings,
                                       uta_db)
    return(alignment_mapper)



def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_json", nargs='*',
    			help="Input file(s), JSON format")
    parser.add_argument("-n", "--track_name", default="Variant Effect Maps",
                        help="Name of the composite track")
    parser.add_argument('-t', "--trackDb", default="trackDb.txt",
                        help="Output trackDb file")
    parser.add_argument('-b', "--bed_dir", help="output bed directory")
    parser.add_argument('-l', "--location_matrix_dir",
                        help="output location matrix directory")
    #parser.add_argument('-w', "--wig_dir", default="wig",
    #                    help="Output wig file directory")
    parser.add_argument('-d', "--debug", type=bool, default=True)
    args = parser.parse_args()
    input_files = []
    for arg in args.input_json:
        input_files += glob(arg)
    args.input_files = input_files
    return(args)

def print_wig_header(output_wig_fp):
    output_wig_fp.write("track type=wiggle_0 name=\"mavedb\"",
                        "description=\"MaveDB Scores\"")

def write_header(track_name, trackdb_fp):
    trackdb_fp.write("track %s\n" % track_name)
    trackdb_fp.write("compositeTrack on\n")
    trackdb_fp.write("\n")
    
    
def add_trackdb_entry(mave_json, track_name, score_set_name, url_base,
                      trackdb_fp):
    """
    Make an entry for this MAVE experiment in trackDb
    """
    trackdb_fp.write("\ttrack %s\n" % score_set_name)
    trackdb_fp.write("\tparent %s on\n" % track_name)
    trackdb_fp.write("\tbigDataUrl %s/%s.bw\n" % (url_base, score_set_name))
    trackdb_fp.write("\tshortLabel %s\n" % mave_json["urn"])
    trackdb_fp.write("\tlongLabel %s\n" % mave_json["title"])
    trackdb_fp.write("\ttype bigWig\n")
    trackdb_fp.write("\tvisibility dense\n")
    trackdb_fp.write("\n")


def get_endpoints(mapped_score):
    """
    Given a mapped score, return its start and end endpoints
    """
    variant = mapped_score["variation"]
    assert(variant["type"] == "Allele")
    location = variant["location"]
    assert(location["type"] == "SequenceLocation")
    interval = location["interval"]
    assert(interval["type"] == "SequenceInterval")
    assert(interval["start"]["type"] == "Number")
    start_pos = interval["start"]["value"]
    assert(interval["end"]["type"] == "Number")
    end_pos = interval["end"]["value"]
    return(start_pos, end_pos)
           

def map_to_genome(alignment_mapper, reference_type, reference_accession,
                  start_pos, end_pos):
    """
    Given a reference type and accession, and start and end positions on that
    accessioned sequence, map the interval to the genome.  Return the chrom,
    genomic start and genomic end
    """
    assert(reference_type == "dna" or reference_type == "protein")
    if reference_type == "dna":
        genomic_start = start_pos
        genonic_end = end_pos
        chrom = reference_accession
    elif reference_type == "protein":
        query_url = ("http://127.0.0.1:8000/cool_seq_tool/alignment_mapper/p_to_g" +
                     "?p_ac=%s&p_start_pos=%d&p_end_pos=%d&residue_mode=inter-residue" +
                     "&target_genome_assembly=GRCh38") \
                     % (reference_accession, start_pos, end_pos)
        response = requests.get(query_url)
        resp = response.json()
        if "g_data" in resp:
            chrom = resp["g_data"]["g_ac"]
            genomic_start = resp["g_data"]["g_start_pos"]
            genomic_end = resp["g_data"]["g_end_pos"]
        else:
            chrom = None
            genomic_start = None
            genomic_end = None
    return(chrom, genomic_start, genomic_end)

    
async def old_map_to_genome(alignment_mapper, reference_type, reference_accession,
                        start_pos, end_pos):
    """
    Given a reference type and accession, and start and end positions on that
    genomic start and genomic end
    """
    assert(reference_type == "dna" or reference_type == "protein")
    if reference_type == "dna":
        genomic_start = start_pos
        genonic_end = end_pos
        chrom = reference_accession
    elif reference_type == "protein":
        resp, w = await alignment_mapper.p_to_g(reference_accession, start_pos,
                                                end_pos,
                                                residue_mode=ResidueMode.INTER_RESIDUE,
                                                target_genome_assembly=Assembly.GRCH38)
        chrom = resp["g_data"]["g_ac"]
        genomic_start = resp["g_data"]["g_start_pos"]
        genomic_end = resp["g_data"]["g_end_pos"]
    return(chrom, genomic_start, genomic_end)


def score_entry_to_label(post_mapped_score):
    (ref_start, ref_end) = get_endpoints(post_mapped_score)
    label = "%d_%d" % (ref_start, ref_end)
    return(label)


def get_genomic_endpoints(alignment_mapper, reference_type, reference_accession,
                          post_mapped_score):
    (ref_start, ref_end) = get_endpoints(post_mapped_score)
    (chrom, genomic_start, genomic_end) = map_to_genome(alignment_mapper,
                                                        reference_type,
                                                        reference_accession,
                                                        ref_start,
                                                        ref_end)
    label = score_entry_to_label(post_mapped_score)
    return(chrom, genomic_start, genomic_end, label)


    
def add_location_to_bed(bed_fp, post_mapped_score, reference_type,
                        reference_accession, alignment_mapper):
    """
    Given a mapped position, add it to the bed file, after translating the
    coordinates to genomic coordinates
    """

    
    
def create_bed(mave_json, score_set_name,
               alignment_mapper, bed_dir):
    """
    Create a bed file that shows the positions scored by genomic coordinates,
    plus a location matrix file that shows the intensity of the score
    """
    bed_filename = "%s/%s.bed" % (bed_dir, score_set_name)
    with open(bed_filename, "w") as bed_fp:
        #bed_fp.write("track type=bed name=\"%s\" description=\"%s\"\n" %
        #             ( score_set_name, mave_json["title"]))
        reference_type = mave_json["mapped_reference_sequence"]["sequence_type"]
        reference_accession = mave_json["mapped_reference_sequence"]["sequence_id"]
        labels_in_bed = {}
        for mapped_score in mave_json["mapped_scores"]:
            (chrom_name, start_pos,
             end_pos, label) = get_genomic_endpoints(alignment_mapper,
                                                     reference_type,
                                                     reference_accession,
                                                     mapped_score["post_mapped"])
            print("working on", chrom_name, start_pos, end_pos, label)
            if not label in labels_in_bed:
                this_chrom_versioned = re.sub("NC_(0)+", "chr", chrom_name)
                this_chrom = re.split("\.", this_chrom_versioned)[0]
                bed_fp.write("%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t0,0,0\n"
                             % (this_chrom, start_pos, end_pos, label,
                                start_pos, end_pos))
                labels_in_bed[label] = 1

        

def gather_unique_alts(mave_json):
    """
    Generate a hash in which each key is an observed alternative allele
    and the value is a numeric index
    """
    unique_alts = {}
    unique_alts_observed = 0
    for mapped_score in mave_json["mapped_scores"]:
        this_alt = mapped_score["post_mapped"]["variation"]["state"]
        this_alt_seq = this_alt["sequence"]
        if not this_alt_seq in unique_alts:
            unique_alts[this_alt_seq] = unique_alts_observed
            unique_alts_observed += 1
    return(unique_alts)


def assemble_scores_per_label(mave_json, alts_scored):
    scores_per_label = {}
    min_score = None
    max_score = None
    for mapped_score in mave_json["mapped_scores"]:
        this_score = mapped_score["score"]
        if min_score is None or this_score < min_score:
            min_score = this_score
        if max_score is None or this_score > max_score:
            max_score = this_score
        locus = mapped_score["post_mapped"]
        locus_label = score_entry_to_label(locus)
        if not locus_label in scores_per_label:
            scores_per_label[locus_label] = [None] * len(alts_scored)
        this_alt = mapped_score["post_mapped"]["variation"]["state"]["sequence"]
        this_alt_index = alts_scored[this_alt]
        scores_per_label[locus_label][this_alt_index] = this_score
    print("scores per label", scores_per_label)
    return(scores_per_label, min_score, max_score)

                
def create_location_matrix(mave_json, score_set_name, location_matrix_dir):
    """
    Create a location matrix file, that describes the scores at each
    position scored
    """
    #
    # Assemble and output a vector of unique "states", or alt residues
    alts_scored = gather_unique_alts(mave_json)
    #
    # Assemble a data structure listing, per locus, the scores observed
    # for each alt at that locus
    (scores_per_label, min_score,
     max_score) = assemble_scores_per_label(mave_json, alts_scored)
    location_matrix_filename = "%s/%s.lm" % (location_matrix_dir,
                                             score_set_name)
    with open(location_matrix_filename, "w") as lm_fp:
        #
        # Write out the unique alts
        lm_fp.write("gene")
        alt_index = 0
        for this_alt in alts_scored.keys():
            lm_fp.write("\t%s" % this_alt)
            assert(alts_scored[this_alt] == alt_index)
            alt_index += 1
        lm_fp.write("\n")
        #
        # For each locus scored, write out the score per alt
        for label in scores_per_label.keys():
            lm_fp.write(label)
            these_scores = scores_per_label[label]
            for index in range(len(alts_scored)):
                this_score = these_scores[index]
                if this_score is None:
                    lm_fp.write("\t0.0")
                else:
                    lm_fp.write("\t%f" % this_score)
            lm_fp.write("\n")
    return(min_score, max_score)
            
    
def mave_to_track(mave_json, track_name, trackdb_fp, alignment_mapper,
                  bed_dir, location_matrix_dir):
    """
    Add WIG content to the output wig from the specified json object
    """
    # Extract the score set name from the URN.  This is the last token of
    # a colon-delimited URN.  This will be the name for instances where
    # reserved characters don't work well, such as filenames
    score_set_name = re.split(":", mave_json["target"]["scoreset"])[-1]
    #add_trackdb_entry(mave_json, track_name, score_set_name, bedDir,
    #                  trackdb_fp)
    create_bed(mave_json, score_set_name, alignment_mapper, bed_dir)
    (min_score,max_score) = create_location_matrix(mave_json, score_set_name,
                                                   location_matrix_dir)
    return(min_score, max_score)
    

def main():
    alignment_mapper = alignment_mapper_setup()
    args = parse_args()
    if not os.path.exists(args.bed_dir):
        os.makedirs(args.bed_dir)
    if not os.path.exists(args.location_matrix_dir):
        os.makedirs(args.location_matrix_dir)
    with open(args.trackDb, "w") as trackdb_fp:
        write_header(args.track_name, trackdb_fp)
        for this_input_file in args.input_files:
            if args.debug:
                print("Working on input file", this_input_file)
                mave_json = json.load(open(this_input_file))
                (min_score, max_score) = mave_to_track(mave_json,
                                                       args.track_name,
                                                       trackdb_fp,
                                                       alignment_mapper,
                                                       args.bed_dir,
                                                       args.location_matrix_dir)
    print(args.input_files, "score_range", min_score, max_score)

if __name__ == "__main__":
    main()
    
