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
import shutil

from bigHeat import runBigHeat


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_json", nargs='*',
    			help="Input file(s), JSON format")
    parser.add_argument("-n", "--track_name", default="Variant_Effect_Maps",
                        help="Name of the composite track")
    parser.add_argument('-t', "--trackDb", default="hub/hg38/trackDb.txt",
                        help="Output trackDb file")
    parser.add_argument("--bed_dir", help="output bed directory")
    parser.add_argument('-c', "--coordinates_dir", default=None,
                        help="Optional output directory for coordinate mappings")
    parser.add_argument('-l', "--location_matrix_dir",
                        help="output location matrix directory")
    parser.add_argument("--hub_dir", default="./hub",
                        help="hub base directory")
    parser.add_argument("--hub_file",
                        default="src/browser_config_files/hub.txt",
                        help="template hub.txt file")
    parser.add_argument("--html", default="src/browser_config_files/description.html",
                        help="Template HTML description file")
    parser.add_argument("--genomes_file",
                        default="src/browser_config_files/genomes.txt",
                        help="template genomes.txt file")
    parser.add_argument("--genome_name", default="hg38",
                        help="Genome name: must match genomes file")
    parser.add_argument('-b', "--bigBed_dir", help="bigBed output directory")
    parser.add_argument('-s', "--chrom_sizes", help="Pathname to the chrom sizes file")
    parser.add_argument('-d', "--debug", type=bool, default=False)
    args = parser.parse_args()
    input_files = []
    for arg in args.input_json:
        input_files += glob(arg)
    args.input_files = input_files
    return(args)


def write_trackDb(track_name, trackdb_fp, html_description_filename):
    trackdb_fp.write("track %s\n" % re.sub('\s', "_", track_name))
    trackdb_fp.write("superTrack on\n")
    trackdb_fp.write("shortLabel Variant_Effect_Maps\n")
    trackdb_fp.write("longLabel This track contains variant effect maps from the MaveDB repository\n")
    trackdb_fp.write("descriptionUrl %s\n" % html_description_filename)
    trackdb_fp.write("\n")
    
def append_to_file(pathname, directive, argument):
    with open(pathname, 'a') as fp:
        fp.write("%s %s\n" % (directive, argument))

        
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
           

def map_to_genome(reference_type, reference_accession,
                  start_pos, end_pos):
    """
    Given a reference type and accession, and start and end positions on that
    accessioned sequence, map the interval to the genome.  Return the chrom,
    genomic start and genomic end
    """
    assert(reference_type == "dna" or reference_type == "protein")
    if reference_type == "dna":
        genomic_start = start_pos
        genomic_end = end_pos
        chrom = reference_accession
    elif reference_type == "protein":
        query_url = ("http://127.0.0.1:8000/cool_seq_tool/alignment_mapper/p_to_g" +
                     "?p_ac=%s&p_start_pos=%d&p_end_pos=%d&residue_mode=inter-residue" +
                     "&target_genome_assembly=GRCh38") \
                     % (reference_accession, start_pos, end_pos)
        response = requests.get(query_url)
        resp = response.json()
        chrom = None
        genomic_start = None
        genomic_end = None
        if "g_data" in resp:
            if resp["g_data"] is not None:
                chrom = resp["g_data"]["g_ac"]
                genomic_start = resp["g_data"]["g_start_pos"]
                genomic_end = resp["g_data"]["g_end_pos"]
    return(chrom, genomic_start, genomic_end)

    


def score_entry_to_label(post_mapped_score):
    (ref_start, ref_end) = get_endpoints(post_mapped_score)
    label = "%d_%d" % (ref_start, ref_end)
    return(label)


def get_genomic_endpoints(reference_type, reference_accession,
                          post_mapped_score):
    (ref_start, ref_end) = get_endpoints(post_mapped_score)
    (chrom, genomic_start, genomic_end) = map_to_genome(reference_type,
                                                        reference_accession,
                                                        ref_start,
                                                        ref_end)
    label = score_entry_to_label(post_mapped_score)
    return(chrom, genomic_start, genomic_end, ref_start, ref_end, label)



    
    
def create_bed(mave_json, score_set_name,
               bed_dir, coordinates_dir, debug=False):
    """
    Create a bed file that shows the positions scored by genomic coordinates,
    plus a location matrix file that shows the intensity of the score
    """
    bed_filename = "%s/%s.bed" % (bed_dir, score_set_name)
    with open(bed_filename, "w") as bed_fp:
        if coordinates_dir is not None:
            coordinates_filename = "%s/%s.tsv" % (coordinates_dir,
                                                  score_set_name)
            coordinates_fp = open(coordinates_filename, "w")
            coordinates_fp.write("Reference_type\tReference_Accession\tReference_Start\tReference_End" \
                                 + "\tChrom\tGehomic_Start\tGenomic_End\n")
        if "mapped_reference_sequence" in mave_json:
            if "sequence_type" in mave_json["mapped_reference_sequence"]:
                reference_type = mave_json["mapped_reference_sequence"]["sequence_type"]
                if "sequence_accessions" in mave_json["mapped_reference_sequence"]:
                    reference_accession = mave_json["mapped_reference_sequence"]["sequence_accessions"][0]
                    assert(len(mave_json["mapped_reference_sequence"]["sequence_accessions"]) == 1)
                else:
                    reference_accession = mave_json["mapped_reference_sequence"]["sequence_id"]
                labels_in_bed = {}
                for mapped_score in mave_json["mapped_scores"]:
                    (chrom_name, genomic_start, genomic_end,
                     reference_start, reference_end, label)  = get_genomic_endpoints(reference_type,
                                                                                     reference_accession,
                                                                                     mapped_score["post_mapped"])
                    if debug:
                        print("working on", chrom_name, genomic_start, genomic_end, label)
                    if chrom_name is not None:
                        if not label in labels_in_bed:
                            this_chrom_versioned = re.sub("NC_(0)+", "chr", chrom_name)
                            this_chrom = re.split("\.", this_chrom_versioned)[0]
                            if this_chrom == "chr23":
                                this_chrom = "chrX"
                            elif this_chrom == "chr24":
                                this_chrom = "chrY"
                            bed_fp.write("%s\t%d\t%d\t%s\t0\t+\t%d\t%d\t0,0,0\n"
                                         % (this_chrom, genomic_start, genomic_end, label,
                                            genomic_start, genomic_end))
                            if coordinates_dir is not None:
                                coordinates_fp.write("%s\t%s\t%d\t%d\t%s\t%d\t%d\n" 
                                                     % (reference_type, reference_accession,
                                                        reference_start, reference_end,
                                                        this_chrom, genomic_start, genomic_end))
                            labels_in_bed[label] = 1
    return(bed_filename)

        

def gather_unique_alts(mave_json):
    """
    Generate a hash in which each key is an observed alternative allele
    and the value is a numeric index
    """
    unique_alts = {}
    unique_alts_observed = 0
    if "mapped_scores" in mave_json:
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
    if "mapped_scores" in mave_json:
        for mapped_score in mave_json["mapped_scores"]:
            this_score = mapped_score["score"]
            if this_score is not None:
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
    location_matrix_filename = "%s/%s.lm" % (location_matrix_dir, score_set_name)
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
    return(min_score, max_score, location_matrix_filename)
            
    
def mave_to_track(mave_json, track_name, trackdb_fp, coordinates_dir,
                  bed_dir, location_matrix_dir, bigBed_dir, chrom_sizes,
                  debug):
    # Extract the score set name from the URN.  This is the last token of
    # a colon-delimited URN.  This will be the name for instances where
    # reserved characters don't work well, such as filenames
    score_set_name = re.split(":", mave_json["urn"])[-1]
    gene_name = re.split("\s", mave_json["targetGene"]["name"])[0]
    shortLabel = gene_name + "_" + score_set_name
    longLabel = mave_json["title"]
    bed_file = create_bed(mave_json, score_set_name, 
                          bed_dir, coordinates_dir, debug)
    (min_score,max_score,
     location_matrix_file) = create_location_matrix(mave_json,
                                                    score_set_name,
                                                    location_matrix_dir)
    my_bigBed_dir = "%s/%s" % (bigBed_dir, score_set_name)
    if not os.path.exists(my_bigBed_dir):
        os.makedirs(my_bigBed_dir)
    runBigHeat(bed_file, location_matrix_file, chrom_sizes, my_bigBed_dir,
               shortLabel=shortLabel, longLabel=longLabel)
    return(my_bigBed_dir)
    
def filename_from_pathname(pathname):
    return(re.split("[/]", pathname)[-1])

def main():
    args = parse_args()
    if not os.path.exists(args.bed_dir):
        os.makedirs(args.bed_dir)
    if not os.path.exists(args.location_matrix_dir):
        os.makedirs(args.location_matrix_dir)
    if args.coordinates_dir is not None:
        if not os.path.exists(args.coordinates_dir):
            os.makedirs(args.coordinates_dir)
    if not os.path.exists(args.hub_dir):
        os.makedirs(args.hub_dir)
    shutil.copy(args.hub_file, args.hub_dir)
    shutil.copy(args.genomes_file, args.hub_dir)
    shutil.copy(args.html, args.hub_dir)
    append_to_file(args.hub_dir + "/" + filename_from_pathname(args.hub_file),
                   "description", filename_from_pathname(args.html))
    genomes_dir = args.hub_dir + "/" + args.genome_name
    if not os.path.exists(genomes_dir):
        os.makedirs(genomes_dir)
    shutil.copy(args.html, genomes_dir)
    bigBed_dir = genomes_dir +"/" + "bigBed"
    if not os.path.exists(bigBed_dir):
        os.makedirs(bigBed_dir)
    trackdb_filename = genomes_dir + "/trackDb.txt"
    with open(trackdb_filename, "w") as trackdb_fp:
        #
        # When writing the trackDb file, extract the filename from the pathname, reflecting
        # that the file has been copied to the directory with 
        write_trackDb(args.track_name, trackdb_fp, re.split("[/]", args.html)[-1])
        for this_input_file in args.input_files:
            if args.debug:
                print("Working on input file", this_input_file)
            mave_json = json.load(open(this_input_file))
            (this_bigBed_dir) = mave_to_track(mave_json,
                                              args.track_name,
                                              trackdb_fp,
                                              args.coordinates_dir,
                                              args.bed_dir,
                                              args.location_matrix_dir,
                                              bigBed_dir,
                                              args.chrom_sizes,
                                              args.debug)
            #
            # If bigBed files were generated successfully, then include this
            # experiment's trackDb file into the larger trackDb.  Else, remove
            # this bigBed directory.
            if (len(glob(this_bigBed_dir + "/*bb")) > 0):
                relative_bigBed_dir = re.sub(genomes_dir, "",
                                              this_bigBed_dir)
                trackdb_fp.write("include .%s/trackDb.ra\n"
                                 % (relative_bigBed_dir))
                if args.debug:
                    print("including score set", this_bigBed_dir)
            else:
                os.rmdir(this_bigBed_dir)
                if args.debug:
                    print("omitting score_set", this_bigBed_dir)

if __name__ == "__main__":
    main()
    
