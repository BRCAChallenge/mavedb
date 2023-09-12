#!/usr/bin/env python

import argparse
import json

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_json", help="Input file(s), JSON format")
    args = parser.parse_args()
    mave_json = json.load(open(args.input_json))
    for mapped_score in mave_json["mapped_scores"]:
        score = mapped_score["score"]
        ref = mapped_score["post_mapped"]["vrs_ref_allele_seq"]
        variation = mapped_score["post_mapped"]["variation"]
        start = variation["location"]["interval"]["start"]["value"]
        end = variation["location"]["interval"]["end"]["value"]
        alt = variation["state"]["sequence"]
        print("%d\t%d\t%s\t%s\t%s" % (start, end, ref, alt, score))

if __name__ == "__main__":
    main()
