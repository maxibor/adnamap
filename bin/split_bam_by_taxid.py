#!/usr/bin/env python
import argparse
from pathlib import Path
from multiprocessing.pool import ThreadPool
from functools import partial
import pysam


def read_bam(bamfile):
    """
    Read a BAM file and return a list of reads.
    """
    bam = pysam.AlignmentFile(bamfile, "rb")
    taxon_levels = ["species", "strain"]
    taxons = dict()
    refids = dict()
    for read in bam.fetch():
        rank = read.get_tag("XR")
        if rank in taxon_levels:
            taxon = read.get_tag("XT")
            if taxon not in taxons:
                taxons[taxon] = {'read': [read], 'name': read.get_tag('XN'), 'ref': {read.reference_name: read.reference_length} }
            else:
                if read.reference_name not in taxons[taxon]['ref']:
                    taxons[taxon]['ref'][read.reference_name] = read.reference_length
                taxons[taxon]['read'].append(read)
        else:
            continue
    return taxons


def write_bam(taxon_dict_item, basename):
    taxid, taxon_dict_entry = taxon_dict_item
    header = {
        "HD": {"VN": "1.0"},
        "SQ": [
            {"LN": taxon_dict_entry['ref'][refname], "SN": refname} for refname in taxon_dict_entry['ref']
        ]
    }
    refid_dict = dict()
    for i, ref in enumerate(taxon_dict_entry['ref'].keys()):
        refid_dict[ref] = i
    with pysam.AlignmentFile(f"{basename}_{taxid}.bam", "wb", header=header) as out:
        for read in taxon_dict_entry['read']:
            read.reference_id = refid_dict[read.reference_name]
            out.write(read)

def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Split BAM file created by SAM2LCA by species TAXID",
        epilog="Example: python split_bam_by_taxid.py myfile.sam2lca.bam",
    )
    parser.add_argument(
        "bam",
        metavar="bam",
        type=Path,
        help="sam2lca bam output file",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        help="number of threads",
    )


    return parser.parse_args(argv)

def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    taxons = read_bam(args.bam)
    write_bam_partial = partial(write_bam, basename=Path(args.bam).stem)
    with ThreadPool(args.threads) as p:
        p.map(write_bam_partial, taxons.items())

if __name__ == "__main__":
    main()
