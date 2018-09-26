import martian
import cPickle as pickle
import pandas as pd
from tools.transcripts import *
import numpy as np
import pysam
from tools.dataOps import *
from tools.misc import *
from collections import *

__MRO__ = '''
stage SPLIT_AND_EXTRACT(
    in   gp        transcripts,
    in   int       lower_size_cutoff,
    in   int       upper_size_cutoff,
    in   bam       possorted_bam,
    in   tsv       cell_barcodes,
    in   string[]  valid_chroms,
    in   string    kit_type,
    out  csv       results,
) split using (
    in  string[]   tx_subset,
    out pickle     subset_results,
)
'''

def split(args):
    # validate
    if args.kit_type not in ["5'", "3'"]:
        martian.exit("Kit type is not one of 5' or 3'.")
    # group by gene
    tx_dict = get_gene_pred_dict(args.transcripts)
    tx_by_name2 = defaultdict(list)
    valid_chroms = set(args.valid_chroms)
    for tx in tx_dict.itervalues():
        if tx.chromosome in valid_chroms:
            tx_by_name2[tx.name2].append(tx)

    singletons = [x[0] for x in tx_by_name2.itervalues() if len(x) == 1 and
                  args.lower_size_cutoff <= len(x[0]) <= args.upper_size_cutoff]

    avg = np.mean([len(x) for x in singletons])
    med = np.median([len(x) for x in singletons])
    martian.log_info('{} singleton genes under consideration (avg size = {} median size = {}'.format(len(singletons),
                                                                                                     avg, med))
    def tx_to_str(singletons):
        for x in singletons:
            yield x.get_gene_pred()
    chunks = [{'tx_subset': list(x), '__mem_gb': 4} for x in grouper(tx_to_str(singletons), 20)]
    return {'chunks': chunks, 'join': {'__mem_gb': 32}}


def main(args, outs):
    in_bam = pysam.Samfile(args.possorted_bam)
    bcs = {x.rstrip() for x in open(args.cell_barcodes)}
    txs = [GenePredTranscript(x) for x in args.tx_subset]
    results = [[tx.name, find_recs(tx, in_bam, bcs)] for tx in txs]
    outs.pickle = martian.make_path('subset_results.pickle')
    with open(outs.pickle, 'w') as outf:
         pickle.dump(results, outf)


def join(args, outs, chunk_defs, chunk_outs):
    outs.coerce_strings()
    results = flatten_list_of_lists([pickle.load(open(chunk.subset_results)) for chunk in chunk_outs])
    results = flatten_list_of_lists(results)
    tx_dict = get_gene_pred_dict(args.transcripts)
    results = parse_results(results, tx_dict, args.kit_type)
    outs.results = martian.make_path('results.csv')
    df = pd.DataFrame(results, columns=['tx_id', 'tx_position', 'read_molecule_fraction'])
    df.to_csv(outs.results)


def find_recs(tx, bam, bcs):
    rec_map = defaultdict(set)
    all_keys = set()
    for e in tx.exon_intervals:
        for rec in bam.fetch(tx.chromosome, e.start, e.stop):
            if rec.has_tag('CB') and rec.has_tag('UB') and rec.get_tag('CB') in bcs:
                key = (rec.get_tag('UB'), rec.get_tag('CB'))
                all_keys.add(key)
                for p in rec.get_reference_positions():
                    tx_p = tx.chromosome_coordinate_to_mrna(p)
                    if tx_p is None:  # may align off end
                        continue
                    rec_map[tx_p].add(key)
    tot_keys = len(all_keys)
    return [[x, 1.0 * len(y) / tot_keys] for x, y in sorted(rec_map.iteritems(), key=lambda x: x[0])]


def parse_results(results, tx_dict, kit_type):
    r = []
    for tx_id, data in pairwise(results):
        if len(data) == 0:
            continue
        xvals, yvals = zip(*data)
        tx = tx_dict[tx_id]
        if len(tx) > 10000:
            continue
        if kit_type == "3'":  # flip around 3' data
            xvals = [len(tx) - x for x in xvals[::-1]]
        for xval, yval in zip(xvals, yvals):
            r.append([tx_id, xval, yval])
    return r
