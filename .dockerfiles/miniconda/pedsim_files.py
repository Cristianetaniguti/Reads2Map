import csv
import argparse
import pandas as pd
from pyfasta import Fasta


def _process_indels_coordinates(record, fasta):
    chrom = record['seq_id']
    if record["type"] == "-":
        fix = 1
        start, end = record['orig_position'] - fix, record['orig_position'] - fix + record['indel_length']
        reference = fasta.sequence({'chr': chrom, 'start': start, 'stop': end})
        alternative = list(reference)[0]
        position = record['orig_position'] - fix
        return [record['seq_id'], position, reference, alternative]
    
    if record["type"] == "+":
        start, end = record['orig_position'], record['orig_position']
        reference = fasta.sequence({'chr': chrom, 'start': start, 'stop': end})
        alternative = reference + record['indel_sequence']
        position = record['orig_position']
        return [record['seq_id'], position, reference, alternative]

    raise ValueError("Unknown value type: %s" % record["type"])

def produce_vcf_like_data(snps, indels, reference):
    f = Fasta(reference)
    ref_bases = indels.apply(lambda row: _process_indels_coordinates(row, f), axis=1)
    final_indels = pd.DataFrame(list(ref_bases), columns=["chr", "pos", "ref", "alt"])
    final_snps = snps[["seq_id", "orig_position", "orig_base", "new_base"]]
    final_snps.columns = ["chr", "pos", "ref", "alt"]
    return pd.concat([final_indels, final_snps], ignore_index=True).sort_values(by=['chr', 'pos'])

def produce_markers_list(variants):
    markers_names = ["M%03d" % i for i in range(1, variants.shape[0] + 1)]
    markers_positions = list(variants.apply(lambda row: row['pos']/1000000 * 4.63, axis=1))
    markers = pd.DataFrame({'marker':markers_names, 'chromosome': variants["chr"].tolist(), 'position': markers_positions})
    return markers

def produce_founders_individuals(markers, variants):
    return pd.DataFrame({
        'marker': markers['marker'].tolist(),
        'P1_1': variants['ref'].tolist(),
        'P1_2': variants['ref'].tolist(),
        'P2_1': variants['alt'].tolist(),
        'P2_2': variants['alt'].tolist(),
        })


def produce_parameters(seed):
    return dict(
        PLOIDY=2,
        MAPFUNCTION='HALDANE',
        MISSING = 'NA',
        CHROMFILE = "inb.chrom",
        POPTYPE = 'F2',
        SEED = seed,
        POPSIZE = 150,
        MAPFILE = "mapfile.map",
        FOUNDERFILE = "founderfile.gen",
        OUTPUT = "sim_inb",
    )

def chromosome_structure(markers):
    chromossomes = markers['chromosome'].tolist()
    markers = markers['position'].tolist()
    return pd.DataFrame({
        'chromosome': [chromossomes[0]],
        'length': [max(markers)],
        'centromere': [max(markers)/2],
        'prefPairing': [0],
        'quadrivalents': [0]
    })


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--indels', required=True, help='Pirs list of indel variants')
    parser.add_argument('--snps', required=True, help='Pirs list of snp variants')
    parser.add_argument('--reference', required=True, help='Fasta sequence of the reference chromosome')
    parser.add_argument('--seed', required=True, help='Seed to be applied on random steps')
    args = parser.parse_args()

    # outputs
    founder_file = "founders.txt"
    map_file = "markers.txt"
    chromosome_file = "chromosome.txt"
    parameters_file = "parameters.txt"
    variants_file = "tot_mks.txt"

    indels = pd.read_csv(args.indels, sep="\t", skiprows=3, 
                         names=['seq_id', 'orig_position', 'new_position', 'type', 'indel_length', 'indel_sequence'], 
                         dtype= {"orig_position": int, "indel_length": int})
    snps = pd.read_csv(args.snps, sep="\t", skiprows=3,  
                       names=['seq_id', 'orig_position', 'new_position', 'orig_base', 'new_base'],
                       dtype= {"orig_position": int, "new_position": int})
    
    variants = produce_vcf_like_data(snps, indels, args.reference)
    if len(list(set(variants['chr'].tolist()))) > 1:
        raise ValueError("We only simulate one chromosome each time")
    variants.reset_index()[["chr", "pos", "ref", "alt"]].to_csv(variants_file, sep="\t", quoting=csv.QUOTE_NONNUMERIC, index=False, header=False)

    markers = produce_markers_list(variants)
    markers[["marker", "chromosome", "position"]].to_csv(map_file, sep="\t", index=False)
    
    founders = produce_founders_individuals(markers, variants)
    founders[["marker", "P1_1", "P1_2", "P2_1", "P2_2"]].to_csv(founder_file, sep="\t", index=False)

    chromosome = chromosome_structure(markers)
    chromosome[["chromosome", "length", "centromere", "prefPairing", "quadrivalents"]].to_csv(chromosome_file, sep="\t", index=False)

    parameters = produce_parameters(args.seed)
    with open(parameters_file, "w+") as out:
        for key, value in parameters.items():
            out.write("%s = %s\n" % (key, value))
