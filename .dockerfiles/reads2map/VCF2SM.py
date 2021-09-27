#!/usr/bin/python

import sys
import os
import argparse
import subprocess
import shlex
from os.path import expanduser
import getopt
from HTMLParser import HTMLParser
import multiprocessing as mp
import re
from collections import namedtuple
import ast
import datetime


# Search for SuperMASSA dependencies (exit early if any missing)
try:
    from copy import deepcopy
    from pprint import pprint
    import cProfile
    import itertools
    import math
    import matplotlib
    import numpy
    import pylab as P
    import random
    import tempfile
except ImportError, e:
    sys.stderr.write('Error caught before starting the script.\n')
    sys.stderr.write('SuperMASSA dependencies missing.\n')
    sys.stderr.write('\t' + str(e) + '\n')
    sys.exit()


# Variables to store filtering results
SM_fail = -9
var_filter_posterior = -1
var_filter_ploidy = -2
var_filter_call_rate = -3
filters = {"posterior":0,
           "ploidy":0,
           "call_rate":0,
           "two_alleles":0,
           "min_depth":0,
           "max_depth":0}


class MyHTMLParser(HTMLParser):
    def __init__(self, geno_count, inference):
        HTMLParser.__init__(self)
        self.in_header_result = False
        self.in_genotype_result = False
        self.inference = inference
        self.allele_doses = ['.'] * geno_count
        self.posterior_pattern = re.compile('Pr\(G = g\*, parameter\* \| D\) = (?P<posterior>[0-9].[0-9]*)')
        self.ploidy_pattern = re.compile('\'ploidy\': (?P<ploidy>[0-9]*)')
        self.parents_pattern = re.compile('\'parameter\': (?P<parents>\(\([0-9]*, [0-9]*\), \([0-9]*, [0-9]*\)\))')

    @property
    def genotypes(self):
        return self.allele_doses

    @property
    def posterior(self):
        return self.posterior

    @property
    def ploidy(self):
        return self.ploidy

    @property
    def parents(self):
        return self.parents

    def handle_data(self, data):
        if data == "Results":
            self.in_header_result = True
        if data == "Predicted genotypes":
            self.in_header_result = False
            self.in_genotype_result = True

        if self.in_header_result:
            post_match = self.posterior_pattern.search(data)
            if post_match:
                self.posterior = float(post_match.group('posterior'))
            ploidy_match = self.ploidy_pattern.search(data)
            if ploidy_match:
                self.ploidy = int(ploidy_match.group('ploidy'))
            if self.inference == "f1":
                parents_match = self.parents_pattern.search(data)
                if parents_match:
                    self.parents = parents_match.group('parents')

        if self.in_genotype_result:
            data_split = data.strip().split("\t")
            if data_split[0] not in ["Predicted genotypes", ""]:
                cur_idx = int(data_split[0])
                self.allele_doses[cur_idx] = eval(data_split[1])


def getrange(range_limits):
    my_split = range_limits.split(':')
    try:
        if len(my_split) == 1:
            res = int(range_limits) + 8
        elif len(my_split) == 2:
            res = range(int(my_split[0]) + 8, int(my_split[1]) + 9) # the hardcoded 8/9 values remove 8 mandatory and the format columns from a vcf file
        else:
            raise ValueError("More than 2 values to expand range from.")
    except ValueError:
        print "Could not expand range limits from '%s'" % range_limits
        raise
    return res


def iter_flatten(iterable, toint = False):
    it = iter(iterable)
    for e in it:
        if isinstance(e, (list, tuple)):
            for f in iter_flatten(e, toint):
                yield f
        else:
            if toint:
                yield int(e)
            else:
                yield e


def vcf_create_record(in_line, allele_depth_arg, use_parents, variant_ploidy, variant_genotypes, variant_parents, idx, par1Idx, par2Idx):
    vals = in_line.split()

    # Get allele depths (do this again here, to avoid passing depths to multiple functions)
    temp_fmt = vals[8].split(':')
    if len(allele_depth_arg) == 1:
        depth_idx = temp_fmt.index(allele_depth_arg[0])
        allele_dp = [get_AD_idx_single(vals[i], depth_idx) for i in idx]
        if use_parents:
            allele_dp_par1 = [get_AD_idx_single(vals[i], depth_idx) for i in par1Idx]
            allele_dp_par2 = [get_AD_idx_single(vals[i], depth_idx) for i in par2Idx]
    else:
        depth_idx = (temp_fmt.index(allele_depth_arg[0]),
                     temp_fmt.index(allele_depth_arg[1]))
        allele_dp = [get_AD_idx_double(vals[i], depth_idx) for i in idx]
        if use_parents:
            allele_dp_par1 = [get_AD_idx_double(vals[i], depth_idx) for i in par1Idx]
            allele_dp_par2 = [get_AD_idx_double(vals[i], depth_idx) for i in par2Idx]

    combined_allele_dp_par1 = [0,0]
    combined_allele_dp_par2 = [0,0]
    if use_parents:
        nPar1 = par1Idx.__len__()
        nPar2 = par2Idx.__len__()
        combined_allele_dp_par1 = [sum([int(allele_dp_par1[i][j]) for i in xrange(0, nPar1)]) for j in xrange(0,2)]
        combined_allele_dp_par2 = [sum([int(allele_dp_par2[i][j]) for i in xrange(0, nPar2)]) for j in xrange(0,2)]

    depth_samples = [i for i in iter_flatten(allele_dp, True)]
    par1_depth = sum(combined_allele_dp_par1)
    par2_depth = sum(combined_allele_dp_par2)
    depth_total = sum(depth_samples) + par1_depth + par2_depth
    total_calls = sum(i != '.' for i in variant_genotypes)
    if variant_parents:
        total_calls += 2

    # Begin modifying input line to include SuperMASSA results
    variant_line = vals[0]

    # No quality score is assigned
    vals[5] = "."

    # For F1 progenies, if there is no depth information for both parents, SuperMASSA arbitrarily assigns genotypes
    # We mark this here for downstream fiters
    if use_parents and par1_depth == 0 and par2_depth == 0:
        vals[6] = "PAR_GT"
    else:
        # All other output variants already passed filtering criteria
        vals[6] = "PASS"

    # INFO fields: number of samples with data, total depth and ploidy estimate
    vals[7] = ("NS=" + str(total_calls) +
               ";DP=" + str(depth_total) +
               ";Pl=" + str(variant_ploidy))

    # FORMAT field: genotype and allele depths
    vals[8] = "GT:AD:DP"

    for i in vals[1:9]:
        variant_line += "\t" + i

    if variant_parents:
        variant_parents = ast.literal_eval(variant_parents)
        par1_call, par2_call = variant_parents
        # CAUTION!!
        # If there is allele depth data for the second parent, but not for the first, SuperMASSA swaps their genotypes!
        # We correct this here
        if par1_depth == 0 and par2_depth > 0:
            par1_call, par2_call = par2_call, par1_call

        par1_geno_temp = ['0'] * int(par1_call[0]) + ['1'] * int(par1_call[1])
        par1_geno = '/'.join(par1_geno_temp)
        par1_allele_depths = ','.join(map(str, combined_allele_dp_par1))
        variant_line += "\t" + ':'.join([par1_geno, par1_allele_depths, str(par1_depth)])

        par2_geno_temp = ['0'] * int(par2_call[0]) + ['1'] * int(par2_call[1])
        par2_geno = '/'.join(par2_geno_temp)
        par2_allele_depths = ','.join(map(str, combined_allele_dp_par2))
        variant_line += "\t" + ':'.join([par2_geno, par2_allele_depths, str(par2_depth)])

    for i in xrange(len(variant_genotypes)):
        if variant_genotypes[i] == '.':
            new_geno = '/'.join(['.'] * variant_ploidy)
        else:
            geno_temp = ['0'] * int(variant_genotypes[i][0]) + ['1'] * int(variant_genotypes[i][1])
            new_geno = '/'.join(geno_temp)
        ind_allele_depths = ','.join(allele_dp[i])
        ind_depth = int(allele_dp[i][0]) + int(allele_dp[i][1])
        variant_line += "\t" + ':'.join([new_geno, ind_allele_depths, str(ind_depth)])
    variant_line += "\n"
    return variant_line


def run_SM(in_list):
    SuperMASSA_cmd, args, use_parents, in_line, geno_count, idx, par1Idx, par2Idx, cur_file = in_list
    cmd_args = [expanduser(e) for e in shlex.split(SuperMASSA_cmd)]

    try:
        p = subprocess.check_output(cmd_args, stderr=subprocess.STDOUT, shell=False)
    except subprocess.CalledProcessError:
        # SuperMASSA did not succeed
        os.remove(cur_file)
        if use_parents:
            cur_par_file = os.path.splitext(cur_file)[0] + "_PARENTS.txt"
            os.remove(cur_par_file)
        # For now, simply ignore this variant site
        return SM_fail
    else:
        os.remove(cur_file)
        if use_parents:
            cur_par_file = os.path.splitext(cur_file)[0] + "_PARENTS.txt"
            os.remove(cur_par_file)

        if p[0:5] == "Error":
            raise RuntimeError("Error when running SuperMASSA. Please check command line arguments.")

        if p[0:16] == "<h1>Results</h1>":
            parser = MyHTMLParser(geno_count, args.inference)
            parser.feed(p)
            parser.close()
            variant_genotypes = parser.genotypes
            variant_posterior = parser.posterior
            variant_ploidy = parser.ploidy
            if args.inference == "f1":
                variant_parents = parser.parents
            else:
                variant_parents = None

            called_genos = sum(i != '.' for i in variant_genotypes)
            frac_called_genos = float(called_genos) / len(variant_genotypes)

            if variant_posterior < args.post:
                return var_filter_posterior
            elif variant_ploidy not in args.ploidy_filter:
                return var_filter_ploidy
            elif frac_called_genos < args.callrate:
                return var_filter_call_rate
            else:
                res = vcf_create_record(in_line, args.allele_depth, use_parents, variant_ploidy, variant_genotypes, variant_parents, idx, par1Idx, par2Idx)
                return res
        else:
            return SM_fail


def check_parents(args, parser):
    use_parents = False
    par1_avail = args.par1_pattern or args.par1_range
    par2_avail = args.par2_pattern or args.par2_range
    if args.inference == "f1":
        if par1_avail and par2_avail:
            use_parents = True
        elif par1_avail or par2_avail:
            parser.error("Please set both parents 1 and 2 or leave both empty")
    return use_parents


def check_AD_field(args, parser):
    temp_split = args.allele_depth.split("/")
    if len(temp_split) == 1 or len(temp_split) == 2:
        args.allele_depth = temp_split
    else:
        parser.error("Please inform only one or two allele depth field(s) (separated by '/')")


def check_SM_args(args, parser):
    if not args.SMscript or not args.inference:
        parser.error("Please provide the location of the SuperMASSA script and the type of inference")
    #######################################
    # CHECK FOR CORRECT SUPERMASSA SCRIPT #
    #######################################

    if args.inference == 'f1':
        ploidy_range = args.ploidy_range.split(':')
        if len(ploidy_range) == 3 and ploidy_range[1] != ploidy_range[0]:
            ploidy_range = [int(e) for e in ploidy_range]
            ploidies = range(ploidy_range[0], ploidy_range[1] + 1, ploidy_range[2])
            if any([e % 2 == 1 for e in ploidies]):
                parser.error("SuperMASSA cannot use odd ploidies for f1 inference")

    if args.ploidy_filter:
        allowed_ploidy_range = args.ploidy_filter.split(':')
    else:
        # if no ploidy filter is set by the user, accept all in tested range
        allowed_ploidy_range = args.ploidy_range.split(':')
    allowed_ploidy_range = [int(e) for e in allowed_ploidy_range]
    if len(allowed_ploidy_range) == 1:
        args.ploidy_filter = allowed_ploidy_range
    elif len(allowed_ploidy_range) == 2:
        args.ploidy_filter = range(allowed_ploidy_range[0], allowed_ploidy_range[1] + 1, 1)
    elif len(allowed_ploidy_range) == 3:
        args.ploidy_filter = range(allowed_ploidy_range[0], allowed_ploidy_range[1] + 1, allowed_ploidy_range[2])
    else:
        parser.error("Cannot expand ploidy range filter.")

    args.SMscript = expanduser(args.SMscript)


def expand_file_names(args, parser):
    in_file_list = []
    out_file_list = []

    if '+' in args.input or '+' in args.output:
        if '+' not in args.input or '+' not in args.output:
            parser.error("If multiple files are supplied, both input and output file paths must contain '+'")
        if args.sF and args.eF:
            in_file_list = [args.input.replace('+', str(f)) for f in xrange(args.sF[0], args.eF[0] + 1)]
            out_file_list = [args.output.replace('+', str(f)) for f in xrange(args.sF[0], args.eF[0] + 1)]
        else:
            parser.error("Please set sF and eF if file paths contain '+'")
    else:
        in_file_list = [args.input]
        out_file_list = [args.output]

    File_names = namedtuple('IOFiles', ['in_files', 'out_files'])
    return File_names(in_file_list, out_file_list)


def create_dirs(out_file_list):
    base_dir = os.path.dirname(os.path.abspath(out_file_list[0]))
    if not os.path.exists(base_dir):
        os.makedirs(base_dir)
    SM_input_dir = base_dir + '/InputFiles/'
    if not os.path.exists(SM_input_dir):
        os.makedirs(SM_input_dir)
    return SM_input_dir


def get_par_indices(args, parser, header, which_parent):
    if which_parent == 1:
        par_pattern = args.par1_pattern
        par_range = args.par1_range
    elif which_parent == 2:
        par_pattern = args.par2_pattern
        par_range = args.par2_range
    else:
        parser.error("Unknown parent")

    par_idx = []
    if par_pattern:
        par_idx = [i for i, v in enumerate(header) if any(p in header[i] for p in par_pattern) and i > 8] # the hardcoded 8 value removes 8 mandatory and the format columns from a vcf file
    elif par_range:
        par_idx = [getrange(i) for i in par_range]
        par_idx = [i for i in iter_flatten(par_idx, False)]
    if par_idx == []:
        parser.error("Could not find any samples matching pattern for parent " + str(which_parent))
    return par_idx


def get_sample_indices(args, parser, header, par1Idx, par2Idx):
    idx = []
    if args.geno_pattern:
        idx = [i for i, v in enumerate(header) if any(p in header[i] for p in args.geno_pattern) and i > 8] # the hardcoded 8 value removes 8 mandatory and the format columns from a vcf file
    elif args.geno_range:
        idx = [getrange(i) for i in args.geno_range]
        idx = [i for i in iter_flatten(idx, False)]
    else:
        # if no pattern or range is provided, use all samples in the VCF file
        idx = [i for i, v in enumerate(header) if i > 8] # the hardcoded 8 value removes 8 mandatory and the format columns from a vcf file
    # remove parent indices (if any) from sample indices
    idx = [i for i in idx if i not in par1Idx and i not in par2Idx]
    if idx == []:
        parser.error("Could not find any samples matching pattern")
    return idx


def get_AD_idx_single(vals, depth_idx):
    temp_split = vals.split(':')
    if temp_split != ['.'] and temp_split[depth_idx] != '.':
        return [e.replace('.', '0') for e in temp_split[depth_idx].split(',')]
    else:
        return ['0','0']


def get_AD_idx_double(vals, depth_idx):
    temp_split = vals.split(':')
    if temp_split != ['.']:
        return [temp_split[depth_idx[0]].replace('.', '0'), temp_split[depth_idx[1]].replace('.', '0')]
    else:
        return ['0','0']


def process_variant(args, parser, n_vars, l, f, idx, SM_input_dir, use_parents, par1Idx, par2Idx, jobs):
    vals = l.split()
    alt_count = len(vals[4].split(','))

    global filters

    if alt_count == 1 and vals[4] != '.' and vals[3] != '.':
        # We currently ignore variants with more than 2 alleles
        temp_fmt = vals[8].split(':')
        if len(args.allele_depth) == 1:
            try:
                depth_idx = temp_fmt.index(args.allele_depth[0])
            except ValueError:
                parser.error(str(args.allele_depth) + " field not present in format")
            allele_depths = [get_AD_idx_single(vals[i], depth_idx) for i in idx]
            if use_parents:
                allele_depths_par1 = [get_AD_idx_single(vals[i], depth_idx) for i in par1Idx]
                allele_depths_par2 = [get_AD_idx_single(vals[i], depth_idx) for i in par2Idx]
        else:
            try:
                depth_idx = (temp_fmt.index(args.allele_depth[0]),
                             temp_fmt.index(args.allele_depth[1]))
            except ValueError:
                parser.error(str(args.allele_depth[0]) + "/" + str(args.allele_depth[1]) + " fields not present in format")
            allele_depths = [get_AD_idx_double(vals[i], depth_idx) for i in idx]
            if use_parents:
                allele_depths_par1 = [get_AD_idx_double(vals[i], depth_idx) for i in par1Idx]
                allele_depths_par2 = [get_AD_idx_double(vals[i], depth_idx) for i in par2Idx]

        # Check for average allele depth
        depth_samples = sum(iter_flatten(allele_depths, True))
        avg_depth = depth_samples/(len(idx) * 1.0)

        if avg_depth < args.minimum_depth:
            filters["min_depth"] += 1
        elif avg_depth > args.maximum_depth:
            filters["max_depth"] += 1
        else:
            if vals[2] != '.':
                marker_name = vals[2]
            else:
                marker_name = vals[0] + '_' + vals[1]
            tempfilename = SM_input_dir + os.path.basename(f) + '_var' + str(n_vars) + '_' + marker_name + '.txt'
            tempFile = file(tempfilename, 'w')
            geno_count = idx.__len__()
            [tempFile.write(str(i) + '\t' + allele_depths[i][0] + '\t' + allele_depths[i][1] + '\n') for i in xrange(geno_count)]
            tempFile.close()

            combined_allele_depths_par1 = [0,0]
            combined_allele_depths_par2 = [0,0]
            if use_parents:
                tempParfilename = SM_input_dir + os.path.basename(f) + '_var' + str(n_vars) + '_' + marker_name + '_PARENTS.txt'
                tempParFile = file(tempParfilename, 'w')
                nPar1 = par1Idx.__len__()
                nPar2 = par2Idx.__len__()
                [tempParFile.write('PAR1\t' + allele_depths_par1[i][0] + '\t' + allele_depths_par1[i][1] + '\n') for i in xrange(0, nPar1)]
                [tempParFile.write('PAR2\t' + allele_depths_par2[i][0] + '\t' + allele_depths_par2[i][1] + '\n') for i in xrange(0, nPar2)]
                tempParFile.close()
                combined_allele_depths_par1 = [sum([int(allele_depths_par1[i][j]) for i in xrange(0, nPar1)]) for j in xrange(0,2)]
                combined_allele_depths_par2 = [sum([int(allele_depths_par2[i][j]) for i in xrange(0, nPar2)]) for j in xrange(0,2)]

            ############################
            # CHECK FOR PYTHON VERSION #
            ############################
            SuperMASSA_cmd = "python " + args.SMscript
            SuperMASSA_cmd += " --file=" + tempfilename
            SuperMASSA_cmd += " --inference=" + args.inference
            SuperMASSA_cmd += " --ploidy_range=" + args.ploidy_range
            SuperMASSA_cmd += " --naive_posterior_reporting_threshold=" + str(args.naive_reporting)
            SuperMASSA_cmd += " --print_genotypes"
            if args.exact:
                SuperMASSA_cmd += " --optimal_inference"
            if args.sigma_range:
                SuperMASSA_cmd += " --sigma_range=" + args.sigma_range
            if use_parents:
                SuperMASSA_cmd += " --f1_parent_data=" + tempParfilename

            job = (SuperMASSA_cmd, args, use_parents, l, geno_count, idx, par1Idx, par2Idx, tempfilename)
            jobs.append(job)
    else:
        filters["two_alleles"] += 1


def write_VCF_metadata(out_file, args, use_parents):
    now = datetime.datetime.now()
    out_file.write('##fileformat=VCFv4.2\n')
    out_file.write('##VCF2SMCommandLine=<ID=VCF2SM,Time="' + now.strftime("%Y-%m-%d %H:%M") + '",CommandLineOptions="' + str(args) + '">\n')
    out_file.write('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with a genotype call (including parents, for F1 progenies)">\n')
    out_file.write('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth (including reads from parents and uncalled genotypes)">\n')
    out_file.write('##INFO=<ID=Pl,Number=1,Type=Integer,Description="Estimated ploidy">\n')
    if use_parents:
        out_file.write('##FILTER=<ID=PAR_GT,Description="No allele depth data for parents (arbitrary genotypes)">\n')
    out_file.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
    out_file.write('##FORMAT=<ID=AD,Number=2,Type=Integer,Description="Allelic depths for the reference and alternate alleles">\n')
    out_file.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')


def process_input_VCF_files(args, parser, SM_input_dir, use_parents, in_file_list, out_file_list):
    global filters
    passed = 0
    failed = 0

    for k in xrange(len(in_file_list)):
        f = in_file_list[k]
        input_file_conn = file(f, 'r')

        out_filename = out_file_list[k]
        out_file = file(out_filename, 'w')

        write_VCF_metadata(out_file, args, use_parents)

        for l in input_file_conn:
            if l[0] == '#':
                if l[1] == '#':
                    # Do nothing: ignore original VCF metadata
                    continue
                elif l[1] != '#':
                    header = l.split()

                    par1Idx = []
                    par2Idx = []
                    if use_parents:
                        par1Idx = get_par_indices(args, parser, header, 1)
                        par2Idx = get_par_indices(args, parser, header, 2)

                    idx = get_sample_indices(args, parser, header, par1Idx, par2Idx)

                    names = [header[i] for i in idx]

                    # Write VCF header
                    out_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT")
                    if args.inference == "f1":
                        out_file.write("\tPARENT1\tPARENT2")
                    [out_file.write("\t" + i) for i in names]
                    out_file.write("\n")

                    break

        jobs = []
        n_vars = 0
        for l in input_file_conn:
            if l[0] != '#':
                n_vars += 1
                process_variant(args, parser, n_vars, l, f, idx, SM_input_dir, use_parents, par1Idx, par2Idx, jobs)

        pool = mp.Pool(processes = args.threads)
        results = pool.map(run_SM, jobs)

        for res in results:
            if res == SM_fail:
                failed += 1
            elif res == var_filter_posterior:
                filters["posterior"] += 1
            elif res == var_filter_ploidy:
                filters["ploidy"] += 1
            elif res == var_filter_call_rate:
                filters["call_rate"] += 1
            elif res != None:
                passed += 1
                out_file.write(str(res))
        out_file.close()

        input_file_conn.close()

    if os.listdir(SM_input_dir) == []:
        os.rmdir(SM_input_dir)

    total_filtered = sum(filters.values())
    print (' RESULTS '.center(80, '-'))
    print "SuperMASSA failed variants: " + str(failed)
    print ('-' * 80)
    print "Variants removed by filtering criteria: " + str(total_filtered)
    print "\tLow average depth: " + str(filters["min_depth"])
    print "\tHigh average depth: " + str(filters["max_depth"])
    print "\tLow call rate: " + str(filters["call_rate"])
    print "\tLow ploidy posterior probability: " + str(filters["posterior"])
    print "\tPloidy not in the specified range: " + str(filters["ploidy"])
    print "\tNon-biallelic variant: " + str(filters["two_alleles"])
    print ('-' * 80)
    print "Variants kept: " + str(passed)
    print ('-' * 80)


def main():
    parser = argparse.ArgumentParser()

    optional = parser._action_groups.pop()
    required = parser.add_argument_group('required named arguments')

    required.add_argument("-i", "--input", type=str, help="Path to VCF input file", required = True)
    required.add_argument("-o", "--output", type=str, help="Path to VCF output file", required = True)
    required.add_argument("-S", "--SMscript", type=str, help="Path to the SuperMASSA script", required = True)
    required.add_argument("-I", "--inference", type=str, help="Inference mode for SuperMASSA", required = True)

    group_samples = optional.add_mutually_exclusive_group()
    group_samples.add_argument("-g", "--geno_pattern", nargs="+", type=str, help="Pattern(s) of genotype IDs to include as samples")
    group_samples.add_argument("-r", "--geno_range", nargs="+", type=str, help="Indices of genotypes to include as samples (list or range in the form 1:150)")
    group_parent1 = optional.add_mutually_exclusive_group()
    group_parent1.add_argument("-1", "--par1_pattern", nargs="+", type=str, help="Pattern(s) of genotype IDs for parent 1")
    group_parent1.add_argument("-k", "--par1_range", nargs="+", type=str, help="Indices of genotypes to include as parent 1 (list or range in the form 1:5)")
    group_parent2 = optional.add_mutually_exclusive_group()
    group_parent2.add_argument("-2", "--par2_pattern", nargs="+", type=str, help="Pattern(s) of genotype IDs for parent 2")
    group_parent2.add_argument("-l", "--par2_range", nargs="+", type=str, help="Indices of genotypes to include as parent 2 (list or range in the form 1:5)")

    optional.add_argument("--sF", nargs=1, type=int, help="Starting file")
    optional.add_argument("--eF", nargs=1, type=int, help="Ending file")
    optional.add_argument("-a", "--allele_depth", type=str, default="AD", help="VCF field(s) to get allele depths from")
    optional.add_argument("-d", "--minimum_depth", type=float, default=0, help="Minimum average depth per sample for variant site to be processed (not including parents, for F1 progenies)")
    optional.add_argument("-D", "--maximum_depth", type=float, default=float("inf"), help="Maximum average depth per sample for variant site to be processed (not including parents, for F1 progenies)")
    optional.add_argument("-e", "--exact", action='store_true', help="Perform exact inference (default FALSE)")
    optional.add_argument("-M", "--ploidy_range", type=str, default = "2:16", help="Ploidy range to test with SuperMASSA")
    optional.add_argument("-V", "--sigma_range", type=str, help="Sigma range for SuperMASSA")
    optional.add_argument("-f", "--ploidy_filter", type=str, help="Range of ploidies to include in the output")
    optional.add_argument("-p", "--post", type=float, default = 0.0, help="Minimum posterior probability to keep variant")
    optional.add_argument("-n", "--naive_reporting", type=float, default = 0.0, help="Naive posterior reporting threshold")
    optional.add_argument("-c", "--callrate", type=float, default = 0.0, help="Minimum call rate to keep variant")
    optional.add_argument("-t", "--threads", type=int, default=1, help="Maximum number of threads to use")

    parser._action_groups.append(optional)
    args = parser.parse_args()

    use_parents = check_parents(args, parser)

    check_AD_field(args, parser)

    check_SM_args(args, parser)

    io_files = expand_file_names(args, parser)
    in_file_list = io_files.in_files
    out_file_list = io_files.out_files
    if not all(os.path.isfile(e) for e in in_file_list):
        parser.error("Input file doest not exist")

    SM_input_dir = create_dirs(out_file_list)

    process_input_VCF_files(args, parser, SM_input_dir, use_parents, in_file_list, out_file_list)


if __name__ == "__main__":
    main()
