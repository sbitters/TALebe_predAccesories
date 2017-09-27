#!/usr/bin/env python3
#  -*- coding: utf-8 -*-
#
# 2017-05-10
# STB
#

import os
import sys
import regex as re
import pandas as pd
from collections import OrderedDict
from optparse import OptionParser

prog_version = "0.2"

MIT_license = """Copyright 2017 Sven T. Bitters (sven.bitters@gmail.com)

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
"""


def input_parse(MIT_license):
    # Parse inputs

    parser = OptionParser()
    parser.add_option("-l", "--license", action="store_true", dest="license_text", default=False,
                      help="show license information and exit")
    parser.add_option("-e", "--EBE", dest="EBE_f_loc", default="",
                      help="String. Mandatory. Path to EBE predictor's output file that will be reformatted in tablemaker.")
    parser.add_option("-p", "--predictor", dest="predictor_identity", default="",
                      help="String. Mandatory. Name of the EBE prediction tool that generated the data supplied in -e."
                           "\nOptions: TALVEZ, TALgetter, TALENT")
    parser.add_option("-c", "--column", dest="target_col", default="",
                      help="String. Mandatory. Exact name of the column in the prediction tool's output file that shall be used.")
    parser.add_option("-o", "--output", dest="output_loc", default="",
                      help="String. Optional. Complete path and filename for TALVEZ_tablemaker's output file.")
    parser.add_option("-s", "--strand", dest="target_strand", default="pos",
                      help="String. Mandatory. TALVEZ looks for EBEs on both strands of the input DNA sequence. Report only"
                           "EBEs on the positive or negative strand or on both? This option is only applicable if the"
                           "EBE predictor outputs not only predictions for one strand."
                           "\nOptions: pos, neg, both"
                           "\nDefault: pos")
    parser.add_option("-r", "--promotor", dest="prom_length", default=0,
                      help="Integer. Optional.")


    (parserargs, opts) = parser.parse_args()

    show_license = parserargs.license_text
    if show_license:
        print(MIT_license)
        sys.exit()

    if not parserargs.EBE_f_loc and not parserargs.target_col:
        sys.exit("\nError. Please specify a path to your EBE predictor's output file and a column therein.\n")

    if not parserargs.EBE_f_loc:
        sys.exit("\nError. Please specify a path to your EBE predictor's output file.\n")

    if not parserargs.target_col:
        sys.exit("\nError. Please specify a column of your EBE predictor's output file.\n")

    if not parserargs.predictor_identity:
        sys.exit("\nError. Please give the name of the EBE predictor that created the file you want to process.\n")


    selector_set = predictor_select(parserargs.predictor_identity)

    return parserargs, selector_set


def predictor_select(ebe_predictor):
    if ebe_predictor == "TALVEZ":
        pos_strand = 4
        pos_start = 6
        start_equalizer = 0
        def_strand = ["+strand", "-strand"]
        pos_seqid = 2
        pos_sequence = 8
    elif ebe_predictor == "TALgetter":
        pos_strand = 7
        pos_start = 2
        start_equalizer = 0
        def_strand = ["+", "-"]
        pos_seqid = 1
        pos_sequence = 4
    elif ebe_predictor == "TALENT":
        pos_strand = 7
        pos_start = 4
        start_equalizer = -2
        def_strand = ["+", "-"]
        pos_seqid = 1
        pos_sequence = 10
    else:
        sys.exit("\nError. Please select an EBE prediction tool.\n")

    selector_set = (ebe_predictor, pos_strand, pos_start, start_equalizer, def_strand, pos_seqid, pos_sequence)

    return selector_set


def read_EBEpred(input_file, selector_set):
    print("\nReading " + selector_set[0] + "'s output file...")

    TAL_ID_list = list()
    EBEpred_lines = list()
    skipfirstline = 0
    EBEpred_header = ""

    with open(input_file, "r") as input_handle:
        for line in input_handle:
            line = line.rstrip()

            splitline = line.split("\t")

            if skipfirstline == 0:
                EBEpred_header = splitline
                skipfirstline = 1
            else:
                EBEpred_lines.append(line)

                if selector_set[0] == "TALVEZ":
                    TAL_ID = re.sub(">", "", splitline[0])
                else:
                    TAL_ID = splitline[0]

                if TAL_ID not in TAL_ID_list:
                    TAL_ID_list.append(TAL_ID)

    return TAL_ID_list, EBEpred_header, EBEpred_lines


def create_empty_list(TAL_ID_list, empty_list_obj):
    empty_list = list()

    for null in TAL_ID_list:
        empty_list.append(empty_list_obj)

    return empty_list


def create_empty_dict(EBEpred_lines, strand, selector_set, TAL_ID_list, prom_length):
    empty_dict = OrderedDict()
    empty_EBE_dict = OrderedDict()

    zero_list = create_empty_list(TAL_ID_list, 0)
    empty_list = create_empty_list(TAL_ID_list, [])

    for line in EBEpred_lines:
        splitline = line.split("\t")

        if selector_set[0] == "TALVEZ":
            SEQ_ID = re.sub(">", "", splitline[2])
        else:
            SEQ_ID = splitline[selector_set[5]]

        gene = SEQ_ID.split(".r")[0]
        gene = gene.split(".g")[0]
        gene = gene.split(".c")[0]
        gene = gene.split(".i")[0]

        EBEstrand = splitline[selector_set[1]]
        EBE_start = int(splitline[selector_set[2]]) + selector_set[3] - prom_length

        if strand == EBEstrand:
            SEQ_ID_key = SEQ_ID + "." + str(EBE_start)
            empty_dict[SEQ_ID_key] = zero_list.copy()
            empty_EBE_dict[gene] = empty_list.copy()
        elif not strand:
            SEQ_ID_key = SEQ_ID + "." + str(EBE_start) + "_" + EBEstrand
            empty_dict[SEQ_ID_key] = zero_list.copy()
            empty_EBE_dict[gene] = empty_list.copy()

    return empty_dict, empty_EBE_dict


def fill_dict(EBEpred_lines, mydict, EBE_dict, header, selected_col, strand, selector_set, TAL_ID_list, prom_length):
    print("Handling data...")

    EBEseq_list = list()
    for null in TAL_ID_list:
        EBEseq_list.append([].copy())

    element_no = header.index(selected_col)

    for line in EBEpred_lines:
        splitline = line.split("\t")

        if selector_set[0] == "TALgetter":
            if splitline[selector_set[6]][0].upper() == "T":
                TAL_ID = splitline[0]
                SEQ_ID = splitline[1]
            else:
                continue

        elif selector_set[0] == "TALVEZ":
            if splitline[selector_set[6]][0].upper() == "T":
                TAL_ID = re.sub(">", "", splitline[0])
                SEQ_ID = re.sub(">", "", splitline[2])
            else:
                continue

        elif selector_set[0] == "TALENT":
            TAL_ID = splitline[0]
            SEQ_ID = splitline[1]

        myindex = TAL_ID_list.index(TAL_ID)

        gene = SEQ_ID.split(".r")[0]
        gene = gene.split(".g")[0]
        gene = gene.split(".c")[0]
        gene = gene.split(".i")[0]

        EBEstrand = splitline[selector_set[1]]
        EBE_start = int(splitline[selector_set[2]]) + selector_set[3] - prom_length
        EBEsequence = splitline[selector_set[6]]

        EBEpred_element = splitline[element_no]

        if selector_set[0] == "TALgetter":
            if selected_col == "Score":
                EBEpred_element = str(round(float(EBEpred_element), 3))

        if strand == EBEstrand:
            SEQ_ID_key = SEQ_ID + "." + str(EBE_start)
            mydict[SEQ_ID_key][myindex] = EBEpred_element

            EBE_dict[gene][myindex] = EBE_dict[gene][myindex] + [EBE_start]
            EBE_dict[gene][myindex] = sorted(list(set(EBE_dict[gene][myindex])))

            EBElist = EBEseq_list[myindex]
            EBElist = handle_EBEseq(EBElist, TAL_ID, SEQ_ID_key, EBEsequence, selector_set)
            EBEseq_list[myindex] = EBElist

        elif not strand:
            SEQ_ID_key = SEQ_ID + "." + str(EBE_start) + "." + EBEstrand
            mydict[SEQ_ID_key][myindex] = EBEpred_element

            EBE_dict[gene][myindex] = EBE_dict[gene][myindex] + [str(EBE_start) + "/" + EBEstrand]

    return mydict, EBE_dict, EBEseq_list


def handle_EBEseq(EBElist, TAL_ID, SEQ_ID_key, EBEsequence, selector_set):
    EBEseq_dict_key = SEQ_ID_key

    if selector_set[0] == "TALENT":
        EBEsequence = EBEsequence.split(";")
        EBEsequence = EBEsequence[1][-1] + EBEsequence[0]

    EBElist.append([EBEseq_dict_key, EBEsequence])

    return EBElist


def save_dict(mydict, TAL_ID_list, output_file):
    record = list()
    for keys, values in mydict.items():
        keys = [keys]
        item = keys + values

        record.append((item))

    df_header = ["Gene_Symbol"] + TAL_ID_list
    df = pd.DataFrame.from_records(record, columns=df_header)

    df.to_csv(output_file, sep="\t", index=False)


def save_seq(TAL_ID_list, EBEseq_list, output_seq, selector_set):
    for ii in range(0, len(TAL_ID_list)):
        TALE = TAL_ID_list[ii]

        fasta_output = output_seq + "_" + TALE + ".fa"
        with open(fasta_output, "w") as fasta_handle:
            for sublist in EBEseq_list[ii]:

                fasta_header = sublist[0]
                sequence = sublist[1]

                fasta_handle.write(">" + selector_set[0] + "_" + TALE + "_" + fasta_header + "\n")
                fasta_handle.write(sequence + "\n")


parser_arguments, selector_set = input_parse(MIT_license)
input_file = parser_arguments.EBE_f_loc
EBEpred_column = parser_arguments.target_col
prom_length = int(parser_arguments.prom_length)

output_file = parser_arguments.output_loc
if not output_file:
    output_file = input_file + "_tablemaker_" + EBEpred_column

output_dir = os.path.dirname(output_file)
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)

output_seq_dir = output_dir + "/" + "EBE_sequences" + "/"
output_seq = output_seq_dir + os.path.basename(output_file) + "_EBEs"
if not os.path.isdir(output_seq_dir):
    os.makedirs(output_seq_dir)

target_strand = parser_arguments.target_strand
if target_strand == "pos":
    target_strand = selector_set[4][0]
elif target_strand == "neg":
    target_strand = selector_set[4][1]
else:
    target_strand = ""

TAL_ID_list, EBEpred_header, EBEpred_lines = read_EBEpred(input_file, selector_set)

mydict, EBE_dict = create_empty_dict(EBEpred_lines, target_strand, selector_set, TAL_ID_list, prom_length)

mydict, EBE_dict, EBEseq_list = fill_dict(EBEpred_lines, mydict, EBE_dict, EBEpred_header, EBEpred_column, target_strand, selector_set, TAL_ID_list, prom_length)

print("Saving results to " + output_file)
save_dict(mydict, TAL_ID_list, output_file)

output_tmp = output_file + "_EBEs_tmp"
save_dict(EBE_dict, TAL_ID_list, output_tmp)

output_multi = output_file + "_EBEs"
print("Saving results to " + output_multi)
empty_regex = re.compile(r"\[\]")
bracket_regex = re.compile(r"\[|'|\]")
with open(output_tmp, "r") as tmp_handle:
    with open(output_multi, "w") as multi_handle:
        for line in tmp_handle:
            line = re.sub(empty_regex, "n.a.", line)
            line = re.sub(bracket_regex, "", line)
            multi_handle.write(line)

os.remove(output_tmp)

save_seq(TAL_ID_list, EBEseq_list, output_seq, selector_set)

print("done!")
