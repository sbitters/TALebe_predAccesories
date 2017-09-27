#!/bin/bash
#
#
# Copyright 2017 Sven T. Bitters (sven.bitters@gmail.com)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of
# this software and associated documentation files (the "Software"), to deal in
# the Software without restriction, including without limitation the rights to
# use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
# the Software, and to permit persons to whom the Software is furnished to do so,
# subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
# FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
# COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
# IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
# CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


now="$(date +'%Y%m%d%H%M%S')"

input_dna="?"
cupstream=2
forwardonly=false
numprocs=1
outfile="~/TALENT/$now-TALENT_result"
weight=0.9
cutoffmult=3.0
input_rvd=""
TALVEZ="false"
help=false

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -input)
    input_dna="$2"
    shift # past argument
    ;;
    -c)
    cupstream="$2"
    shift # past argument
    ;;
    -f)
    forwardonly="$2"
    shift # past argument
    ;;
    -n)
    numprocs="$2"
    shift # past argument
    ;;
    -o)
    outfile="$2"
    shift # past argument
    ;;
    -w)
    weight="$2"
    shift # past argument
    ;;
    -x)
    cutoffmult="$2"
    shift # past argument
    ;;
    -rvd)
    input_rvd="$2"
    shift # past argument
    ;;
    -format)
    TALVEZ="$2"
    shift # past argument
    ;;
    *)
		help=true
    ;;
esac
shift # past argument or value
done
# based on code by Bruno Bronosky:
# http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash/14203146#14203146

output_dir=$(dirname $outfile)
output_dir="$output_dir/"
outfile_orig="$outfile"

if [ "$forwardonly" = true ]; then
	fo="--forwardonly"
else
	fo=""
fi

if [ $input_dna == "?" ]; then
	echo ""
	echo "Usage: bash runTALENT.sh [options]"
	echo "Options:"
	echo "-input	| String. Complete path to a (multi-)FASTA file containing the DNA sequences to be searched for EBEs."
	echo "-c	| Integer. Sets the allowed upstream bases; 0 for T only, 1 for C only, 2 for either. Default: 2"
	echo "-f	| String. Only search the forward strand of the sequence? Options: true, false. Default: false"
	echo "-n	| Integer. The number of processors to use. Default: 1"
	echo "-o	| String. Template filename to which output will be written. Both a tab-delimited file and gff3 file will be produced. Default: $outfile.gff3"
	echo "-w	| Float. User-defined weight. Default: 0.9"
	echo "-x	| Float. Multiple of best score at which potential sites will be filtered. Default: 3.0"
	echo "-rvd	| String. Complete path of a file containing TALE RVDs - either in TALVEZ/TALgetter format:"
	echo "	>AvrBs4<tab>NI-NG-NI-NI-NG-NG-NI-NS-NG-NI-NS-NG-HD-HD-NS-HD-NG-NG"
	echo "	or in TALENT format:"
	echo "	AvrBs4<tab>NI NG NI NI NG NG NI NS NG NI NS NG HD HD NS HD NG NG"
	echo "	Note: replace <tab> with an actual tab!"
	echo "-format	| String. Is the -rvd in TALVEZ/TALgetter format? Options: true, false. Default: false"
	echo ""
	echo ">>> See also below for TALENT's own -h info! However, this information overrides any information given by TALENT itself."
	echo ""
	echo ""
	talesf "-h"
else

	if [ ! -d $output_dir ]; then
		mkdir -p $output_dir;
	fi
	cd $output_dir

	touch "$outfile_orig"

	mapfile rvdArray < $input_rvd

	firstline=true
	for rvd_line in "${rvdArray[@]}"; do

		if [ "$TALVEZ" = "true" ]; then
			rvd="$(cut -d' ' -f2 <<<$rvd_line)"
			rvd=${rvd//[-]/ }

			TALE_name="$(cut -d' ' -f1 <<<$rvd_line)"
			TALE_name="$(cut -d'>' -f2 <<<$TALE_name)"

		else
			rvd=$rvd_line
		fi
		outfile="$output_dir$now-$TALE_name-prim-tmp"
		talesf "-c "$cupstream $fo "-n "$numprocs "-o"$outfile "-w "$weight "-x "$cutoffmult $input_dna "$rvd"

		touch "$output_dir$now-$TALE_name-tmp"
		while read line
		do
			if [[ "$line" == "#"* ]]; then
				if [ $firstline = true ]; then
					header="TALE GeneID Source Type Start End Score Strand Qual Target_RVD Target_Sequence Match_Sequence"
					echo -e $header >> "$output_dir$now-$TALE_name-tmp"
					firstline=false
				fi
			elif [[ "$line" == "" ]]; then
				echo ""
			else
				echo -e $TALE_name" "$line >> "$output_dir$now-$TALE_name-tmp"
			fi
		done < "$outfile.gff3"

		cat "$output_dir$now-$TALE_name-tmp" >> "$outfile_orig-tmp"
		rm "$output_dir$now-$TALE_name-tmp"
		rm "$outfile.gff3"
		rm "$outfile.txt"

	done

awk -v OFS="\t" '$1=$1' "$outfile_orig-tmp" > "$outfile_orig.gff3"
rm "$outfile_orig-tmp"
rm "$outfile_orig"

fi

echo "Done."

