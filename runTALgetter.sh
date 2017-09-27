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
output_full="~/TALgetter/$now-TALgetter_output.tsv"
uo=0
do=0
input_rvd=""
TALVEZ="false"
pval="COARSE"
pthresh=1.0
top=100
strand="false"

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -input)
    input_dna="$2"
    shift # past argument
    ;;
    -o)
    output_full="$2"
    shift # past argument
    ;;
    -uo)
    uo="$2"
    shift # past argument
    ;;
    -do)
    do="$2"
    shift # past argument
    ;;
    -rvd)
    input_rvd="$2"
    shift # past argument
    ;;
    -f)
    TALVEZ="$2"
    shift # past argument
    ;;
    -pval)
    pval="$2"
    shift # past argument
    ;;
    -pthresh)
    pthresh="$2"
    shift # past argument
    ;;
    -t)
    top="$2"
    shift # past argument
    ;;
    -strand)
    strand="$2"
    shift # past argument
    ;;
    *)
            # unknown option
    ;;
esac
shift # past argument or value
done
# based on code by Bruno Bronosky:
# http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash/14203146#14203146

filename=$(basename $output_full)
output=$(dirname $output_full)
output=$output"/"


if [ $input_dna == "?" ]; then
	echo ""
	echo "Usage: bash runTALgetter.sh [options]"
	echo "Options:"
	echo "-input	| String. Complete path to a (multi-)FASTA file containing the DNA sequences to be searched for EBEs."
	echo "-uo	| Integer. Upstream offset. Number of positions ignored at 5' end of each sequence. Default: 0"
	echo "-do	| Integer. Downstream offset. Number of positions ignored at 3' end of each sequence. Default: 0"
	echo "-model	| String. Model type. TALgetter is the default model that uses individual binding specificities for each RVD. TALgetter13 uses binding specificities that only depend on amino acid 13, i.e., the second amino acid of the repat.While TALgetter is recommended in most cases, the use of TALgetter13 may be beneficial if you search for target sites of TAL effector with many rare RVDs, for instance YG, HH, or S*. Options: TALgetter, TALgetter13. Default: TALgetter"
	echo "-pval	| String. PVals. Computation of p-Values. Options: NONE, COARSE, FINE. Default: COARSE"
	echo "-pthresh	| Float. p-Value. Filter the reported hits by a maximum p-Value. A value of 0 or 1 switches off the filter. Options: real number between 0.0 and 1.0, both included. Default: 1.0"
	echo "-top	| Integer. Maximum number of target sites. Limits the total number of reported target sites in all input sequences, ranked by their score. Options: integer number between 1 and 10000, both included. Default: 100"
	echo "-strand	| String. Both strands. Search both strands of the input sequence for target sites. Options: true, false. Default: false"
	echo "-train	| String. OPTIONAL! Training data. The input data to use for training the model. Supply the path to an annotated FASTA file."
	echo "-o	| String. Full path for the output file. Default: $output$filename"
	echo "-rvd	| String. Complete path of a file containing TALE RVDs in TALVEZ/TALgetter format:"
	echo "	>AvrBs4<tab>NI-NG-NI-NI-NG-NG-NI-NS-NG-NI-NS-NG-HD-HD-NS-HD-NG-NG"
	echo "	Note: replace <tab> with an actual tab!"
	echo "-format	| String. Is the -rvd in TALVEZ/TALgetter format? Options: true, false. Default: false"
	echo ""
	echo ">>> See also below for TALgetter's own help! However, this information overrides any information given by TALgetter itself."
	echo ""
	echo ""
	TALgetter
else

	if [ ! -d $output ]; then
		mkdir -p $output;
	fi
	cd $output

	mapfile rvdArray < $input_rvd

	rm -f $filename
	touch $filename

	firstline=true
	for rvd_line in "${rvdArray[@]}"; do

		if [ "$TALVEZ" = "true" ]; then
			rvd="$(cut -d' ' -f2 <<<$rvd_line)"
			TALE_name="$(cut -d' ' -f1 <<<$rvd_line)"
			TALE_name="$(cut -d'>' -f2 <<<$TALE_name)"

		else
			rvd=$rvd_line
		fi

		echo $rvd

		TALgetter "input="$input_dna "uo="$uo "do="$do "rvd="$rvd "pval="$pval "pthresh="$pthresh "top="$top "strand="$strand > $now"_TALgetter_tmp_"$TALE_name

		touch $now"_TALgetter_tmp_"$TALE_name"_new"
		while read line
		do
			if [[ "$line" == "#"* ]]; then
				if [ $firstline = true ]; then
					header="TALE GeneID Position Distance_to_end Sequence Matches Score Strand p-value E-value"
					echo -e $header >> $now"_TALgetter_tmp_"$TALE_name"_new"
					firstline=false
				fi
			elif [[ "$line" == "" ]]; then
				echo ""
			else
				echo -e $TALE_name" "$line >> $now"_TALgetter_tmp_"$TALE_name"_new"
			fi
		done < $now"_TALgetter_tmp_"$TALE_name

		cat $now"_TALgetter_tmp_"$TALE_name"_new" >> tmp_$filename
		rm $now"_TALgetter_tmp_"$TALE_name
		rm $now"_TALgetter_tmp_"$TALE_name"_new"

	done

awk -v OFS="\t" '$1=$1' "tmp_$filename" > $filename
rm "tmp_$filename"

fi

echo "Done."

