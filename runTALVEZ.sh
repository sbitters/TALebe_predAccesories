#!/bin/bash

# Wrapper script for TALVEZ by Alvaro L Perez-Quintero
#
# Pérez-Quintero AL, Rodriguez-R LM, Dereeper A, López C, Koebnik R, et al. 
# (2013) An Improved Method for TAL Effectors DNA-Binding Sites Prediction 
# Reveals Functional Convergence in TAL Repertoires of Xanthomonas oryzae 
# Strains. PLOS ONE 8(7): e68464. https://doi.org/10.1371/journal.pone.0068464
#
# The variables $talfile and $fastafile must contain absolute paths to their 
# respective TAL_file or FASTA file. All other all other variables containing 
# the path to a file must be given relative to the home directory of 
# TALVEZ_X.X.pl (e.g. ~/bin/TALVEZ_3.2/). The $output variable is
# merely a prefix for other files and must not contain any path info. The 
# default parameters of TALVEZ are retained in this script - if you do 
# not specify your own value, the default value will be handed to TALVEZ.
#
# Output files wil be located in the same directory as your TAL_file.
#
# In order to use this script, you must add TALVEZ_X.X.pl and its directory to 
# your path variable!
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

TALVEZ_loc=$(which TALVEZ)
TALVEZ_file=$(basename $TALVEZ_loc)
TALVEZ_dir=$(dirname $TALVEZ_loc)

# Default values as per TALVEZ_3.2
t_arg=0
a_arg=6
l_arg=0
v_arg=0.0001
e_arg="mat1"
z_arg="mat2"
p_arg="tmp"
talfile=""
fastafile=""
output="output"

while [[ $# -gt 1 ]]
do
key="$1"

case $key in
    -t)
    t_arg="$2"
    shift
    ;;
    -a)
    a_arg="$2"
    shift
    ;;
    -l)
    l_arg="$2"
    shift
    ;;
    -v)
    v_arg="$2"
    shift
    ;;
    -e)
    e_arg="$2"
    shift
    ;;
    -z)
    z_arg="$2"
    shift
    ;;
    -p)
    p_arg="$2"
    shift
    ;;
    --TALfile)
    talfile="$2"
    shift
    ;;
    --FASTAfile)
    fastafile="$2"
    shift
    ;;
    -o|--output)
    output="$2"
    shift
    ;;
    *)
            # unknown option
    ;;
esac
shift
done
# based on code by Bruno Bronosky:
# http://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash/14203146#14203146


if [ $1 == "-h" ] || [ $1 == "--help" ] || [ -z "${talfile}" ] || [ -z "${fastafile}" ]; then
	perl TALVEZ -h
else

	maindir=$(dirname $talfile)
	maindir="$maindir"

	talfile=$(basename $talfile)
	fastafile=$(basename $fastafile)

	logfile=$maindir/$now'_TALVEZ_log.txt'

	{
	if [ "$maindir" != "$TALVEZ_dir" ]; then
		cp -af $TALVEZ_dir/$e_arg $maindir/$e_arg
		cp -af $TALVEZ_dir/$z_arg $maindir/$z_arg
		cp -af $TALVEZ_dir/simplescancode $maindir/simplescancode
		cp -af $TALVEZ_loc $maindir/$TALVEZ_file

		chmod 777 $maindir/$e_arg
		chmod 777 $maindir/$z_arg
		chmod -R 777 $maindir/simplescancode
		chmod 777 $maindir/$TALVEZ_file
	fi

	cd $maindir
	echo "perl TALVEZ.pl -t $t_arg -a $a_arg -l $l_arg -v $v_arg -e $e_arg -z $z_arg -p $p_arg $talfile $fastafile $output"
	perl TALVEZ.pl -t $t_arg -a $a_arg -l $l_arg -v $v_arg -e $e_arg -z $z_arg -p $p_arg $talfile $fastafile $output

	if [ "$maindir" != "$TALVEZ_dir" ]; then
		rm $maindir/$e_arg
		rm $maindir/$z_arg
		rm -r $maindir/simplescancode
		rm -r $maindir/$p_arg
		rm $maindir/$TALVEZ_file
	fi

	echo "Output files are located in $maindir"
	echo "Done."
	} 2>&1 | tee -a $logfile

fi


#last edited by STB, 2017-04-27
