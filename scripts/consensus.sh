if [ -z  "$1" ] || [ -z "$2" ]; 
then 
    echo "Usage: $(basename $0) <ngsdir> <reference>";
    exit 1
fi;
makefile=/media/VD_Research/Admin/PBS/Software/bioframework/scripts/Makefile

# make -f $(dirname $0)/Makefile NGSDIR=$1
make -f $makefile NGSDIR=$1 REF=$2
