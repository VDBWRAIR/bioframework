if [ -z  "$1" ]; 
then 
    echo "Usage: $(basename $0) <ngsdir>";
    exit 1
fi;
makefile=/media/VD_Research/Admin/PBS/Software/bioframework/scripts/Makefile

# make -f $(dirname $0)/Makefile NGSDIR=$1
make -f $makefile NGSDIR=$1
