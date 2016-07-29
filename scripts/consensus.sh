if [ -z  "$1" ]; 
then 
    echo "Usage: $(basename $0) <ngsdir>";
    exit 1
fi;

make -f $(dirname $0)/Makefile NGSDIR=$1
