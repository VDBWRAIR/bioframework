if [ -e  $boo ]; 
then 
    echo "Usage: $(basename $0) <ngsdir>";
    exit 1
fi;

make NGSDIR=$1
