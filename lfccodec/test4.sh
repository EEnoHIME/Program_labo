#! /bin/sh

#for image in 108 109 110 111 112 113
for image in 4
do

DIR="/rda2/DATABASE/TMW"
DIR2="/rda1/users/miyaza/lfc_pgm"
DIR3="/rda1/users/miyaza/lfc_pgm/Buildings_8bit"
FLAG="-o -f"$@""

## IMAGE ##
 if   [ "$image" = 1 ]; then IMG=$DIR2/output0007.pgm;     NAME="image0007"
 elif [ "$image" = 2 ]; then IMG=$DIR2/output0017.pgm;     NAME="image0017"
 elif [ "$image" = 3 ]; then IMG=$DIR2/output0019.pgm;     NAME="image0019"
 elif [ "$image" = 4 ]; then IMG=$DIR2/output0027.pgm;     NAME="image0027"
 
 elif [ "$image" = 5 ]; then IMG=$DIR3/Black_Fence_LFR_RawImg.pgm;     NAME="Black_Fence"
 elif [ "$image" = 6 ]; then IMG=$DIR3/Palais_du_Luxembourg_LFR_RawImg.pgm;       NAME="Palais_du_Luxembourg"
 elif [ "$image" = 7 ]; then IMG=$DIR3/Pillars_LFR_RawImg.pgm;       NAME="Pillars"
 elif [ "$image" = 8 ]; then IMG=$DIR3/Red_&_White_Building_LFR_RawImg.pgm;     NAME="Red_&_White_Building"
 elif [ "$image" = 9 ]; then IMG=$DIR2/Rolex_Learning_Center_LFR_RawImg.pgm;         NAME="Rolex_Learning_Center"
 
 elif [ "$image" = 10 ]; then IMG=$DIR/lennagrey.pgm;   NAME="lennagrey"
 elif [ "$image" = 11 ]; then IMG=$DIR/noisesquare.pgm; NAME="noisesquare"
 elif [ "$image" = 12 ]; then IMG=$DIR/peppers.pgm;     NAME="peppers"
 elif [ "$image" = 13 ]; then IMG=$DIR/shapes.pgm;      NAME="shapes"


 fi

LOG="LOG/$NAME.txt"
MRP="mrpfile/$NAME.mrp"
PGM="Decoded_File/$NAME.pgm"

echo $FLAG > $LOG
echo `hostname` >> $LOG
	whoami | tee -a $LOG
	date | tee -a $LOG
	nice ./encmrp $FLAG $IMG $MRP | tee -a $LOG
	nice ./decmrp $MRP $PGM | tee -a $LOG
	if cmp -s $IMG $PGM;
	then
		echo "OK!" | tee -a $LOG
    rm -f $PGM
	else
		echo "ERROR!" | tee -a $LOG
    rm -f $PGM
		exit
	fi
done
