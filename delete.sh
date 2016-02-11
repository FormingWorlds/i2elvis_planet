#!/bin/sh

#get simulation name
SIMNAME=$(pwd | cut -d "/" -f 7)
COUNT=$(ls $SIMNAME*.prn | wc -l)
FILELIST=`ls $SIMNAME*.prn`
LENGTH=`expr length $SIMNAME + 1`

#get latest file
#echo `ls $SIMNAME*.prn | cut -d "." -f 1 | cut -b $LENGTH- | sort -n | tail -n 1`
LATEST=`ls $SIMNAME*.prn | cut -d "." -f 1 | cut -b $LENGTH- | sort -n | tail -n 1`
LATEST=${SIMNAME}${LATEST}".prn"

#delete all *0.prn, *2-9.prn BUT NOT latest file
for FILENAME in $FILELIST ; do
	if [[ "$FILENAME" != [A-Za-z0-9]*1.prn ]] ; then
		if [ "$FILENAME" != "$LATEST" ] ; then
			rm ${FILENAME}
		fi
	fi
done

rm lsf*


