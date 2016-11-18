#!/bin/sh

#get list of all data dumps in folder
FILELIST=$(ls *.prn)

#get latest file
LATEST=$(ls *.prn | cut -d "." -f 1 | cut -b 4- | sort -n | tail -n 1)
LATEST=dd_${LATEST}.prn

#delete all *0.prn, *2-9.prn BUT NOT latest file
for FILENAME in $FILELIST ; do
        if [[ "$FILENAME" != [A-Za-z0-9].prn ]] ; then
                if [ "$FILENAME" != "$LATEST" ] ; then
                        rm ${FILENAME}
                fi
        fi
done
#rm lsf*
