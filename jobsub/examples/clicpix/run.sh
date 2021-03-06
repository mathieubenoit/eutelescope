#!/bin/sh

first="4395"
last="4410"
RUNLIST="runlist-150.csv"

#first="33"
#last="33"
#RUNLIST="runlist-20.csv"

modus="straight"
#modus="daf"
#modus="gbl"

DRY=

for i in `seq $first $last`; do

 jobsub.py $DRY -c configTP3Math.cfg  converter  $i
 jobsub.py  $DRY -c configTP3Math.cfg  clusearch $i
 #jobsub.py  $DRY -c configTP3Math.cfg -csv $RUNLIST filter $i

if [[ $modus == "straight" ]]; then
 jobsub.py $DRY -c configTP3Math.cfg  hitmaker   $i
# alignment using straight line assumption
 jobsub.py $DRY -c configTP3Math.cfg  align      $i
# fitter using broken line implementation by F.Zarnezki
 jobsub.py $DRY -c configTP3Math.cfg  fitter $i

elif [[ $modus == "daf" ]]; then
 jobsub.py $DRY -c configTP3Math.cfg -csv $RUNLIST hitmaker   $i
 jobsub.py $DRY -c configTP3Math.cfg -csv $RUNLIST aligndaf   $i
 jobsub.py $DRY -c configTP3Math.cfg -csv $RUNLIST trackdaf   $i

elif [[ $modus == "gbl" ]]; then
# FOR NEW gbl ONE NEEDS TO GET HITS IN LOCAL COORDINATE SYSTEM:
 jobsub.py  $DRY -c configTP3Math.cfg -csv $RUNLIST hitlocal $i
# Exhautsive and Helix track search results are not identical - to be investigated (perhaps trivially explained)
# jobsub.py  $DRY -c configTP3Math.cfg -csv $RUNLIST tracksearchExh $i
 jobsub.py  $DRY -c configTP3Math.cfg -csv $RUNLIST tracksearchHelix $i
# echo jobsub.py  $DRY -c configTP3Math.cfg -csv $RUNLIST aligngbl $i
#
 bash alignment.sh -r ${i} -l ${RUNLIST}
#### jobsub.py  $DRY -c configTP3Math.cfg -csv $RUNLIST aligngbl $i
#
# 10 iterations to aligne 6 planes (6D)
 jobsub.py  $DRY -c configTP3Math.cfg -o GearAlignedFile="gear-${i}-010.xml" -csv $RUNLIST trackgbl $i

fi
done

#
