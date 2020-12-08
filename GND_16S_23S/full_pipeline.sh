#! /bin/bash

# run TreeN93
bash run_treeN93.sh 

# ani binning
python ani_binning.py
sed -i '' -e "s/\[\&R\] //g" *nwk

# compute support for the samples
for x in 16S*nwk; do
    outfile=`echo $x | sed -e "s/16S/16S_support/g"`
    while read t16; do
        bs=`mktemp`
        tr=`mktemp`
        echo $t16 > $tr
        nw_prune -v RAxML_bootstrap.T16S-SINA-1000-bootcutoff `nw_labels -I $tr` > $bs
        nw_support -p $tr $bs >> $outfile
        rm $bs $tr
    done <$x
done    

for x in 23S*nwk; do
    outfile=`echo $x | sed -e "s/23S/23S_support/g"`
    while read t23; do
        bs=`mktemp`
        tr=`mktemp`
        echo $t23 > $tr
        nw_prune -v RAxML_bootstrap.T23S-SINA-1000-bootcutoff `nw_labels -I $tr` > $bs
        nw_support -p $tr $bs >> $outfile
        rm $bs $tr
    done <$x
done    

# collapse bootstrap support and compute qdist
for x in 16S_support*nwk; do
    y=`echo $x | sed -e "s/16S/23S/g"`
    outfile75=`echo $x | sed -e "s/16S/qdist75_weighted/g" -e "s/\.nwk/\.txt/g"`
    tr1=`mktemp`
    tr2=`mktemp`

    tempout1=`mktemp`
    tempout2=`mktemp`
    tempout3=`mktemp`

    nw_ed $x 'i & b < 75' o > $tr1
    nw_ed $y 'i & b < 75' o  > $tr2
    pairs_quartet_dist -v $tr1 $tr1 | awk '{print $7;}' > $tempout1
    pairs_quartet_dist -v $tr2 $tr2 | awk '{print $7;}' > $tempout2
    pairs_quartet_dist -v $tr1 $tr2 | awk '{print $2,$5,$7;}' > $tempout3
    paste $tempout1 $tempout2 $tempout3 | awk '{print 1-($4+($1+$2-$5)/3)/$3,$4,($1+$2-$5)/3,$3+$5-$1-$2,$5;}' > $outfile75
    
    rm $tempout1 $tempout2 $tempout3
done    

grep . pdist_* | sed -e "s/pdist_//g" -e "s/\.txt:/ /g" > pdist.dat

# combine data
for treefile in 16S_[0-9]*nwk; do
    range=`echo $treefile | sed -e "s/16S_//g" -e "s/\.nwk//g"`
    pdistfile=`echo $treefile | sed -e "s/16S/pdist/g" -e "s/nwk/txt/g"`
    qdistfile=`echo $treefile | sed -e "s/16S/qdist75_weighted_support/g" -e "s/nwk/txt/g"`
    temp1=`mktemp`
    temp2=`mktemp`
    temp3=`mktemp`
# leaf labels
    while read tree; do 
        echo $tree | nw_labels -I -t - 
    done < $treefile > $temp1
# pdist
    i=45
    while [ $i -le 2250 ]; do 
        head -n$i $pdistfile | tail -n45 | numlist -med; i=$((i+45))
    done > $temp2
    i=45
    while [ $i -le 2250 ]; do 
        head -n$i $pdistfile | tail -n45 | numlist -avg; i=$((i+45))
    done > $temp3
# combine
    paste $temp1 $temp2 $temp3 $qdistfile | sed -e "s/^/$range /g" | tr "\t" " "
    rm $temp1 $temp2 $temp3
done > combined_data.txt   
