old=20000
for ii in `seq 0 7`; do 
    new=$(echo "20000 + $ii*1024" | bc -l)
    echo "Replacing $old with $new"
    sed -i.bak "s/$old/$new/g" mk_mdls.sh
    #cat mk_mdls.sh 
    sbatch mk_mdls.sh 
    old=$new
done 

sed -i "s/$old/20000/g" mk_mdls.sh
