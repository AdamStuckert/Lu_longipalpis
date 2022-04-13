#/bin/sh
# USAGE: ./phased_iqtree.sh -d ${SCR}/longipalpis_ms/phased
while getopts d: option
do
case "${option}"
in
d) DIR=${OPTARG};;
esac
done

## $1 == directory with phased phylip files

#cd $DIR
mkdir ${DIR}/iqtree

CURDIR=$(pwd)
mkdir phased_iqtree
#cd phased_iqtree

cd $DIR

for phyl in $(ls *phylp)
do
win=$(basename $phyl | sed "s/.phylp//g")

cat << EOF > ${CURDIR}/phased_iqtree/Tree_${win}.job
#!/bin/bash
#SBATCH -p general
#sbatch -J Tree
#SBATCH --cpus-per-task=4
#SBATCH -t 5:00:00
#SBATCH --mem=4g
#SBATCH -o trees_%j.out


ml iqtree/1.6.12

cd ${DIR}
iqtree -redo -pre iqtree/${win} -st DNA -nt 4 -ntmax 4 -bb 1000 -s ${phyl}

EOF

sbatch ${CURDIR}/phased_iqtree/Tree_${win}.job
done
