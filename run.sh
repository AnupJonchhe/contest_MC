if [ $# -ne 5 ]
then
  echo Usage ./run.sh upfile downfile outputfile_on_ssd1 indexfilesdirectory_on_ssd1 indexfilesdirectory_on_ssd2
  exit 1
fi

upfile=$1
downfile=$2
outfile=$3
path=$4
path2=$5

./prog_avx_final --path $path --path2 $path2 --out $outfile --up $upfile --down $downfile --num-threads 16 --test
