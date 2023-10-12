for ((i=3;i<=31;i+=2))
do
    echo "$i"
    time ./UpPipe build -i ~/data23/Sac_"$i"mer_60 -f ~/data/Saccharomyces_cerevisiae/Saccharomyces_cerevisiae.R64-1-1.cdna.all.fa -k "$i" -d 60
done

for ((i=31;i>=3;i-=2))
do
    echo "$i"
    time ./UpPipe alignment -i ~/data23/Sac_"$i"mer_60 -o out  -f ~/data/experiment/RNA_read/100K.fastq -r 40
done