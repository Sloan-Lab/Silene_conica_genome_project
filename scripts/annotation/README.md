# Round 1 MAKER
###### (see https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2 for more information)

##### We've run
`bash ./round1_run_maker.sh 2>&1 | tee round1_run_maker.log`

##### Now let's process the outputs of the first round of MAKER
```
cd sconica_rnd1.maker.output
gff3_merge -s -d sconica_rnd1_master_datastore_index.log > sconica_rnd1.all.maker.gff
fasta_merge -d sconica_rnd1_master_datastore_index.log
# GFF w/o the sequences
gff3_merge -n -s -d sconica_rnd1_master_datastore_index.log > sconica_rnd1.all.maker.noseq.gff
```
##### Also get a summary (summary.sh originates from the falling repo https://github.com/ISUgenomics/common_scripts.git)
`less sconica_rnd1.all.maker.gff |awk '$3=="mRNA"' |grep "mRNA-1" |awk '{print $5-$4}' | summary.sh`


##### Now let's train SNAP
```
mkdir snap
mkdir snap/round1
cd snap/round1
# export 'confident' gene models from MAKER and rename to something meaningful
maker2zff -x 0.25 -l 50 -d ../../sconica_rnd1.maker.output/sconica_rnd1_master_datastore_index.log
rename 's/genome/sconica_rnd1.zff.length50_aed0.25/g' *
# gather some stats and validate
fathom sconica_rnd1.zff.length50_aed0.25.ann sconica_rnd1.zff.length50_aed0.25.dna -gene-stats > gene-stats.log 2>&1
fathom sconica_rnd1.zff.length50_aed0.25.ann sconica_rnd1.zff.length50_aed0.25.dna -validate > validate.log 2>&1
# collect the training sequences and annotations, plus 1000 surrounding bp for training
fathom sconica_rnd1.zff.length50_aed0.25.ann sconica_rnd1.zff.length50_aed0.25.dna -categorize 1000 > categorize.log 2>&1
fathom uni.ann uni.dna -export 1000 -plus > uni-plus.log 2>&1
# create the training parameters
mkdir params
cd params
forge ../export.ann ../export.dna > ../forge.log 2>&1
cd ..
# assembly the HMM
hmm-assembler.pl sconica_rnd1.zff.length50_aed0.25 params > sconica_rnd1.zff.length50_aed0.25.hmm
```

##### Now let's train AUGUSTUS
###### cut out the genes 
```
awk -v OFS="\t" '{ if ($3 == "mRNA") print $1, $4, $5 }' ../../sconica_rnd1.maker.output/sconica_rnd1.all.maker.noseq.gff | \
  awk -v OFS="\t" '{ if ($2 < 1000) print $1, "0", $3+1000; else print $1, $2-1000, $3+1000 }' | \
  bedtools getfasta -fi ../../s_conica_v.0.2.racon6.fasta -bed - -fo sconica_rnd1.all.maker.transcripts1000.fasta
```
###### run BUSCO to train AUGUSTUS
`run_BUSCO.py -i sconica_rnd1.all.maker.transcripts1000.fasta -o sconica_rnd1_maker -l ~/bioinformatics/busco_db/eudicotyledons_odb10 -m genome -c 56 --long -z --augustus_parameters='--progress=true'`



