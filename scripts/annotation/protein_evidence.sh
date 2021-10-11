# record for generating protein hints needed to annotation Silene conica genome

# First pull the data from google drive
gdown https://drive.google.com/file/d/1J9d9iIaG3eJgg6zWwAIybnhdQKjkwqPG/view?usp=sharing

# Now we need to collapse isoforms for datasets we don't have a genome for. We do this with cd-hit, following the steps described in 
# cDNA-Cupcake: https://github.com/Magdoll/cDNA_Cupcake

# First, rename the headers in fasta since cd-hit-est doesn't like the isoseq headers

awk '/^>/{print ">" ++i; next}{print}' agrostemma.clustered.hq.fasta > agrostemma.clustered.hq_mod.fasta

# now cd-hit
cd-hit-est -i agrostemma.clustered.hq_mod.fasta -o agrostemma.clustered.hq_mod.collapse.fasta -c 0.99 -T 6 -G 0 -aL 0.90 -AL 100 -aS 0.99 -AS 30 

# this process is repeated for each of the non-conica iso-seq datasets

# Now let's collapse isoforms for S. conica where we do have a genome to map to

# first we map the reads with minimap2

minimap2 -t 55 -ax splice -uf --secondary=no -C5 -O6,24 -B4 s_conica_v.0.2.racon6.fasta conica.clustered.hq.fasta > hq_isoforms.fasta.sam 2> hq_isoforms.fasta.sam.log

# now we use cDNA-Cupcake scripts to collapse isoforms
sort -k 3,3 -k 4,4n hq_isoforms.fastq.sam > hq_isoforms.fastq.sorted.sam
collapse_isoforms_by_sam.py --input conica.clustered.hq.fasta \
-s hq_isoforms.fasta.sorted.sam --dun-merge-5-shorter -o sconica.transcripts

# we now have our target transcripts in the form of sconica.transcripts.collapsed.rep.fa

# Now we need to identify the coding portions of all of our transcripts and also translate them
# We'll use the Transdecoder pipeline to do this. Assume we've setup both the swiss_prot and pfam databases for the
# following to work. The basic outline is found here: https://github.com/TransDecoder/TransDecoder/wiki

TransDecoder.LongOrfs -t agrostemma.clustered.hq_mod.collapse.fasta

# Now we'll use diamond (v.2.0.11) to generate blast evidence against swissprot
diamond blastp --query agrostemma.clustered.hq_mod.collapse.fasta --max-target-seqs 1 --ultra-sensitive --threads 55 --db uniprot_sprot.dmnd --evalue 1e-5 -out --outfmt 6 --out diamond.agrostemma.outfmt6

# Now let's use hmmer to generate matches against pfam. This is super slow so I'm going to use parallels, see script hmmscan.sh, though the basic form is:
hmmscan --cpu 1 --domtblout ${sample}.pfam.domtblout ~/bioinformatics/pfam/Pfam-A.hmm ${sample}

# let's split up the longest_orfs.pep into 100 smaller fasta files; pyfasta can be pulled using bioconda
pyfasta split -n 100 longest_orfs.pep

# and we call it as follows:
ls longest_orfs.pep.* | parallel -j55 -k bash hmscan.sh {}

# because the resultant files have a header we need to cut this and then cat these files; see trim_header.sh and use the same as hmmscan.sh
ls *.pfam.domtblout | parallel -j12 -k bash trim_headers.sh {}
cat *.out > compile.pfam.domtblout
head -n 3 longest_orfs.pep.000.pfam.domtblout | cat - compile.pfam.domtblout > pfam.agrostemma.domtblout

# now we can finish up the Transdecoder pipeline ande get our target proteins
TransDecoder.Predict -t agrostemma.clustered.hq_mod.collapse.fasta --retain_pfam_hits pfam.agrostemma.domtblout \
 --retain_blastp_hits diamond.agrostemma.outfmt6

# Because of the original renaming we did above we need to make sure our different protein files have unique names, so we do:
awk '/^>/{print ">agrostemma" ++i; next}{print}' agrostemma.clustered.hq_mod.collapse.fasta.transdecoder.pep \
> agrostemma.clustered.hq_mod.collapse.fasta.transdecoder.rename.pep

# Now we just cat things up and get ready to run MAKER2!
cat *transdecoder.rename.pep > proteins.fasta
