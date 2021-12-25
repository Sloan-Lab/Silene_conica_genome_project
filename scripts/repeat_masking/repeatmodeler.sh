# environment created with conda, see conda.lst for necessary packages and versions.

# Build RepeatModeler database
BuildDatabase -name s_conica s_conica_v.0.2.racon6.fasta

# Run RepeatModeler to discover repeat families
RepeatModeler -pa 28 -LTRStruct -database s_conica

# I just want to use the TEs for masking, avoid the Unknowns - some of those might be genes
grep "DNA\|RC\|LTR\|LINE\|SINE" consensi.fa.classified | sed 's/>//g' | awk '{print $1}' > conica_TEs

# Let's linearize fasta sequences (see fastx-toolkit, http://hannonlab.cshl.edu/fastx_toolkit/commandline.html)
fasta_formatter -i consensi.fa.classified -o conica_assembly-families.fa -w 0
#while read line; grep -A1 $line conica_assembly-families.fa >> conica_TEs.fa; done < conica_TEs

# Let's extract the target TE families we want masked (see seqtk, https://github.com/lh3/seqtk)
seqtk subseq conica_assembly-families.fa conica_TEs > conica_TEs.fa

# Now we can run RepeatModeler on the assembly and generate gff of repeats
RepeatMasker -pa 55 -no_is -norna -a -gff -lib conica_TEs.fa -dir RMask_conica_TEs_noIS s_conica_v.0.2.racon6.fasta

# Generate the gff for complex repeats to use in MAKER
rmOutToGFF3.pl s_conica_v.0.2.racon6.fasta.out > s_conica_v.0.2.racon6.fasta.out.gff3

# isolate complex repeats
grep -v -e "Satellite" -e ")n" -e "-rich" s_conica_v.0.2.racon6.fasta.out.gff3 > s_conica_v.0.2.racon6.fasta.out.complex.gff3

# reformat to work with MAKER
cat s_conica_v.0.2.racon6.fasta.out.complex.gff3 | \
perl -ane '$id; if(!/^\#/){@F = split(/\t/, $_); chomp $F[-1];$id++; $F[-1] .= "\;ID=$id"; $_ = join("\t", @F)."\n"} print $_' \ 
> s_conica_v.0.2.racon6.fasta.out.complex.reformat.gff3
