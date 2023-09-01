#!/usr/bin/perl

use strict;
use warnings;
use Bio::SearchIO; 
use sloan;
use List::Util qw(min max);

my $usage = "\nUSAGE: $0 blastFile outputName Evalue_cutoff\n\n";

my $blastFile = shift or die ($usage);
my $output = shift or die ($usage);
my $evalue_cutoff = shift or die ($usage);

my $FH1 = open_output ("$output\.hsps.queries.txt");
my $FH2 = open_output ("$output\.hsps.hits.txt");


print $FH1 "SeqName\tSeqLength\tHitStart\tHitEnd\tHitLength\tPercentID\tE-value\n";
print $FH2 "SeqName\tSeqLength\tHitStart\tHitEnd\tHitLength\tPercentID\tE-value\n";


my $SearchIO_obj = new Bio::SearchIO(-format => 'blast', -file   => $blastFile);


while (my $result_obj = $SearchIO_obj->next_result){
	my $query_name = $result_obj->query_name;
	my $query_len = $result_obj->query_length;
	while (my $hit_obj = $result_obj->next_hit){
		my $hit_name = $hit_obj->name;
		my $hit_len = $hit_obj->length;
  		while (my $hsp_obj = $hit_obj->next_hsp){
			my $hsp_evalue = $hsp_obj->evalue;
			$hsp_evalue <= $evalue_cutoff or last;
			my $percentID = 100 * $hsp_obj->frac_identical;
			my $query_start = $hsp_obj->start('query');
			my $query_end = $hsp_obj->end('query');
			my $hit_start = $hsp_obj->start('hit');
			my $hit_end = $hsp_obj->end('hit');

			print $FH1 "$query_name\t$query_len\t", min ($query_start, $query_end), "\t", max ($query_start, $query_end), "\t", max ($query_start, $query_end) - min ($query_start, $query_end)	+ 1, "\t$percentID\t$hsp_evalue\n";		
			print $FH2 "$hit_name\t$hit_len\t", min ($hit_start, $hit_end), "\t", max ($hit_start, $hit_end), "\t", max ($hit_start, $hit_end) - min ($hit_start, $hit_end)	+ 1, "\t$percentID\t$hsp_evalue\n";		
		}
	}
}
	



