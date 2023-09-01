#!/usr/bin/perl

use strict;
use warnings;
use Bio::SearchIO; 
use sloan;

my $usage = "\nUSAGE: $0 blastFile outputName Evalue_cutoff\n\n";

my $blastFile = shift or die ($usage);
my $output = shift or die ($usage);
my $evalue_cutoff = shift or die ($usage);

my $FH1 = open_output ("$output\.queries.txt");
my $FH2 = open_output ("$output\.hits.txt");


print $FH1 "SeqName\tLength\tBlastCoverage\n";
print $FH2 "SeqName\tLength\tBlastCoverage\n";


my $SearchIO_obj = new Bio::SearchIO(-format => 'blast', -file   => $blastFile);

my %repHash;

my %query_lengths;
my %hit_lengths;
my %query_cov;
my %hit_cov;
my %query_sitesHoH;
my %hit_sitesHoH;



while (my $result_obj = $SearchIO_obj->next_result){
	my $query_name = $result_obj->query_name;
	my $query_len = $result_obj->query_length;
	my $any_hits = 0;
	$query_lengths{$query_name} = $query_len;
	while (my $hit_obj = $result_obj->next_hit){
		my $hit_evalue = $hit_obj->significance;
		$hit_evalue <= $evalue_cutoff and $any_hits = 1;
		my $hit_name = $hit_obj->name;
		my $hit_len = $hit_obj->length;
		$hit_lengths{$hit_name} = $hit_len;
  		while (my $hsp_obj = $hit_obj->next_hsp){
			my $hsp_evalue = $hsp_obj->evalue;
			$hsp_evalue <= $evalue_cutoff or next;
			my $query_start = $hsp_obj->start('query');
			my $query_end = $hsp_obj->end('query');
			my $hit_start = $hsp_obj->start('hit');
			my $hit_end = $hsp_obj->end('hit');
			for (my $i = $query_start; $i <= $query_end; ++$i){
				unless (exists ($query_sitesHoH{$query_name}->{$i})){
					++$query_cov{$query_name};
					$query_sitesHoH{$query_name}->{$i} = 1;
				}
			}
			for (my $i = $hit_start; $i <= $hit_end; ++$i){
				unless (exists ($hit_sitesHoH{$hit_name}->{$i})){
					++$hit_cov{$hit_name};
					$hit_sitesHoH{$hit_name}->{$i} = 1;
				}
			}
		}
	}
	$any_hits or $query_cov{$query_name} = 0;
}
	

foreach (sort keys %query_lengths){
	print $FH1 "$_\t$query_lengths{$_}\t$query_cov{$_}\n";
}

foreach (sort keys %hit_lengths){
	print $FH2 "$_\t$hit_lengths{$_}\t$hit_cov{$_}\n";
}



