#!/usr/bin/env perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: $0 bed_file ChromosomeSummary  window_size step_size\n\n";

my $bedfile = shift or die ($usage);
my $chromFile = shift or die ($usage);
my $window_size = shift or die ($usage);
my $step_size = shift or die ($usage);


my %chromHoA;

my @bed_lines = file_to_array($bedfile);
shift @bed_lines;

foreach (@bed_lines){
	chomp $_;
	my @sl = split (/\t/, $_);
	$chromHoA{$sl[0]}[$sl[6]] = 1;
}

my @chromLines = file_to_array($chromFile);
shift @chromLines;

print "Chromosome\tWindowPosition\tGenes_per_kb\n";

foreach (@chromLines){
	chomp $_;
	my @sl = split (/\t/, $_);
	my $length = $sl[2];
	
	
	my $window_start = 1;
	
	while ($window_start + $window_size - 1 <= $length){
		my $gene_count = 0;
		for (my $i = $window_start; $i < $window_start + $window_size; ++$i){
			if (exists ($chromHoA{$sl[0]}[$i])){
				++$gene_count;
			}
		}
		
		print "$sl[0]\t", $window_start + $window_size/2 - 1, "\t", 1000 * $gene_count/$window_size, "\n";
		
		$window_start += $step_size;
	}
	
}