#!/usr/bin/env perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: $0 blastSummary ChromosomeSummary E-value_cutoff window_size step_size\n\n";

my $blastFile = shift or die ($usage);
my $chromFile = shift or die ($usage);
my $evalueCutoff = shift or die ($usage);
my $window_size = shift or die ($usage);
my $step_size = shift or die ($usage);


my %chromHoA;

my @blast_lines = file_to_array($blastFile);
shift @blast_lines;

foreach (@blast_lines){
	chomp $_;
	my @sl = split (/\t/, $_);
	for (my $i = $sl[2]; $i <= $sl[3]; ++$i){
		$chromHoA{$sl[0]}[$i] = 1;
	}
}

my @chromLines = file_to_array($chromFile);
shift @chromLines;

print "Chromosome\tWindowPosition\tCoverage\n";

foreach (@chromLines){
	chomp $_;
	my @sl = split (/\t/, $_);
	my $length = $sl[2];
	
	
	my $window_start = 1;
	
	while ($window_start + $window_size - 1 <= $length){
		my $coverage = 0;
		for (my $i = $window_start; $i < $window_start + $window_size; ++$i){
			if (exists ($chromHoA{$sl[0]}[$i])){
				++$coverage;
			}
		}
		
		print "$sl[0]\t", $window_start + $window_size/2 - 1, "\t", $coverage/$window_size, "\n";
		
		$window_start += $step_size;
	}
	
}