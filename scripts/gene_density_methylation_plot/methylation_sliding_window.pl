#!/usr/bin/env perl

use strict;
use warnings;
use sloan;

my $usage = "\nUSAGE: $0 bed_file ChromosomeSummary  window_size step_size\n\n";

my $bedfile = shift or die ($usage);
my $chromFile = shift or die ($usage);
my $window_size = shift or die ($usage);
my $step_size = shift or die ($usage);


my %chromHoA_meth;
my %posHoA;

my @bed_lines = file_to_array($bedfile);

foreach (@bed_lines){
	chomp $_;
	my @sl = split (/\t/, $_);
	$posHoA{$sl[0]}[$sl[1]]= 1;
	$chromHoA_meth{$sl[0]}[$sl[1]] = $sl[10];
}

my @chromLines = file_to_array($chromFile);
shift @chromLines;

print "Chromosome\tWindowPosition\tMethProp\n";

foreach (@chromLines){
	chomp $_;
	my @sl = split (/\t/, $_);
	my $length = $sl[2];
	
	
	my $window_start = 1;
	
	while ($window_start + $window_size - 1 <= $length){
		my $site_count = 0;
		my $methsum = 0;
		for (my $i = $window_start; $i < $window_start + $window_size; ++$i){
			if (exists ($posHoA{$sl[0]}[$i])){
				++$site_count;
				$methsum += $chromHoA_meth{$sl[0]}[$i];
			}
		}
		
		if ($site_count){
			print "$sl[0]\t", $window_start + $window_size/2 - 1, "\t", $methsum/$site_count, "\n";		
		}
		$window_start += $step_size;
	}
}