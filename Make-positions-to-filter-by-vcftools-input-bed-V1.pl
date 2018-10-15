#!/usr/bin/perl
#This script gets a bed file as input and reporduce these positions in one entry line
#For example
#chr1	5	8 (bed input)
#chr	6
#chr	7
#chr	8
#Nima Rafati 2014-02-13
my $usage="./*.pl file.bed\n";
open(inF1,$ARGV[0]) || die print $usage;
while(<inF1>)
{
	chomp($_);
	@lineArr=split("\t",$_);
	$chr=$lineArr[0];
	$start=$lineArr[1];
	$end=$lineArr[2];
	for($start;$start<$end+1;$start++)
	{
		print "$chr\t$start\n";
	}
}
close inF1;

