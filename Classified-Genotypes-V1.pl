#!/usr/bin/perl -l
#This script extracts the count of possible genotypes extracted from GATK by VariantsToTable into different files
#Nima Rafati 20140218 (nima.rafati@imbim.uu.se)
#use warnings;
my $fileGT=$ARGV[0];
my @lineArr=();
my $usage = "./*.pl GT.file\n";
my $cntr=0;
my $switchRead=1;
my $ref="";
my %gtHash=();
my @alleleArr=();
my $alleleHash=();
my $ref="";
my $alt="";
my @outArr=();
my $out="";
open(inF1,$fileGT) || die print $usage;
while(<inF1>)
{
	$cntr++;
	@lineArr=split("\t",$_);
	$chr=$lineArr[0];
	$pos=$lineArr[1];
	$ref=$lineArr[2];
	$alt=$lineArr[3];
	$switchRead=1;
	if ($cntr>=1 && $lineArr[0] ne "CHROM")
	{
		chomp($_);
		@lineArr=split("\t",$_);
		for(my $i=4;$i<scalar(@lineArr);$i++)
		{
			@alleleArr=split("\/",$lineArr[$i]);
			#Hom Ref
			if ($alleleArr[0] eq $ref && $alleleArr[1] eq $ref)
			{
				$out.="\tHom.Ref";
			}
			#Het reciprocal (A/T or T/A)
			elsif($alleleArr[0] eq $ref && $alleleArr[1] eq $alt)
			{
				$out.="\tHet";
			}
			elsif ($alleleArr[0] eq $alt && $alleleArr[1] eq $ref)
                        {
				$out.="\tHet";
			}
			#Hom Alt
                        elsif ($alleleArr[0] eq $alt && $alleleArr[1] eq $alt)
                        {
				$out.="\tHom.Alt";
			}
			#No call (Deletion? =./.)
			elsif($alleleArr[0] eq "." && $alleleArr[1] eq ".")
			{
				$out.="\tDel";
			}

		}
		print "$chr\t$pos\t$ref\t$alt$out";
		$out="";
	}
}
close inF1;
