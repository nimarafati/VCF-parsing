#!/usr/bin/perl
#This script extract genotypes and generate haplotypes
#Nima Rafati 20140207
use strict;
use Getopt::Long;
#my $vcfFile=$vcf_file;
my $startSwitch=0;
my @lineArr=();
my $chr="";
my $position=0;
my $ref="";
my $var="";
my $qual="";
my $numSamples=0;
my $format="";
my @sampleArr=();
my $genotype="";
my @formatArr=();
my @gtArr=();
my $cntr=0;
my $l=0;
my $cntr=0;
my $gtPos=0;
my $g1="Hap1";
my $g2="Hap2";
my $size=0;
my %hapHash=();
my $s=0;
my $cntrHap=0;
my $helpFlag="";
my $formatFlag="";
my $vcfFile="";
my $phasedFlag="F";
my $fasta=1;
my $isPhase="F";
my $hap1="";
my $hap2="";
my $countme=0;
my $switchMiss="F";
my $sampleDP=10;
my $dpPos="NA";
my $refCnt=0;my $altVCnt=0;my $totalDP=0;
my $goToMiss=0;
my $fasta_tandem="";
my $usage="./extract-GT-V1.pl -format Fasta/Tab -vcf SNP.vcf
Nima Rafati 20140209 (nima.rafati\@imbim.uu.se)
-format tab/fasta 
		the defualt output is fasta file for each haplotype of each individuals:
		>ind1-1
		hap-1
		>ind1-2
		hap-2
		and the tab format makes a matrix for haplotypes:
		ind1-1	h	a	p	-	1
		ind1-2	h	a	p	-	2
-vcf	input vcf file
-phased	T/F (defualt is F meaning the vcf is not phased)
-missingGT	T/F (defualt is F meaning does not report missing GT otherwise\nindividuals/populations with missing GT will be represented by -1)
-sampleDP	10 (set depth value to read each position default is 10)
-fasta_tandem  T/F This option will concatenate the genotypes and each individual will have only one fasta sequence. PLEASE MAKE SURE TO USE IT TOGETHER WITH \"-format fasta\" (default F)\n";
&GetOptions ('h'=>\$helpFlag,
			 'format=s' =>\$formatFlag,
			 'vcf=s' =>\$vcfFile,
			 'phased=s' =>\$phasedFlag,
			 'missingGT=s' =>\$switchMiss,
			 'sampleDP=s' =>\$sampleDP,
			'fasta_tandem=s' =>\$fasta_tandem);
if($helpFlag)
{
	die print $usage;
}

unless($vcfFile)
{
	die print $usage;
}

if ($formatFlag eq "" || $formatFlag eq "fasta")
{
	$fasta=1;
}
elsif($formatFlag eq "tab")
{
	$fasta=0;
}
if($phasedFlag eq "T")
{
	my $isPhase="T";	
}
else
{
	my $isPhase="F";
}
if($switchMiss eq "T")
{
	$switchMiss=1;
}
else
{
	$switchMiss=0;
}
if($sampleDP eq "")
{
	$sampleDP=10;
}
if($fasta_tandem eq "")
{
	my $fasta_tandem="F";
}
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	CG1	CG2	CG3	CG4	CG5	CG6
#Consensus-4	39	.	T	C	.	PASS	AR2=0;DR2=0;AF=0.267	GT:DS:GP	0|0:0.533:0.541,0.385,0.074	0|0:0.533:0.541,0.385,0.074	0|1:0.533:0.541,0.385,0.074	0|0:0.533:0.541,0.385,0.074	0|1:0.533:0.541,0.385,0.074	0|0:0.533:0.541,0.385,0.074
open(outF1,">just-called-position.bed");
open(inF1,$vcfFile) || die "$usage.\n";
while(<inF1>)
{
	$countme++;
	chomp($_);
	@lineArr=split("\t",$_);
	my $size=scalar(@lineArr);
	if($startSwitch == 1)
	{
		$chr=$lineArr[0];
		$position=$lineArr[1];
		$ref=$lineArr[3];
		$var=$lineArr[4];
		$qual=$lineArr[5];
		$format=$lineArr[8];
		##Get the GT position in format column
		@gtArr=split(":",$lineArr[8]);
		foreach $l (@gtArr)
		{
			if($l eq "GT") 
			{
				$gtPos=$cntr;
			}
			if($l eq "AD")
			{
				$dpPos=$cntr;
			}
			$cntr++;
		}
		##Retrieve the genotype based on the identified GT position in previous for loop
		#Var allele will be always in second haplotype
		for ($s=9;$s<$size;$s++)
		{
#			print "$position\t$lineArr[$s]";<STDIN>; 
						
			if ($lineArr[$s] ne "./.")
			{
			#	print $lineArr[$s],"\t",$s;<STDIN>;
				print outF1 "$chr\t",$position-1,"\t$position\t$chr.$position\n";
				@formatArr=split(":",$lineArr[$s]);
#				print "$lineArr[$s]\t$phasedFlag";<STDIN>;
				#Get genotypes
				if ($phasedFlag eq "T")
				{
					$formatArr[$gtPos]=~ m/(.*)\|(.*)/;
					if($1==0){$hap1=$ref;}elsif($1==1){$hap1=$var;}elsif($1 eq "."){$hap1="-";}
					if($2==0){$hap2=$ref;}elsif($2==1){$hap2=$var;}elsif($2 eq "."){$hap2="-";}
#					$hap1=$1; $hap2=$2;
#					print "HERE:$hap1\t$hap2";<STDIN>;
				}
				else
				{
					$formatArr[$gtPos]=~ m/(.*)\/(.*)/;
					if($1==0){$hap1=$ref;}elsif($1==1){$hap1=$var;}elsif($1 eq "."){$hap1="-";}
					if($2==0){$hap2=$ref;}elsif($2==1){$hap2=$var;}elsif($2 eq "."){$hap2="-";}
#				$hap1=$1; $hap2=$2;
#					print "$countme $hap1\t$hap2";<STDIN>;
				}
				#Get depth
#				print $dpPos;<STDIN>;
				if (length($formatArr[$dpPos])<=1)
				{
					$goToMiss=1;
				}
				else
				{
					$formatArr[$dpPos]=~ m/(.*),(.*)/;
					my $refCnt=$1;
					my $altVCnt=$2;
					my $totalDP=$1+$2;
					if($totalDP>=$sampleDP/2)
					{
						$goToMiss=0;
					}
					else
					{
						$goToMiss=1;
					}
				}
				if($dpPos eq "NA")
				{
					$goToMiss=0;
				}
				#Put haplotypes in hash
				if ($fasta==1 && $goToMiss==0)
				{
						$hapHash{$sampleArr[$s-9]}{$g1}=$hapHash{$sampleArr[$s-9]}{$g1}."$hap1";
						$hapHash{$sampleArr[$s-9]}{$g2}=$hapHash{$sampleArr[$s-9]}{$g2}."$hap2";
				}
				elsif($fasta==0 && $goToMiss==0)
				{
						$hapHash{$sampleArr[$s-9]}{$g1}=$hapHash{$sampleArr[$s-9]}{$g1}."$hap1\t";
						$hapHash{$sampleArr[$s-9]}{$g2}=$hapHash{$sampleArr[$s-9]}{$g2}."$hap2\t";
				}
				if($goToMiss==1)
				{
	                               if ($fasta==1)
        	                        {       
                	                        $hapHash{$sampleArr[$s-9]}{$g1}=$hapHash{$sampleArr[$s-9]}{$g1}."-";
                        	                $hapHash{$sampleArr[$s-9]}{$g2}=$hapHash{$sampleArr[$s-9]}{$g2}."-";
                                	}               
	                                else    
        	                        {       
                	                        $hapHash{$sampleArr[$s-9]}{$g1}=$hapHash{$sampleArr[$s-9]}{$g1}."-\t";
                        	                $hapHash{$sampleArr[$s-9]}{$g2}=$hapHash{$sampleArr[$s-9]}{$g2}."-\t";
                                	}
					$goToMiss=0;
				}
#print "GTGTGTGTGT\n$_\n$hapHash{$sampleArr[$s-9]}{$g1}\n$hapHash{$sampleArr[$s-9]}{$g2}\n";<STDIN>;
			}
			elsif($lineArr[$s] eq "./." && $switchMiss==1)
			{
				print outF1 "$chr\t",$position-1,"\t$position\t$chr.$position\n";
                                if ($fasta==1)
                                {
					$hapHash{$sampleArr[$s-9]}{$g1}=$hapHash{$sampleArr[$s-9]}{$g1}."-";
                                        $hapHash{$sampleArr[$s-9]}{$g2}=$hapHash{$sampleArr[$s-9]}{$g2}."-";
                                }
                                else
                                {
                                        $hapHash{$sampleArr[$s-9]}{$g1}=$hapHash{$sampleArr[$s-9]}{$g1}."-\t";
                                        $hapHash{$sampleArr[$s-9]}{$g2}=$hapHash{$sampleArr[$s-9]}{$g2}."-\t";
                                }
			}
		}
		$cntr=0;
	}
	if ($lineArr[0] eq "#CHROM")
	{
		$numSamples=$size-9;
		##Get sample IDs
		for (my $i=9;$i<$size;$i++)
		{	
			if($lineArr[$i] ne "")
			{
				push (@sampleArr,$lineArr[$i]);
				$cntrHap++;
			}
		}
		$startSwitch=1;
	}	
}
close inF1;
close outF1;
$cntrHap=0;
#Printing haplotypes
foreach my $k (sort keys %hapHash)
{
	foreach my $g (sort keys %{$hapHash{$k}})
	{
		$cntrHap++;
		if($fasta == 1 && $fasta_tandem eq "F")
		{
			print ">$k-$cntrHap\n$hapHash{$k}{$g}\n";
		}
		elsif($fasta == 0)
		{
			print "$k-$cntrHap\t$hapHash{$k}{$g}\n";
		}
	}
	if($fasta == 1 && $fasta_tandem == "T")
	{
		print ">$k\n$hapHash{$k}{$g1}$hapHash{$k}{$g2}\n";
	}
	$cntrHap=0;
}
