#!/usr/bin/perl

# NOTES
# version 2, Feb2018 (original written in Dec 2012)
# A merge between the scripts findPosOnfAlb15Chrom_BACKtranslation.pl and
# findPosOnfAlb15Chrom_BACKtranslationNCBI.pl - this version can handle both
# NCBI chromosomes and any kind of custom made fAlb15 scaffolds. 

my $usage = "
# # # # # #
# translateChromToScaff.pl
# version 2, Feb 2018
# Author: Linnéa Smeds
#   
# NOTE 1: If your chromsomes are fAlb15, you get fAlb15 scaffolds, and if they 
# are FicAlb1.5, you get FicAlb1.5 scaffolds. To translate between scaffold 
# versions, use translateScaffToScaff.pl
#
# NOTE 2: This version has a flag for changing the linkage map version. For fAlb15
# the default is 20140121, for FicAlb1.5 it's 20130221.
#
# NOTE 3: Please note that the script cannot handle cases where start and stop
# not both are either exclusive or inclusive (like BED-FORMAT where start is 
# inclusive and stop is not). This is because start and stop can switch place
# in case of minus-oriented scaffolds. If you want to translate BED files, 
# please add 1 to all your start positions before translating, then rearrange
# your output so that start is always smaller than stop, and finally subtract
# 1 from all start positions again. 
# 
# ===========================================================================
# DESCRIPTION:
# Takes a list with positions on the concatenated chromosome, and finds the 
# location on the original scaffold given a certain gap size between 
# scaffolds.
# For each position column, the script returns two columns - one for the scaffold
# and one for the new position.
# When translating from NCBI assembly FicAlb1.5, gapsize, version and level are
# fixed and any attempt to set these flags will be overruled.
# 
# USAGE:
# perl translateChromToScaff.pl -in=file [ -name=fAlb15 -posCol=\"2,3\" -chrCol=1, 
#	-gap=5000 -level=all|strict -vers=20130221 -out=newfile -h]
#
# PARAMETERS [Default in brackets]
-in=file\t Any kind of table with chromosomes and positions
-name=string\t Name of assembly (fAlb15 or FicAlb1.5). [fAlb15]
-posCol=\"n,n\"\t Column or columns with positions to translate. separated by
\t\t a comma \",\". [2]
-chrCol=n\t The column that contains the chromosome names. [1]
-gap=n\t\t How many Ns that are inserted between the scaffolds [5000]
-level=all|strict What scaffolds were used for the chromosomes, strict means
\t\t  only scaffolds where we are absolutely sure or position/
\t\t  direction, all means all possible scaffolds.
-version=n\t Linkage map version 20130221 or 20140121 [20140121]
-out=file\t Name of output file. [add suffix \"chrompos\" to infile]
-h\t\t Print this help message.
#============================================================================\n";

use strict;
use warnings;
use Getopt::Long;

my ($in,$name,$posCol,$chrCol,$level,$version,$out,$help);
my $gapsize = 5000;	#Default value (updated this 20130227 because gap=0 was
			#interpreted as "undefined" when testing unless($gapsize)
my $liftOverDir="/proj/sllstore2017033/repos/assembly/liftOver";	# Updated this 20210705 after repos move
GetOptions(
  	"in=s" => \$in,
	"name=s" => \$name,
   	"posCol=s" => \$posCol,
  	"chrCol=s" => \$chrCol,
	"gap=i" => \$gapsize,
	"level=s" => \$level,
	"version=s" => \$version,
	"out=s" => \$out,
	"h" => \$help);

#-------------------------------------------------------------------------------
#Checking input
print "\n";

if($help) {
	die $usage . "\n";
}
unless($in) {
	die "ERROR: No infile was given! \n" . $usage ."\n";
}
unless($posCol) {
	$posCol = 2;
}	
if($chrCol) {
	$chrCol=$chrCol-1;
}
else {
	$chrCol = 0;
}

if($level) {
	unless($level eq "strict" || $level eq "all"){
		die "ERROR: Flag -level must be either strict or all! \n".$usage."\n";
	}
}else {
	$level="strict";
}

if($version) {
	unless($version==20130221 || $version==20140121) {
		die "ERROR: Unknown linkage map version $version, possible versions are 20130221 or 20140121!\n Note that version flag is fixed set for translations from NCBI version FicAlb1.5";
	}
}
else {
	$version=20140121;
}

my $choice;
if($name) {
	if($name eq "fAlb15") {
		$choice="$liftOverDir/$name.chrom.$level.$version.txt";
	}
	elsif($name eq "FicAlb_1.5" || $name eq "NCBI" || $name eq "FicAlb1.5") {
		if($name eq "NCBI" || $name eq "FicAlb_1.5") {
			print "Assembly name is $name - I assume you mean \"FicAlb1.5\"!\n";
			$name="FicAlb1.5";
		}
		$choice="$liftOverDir/$name.chrom.txt";
		print "NOTE that FicAlb1.5 chromosome assembly has a fixed set of \nscaffolds and flags -gap, -level and version cannot be used!\n";
		$gapsize=5000;
		$version=20130221;
	}
	else {
		die "ERROR: unrecognized assembly name $name! Please use \"fAlb15\" or \"FicAlb1.5\"!\n";
	}
}
else {
	$name="fAlb15";
	$choice="$liftOverDir/$name.chrom.$level.$version.txt";
}


unless($out) {
	$out= $in.".scaffpos";
}
unless(-e $in) {
	die "ERROR: File $in doesn't exist!\n";
}

if($in =~ m/\.bed$/) {
	print "\nWARNING! The script cannot handle bed format properly - see usage!\n";
}


my $readableChrCol=$chrCol+1;

print "--------------------------------------------\n";
print "Translating $name chromosomes to $name scaffolds...\n";
print "Will look for chromosome name in column $readableChrCol and positions in column[s] $posCol.\n";
print "Parameters are gapsize=$gapsize, level=$level, version=$version.\n";
print "Translation table: $choice\n";
print "Output is printed to file $out.\n";
print "--------------------------------------------\n";



# ----------------------------------------------------------------------------------------
# 


# Make hash of the scaffolds
my %scaffs = ();
open(IN, $choice);
my $prev = "";
my $tmpstart=0;
while(<IN>) {
	my  @tab = split(/\s+/, $_); # Columns are chr, scaf, len, dir, [type, color, link]

	$scaffs{$tab[0]}{$tab[1]}{'len'}=$tab[2];
	$scaffs{$tab[0]}{$tab[1]}{'dir'}=$tab[3];
	
	if($prev ne $tab[0]) {
		$tmpstart = 0;
	}
	$scaffs{$tab[0]}{$tab[1]}{'offset'}=$tmpstart;

	$prev=$tab[0];
	$tmpstart+=$tab[2];
	$tmpstart+=$gapsize;
}


# Open outfile 
open(OUT, ">$out");


# Go through the infile, change the name and position columns
# and print the line with the new scaffold name and positions 
open(IN, $in);
my $warningCnt=0;
my $cnt=0;
while(<IN>) {
	chomp($_);
	my @tab = split(/\s+/, $_);

	my @col = split(/,/, $posCol);
	my $chrom = $tab[$chrCol];
	my $printFlag = "on";
	
#	print "STDIN: chrom is $chrom\n";
	my $warnFlag="off";
	foreach my $c (@col) {
		my $c=$c-1;	#index is real col_no-1
		
		my ($scaff, $pos) = ("-","-");

		foreach my $key (sort {$scaffs{$chrom}{$b}{'offset'} <=> $scaffs{$chrom}{$a}{'offset'}} keys %{$scaffs{$chrom}}) {
#			print "STDIN: in the for loop, keys is $key and offset is ".$scaffs{$chrom}{$key}{'offset'}."\n";
			if($scaffs{$chrom}{$key}{'offset'}<$tab[$c]){
				$scaff=$key;
				$pos = "";
				if($scaffs{$chrom}{$key}{'dir'} eq "+") {
					$pos = $tab[$c]-$scaffs{$chrom}{$key}{'offset'};
				}
				else {
					$pos = $scaffs{$chrom}{$key}{'len'}-($tab[$c]-$scaffs{$chrom}{$key}{'offset'})+1;
				}
				last;
			}
		}
		if($scaff eq "-") {
			print STDERR "WARNING: Can't find $chrom in the translation file - are you using the correct version?\n";
			$tab[$c]="UNKNOWN\tUNKNOWN";
			$warnFlag="on";
		}
		else {
			if($pos<1) {
				print STDERR "WARNING: ".$chrom.":".$tab[$c]." is outside the $scaff range! Changing $pos to 1\n";
				$pos=1;
				$warnFlag="on";
			}
			if($pos>$scaffs{$chrom}{$scaff}{'len'}) {
				print STDERR "WARNING: ".$chrom.":".$tab[$c]." is outside the $scaff range! Changing $pos to ".$scaffs{$chrom}{$scaff}{'len'}."\n";
				$pos=$scaffs{$chrom}{$scaff}{'len'};
				$warnFlag="on";			
			}
			$tab[$c]=$scaff."\t".$pos;
		}
	}
	if($warnFlag eq "on") {
		$warningCnt++;
	}
	print OUT join("\t", @tab) . "\n";
	$cnt++;
}
close(IN);


print "--------------------------------------------\n";
print "Done!\n";
print "Processed $cnt lines\n";
my $frac=$warningCnt/$cnt;
if($frac>0.1) {
	my $perc=int($frac*100);
	print "You recieved warnings for $warningCnt lines ($perc% of the input file)\n";
	print "A warning means a position lies outside of the scaffold range\n";
	print "(for example in gaps) but a high fraction of warnings might suggest\n";
	print "that the input isn't in the given format (=wrong assembly version).\n";
	print "Please make sure you are using the right settings!\n"; 
}
