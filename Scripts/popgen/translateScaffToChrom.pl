#!/usr/bin/perl

# NOTES
# This is a merge between the scripts findPosOnfAlb15Chrom.pl and findPosOnFicAlb1.5ChromNCBI.pl
# Version 1.4, Apr 2015 (original written in Dec 2012)
# NEW v1.4: Fixed that any comment lines (starting with \"#\") is printed unchanged.
# NEW v1.3 (Sep 2013):  Added options to include unassigned scaffolds, and 
# choose input format gtf (set the position columns automatically)
# NEW v1.2: BUGFIXED! Before 12 Sept 2013, all positions on reverse scaffolds
# was one base too small (pos 1000 should be 1001, etc).    
# NEW v1.1 (Feb 2013): THIS VERSION IS UPDATED TO INCLUDE THE CORRECT SCAFFOLD N00498 
# INSTEAD OF THE WRONG N00546.


my $usage = "
# # # # # #
# translateScaffToChrom.pl
# Version 2, Feb 2018
# Author: Linnéa Smeds 
#

# NOTE 1: If your scaffolds are fAlb15, you get fAlb15 chromosomes, and if they 
# are FicAlb1.5, you get FicAlb1.5 (NCBI) chromosomes. To translate between  
# scaffold versions, use translateScaffToScaff.pl
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
# 1 from all start positions again. Also the orientation column has to be 
# fized manually, see below!

# NOTE 4: By default, the script cannot handle orientation columns properly
# (so any + or - will be left unchanged even if the scaffold lies reverse on
# the chromosome). Using flag format=gtf makes it work for gtf files!!

# ===========================================================================
# Takes a list with positions and finds the location on the concatenated 
# chromosome, given a certain gap size between chromosomes (and what scaffold
# that should be included - all or only the anchored or linked ones).
# 
# USAGE:
# perl findPosOnfAlb15Chrom.pl -in=file [-posCol=\"2,3\" -scafCol=1, -gap=5000 
#	-level=all|strict -out=newfile -h]

# All input flags: [default in brackets]
-in=file\t Any kind of table with scaffolds and positions (in version fAlb15!)
-name=string\t Name of assembly (fAlb15 or FicAlb1.5). [fAlb15]
-posCol=\"n,n\"\t Column or columns with positions to translate. separated by
\t\t a comma \",\". [2]
-scafCol=n\t The column that contains the scaffold names. [first column]
-gap=n\t\t How many Ns that are counted between the scaffolds [5000]
-level=all|strict What scaffolds to use for the chromosomes, strict means
\t\t  only scaffolds where we are absolutely sure or position/direction, all
\t\t  means all possible scaffolds.
-un\t\t If set, regions in unassigned scaffolds will be printed to output,
\t\t unchanged. [not used]
-format=gtf\t When used, scaffold and position columns are chosen 
\t\t automatically according to the format standard, and the strand (col 7)
\t\t is reversed if the scaffold has a negative orientation. [not used]
-version=n\t Linkage map version 20130221 or 20140121 [20140121, but 20130221 
\t\t  for FicAlb1.5]
-out=file\t Name of output file. [add suffix \"chrompos\" to infile]
-h\t\t Print this help message.
#============================================================================\n";

use strict;
use warnings;
use Getopt::Long;

my ($in,$name,$posCol,$scafCol,$level,$unassigned,$format,$version,$out,$help);
my $gapsize = 5000;	#Default value (updated this 20130227 because gap=0 was
				#interpreted as "undefined" when testing unless($gapsize)
my $liftOverDir="/proj/sllstore2017033/repos/assembly/liftOver";	# Updated this 20210705 after repos move
GetOptions(
  	"in=s" => \$in,
 	"name=s" => \$name,
  	"posCol=s" => \$posCol,
  	"scafCol=s" => \$scafCol,
	"gap=i" => \$gapsize,
	"level=s" => \$level,
	"un" => \$unassigned,
	"format=s" => \$format,
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
if($posCol) {
	if($format) {
		if($format eq "gtf"  && $posCol ne "4,5") {
			die "ERROR: Parameter conflict: -gtf is not compatible with -posCol=$posCol.\n\tLeave -posCol out and it will be set to \"4,5\"\n";
		}
		# BED IS NOT WORKING PROPERLY
	#	if($format eq "bed" && $posCol ne "2,3") {
	#		die "ERROR: Parameter conflict: -bed is not compatible with -posCol=$posCol.\n\tLeave -posCol out and it will be set to \"2,3\"\n";
	#	}
	}
}
else {
	$posCol = 2;
	
}	
if($scafCol) {
	if($format) {
		if($format eq "gtf" && $scafCol != 1) {
			die "ERROR: Parameter conflict: -gtf is not compatible with -scafCol=$scafCol.\n\tLeave -scafCol out and it will be set to 1\n";
		}
		#  BED IS NOT WORKING PROPERLY
	#	if($format eq "bed" && $scafCol != 1) {
	#		die "ERROR: Parameter conflict: -bed is not compatible with -scafCol=$scafCol.\n\tLeave -scafCol out and it will be set to 1\n";
	#	}
	}
	$scafCol=$scafCol-1;
	
}
else {
	$scafCol = 0;
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
		die "ERROR: Unknown linkage map version $version, possible versions are 20130221 or 20140121!\n Note that version flag is fixed for translations from NCBI version FicAlb1.5";
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


if($unassigned){
	$unassigned="yes";
}
else {
	$unassigned="no";
}
if($format) {
	if($format eq "gtf") {
		$posCol = "4,5";
	}
	#  BED IS NOT WORKING PROPERLY
	#elsif($format eq "bed") {
	#	$posCol = "2,3";
	#}
	else {
		die "ERROR: Unknown format $format. Known formats are gtf.\n".$usage."\n";
	}
}
else {
	$format="no";
}		 

unless($out) {
	$out= $in.".chrompos";
}
unless(-e $in) {
	die "ERROR: File $in doesn't exist!\n";
}

my $readableCol=$scafCol+1;
print "--------------------------------------------\n";
print "Translating $name scaffolds to $name chromosomes...\n";
print "Will look for scaffold name in column $readableCol and positions in column[s] $posCol.\n";
print "Parameters are gapsize=$gapsize, level=$level, version=$version.\n";
if($unassigned eq "yes") {
	print "Unassigned scaffolds will be printed unchanged\n";
}
unless($format eq "no") {
	print "Format flag is set to $format\n";
}
print "Translation table: $choice\n";
print "Output is printed to file $out.\n";
print "--------------------------------------------\n";



# ----------------------------------------------------------------------------------------
# Make hash of the scaffolds
my %scaffs = ();
open(IN, $choice);
my $prev = "";
my $tmpstart=0;
while(<IN>) {
	my  @tab = split(/\s+/, $_); # Columns are chr, scaf, len, dir, [type, color, link]

	$scaffs{$tab[1]}{'chr'}=$tab[0];
	$scaffs{$tab[1]}{'len'}=$tab[2];
	$scaffs{$tab[1]}{'dir'}=$tab[3];
	
	if($prev ne $tab[0]) {
		$tmpstart = 0;
	}
	$scaffs{$tab[1]}{'offset'}=$tmpstart;

	$prev=$tab[0];
	$tmpstart+=$tab[2];
	$tmpstart+=$gapsize;
}


# Open outfile 
open(OUT, ">$out");


# Go through the infile, change the name and position columns
# and print the line with the new scaffold name and positions 
open(IN, $in);
my $cnt=0;
while(<IN>) {
	if(/^#/){		#ADDED 2015-04-10
		print OUT;
	}
	else {
		chomp($_);
		my @tab = split(/\t+/, $_);

		my @col = split(/,/, $posCol);
		my $scaff = $tab[$scafCol];
		my $badFlag = "off";

		if(defined $scaffs{$scaff}) {

			my $chrom = $scaffs{$scaff}{'chr'};

			foreach my $c (@col) {
				my $c=$c-1;	#index is real col_no-1

				# CALCULATE THE NEW POSITION BASED ON THE OFFSET AND THE DIRECTION!
				my $newpos;
				if($tab[$c]>$scaffs{$scaff}{'len'} || $tab[$c]<0) {	# first check that scaffold is within scaffold range!
					print STDERR "WARNING: ".$tab[$c]." lies outside the $scaff range!! Skip line..\n";
					$badFlag="on";
				}
				else {
					if($scaffs{$scaff}{'dir'} eq "+") {
						$newpos = $scaffs{$scaff}{'offset'}+$tab[$c];		
					}
					elsif($scaffs{$scaff}{'dir'} eq "-") {
						$newpos = $scaffs{$scaff}{'offset'}+$scaffs{$scaff}{'len'}-$tab[$c]+1;		#BUGFIXED 20130912 (added +1)	
					}
					else {
						print "WARNING: $scaff has no assigned direction, assuming +!\n";
						$newpos = $scaffs{$scaff}{'offset'}+$tab[$c];		
					}	
					$tab[$c]=$newpos;
				}

			}
			unless ($badFlag eq "on") {
				$tab[$scafCol]=$scaffs{$scaff}{'chr'};
				# Switch start and stop on negative scaffolds for gtf files
				if($format eq "gtf" && $scaffs{$scaff}{'dir'} eq "-") {
					#print "DEBUG: format is gtf and scaffold orientation is negative!!\n";
					my $tempstart=$tab[4];	
					$tab[4]=$tab[3];
					$tab[3]=$tempstart;
					if($tab[6] eq "+") {
						$tab[6]="-";
					}
					elsif($tab[6] eq "-") {
						$tab[6]="+";
					}
				}
				print OUT join("\t", @tab) . "\n";
			}
				
		}
		else {
			if($unassigned eq "yes") {
				print OUT join("\t", @tab) . "\n";	
			}
			else {
				print $scaff." is not anchored to a chromosome, skipping line...\n";	
			}	
		}
	}
	$cnt++;
}
close(IN);

print "--------------------------------------------\n";
print "Done!\n";
print "Processed $cnt lines\n";
