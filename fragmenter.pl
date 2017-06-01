#!/usr/bin/perl
#written by W. Jeck 2007; modified by E. Earley 2011
#UNC at Chapel Hill, NC; Dept. Biology; Corbin Jones Lab

#Takes a genome file as the first argument (FASTA format) and outputs
#a fasta with many reads (file is the second argument). You also specify 
#read length, coverage, and error rate desired.

#example call
#perl fragmenter.pl reference.fa output.fa read-length coverage artificial-error-frequency


$genomeFile = $ARGV[0];
$outputFile = $ARGV[1];
$readLength = $ARGV[2];
$coverage = $ARGV[3];
$singleBaseErrorFrequency = $ARGV[4];

$singleBaseErrorFrequency = $singleBaseErrorFrequency * 1.3333;

unless (@ARGV) {
	print "Genome Out Readlen coverage\n";
	exit;
}

srand(time * 1000);

###reads in the genome

open(FILE, "$genomeFile") || die ("error opening genome file $!\n");
$junkline = <FILE>;

$genome = "";
while ($line = <FILE>) {
	chomp($line);
	$genome .= $line;
}
$genome =~ s/\r//g;

$size = length($genome);
$readspergenome = $size/$readLength;
$numReads = $coverage * $readspergenome;

print "$numReads fragments for $coverage X coverage\n$size nucleotides in genome\n";

#open output file
open(OUT, ">$outputFile");

for ($i = 1; $i <= $numReads; $i++) {
	print "$i\n";
	$startLocation = rand();
	$direction = rand();

	$startLocation = $startLocation * length($genome);

	if ($direction < .5) {
		#get the reverse compliment of the sequence
		$startLocation -= $readLength;
		if ($startLocation < 0) {
			$startLocation = 0;
		}
		$thisRead = substr($genome, $startLocation, $readLength);
		
		#reverse compliment $thisRead
		
		$thisRead = reverse($thisRead);
		$thisRead = uc($thisRead);
		$thisRead =~ tr/ACTG/tgac/;
		$thisRead = uc($thisRead);
		
	} else {
		#get forward strand
		$thisRead = substr($genome, $startLocation, $readLength);	
	}
	
	$lengthThisRead = length($thisRead);
	
	$doneRead = "";
	for ($j = 0; $j < $lengthThisRead; $j++) {
		$changeBaseRand = rand();
		if ($changeBaseRand < $singleBaseErrorFrequency) {
			#change the base
			$whatBase = rand();
			if ($whatBase < .25) {
				$nowBase = "A";
			} elsif ($whatBase < .50) {
				$nowBase = "T";
			} elsif ($whatBase < .75) {
				$nowBase = "G";
			} else {
				$nowBase = "C";
			}
		} else {
			$nowBase = substr($thisRead, $j, 1);
		}
		$doneRead .= $nowBase;
	}
	
	
	##If you want to insert errors, just randomly act on $thisRead at this point
	#
	##
	$holder = $doneRead;
	if ($doneREAD =~ s/[^ATCGNatcgn]+?//ig) {
		
	} else {
		print OUT ">READ_$i\n$holder\n";
	}
}
