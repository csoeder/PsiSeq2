#!/usr/bin/perl
#12-14-2010	Eric J Earley	UNC @ Chapel Hill	C. Jones Lab

#VERSION 3 notes
#tracks the size of A block and its coordinates

#VERSION 2 notes
#tracks the proportion of As and Bs by position:
#output: pos: X \t #As	\t	#Bs
#


#VERSION 1.0 notes
#simulate introgression of haplotype blocks (size and frequency) in a backcross of X generations with crossover and no selection
#	constant population size
#		randomly pick X individuals each generation (X = population size)
#		no bias for 1 block frequency
#	1 crossover each generation for each chromosome
#	

$POP_SIZE= $ARGV[0];
$CHROM_SIZE = $ARGV[1];
$GEN = $ARGV[2];
$OUT = $ARGV[3];
$BLOCKOUT = $ARGV[4];
$QUALITY = $ARGV[5];
open (BLOCKOUT, ">$BLOCKOUT");
open(OUT, ">$OUT");
open(QUALITY, ">$QUALITY");
unless (@ARGV){
	print "input: <population size><chromosome length><number of generations><position.out><blocks.out><quality.test>\n\n";
	exit;
}

print OUT "For each position along chromosome, how often do we see an A or B block after X generations?\n";
print BLOCKOUT "For each individual, where does the A block start, and how long is it?\npop size: $POP_SIZE; chrom length: $CHROM_SIZE; generations: $GEN\n\n\nIndividual\tA_block_start\tA_block_length\n";
print QUALITY "Does simulation meet expected half-ing of population wide A identity?\n\ngeneration\tA_block size\tobs\texp\n";



#create template AAAAAAAAAAAAAA chromosomes
until ($chrom_counter == $CHROM_SIZE) {
	$chrom .= "A";
	$chrom_counter ++;
}
$chrom_counter = ();

# create recurring backcross parent BBBBBBBBBBBBBBBB
until ($chrom_counter == $CHROM_SIZE) {
	$bc .= "B";
	$chrom_counter++;
}
$chrom_counter = ();

#### CREATE PARENTAL POPULATION OF  AAAAAAAAAA CHROMOSOMES
until ($pop_counter == $POP_SIZE) {
	push(@pop_array, $chrom);
	$pop_counter ++;
}
print  OUT "population size: $POP_SIZE\nchrom length: $CHROM_SIZE\ngenerations: $GEN\n";
print OUT "pos\tA\tB\n";
#print BLOCKOUT "population size: $POP_SIZE\nchrom length: $CHROM_SIZE\ngenerations: $GEN\n";
#print BLOCKOUT "individual\tA_block_start\tlength\n";

######## CROSSOVER
#### 2 daughters made from each parent (@pop_array)
$gen_counter = 0;
until ($gen_counter == $GEN) {
	$gen_counter++;	

	print "generation $gen_counter\n";



	foreach (@pop_array) {
		$random = rand($CHROM_SIZE);
		$pos = sprintf("%.f", $random);
			#EG: pos=1
			#      0123456789
			#       P
			#D1 -> ABBBBBBBBB...
			#D2 -> BAAAAAAAAA...
		#5'
		$daughter1 = substr($_,0,$pos);	#  A[POS->]_______
		$bc_piece = substr($bc, $pos);	#  _[POS->]BBBBB	
		$daughter1 .= $bc_piece;	#pasted together
#print BLOCKOUT "daughter1: $daughter1\n";
	
		#3'
		$daughter2 = substr($_,$pos);	#take 3' chunk of parent to make daughter2
		$bc_piece = substr($bc,0,$pos);	#5' BBBBBB chunk
		$daughter2 = $bc_piece.$daughter2;	#pasted
#print BLOCKOUT "daughter2: $daughter2\n";
		#length of daughters verified
#print BLOCKOUT "pos: $pos\n";
		#pick random daughter	
		$pick = (rand());
#print "pick: $pick\t";
		if ($pick < 0.5) {
			push(@daughter_array,$daughter1);
		}
		if ($pick > 0.5) {
			push (@daughter_array,$daughter2);
		}
	}

	@pop_array = @daughter_array;
	
	@daughter_array = ();

############	go position by position within chrom and count # of As and Bs
	$zero_count = 0;
	$one_count = 0;
	$st_count = 0;
	$individual = 0;
	$A_start = "-";



	foreach(@pop_array) {
		$haplotype = $_;
#print BLOCKOUT $haplotype."\n";
		$individual ++;
		$len = length ($haplotype);
		$prev_base = "B";	#in case 1st position is an A
		until ($st_count == $len) {	#begin counting blocks here
			$base = substr($haplotype, $st_count,1);
			if ($base eq "B") {
				$zero_count++;
				$prev_base = "B";
				if ($gen_counter == $GEN) {
					$pos_hash{$st_count} .= "B";
				}
			}
			if ($base eq "A") {
				$one_count++;
				if ($gen_counter == $GEN) {
					$pos_hash{$st_count} .= "A";
					if ($prev_base eq "B") {
						$A_start = $st_count;
					}
				}
				$prev_base = "A";
			}
			

			$st_count++;
		}

		if ($gen_counter == $GEN) {
			print BLOCKOUT "$individual\t$A_start\t$one_count\n";
#			print BLOCKOUT $haplotype."\n";
#			
#			$sum_one_count += $one_count;	#gather cumulative block sizes for averaging later
			
		}
		$gen_A_total += $one_count;
	

		$zero_count = 0;
		$one_count = 0;
		$st_count = 0;
		$one_block = 0;
		%block_length_hash = ();
		$A_start = "-";
		$one_block_length = 0;
	
	}
#	$average_block_size = $sum_one_count / $pop_size;
	
	
$gen_A_prop = $gen_A_total / (($POP_SIZE * $CHROM_SIZE)*2); #x2 because everyone has 1 whole B chrom
$expected = (0.5)**($gen_counter + 1);
print QUALITY "$gen_counter\t$gen_A_total\tobs=$gen_A_prop\texp=$expected\n";
$gen_A_total = ();
$gen_A_prop = ();
	
	
	
	foreach $key (sort {$a<=>$b} keys %pos_hash) {
		$A = $pos_hash{$key} =~ s/A/A/g;
		if ($A !~ /\d/) {
			$A = 0;
		}
		$B = $pos_hash{$key} =~ s/B/B/g;
		if ($B !~ /\d/) {
			$B = 0;
		}
		print OUT "$key\t$A\t$B\n";
	}
	
}

