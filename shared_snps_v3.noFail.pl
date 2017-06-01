#!/usr/bin/perl


unless (@ARGV) {
	print "\ninput: <hybird_ag_parent.pileup> <parent_ag_parent.pileup> <out.txt>\n\n";
	exit;
}


open(HYBRID, $ARGV[0]) || die "can't find file $ARGV[0]\n";
open(OUT, ">$ARGV[2]");

###get number of contigs (chrom arms, etc.)
print "checking number of contigs...\n";
while ($throw = <HYBRID>) {
	$throw =~ /^(.+?)\t\d/;
	$contig = $1;
	if ( exists ($contig_hash{$contig}) ) {
	} else {
		$contig_hash{$contig} = "";
		print $contig."\n";
	}
	
}
close HYBRID;

foreach $key (keys %contig_hash) {	#$key is an array [0]=key1 [1]=key2
	%hashish = ();
	$mini_counter = 0;	
	$sum_h_cvg = 0;
	$hyb_mm_count = 0;
	$sum_p_cvg = 0;
	$par_mm_count = 0;	


	open(HYBRID, $ARGV[0]);
	###load hybrid hash
	print "###########\nloading $key hybrid hash...\n";
	while ($hybrid_line = <HYBRID>) {
		chomp $hybrid_line;
		if ($hybrid_line =~ /^$key/) {

			@hybrid_line_array = split('\t', $hybrid_line);
			$contig = $hybrid_line_array[0];
			$position = $hybrid_line_array[1];
			$ref_base = $hybrid_line_array[2];
			$cvg = $hybrid_line_array[3];	
			$read_base = $hybrid_line_array[4];	
			$quality = $hybrid_line_array[5];
			$sum_h_cvg += $cvg;
			$mini_counter++;		
			if (($cvg > 1) && ($read_base !~ /[+-]/) && ($read_base =~ /[atcg]/) && ($ref_base !~ /N/) ) {
				&snp;
				if ( length($read_base) > 1 ) {
					$hyb_ref_base = $ref_base;
					$hyb_read_base = $read_base;
					$hyb_cvg = $cvg;
					$hyb_mismatch = $mismatch;
					if ($mismatch =~ /[atcg]/i) {
						$hyb_mm_count++;
					}$hyb_mm_count++;
					if ($hyb_mismatch ne "null") {
						$hashish{$contig}{$position} = "$hyb_ref_base\t$hyb_read_base\t$hyb_cvg\t$hyb_mismatch";
					}		
				}
			}	
			if (($position % 1000000) == 0) {
				print "$key\t$position\n";
			}			
		}
		$hyb_ref_base = ();
		$hyb_read_base = "null";
		$hyb_cvg = ();
		$hyb_mismatch = "null";					

	}
	print "number of hybrid mismatches: $hyb_mm_count out of $mini_counter positions\n";
	if ($hyb_mm_count < 1) {
		print "\tare there only Ns in the reference base column in pileup?\n";
	}
	print "avg cvg hybrid: ".($sum_h_cvg / $mini_counter)."\n";
	$mini_counter = 0;


	close HYBRID;		# closed hybrid file, but still in foreach
##### Parent	
	print "correlating to parent mismatches...\n";
	open(PARENT, $ARGV[1]) || die "can't find file $ARGV[1]\n";
	while ($parent_line = <PARENT>) {
		chomp $parent_line;
		if ($parent_line =~ /^$key/) {		
			@parent_line_array = split('\t', $parent_line);
			$contig = $parent_line_array[0];
			$position = $parent_line_array[1];
			$ref_base = $parent_line_array[2];
			$cvg = $parent_line_array[3];	
			$read_base = $parent_line_array[4];	
			$quality = $parent_line_array[5];		
			$sum_p_cvg += $cvg;
			$mini_counter++;		
			if (($position % 1000000) == 0) {
				print "$key\t$position\n";
			}		
		
			if ( ($cvg > 1) && ($read_base !~ /[+-]/) && ($read_base =~ /[atcg]/i) && ($ref_base !~ /N/)) {
			
				&snp;
				if (length($read_base) > 1) {
					if ($mismatch =~ /[atcg]/i) {
						$par_mm_count++;
					}
					$p_ref_base = $ref_base;
					$p_read_base = $read_base;
					$p_cvg = $cvg;
					$p_mismatch = $mismatch;	
					if ($p_mismatch =~ /[ACTG]/) {
						if (exists ($hashish{$contig}{$position}) ) {
							
							@hyb = split('\t', $hashish{$contig}{$position});
							$hyb_ref_base = $hyb[0];
							$hyb_mismatch = $hyb[3];
			
							if ( ($p_ref_base eq $hyb_ref_base) && ($p_mismatch eq $hyb_mismatch) ) {
								print OUT "$contig\t$position\t1\n";
							} else {
								print OUT "$contig\t$position\t0\n";					
							}
						} else {
								print OUT "$contig\t$position\t0\n";
	
						}		
					}
				}
					
			}
		
		}
	}
	close PARENT;
	
	print "number of parent mismatches: $par_mm_count out of $mini_counter positions\n";
#	print "avg cvg parent: ".($sum_p_cvg / $mini_counter)."\n";
	$mini_counter = 0;

}

# notes:  we are going contig by contig, finding a variant in the hybirdt pileup and storing it in a hash. 
# Next we find that same contig in the parent vs other parent pileup to see if that snp is legit
# info is then printed to tthe outfile print OUT "$contig\t$position\t0\n";		




















#snp call. Plurality wins. # cdj says we should make this tweakable -- 
sub snp {
	$match = 0;
#	$read_base =~ s/[\^\$\$fn0-9\!\*]//gi;
	if ($read_base =~ /(\.|\,)/) {
		$match = $read_base =~ s/(\.|\,)//gi;
	}	
#	$match = $1;
	$read_base =~ s/[^agtc]//gi;	
	$mismatch_cvg = length($read_base);
	$counter = 0;
	$a = 0; $t = 0; $c = 0; $g = 0;
	if ( length($read_base) > 1) {
		until ($counter == length($read_base)) {
			$vote = substr($read_base, $counter, 1);
			if ( ($vote eq "A") || ($vote eq "a") ) {
				$a++;
			}
			if ( ($vote eq "T") || ($vote eq "t") ) {
				$t++;
			}
			if ( ($vote eq "C") || ($vote eq "c") ) {
				$c++;
			}
			if ( ($vote eq "G") || ($vote eq "g") ) {
				$g++;
			}

			$counter++;
		}
		#find winner
		$mismatch = "null";
		if ( ($a > $c) && ($a > $t) && ($a > $g) && ($a > $match) ) {
			$Aratio = $a / ($mismatch_cvg + $match);
			$mismatch = "A";
		}
		if ( ($t > $a) && ($t > $c) && ($t > $g) && ($t > $match) ) {
			$Tratio = $t / ($mismatch_cvg + $match);
			$mismatch = "T";
		}
		if ( ($c > $a) && ($c > $t) && ($c > $g) && ($c > $match) ) {
			$Cratio = $c / ($mismatch_cvg + $match);
			$mismatch = "C";
		}
		if ( ($g > $a) && ($g > $t) && ($g > $c) && ($g > $match) ) {
			$Gratio = $g / ($mismatch_cvg + $match);
			$mismatch = "G";
		}	
	} else {
		$mismatch = "null";
	}
}

