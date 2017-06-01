#!/usr/bin/perl
#written by EJ Earley (2011)
#published within:
#	Earley, EJ, & CD Jones. 2011. Genetics
#see MANUAL.pdf for details




#### PSIseq directory #############################
$PSIseq_dir = "~/PSIseq";  # change the content within quotes to be your actual PSIseq directory
################################################### eg: "/Applications/PSIseq/" if you put PSIseq in the Applications folder




$clear = `clear`;
print $clear;
print "\n\nPSIseq pipeline\n";
print "\tPhenotype-based Selection & Introgression with whole-genome reSEQuencing\n";
print "\tSee Manual.pdf for more details\n\n\n";
print "\t1. Correlate mismatches (experimental vs. species) from pileups: shared_snps_v3.pl\n";
print "\t2. Binning: bos.pl\n";
print "\t3. Binomial Testing: stats.R\n";
print "\t4. Figures: line_graph.R\n";
print "choice: ";
$choice = <STDIN>;

#$dir_name = 'run_' . get_timestamp(); #from anonymous author
#sub get_timestamp {
#	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
#	$year=$year+1900;
#	return $year . sprintf('%02u%02u_%02u%02u%02u', $mon, $mday, $hour, $min, $sec);
#}	#end
#$mkdir = "mkdir $dir_name";
#`$mkdir`;

###################### 1 ########################################
###################### 1 ########################################
###################### 1 ########################################
#shared_snps_v3.pl
if ($choice == 1) {	
	print "chose $choice\n";
	print "### shared SNPs ###\n\n";
	print "hybrid vs parentX pileup file: ";
	$hybrid = <STDIN>;
	chomp $hybrid;
	print "parentY vs parentX pileup file: ";
	$species = <STDIN>;
	chomp $species;
	print "output name (default \"shared_hyb.txt\"): ";
	$shared_hyb = <STDIN>;
	chomp $shared_hyb;
	if ($shared_hyb =~ /./) {	
	} else {
		$shared_hyb = "shared_hyb.txt";
	}
	open(NOTES, ">Run_Notes.txt");
	print NOTES "input hybrid pileup: $hybrid\n";
	print NOTES "input parent pileup: $species\n";
	print NOTES "shared_snps: $shared_hyb\n";
	####correlate	
	&correlate;	#print "perl $PSIseq_dir/shared_snps_v3.pl hyb.mpileup parent.mpileup shared_hyb.txt\n";
	goto BOS;
}

###################### 2 ########################################
###################### 2 ########################################
###################### 2 ########################################
#bos.pl

if ($choice == 2) {	
	print "chose $choice\n";
	print "### Binning ###\n\n";
	print "hybrid vs parent shared SNPs file (shared_snps_v3.pl output): ";
	$shared_hyb = <STDIN>;
	chomp $shared_hyb;
	print "binned output (default \"bos_hyb.txt\"): ";
	$bos_hyb = <STDIN>;
	chomp $bos_hyb;	
	BOS:
	print NOTES "binned SNPs: $bos_hyb\n";	
	&bin;	#print "perl $PSIseq_dir/bos.pl shared_hyb.txt 1000 2000 1000 bos_hyb.txt\n";	
	goto STATS;

}


###################### 3 ########################################
###################### 3 ########################################
###################### 3 ########################################
#stats.R

if ($choice == 3) {	#Binomial Testing
	print "chose $choice\n";
	print "### Binomial Testing ###\n\n";
	print "binned hybrid vs parent file (bos.pl output): ";
	$bos_hyb = <STDIN>;
	chomp $bos_hyb;

	STATS:
	if ($bos_hyb =~ /./) {
	
	} else {
		$bos_hyb = "bos_hyb.txt";
	}
	&binomial;	
	print NOTES "Binomial testing: *.test\n";
	goto FIGURES;
}



###################### 4 ########################################
###################### 4 ########################################
###################### 4 ########################################
#line_graph.R


if ($choice == 4) {	#Figures	
	FIGURES:
	@dir = <*>;
	foreach $file (@dir) {

		if ($file =~ /\.test$/) {
			&figures;
		}
		
	}
	
	@dir = ();
	exit;
}




####################SUBROUTINES#####################
####################SUBROUTINES#####################
####################SUBROUTINES#####################




sub bin {
	print "\n### Bining ###\n";	
	$bos = "perl $PSIseq_dir/bos.pl $shared_hyb 1000 2000 1000 bos_hyb.txt";
	print "$bos\n";
	`$bos`;
}




sub binomial {
	$contig = ();
	$head = ();
	$line = ();
	$key = ();
	%hashish = ();


	print "\n### Binomial ###\n";	
	open(BOS, "$bos_hyb");
	$head = <BOS>;
	chomp $head;
	while ($line = <BOS>) {
		chomp $line;
		$line =~/^([\d\w]+)\s/;
		$contig = $1;
		push(@{$hashish{$contig}},$line);
	}
	foreach $key (keys %hashish) {
		open(CONTIG, ">$key");	#writes to file "$contig" (eg: 2L, 2R, etc.) the chunk of bos_hyb.txt containing that contig name
		print CONTIG $head."\n";
		foreach (@{$hashish{$key}}) {
			print CONTIG $_."\n";
		}
		close CONTIG;
		open(STATS, "stats.R");
		open(NEWSTATS, ">$key".".R");
		while ($line = <STATS>) {
			chomp $line;
			if ($line =~ /^data\<\-read\.table\(\"/) {
				$line = "data<-read.table(\"".$key."\", header=T)";
			}
			if ($line =~ /^write\.table/) {
				$line = "write.table(final, file=\"".$key.".test\", sep=\"\t\", col.names=TRUE, row.names=FALSE)\n";
			}

			print NEWSTATS $line."\n";
		}

		$stats = "R CMD BATCH $key".".R $key".".Rout";
		print "$stats\n";
		`$stats`;
		close STATS;		
		unlink ("$key".".R");		#eg: 2L.R
	}
}




sub figures {
	print "### Figures ###\n";
	open(LINE_GRAPH, "line_graph.R"); 
	open(R, ">$file\.R");
	while ($line = <LINE_GRAPH>) {
		chomp $line;
		if ($line =~ /data\s\=\sread\.table/) {
			$line = "data = read.table('".$file."', header=T)";
		}
		if ($line =~ /pdf\(/) {
			$line = "pdf(file=\"".$file.".pdf\", height = 5, width = 12)";
		}
		if ($line =~ /	xlab\=\"\"\,/) {
#			$test =~ s/\.test//g;
			$line = "	xlab=\"".$file."\",type=\"l\",frame.plot=T, lwd=2, yaxt=\"n\", ylab=\"\", main=\"".$file."\")";
		}
		
		print R $line."\n";
	
	}
	close R;
	
	$r = "R CMD BATCH ".$file.".R ".$file.".Rout";
	print $r."\n";
	`$r`;


}





























































sub correlate {
	print "### Correlate SNPs ###\n";
	$shared_snps = "perl shared_snps_v3.pl ".$hybrid." ".$species." $shared_hyb";
	print "$shared_snps\n";
#	`$shared_snps`;




	open(HYBRID, $hybrid) || die "can't find file $hybrid\n";
	open(OUT, ">$shared_hyb");
	
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
	
	foreach $key (keys %contig_hash) {	
		 %hashish = ();
		 $mini_counter = 0;	
		 $sum_h_cvg = 0;
		 $hyb_mm_count = 0;
		 $sum_p_cvg = 0;
		 $par_mm_count = 0;	
	
	
		open(HYBRID, $hybrid);
		###load hybrid hash
		print "loading $key hybrid hash...\n";
		while ( $hybrid_line = <HYBRID>) {
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
						}
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
	
	
		close HYBRID;		
	##### Parent	
		print "correlating to parent mismatches...\n";
		open(PARENT, $species) || die "can't find file $species\n";
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
		print "avg cvg parent: ".($sum_p_cvg / $mini_counter)."\n";
		$mini_counter = 0;
	
	}
	
	#snp call. Plurality wins.
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

}