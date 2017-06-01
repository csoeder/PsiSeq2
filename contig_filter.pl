#!/usr/bin/perl
#written by Eric J Earley, 2011
#dirty script to filter reference sequence by contig size (eg: Dsec-r1.3.fa has 14k+)
#fasta format only



$file = $ARGV[0];
$size = $ARGV[1];
unless (@ARGV){
        print "Grabs only contigs greater than X\ninput: <in.fa><size> <out>\n";
        exit;
}
open (FILE, $file) || die "oy vey $file\n";
open(OUT, $ARGV[2]);
while ($line =<FILE>) {
	chomp $line; 
	if ($line =~ /^[ATCGN]/) {
                $bglen+= length($line);
        }
	if ($line =~ /^>/) {
		$bgcontig++;
	}

}
close FILE;

open(FILE, $file);
until (eof FILE) {
        SPOT2:
	$line = <FILE>;
        chomp $line;
        if ($line =~ /length=(\d+)/) {
		$length = $1;
		if ($length > $size) {
			SPOT1:
			$smcontig++;
        	        print OUT $line."\n";
                	$line = <FILE>;
			chomp $line;
			$seq .= $line;
			while ($line !~ m/^>/) {
				$line = <FILE>;
				chomp $line;
				if ($line =~ /length=(\d+)/) {
					$length = $1;
					print OUT $seq."\n";
					$smlen+= length($seq);
					$seq = ();
					if ($length > $size) {
						goto SPOT1;
					}else {
						goto SPOT2;
					}
					
				} else {
					$seq .= $line;
                   		}    
			 
			}
                }
	
        }

}
print "size pre-filter $bglen\n";
print "contigs pre-filter $bgcontig\n";
print "size post-filter $smlen\n";
print "contigs post-filter $smcontig\n";
