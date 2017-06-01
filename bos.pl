#!/usr/bin/perl
#Eric J Earley April 2011

#	bos = binning, overlapping, sliding

######### input #############
######### input #############
# results of shared_snps.pl
#format, where header=chromosome arm (or whatever fastA header you have)
#1 = confirmed hybrid mismatch existing within SNDspace
#0 = no hybrid mismatch found, but parental SND exists
#header	1	pos
#header	1	pos
#header	0	pos
#header	0	pos
#header	0	pos
#header	0	pos
#header	0	pos
#header	1	pos
#header	0	pos
#header	1	pos



######## output ##############
######## output ##############
#Contig	BinSize	SlideRate	BinID	Position	OneRatio


#gives ratio: 1s/(1s+0s) for user defined bin sizes and user defined sliding rate



$FILE = $ARGV[0];
$MIN_SIZE = $ARGV[1];
$MAX_SIZE = $ARGV[2];
$RATE = $ARGV[3];
$OUT = $ARGV[4];


open(OUT, ">$OUT");
unless (@ARGV) {
	print "\ninput: <shared_snps.out> <min bin size> <max bin size> <bin growth rate> <out>\n";
	print "If you just want one bin size: <min=X> <max=X+Y> <rate=Y>\n";
	exit;
}

open (FILE, $FILE) || die "can't find $FILE\n";
while (my $line = <FILE>) {
	$line_count++;
	$line=~/^(.+)\s(\d+)\s([10])/;
	$contig = $1;
	$position = $2;
	$elem = $3;
	#contig	1	pos
	push(@{$hashish{$contig}},"$elem\t$position");	#a hash of arrays because we only care about number of SNPs (array index) not position
	$num_contig++;
}
close FILE;



$file_length = $line_count;
$line_count = ();

print "file length = ".$file_length."\n";
print "number of contigs = $num_contig\n";
print OUT "Contig\tBinSize\tSlideRate\tBinID\tPosition\tOneRatio\n";

foreach $key (sort {$a cmp $b} keys %hashish) {

	@contig_array = @{ $hashish{$key} };
	$array_length = @contig_array;

	
	$bin_size = $MIN_SIZE;
	until ($bin_size == $MAX_SIZE) {

		$slide_rate = 0.1 * $bin_size;
		$NEW_START = 1;

		#when end of bin surpasses array length, end loop.
		while ($end_of_bin < $array_length) {

			#populate 1s & 0s within bin
			until ($add_to_bin_counter == $bin_size) {

				$value = substr($hashish{$key}[$NEW_START + $add_to_bin_counter], 0, 1);
#print $value."\n";				
				push(@bin_array, $value);
				$add_to_bin_counter++;
			}
			$bin_num++;
			$size = @bin_array;
			$add_to_bin_counter = ();
			$end_of_bin = $NEW_START + $bin_size;

			foreach $value (@bin_array) {
				@val_a = split('\t', $value);
				$val = $val_a[0];
				$position = $val_a[1];
				if ($val == 1) {
					$ones++;
				}
				if ($val == 0) {
					$zeros++;
				}
			}
			$one_ratio = $ones / ($ones + $zeros);
			@holder_a = split('\t',$hashish{$key}[$NEW_START]);
			$position = $holder_a[1];
			print "$key\t$bin_size\t$slide_rate\t$bin_num\t$position\t$one_ratio\n";
			print OUT "$key\t$bin_size\t$slide_rate\t$bin_num\t$position\t$one_ratio\n";

			$NEW_START = $NEW_START + (0.8 * $bin_size);	#increase bin start site by user-defined sliding rate



			$ones = 0;
			$zeros = 0;
			$one_ratio = 0;
			@bin_array = ();
			close FILE;	

		
	

		}
		$end_of_bin = ();
		$bin_num = ();
		$bin_size = $bin_size + $RATE;
	
	
	}


}