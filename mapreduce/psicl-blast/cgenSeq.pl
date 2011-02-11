#!/usr/bin/perl

$n = $ARGV[0];
$m = $ARGV[1];

for ($j = 0; $j < $m; $j++) {
	print ">lcl|$j\n";
	for ($i = 0; $i < $n; $i++) {
		$r = rand();
		if ($r < 0.25 ) {
			print "A";
		} elsif ($r < 0.5) {
			print "T";
		} elsif ($r < 0.75) {
			print "C";
		} else {
			print "G";
		}
	}
	print "\n";
}
