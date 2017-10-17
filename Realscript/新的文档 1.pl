#!/usr/bin/perl
use strict;
use warnings;
open (IN, <STDIN>);
foreach  (@seq) {
	if ($seq[2*$n] ~= $seq[2$n-1]) {
		$n+ = 1;
	}
}


