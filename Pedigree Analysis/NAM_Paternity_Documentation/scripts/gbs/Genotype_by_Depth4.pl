#!/usr/bin/perl -w
# Modified from Sanzhen Liu
# 12/19/2016

use strict;
use warnings;

# open the file:
open(IN, $ARGV[0]) || die;  #open filehandle
while (<IN>) { #while file handle is open read line by line
	if(/\#\#/){ #if the line starts with ##
        print #print the line by default $_
	}
	if (!/\#\#/) { #if the line starts with #Header row
		my %fd = (); # initiate fd, format-depristo
        my @t = split;  #split by defualt
   		if (/^\#CHROM/) {
            print "#CHROM\t"; #print to file
			print join("\t", @t[1,2,3,4,5,6,7,8]);
			for (my $i=9; $i<= $#t; $i++) {
				printf("\t%s", $t[$i]);
			}
			print "\n";
			next;
		}
		
		print join("\t", @t[0,1,2,3,4,5,6,7,8]);

		 my @format = split(/:/, $t[8]); # format column
		for (my $k = 9; $k <= $#t; $k++) { #loops through each SNP sample combo

			my @depristo = split(/:/, $t[$k]); #array to hold stuff alleles are in 0
			if($depristo[0] =~ m#(0/0|1/1)#){ #if allele is homozygous check allele count
			my $check = $depristo[2];
				if($check < 4){ #if less than 4 reads change to missing
				$depristo[0] = "./.";
				}
				
			}
            my $out = join (":", @depristo);
   			print("\t$out");         

 			}
		print "\n";
	}
}

close IN;
