#!/usr/bin/perl -w
# Author: C. Bernet 2005

#use strict;
use diagnostics;
use Getopt::Long;

my $runmin = -1;
my $runmax = 999999999;
my $regexp = "";
GetOptions ('min=s' => \$runmin,
	    'max=s' => \$runmax, 
	    'regexp=s'=> \$regexp );

my $masterfile = shift;
my $castordir = shift;


print "input files  :\n";

my @rootfiles;

my @tmpfiles = `nsls $castordir`;
foreach my $file (@tmpfiles) {
    if( $file =~ /\.root$/ ) {
	push(@rootfiles, $file); 
	print "\t\t$castordir/$file";
    }
}

print "master       : $masterfile\n";


foreach my $rootfile (@rootfiles) {
    chomp $rootfile;
    print "processing : $rootfile\n";

    if($rootfile =~/(\d+).root/) {
	my $run = $1;
	if( ($run > $runmax) || ($run<$runmin) ) {
	    print "skipping $run > $runmax or $run < $runmin  \n";
	    next;
	}
    }
    
    if($rootfile !~ /$regexp/) {
	print "skipping, does not match $regexp\n";
	next;
    }

    my $outfile = "auto_$masterfile";

    open(IN, "<$masterfile");
    open(OUT, ">$outfile");
    while ( <IN> ) {
	my $line = $_;
	    
	if($line =~ /untracked vstring fileNames/) {
	    print OUT "untracked vstring fileNames = {'rfio:$castordir/$rootfile'}\n";
	}
	elsif($line =~ /untracked string fileName/) {
	    print OUT "untracked string fileName = \"b_$rootfile\" \n";
	}
	else {
	    print OUT "$line";
	}
    }
    close(IN);
    close(OUT);

    `cmsRun $outfile`;
}

