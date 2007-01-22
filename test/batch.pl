#!/usr/bin/perl -w
# Author: C. Bernet 2005

#use strict;
use diagnostics;
use Getopt::Long;
use File::Basename;


my $sfiles=0;
my $filePattern=0;


my $masterfile="pfRootEvent.opt";

my $doNotSubmit=0;
my $help=0;

my $bsub=0;

my $pwd = `pwd`;
chomp $pwd;


my $date = `date +%d%b%Y_%H%M%S`;
my $scandir = "BatchOut_$date";
chomp $scandir;
`mkdir $scandir`;
`echo "scan.pl @ARGV" > scan.pl.log`;
`mv scan.pl.log $scandir`;

GetOptions (
	    'files=s' => \$sfiles,
	    'pattern=s' => \$filePattern, 
	    'master=s' => \$masterfile, 
	    'h' => \$help,
	    'n' => \$doNotSubmit, 
	    'b=s' => \$bsub );

if($help) {
    print "usage : batch.pl -files=\"*.root\" -master myconfig.cfg [-p <regexp>] [-n] [-b \"bsub -q 8nm\"]\n";
    print " -p <regexp>  : specify a regexp to filter input files\n"; 
    print " -n           : do not proceed\n"; 
    print " -b <>        : run on the batch system (LSF)\n";
    exit(1);
}

if($doNotSubmit) {
    print "will do nothing... \n";
}

print "master : $masterfile\n";
print "======== $sfiles ======= \n";


my $ls = "ls";
my $basedir = "$sfiles"; # contains directory path


if( $sfiles =~ /^\/castor/) {
#    print "source files are on castor\n";
    $ls = "nsls";
    $basedir = "$sfiles";
}  

my @rootfiles; # will contain base filenames (no path)

my @tmpfiles = `$ls $sfiles`;
foreach my $file (@tmpfiles) {
    if( $file =~ /\.root$/ ) {
	if( ( $filePattern && 
	      $file =~ /$filePattern/ ) ||
	    !$filePattern ) {
	    push(@rootfiles, basename($file) );
		
	}
    }
}


my $outdir = "$pwd/$scandir";    
my $outcastordir = "$sfiles/$scandir";    
`nsmkdir $outcastordir`;


foreach my $rootfile (@rootfiles) {
    chomp $rootfile;
    my $fullfilename = "$basedir/$rootfile";
    print "processing :  $fullfilename\n";
    
#    if($svalues == 0) {
    my $optfile = "$outdir/auto_$rootfile.cfg";
    my $outrootfile = "$outcastordir/out_$rootfile";
    
    open(IN, "<$masterfile");
    open(OUT, ">$optfile");
    while ( <IN> ) {
	my $line = $_;
	if($line =~ /^\s*\#/) { # commented line
	    print OUT "$line";	    
	    next;
	}
	if($line =~ /untracked vstring fileNames/) {
	    print OUT "untracked vstring fileNames = {'rfio:$fullfilename'}\n";
	}
	elsif($line =~ /untracked string fileName/) {
	    print OUT "untracked string fileName = \"rfio:$outrootfile\" \n";
	}
	else {
	    print OUT "$line";
	}
    }
    close(IN);
    close(OUT);
    

    if(! $doNotSubmit) {
#	    my $outrootfile = "$scandir/out_$rootfile";
	if( !$bsub ) { # standard execution
	    `nice -n 5 cmsRun $optfile`;
#		`mv out.root $outrootfile`;
	}
	else { # batch execution 
	    # print "batch \n";
	    `$bsub  cmsRun $optfile`;
	    # `mv out.root $outrootfile`;		
	}
	#`cp $macro $scandir`;
	#`cp $optfile $scandir`;	    
    }
    
}
