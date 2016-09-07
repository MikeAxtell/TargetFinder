#!/usr/bin/perl -w
use strict;
use Getopt::Std;

my $help = "
Shuffler.pl : Prediction of plant miRNA target sites, including control of false discovery rates, based on shuffled controls

Usage: Shuffler.pl -s <sequence> -d <transcriptome.fasta>

Options:
 -q <string> : Query name (deafults to 'query').
 -c <float> : Max score to test for (default = 8).
 -t <int> : # of threads for Smith-Waterman search (default = 1).
 -n <int> : Number of shuffle permutations used to derived median value (default = 10).
 -f <float> : False discovery rate .. alignments estimated to have FDR at or below this (default = 0.25).
 -r       : Also search the reverse strand for targets.
 -h       : Print this message and quit.

Dependencies (in PATH .. see notes):
 targetfinder.pl  : From github.com/carringtonlab/TargetFinder  .. edit code to match up to your version of ssearch from FASTA
 ssearch36 : Smith-Waterman aligner from FASTA package .. http://fasta.bioch.virginia.edu/fasta_www2/fasta_down.shtml
 uShuffle : Compiled C binary from http://digital.cs.usu.edu/~mjiang/ushuffle/ (rename it to uShuffle)
";

my %opt = ();
# set defaults .. those not set here are set by targetfinder.pl
$opt{'c'} = 8;
$opt{'n'} = 10;

getopts('d:s:q:c:t:n:rh', \%opt);

# give help if requested
if($opt{'h'}) {
    die "$help\n";
}

# check for uShuffle and targetfinder.pl
open(US, "which uShuffle |");
my $us_check = <US>;
close US;
if(!$us_check) {
    die "Required program uShuffle not found in \$PATH\n\n$help\n\n";
}

open(TF, "which targetfinder.pl |");
my $tf_check = <TF>;
close TF;
if(!$tf_check) {
    die "Required program targetfinder.pl not found in \$PATH\n\n$help\n\n";
}

# check for required options -s and -d
unless(($opt{'s'}) and ($opt{'d'})) {
    die "Both options -s and -d are required.\n\n$help\n\n";
}

# build a generic command chunk
my $c_chunk = "-c $opt{'c'} -d $opt{'d'} -p table";
if($opt{'t'}) {
    $c_chunk .= " -t $opt{'t'}";
}
if($opt{'r'}) {
    $c_chunk .= " -r";
}
if($opt{'q'}) {
    $c_chunk .= " -q $opt{'q'}";
}

print STDERR "Searching against input query $opt{'s'}\n";
# Call the real query, store all results
open(TFR, "targetfinder.pl -s $opt{'s'} $c_chunk |");
my @real_hits = <TFR>;
close TFR;

# Store counts of true hits
my %true_simple = ();
foreach my $real_hit (@real_hits) {
    my @rh_fields = split ("\t", $real_hit);
    ++$true_simple{$rh_fields[5]};
}

# Tally up the true hits - cumulative
my %true_cumulative = ();
my $tc = 0;
for(my $i = 0; $i <= $opt{'c'}; $i += 0.5) {
    if(exists($true_simple{$i})) {
	$tc += $true_simple{$i};
    }
    $true_cumulative{$i} = $tc;
}

# prepare for shuffling
my %notOK = ($opt{'s'} => 1);
my %shuf_arrays = ();

# shuffle loop
for(my $i = 1; $i <= $opt{'n'}; ++$i) {
    my $shufseq = get_shuf(\$opt{'s'}, \%notOK);
    
    # test
    #print "shufseq with no newline after it is $shufseq";
    #exit;
    
    $notOK{$shufseq} = 1;
    my %shuf_single = ();
    for(my $j = 0; $j <= $opt{'c'}; $j += 0.5) {
	$shuf_single{$j} = 0;
    }
    print STDERR "Searching against shuffled permutation number $i $shufseq\n";
    open(TF, "targetfinder.pl -s $shufseq $c_chunk |");
    while (<TF>) {
	chomp;
	my @tf_f = split ("\t", $_);
	++$shuf_single{$tf_f[5]};
    }
    close TF;
    for(my $j = 0; $j <= $opt{'c'}; $j += 0.5) {
	push(@{$shuf_arrays{$i}}, $shuf_single{$i});
    }
}

# calculations for shuffles .. medians, cumulative median, and inverse-cumulative median
my %shuf_median_simple = ();
for(my $k = 0; $k <= $opt{'c'}; $k += 0.5) {
    my $med = median(@{$shuf_arrays{$k}});
    $shuf_median_simple{$k} = $med;
}
my %shuf_median_cumulative = ();
my $smc = 0;
for(my $k = 0; $k <= $opt{'c'}; $k += 0.5) {
    $smc += $shuf_median_simple{$k};
    $shuf_median_cumulative{$k} = $smc;
}
my %shuf_median_inverse_cumulative = ();
my $sic;
for(my $k = 0; $k <= $opt{'c'}; $k += 0.5) {
    $sic = $shuf_median_cumulative{$opt{'c'}} - $shuf_median_cumulative{$k};
    $shuf_median_inverse_cumulative{$k} = $sic;
}

# test
print "Score\tTrueSimple\tTrueCumulative\tShufSimple\tShufCumulative\tShufInverseCumulative\n";
for(my $k = 0; $k <= $opt{'c'}; $k += 0.5) {
    print "$k\t";
    print "$true_simple{$k}\t";
    print "$true_cumulative{$k}\t";
    print "$shuf_median_simple{$k}\t";
    print "$shuf_median_cumulative{$k}\t";
    print "$shuf_median_inverse_cumulative{$k}\n";
}


sub median {
    my(@data) = @_;
    my @sorted = sort {$a <=> $b} @data;
    my $length = scalar @sorted;
    my $median;
    my $i = int($length / 2);
    if($length % 2) {
	# odd
	$median = int($sorted[$i]);
    } else {
	# even
	$median = int(($sorted[$i - 1] + $sorted[$i]) / 2);
    }
    return $median;
}
	


sub get_shuf {
    my($input, $hash) = @_;
    my $attempts = 0;
    my $ok = 0;
    my $output;
    until($ok) {
	if($attempts >= 10000) {
	    die "\nFAILURE: Sorry, couldn't find enough valid shuffles of $$input ... is the sequence low complexity?\n";
	}
	my $seed = int(rand(10000));
	open(SHUF, "uShuffle -s $$input -n 1 -k 2 -seed $seed |");
	$output = <SHUF>;
	chomp $output;
	close SHUF;
	unless(exists($$hash{$output})) {
	    $ok = 1;
	}
	++$attempts;
    }
    return $output;
}

    
    


