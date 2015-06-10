#!/usr/bin/env perl
# EPN, Thu Jun 4 14:45:12 2015
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

# Summarize many 'virus report' files output from dnaorg_virus_report.pl.

my $usage  = "\ndnaorg_summarize_virus_reports.pl\n";
$usage .= "\t<file with list of dnaorg_virus_report.pl output files to summarize>\n";
$usage .= "\n";
#$usage .= " BASIC OPTIONS:\n";
#$usage .= "  -F4     <f>: set threshold for anomaly 4 to <f> (gene count/strand order differs from >= <f> fraction of genomes)                       [df: $d#efault_F4]\n";
#$usage .= "\n";

#my $i; # a counter
#&GetOptions( "F4"     => \$F4, 
#             "F5and6" => \$F5and6, 
#             "F7"     => \$F7, 
#             "F8a"    => \$F8a,
#             "F8b"    => \$F8b);

if(scalar(@ARGV) != 1) { die $usage; }
my ($listfile) = (@ARGV);

my @report_A = ();
my @name_A   = ();
my $aidx;
my $ct;
my $name;
my $outline;
my $i;
my $w_name = 0;
my $nreports;

open(IN, $listfile) || die "ERROR unable to open $listfile";
while(my $line = <IN>) { 
  if($line =~ m/\w/) { 
    chomp $line;
    if(! -s $line) { die "ERROR $line file does not exist or is empty"; }
    push(@report_A, $line);
    $nreports++;

    if($line !~ m/\.report$/) { die "ERROR $line file does not end in .report"; }
    $name = $line;
    $name =~ s/\.report$//; # remove trailing .report
    $name =~ s/^.+\///;     # remove everything up to final '/'
    push(@name_A, $name);
    if(length($name) > $w_name) { $w_name = length($name); }
  }
}
close(IN);

my $min_aidx = 1000;
my $max_aidx = 0;
my @outline_A = (); # array of the output lines
my $ridx = 0;

for($i = 0; $i < $nreports; $i++) { 
  my $report = $report_A[$i];
  my $name   = $name_A[$i];

  my $ngenomes   = undef;
  my $nclasses   = undef;
  my $ngenomes_a = undef;

  my @act_A = ();
  $act_A[0] = 0;
  open(IN, $report) || die "ERROR unable to open $report"; 
  while(my $line = <IN>) { 

    if($line =~ /^\# Number-of-classes:\s+(\d+)/) { 
## Number-of-classes: 3
      $nclasses =$1;
    }
    elsif($line =~ /^(\d+)\s+(\S+)/) { 
##anomaly-#   count  description
##---------  ------  -----------
#1                0  genomes have 0 CDS annotations.
      ($aidx, $ct) = ($1, $2);
      if($aidx < $min_aidx) { $min_aidx = $aidx; }
      if($aidx > $max_aidx) { $max_aidx = $aidx; }
      $act_A[$aidx] = ($ct eq "N/A") ? 0 : $ct;
    }
    elsif($line =~ /^total\s+\d+\s+anomalies exist in (\d+) genomes of (\d+)/) { 
#total            2  anomalies exist in 2 genomes of 427 (0.0047)
      ($ngenomes_a, $ngenomes) = ($1, $2);
    }
  }
  close(IN);
  if(! defined $ngenomes)   { die "ERROR didn't read number of genomes in $report"; }
  if(! defined $ngenomes_a) { die "ERROR didn't read number of genomes with anomalies in $report"; }
  if(! defined $nclasses)   { die "ERROR didn't read number of classes in $report"; }

  $outline = sprintf("%-*s  %8d  %8d  %9d  %9.4f", 
                        $w_name, $name, $ngenomes, $nclasses, $ngenomes_a, ($ngenomes_a / $ngenomes));
  for($aidx = $min_aidx; $aidx < $max_aidx; $aidx++) { 
    $outline .= sprintf("  %5d (%6s)", $act_A[$aidx], 
                        ($act_A[$aidx] == 0) ? "-" : sprintf("%6.4f", ($act_A[$aidx] / $ngenomes)));
  }
  $outline .= "\n";
  push(@outline_A, $outline);
}

# output column headers:
my $name_dashes = "";
for (my $i = 0; $i < $w_name; $i++) { $name_dashes .= "-"; }

printf("%-*s  %8s  %8s  %9s  %9s", 
       $w_name, "", "", "", "", "#genomes"); 
for($aidx = $min_aidx; $aidx < $max_aidx; $aidx++) { 
  printf("  %14s", sprintf("anomaly-#%d", $aidx));
}
printf("\n");

printf("%-*s  %8s  %8s  %9s  %9s", 
       $w_name, "name", "#genomes", "fraction", "#classes", "w/anomaly");
for($aidx = $min_aidx; $aidx < $max_aidx; $aidx++) { 
  printf("  %5s %8s", "count", "fraction");
}
printf("\n");

printf("%-*s  %8s  %8s  %9s  %9s", 
       $w_name, $name_dashes, "--------", "--------", "---------", "---------"); 
for($aidx = $min_aidx; $aidx < $max_aidx; $aidx++) { 
  printf("  ----- --------");
}
printf("\n");

foreach $outline (@outline_A) { 
  print $outline;
}
