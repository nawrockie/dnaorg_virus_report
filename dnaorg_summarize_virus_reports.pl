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
my $desc;
my $nsngl;
my $w_name = length("average-per-species");
my $nreports;
my @explanation_A = (); # array of 'explanation' lines
my $stored_explanations = 0;

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
    $name =~ s/\..+$//;     # remove everything after and including final '.'
    push(@name_A, $name);
    if(length($name) > $w_name) { $w_name = length($name); }
  }
}
close(IN);

my $min_aidx = 1000;
my $max_aidx = 0;
my @outline_A = (); # array of the output lines
my @desc_A = ();
my $ridx = 0;

my $tot_ngenomes = 0;
my $tot_nclasses = 0;
my $tot_ngenomes_a = 0;
my $sum_fract_genomes_a = 0.;
my @tot_act_A      = ();
my @tot_act_sngl_A = ();

$tot_act_A[0] = 0;
$tot_act_sngl_A[0] = 0;

for($i = 0; $i < $nreports; $i++) { 
  my $report = $report_A[$i];
  my $name   = $name_A[$i];

  my $ngenomes   = undef;
  my $nclasses   = undef;
  my $ngenomes_a = undef;

  my @act_A      = ();
  my @act_sngl_A = ();
  $act_A[0]      = 0;
  $act_sngl_A[0] = 0;
  $desc_A[0]     = "";

  open(IN, $report) || die "ERROR unable to open $report"; 
  while(my $line = <IN>) { 
    chomp $line;
    if($line =~ /^\# Number-of-classes:\s+(\d+)/) { 
## Number-of-classes: 3
      $nclasses =$1;
    }
    elsif($line =~ /^(\d+)\s+(\S+)\s+(.+)$/) { 
##               count
##           -------------
##anomaly-#   total (sngl)  description
##---------  -------------  -----------
#1               33    (0)  genomes have 0 CDS annotations.
#5              N/A         genomes have a CDS on the positive strand, when >= 0.990 fraction of all genomes do not
      ($aidx, $ct, $desc) = ($1, $2, $3);
      if($aidx < $min_aidx) { $min_aidx = $aidx; }
      if($aidx > $max_aidx) { $max_aidx = $aidx; }
      $nsngl = 0;
      if($desc =~ /^\s*\((\d+)\)\s*/) { 
        $nsngl = $1;
        $desc =~ s/^\s*\(\d+\)\s*//;
      }
      $act_A[$aidx]      = ($ct eq "N/A") ? 0 : $ct;
      $tot_act_A[$aidx] += $act_A[$aidx];
      $act_sngl_A[$aidx]     += $nsngl;
      $tot_act_sngl_A[$aidx] += $nsngl;
        
      $desc_A[$aidx] = $desc;
    }
    elsif($line =~ /^total\s+\d+\s+anomalies exist in (\d+) genomes of (\d+)/) { 
#total            2  anomalies exist in 2 genomes of 427 (0.0047)
      ($ngenomes_a, $ngenomes) = ($1, $2);
    }
    elsif($line =~ /^# Explanation of class definition/) { 
      if(! $stored_explanations) { 
        push(@explanation_A, $line . "\n");
        while($line = <IN>) { 
          push(@explanation_A, $line);
        }
        $stored_explanations = 1;
      }
    }
  }
  close(IN);
  if(! defined $ngenomes)   { die "ERROR didn't read number of genomes in $report"; }
  if(! defined $ngenomes_a) { die "ERROR didn't read number of genomes with anomalies in $report"; }
  if(! defined $nclasses)   { die "ERROR didn't read number of classes in $report"; }

  $tot_ngenomes   += $ngenomes;
  $tot_ngenomes_a += $ngenomes_a;
  $tot_nclasses   += $nclasses;
  $sum_fract_genomes_a += ($ngenomes_a / $ngenomes);

  $outline = sprintf("%-*s  %8d  %8d  %9d  %9.4f |", 
                        $w_name, $name, $nclasses, $ngenomes, $ngenomes_a, ($ngenomes_a / $ngenomes));
  for($aidx = $min_aidx; $aidx <= $max_aidx; $aidx++) { 
    #$outline .= sprintf("  %5d (%6s)", $act_A[$aidx], 
    #                                   ($act_A[$aidx] == 0) ? "-" : sprintf("%6.4f", ($act_A[$aidx] / $ngenomes)));
    $outline .= sprintf(" %5d %4s |", $act_A[$aidx], $act_sngl_A[$aidx] > 0 ? sprintf("(%d)", $act_sngl_A[$aidx]) : "");
  }
  $outline .= "\n";
  push(@outline_A, $outline);
}

# output column headers:
my $name_dashes = "";
for (my $i = 0; $i < $w_name-1; $i++) { $name_dashes .= "-"; }

my $w_count = ($max_aidx - $min_aidx + 1) * 13 + 1;
my $count_dashes = "";
my $count_string = "";
for (my $i = 0; $i < $w_count; $i++) { $count_dashes .= "-"; }
for (my $i = 0; $i < ($w_count/2) - 3; $i++) { $count_string .= " "; }

printf("#%-*s  %8s  %8s  %9s  %9s %s\n", $w_name-1, "", "", "", "", "", $count_string . "counts");
printf("#%-*s  %8s  %8s  %9s  %9s %s\n", $w_name-1, "", "", "", "", "", $count_dashes);


printf("#%-*s  %8s  %8s  %9s  %9s |", 
       $w_name-1, "", "", "", "#genomes", "fraction"); 
for($aidx = $min_aidx; $aidx <= $max_aidx; $aidx++) { 
  printf(" %10s |", sprintf("anomaly#%d", $aidx));
}
printf("\n");

printf("#%-*s  %8s  %8s  %9s  %9s |", 
       $w_name-1, "name", "#classes", "#genomes", "w/anomaly", "w/anomaly");
for($aidx = $min_aidx; $aidx <= $max_aidx; $aidx++) { 
  printf(" %10s |", "   ct sngl");
}
printf("\n");

printf("#%-*s  %8s  %8s  %9s  %9s |", 
       $w_name-1, $name_dashes, "--------", "--------", "---------", "---------"); 
for($aidx = $min_aidx; $aidx <= $max_aidx; $aidx++) { 
  printf(" ----- ---- |");
}
printf("\n");

foreach $outline (@outline_A) { 
  print $outline;
}

printf("#%-*s  %8s  %8s  %9s  %9s |", 
       $w_name-1, "", "", "", "", "");
for($aidx = $min_aidx; $aidx <= $max_aidx; $aidx++) { 
  printf(" %10s |", "");
}
printf("\n");

# total line
printf("%-*s  %8d  %8d  %9d  %9.4f |", 
       $w_name, "total", $tot_nclasses, $tot_ngenomes, $tot_ngenomes_a, $tot_ngenomes_a / $tot_ngenomes);
for($aidx = $min_aidx; $aidx <= $max_aidx; $aidx++) { 
  printf(" %5d %4s |", $tot_act_A[$aidx], $tot_act_sngl_A[$aidx] > 0 ? sprintf("(%d)", $tot_act_sngl_A[$aidx]) : "");
}
printf("\n");

# average line
printf("%-*s  %8.1f  %8.1f  %9.1f  %9.4f |", 
       $w_name, "average-per-species", $tot_nclasses / $nreports, $tot_ngenomes / $nreports, $tot_ngenomes_a / $nreports, $sum_fract_genomes_a / $nreports);
for($aidx = $min_aidx; $aidx <= $max_aidx; $aidx++) { 
  printf(" %5.1f %4s |", $tot_act_A[$aidx] / $nreports, "");
}
printf("\n");

printf("#%-*s  %8s  %8s  %9s  %9s %s\n", $w_name-1, "", "", "", "", "", $count_dashes);

printf("#\n");
printf("#\n");
for($aidx = $min_aidx; $aidx <= $max_aidx; $aidx++) { 
  printf("# anomaly #%d  %s\n", $aidx, $desc_A[$aidx]);
}

printf("#\n");
printf("#\n");
foreach my $line (@explanation_A) {
  print $line;
}
