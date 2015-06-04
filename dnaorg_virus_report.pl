#!/usr/bin/env perl
# EPN, Thu Jun  4 14:45:12 2015
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

# examine an output file from dnaorg_compare_genomes and find and report 'anomalies'

my $usage  = "\ndnaorg_virus_report.pl\n";
$usage .= "\t<output from dnaorg_compare_genomes.pl>\n";
$usage .= "\n";

#my $outdir = undef;
#&GetOptions( "d=s" => \$outdir);

if(scalar(@ARGV) != 1) { die $usage; }
my ($infile) = (@ARGV);


# variables that let us know where we are in the file
my $in_genome_section    = 0;
my $in_summary_section   = 0;
my $past_summary_section = 0;

# Do 2 passes over the file. On the first pass collect stats (CDS lengths, etc.)
# And on the second pass detect anomalies based on those stats.
# 
for(my $p = 1; $p <= 2; $p++) { 
  if($p == 2) { # convert our total lengths into averages
    # TODO: write code to calc avgs
  }

  open(IN, $infile) || die "ERROR unable to open $infile on pass $p";
  while(my $line = <IN>) { 
    if($line =~ m/\w/) { 
      # potentially update variables which tell us where we are
      if((! $in_genome_section) && (! $in_summary_section)) { 
        if($line !~ m/^\#/) { 
          $in_genome_section = 1;
        }
      }
      elsif($in_genome_section) { 
        if($line =~ /^\# Number-of-classes\:\s+(\d+)/) { 
          $in_genome_section   = 0;
          $in_summary_section  = 1;
        }
      }
      elsif($in_summary_section) { 
        if($line =~ m/^# Fetching/) { 
          $in_summary_section   = 0;
          $past_summary_section = 1;
        }
      }

      # sanity check, if we're at the line with column headings for the genome section,
      # make sure we have the fields we think we do:
      my $nexp     = 10; # number of expected fields in genome section
      my $npregene =  9; # number of expected fields before the gene columns in the genome section
      if($line =~ m/^\#accn/) { 
        # this is what we expect:
        #accn       #cds   #pos   #neg  #both  #unkn  strand-string  cls  tot-len       g1     g2     g3     g4     g5     g6     g7     g8
        chomp $line;
        my @elA  = split(/\s+/, $line);
        my @expA = ("#accn", "#cds", "#pos", "#neg", "#both", "#unkn", "strand-string", "cls", "tot-len", "g1");
        my $nexp = 10;
        if(scalar(@elA) < $nexp) { die "ERROR fewer than expected columns in genome section"; }
        for(my $i = 1; $i <= $nexp; $i++) { 
          if($elA[($i-1)] ne $expA[($i-1)]) { die "ERROR genome section header line, token $i is unexpected: %s ne %s", $elA[($i-1)], $expA[($i-1)]; } 
        }
      }    

      # based on where we are, look for anomalies
      if($in_genome_section) { 
        # we already checked that the column headers should be 
        #accn       #cds   #pos   #neg  #both  #unkn  strand-string  cls  tot-len       g1     g2     g3     g4     g5     g6     g7     g8
        #NC_003977      7      7      0      0      0  +++++++          1     3215     2532   1203    846    681    465    639    552
        my @elA  = split(/\s+/, $line);
        my $nel  = scalar(@elA);
        my $accn       = $elA[0];
        my $ncds       = $elA[1];
        my $npos       = $elA[2];
        my $nneg       = $elA[3];
        my $nboth      = $elA[4];
        my $nunkn      = $elA[5];
        my $strand_str = $elA[6];
        my $cls        = $elA[7];
        my $totlen     = $elA[8];
        my @cdslen_A   = ();
        my $ncds2      = $nel - $npregene;
        if($ncds == 0) { # genomes with 0 genes have no 'strand-string' token, fix that:
          $strand_str = "";
          $cls        = $elA[6];
          $totlen     = $elA[7];
          $ncds2++;
        }
        if($ncds ne $ncds2) { die "ERROR discrepancy between stated and inferred number of CDS: $ncds ne $ncds2"; }
        for(my $i = 0; $i < $ncds; $i++) { 
          $cdslen_A[$i] = $elA[($npregene+$i)];
        }

        # if we're in the first pass, gather statistics
        if($p == 1) { 
        }
        
        # if we're in the second pass, look for anomalies
        else { # we're 
          my $a1 = check_a1($ncds);
          my $a2 = check_a2($nboth);
          my $a3 = check_a3($nunkn);

          if($a1 || $a2 || $a3) { 
            printf("A $a1 $a2 $a3 $line"); 
          }
        }
      }
      elsif($in_summary_section) { 
        # do nothing, yet
      }
    }
  }
  close(IN);
}

#############
# SUBROUTINES
#############
# Subroutine: check_a1()       
# Synopsis:   Check for anomaly 1: 0 CDS.
# Args:       $ncds:        number of CDS
#
# Returns:    '1' if 0 cds exist, else '0'

sub check_a1 {
  my $sub_name = "check_a1()";
  my $nargs_exp = 1;

  my ($ncds) = (@_);
  
  return ($ncds == 0) ? '1' : '0';
}

# Subroutine: check_a2()       
# Synopsis:   Check for anomaly 2: >0 CDS on 'both' strands.
# Args:       $nboth: number of CDS on both strands
#
# Returns:    '1' if >0 cds on both strands exist, else '0'

sub check_a2 {
  my $sub_name = "check_a2()";
  my $nargs_exp = 1;

  my ($nboth) = (@_);
  
  return ($nboth > 0) ? '1' : '0';
}

# Subroutine: check_a3()       
# Synopsis:   Check for anomaly 3: >0 CDS on 'unknown' strand.
# Args:       $nunkn: number of CDS on unknown strad
#
# Returns:    '1' if >0 cds on unkwown strand exist, else '0'

sub check_a3 {
  my $sub_name = "check_a3()";
  my $nargs_exp = 1;

  my ($nunkn) = (@_);
  
  return ($nunkn > 0) ? '1' : '0';
}
