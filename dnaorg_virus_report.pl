#!/usr/bin/env perl
# EPN, Thu Jun 4 14:45:12 2015
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);

# examine an output file from dnaorg_compare_genomes and find and report 'anomalies'

my $executable     = $0;
my $F4             = undef;
my $F5and6         = undef;
my $F7             = undef;
my $F8             = undef;
my $default_F4     = 0.99;
my $default_F5and6 = 0.99;
my $default_F7     = 0.25;
my $default_F8     = 0.25;

my $usage  = "\ndnaorg_virus_report.pl\n";
$usage .= "\t<output from dnaorg_compare_genomes.pl>\n";
$usage .= "\n";
$usage .= " BASIC OPTIONS:\n";
$usage .= "  -F4     <f>: set threshold for anomaly 4 to <f> (gene count/strand order differs from >= <f> fraction of genomes)            [df: $default_F4]\n";
$usage .= "  -F5and6 <f>: set threshold for anomaly 5 and 6 to <f> (class has >= 1 gene on a strand which <f> fraction of genomes do not) [df: $default_F5and6]\n";
$usage .= "  -F7     <f>: set threshold for anomaly 7 to <f> (genome length deviates by more than <f> fraction from mean)                 [df: $default_F7]\n";
$usage .= "  -F8     <f>: set threshold for anomaly 8 to <f> (>= 1gene length deviates by more than <f> fraction from class mean)         [df: $default_F8]\n";
$usage .= "\n";

my $i; # a counter
&GetOptions( "F4"     => \$F4, 
             "F5and6" => \$F5and6, 
             "F7"     => \$F7, 
             "F8"     => \$F8);

if(scalar(@ARGV) != 1) { die $usage; }
my ($infile) = (@ARGV);

# store options used, so we can output them 
my $opts_used_short = "";
my $opts_used_long  = "";
if(defined $F4) { 
  $opts_used_short .= "-F4 $F4";
  $opts_used_long  .= "# option:  setting F4 anomaly 4 threshold to $F4 [-F4]\n";
}
if(defined $F5and6) { 
  $opts_used_short .= "-F5and6 $F5and6";
  $opts_used_long  .= "# option:  setting F5and6 anomaly 5 and 6 threshold to $F5and6 [-F5and6]\n";
}
if(defined $F7) { 
  $opts_used_short .= "-F7 $F7";
  $opts_used_long  .= "# option:  setting F7 anomaly 7 threshold to $F7 [-F7]\n";
}
if(defined $F8) { 
  $opts_used_short .= "-F8 $F8";
  $opts_used_long  .= "# option:  setting F8 anomaly 8 threshold to $F8 [-F8]\n";
}

# check for incompatible option values/combinations:
if(defined $F4 && ($F4 < 0. || $F4 > 1.0)) { 
  die "ERROR with -F4 <f>, <f> must be a number between 0 and 1."; 
}
if(defined $F5and6 && ($F5and6 < 0. || $F5and6 > 1.0)) { 
  die "ERROR with -F5and6 <f>, <f> must be a number between 0 and 1."; 
}
if(defined $F7 && ($F7 < 0. || $F7 > 1.0)) { 
  die "ERROR with -F7 <f>, <f> must be a number between 0 and 1."; 
}
if(defined $F8 && ($F8 < 0. || $F8 > 1.0)) { 
  die "ERROR with -F7 <f>, <f> must be a number between 0 and 1."; 
}

# set defaults for variables not set by the user via options
if(! defined $F4)     { $F4     = $default_F4;     }
if(! defined $F5and6) { $F5and6 = $default_F5and6; }
if(! defined $F7)     { $F7     = $default_F7;     }
if(! defined $F8)     { $F8     = $default_F8;     }

# define anomaly descriptions:
my %desc_H = ();
$desc_H{"1"} = "0 CDS annotations.";
$desc_H{"2"} = "at least one gene with coding sequence on both strands";
$desc_H{"3"} = "at least one gene on an unknown strand";
$desc_H{"4"} = sprintf("a CDS order/strand string that differs from >= %.3f fraction of genomes", $F4); 
$desc_H{"5"} = sprintf("a CDS on the positive strand, and >= %.3f fraction of all genomes do not", $F5and6);
$desc_H{"6"} = sprintf("a CDS on the negative strand, and >= %.3f fraction of all genomes do not", $F5and6);
$desc_H{"7"} = sprintf("total length that deviates by more than %.3f fraction from the mean", $F7);
$desc_H{"8"} = sprintf("CDS length that deviates by more than %.3f fraction from the mean for that class", $F7);
$desc_H{"9"} = sprintf("shifting gene order minimizes CDS length deviation from mean");

my $na_types = 9;
my @na_A  = (); #  na_A[$i]: number of genomes with $i anomalies
my @act_A = (); # act_A[$i]: number of genomes that have anomaly $i
for($i = 0; $i <= $na_types; $i++) { 
  $na_A[$i]  = 0;
  $act_A[$i] = 0;
}

###############
# Preliminaries
###############
# output banner
my $script_name = "dnaorg_virus_report.pl";
my $script_desc = "Report anomalies for a set of virus genomes of the same species";
print ("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
print ("# $script_name: $script_desc\n");
print ("# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n");
print ("# command: $executable $opts_used_short $infile\n");
printf("# date:    %s\n", scalar localtime());

if($opts_used_long ne "") { 
  print $opts_used_long;
}

##############################################
# Process the dnaorg_compare_genomes.pl output
##############################################
# variables that let us know where we are in the file
my $in_genome_section    = 0;
my $in_summary_section   = 0;
my $past_summary_section = 0;

my %cls_ngenome_H = (); # key: class name (usually integer), value number of genomes in class
my %cls_haspos_H  = (); # key: class name (usually integer), value '1' if class has >0 genes on positive strand
my %cls_hasneg_H  = (); # key: class name (usually integer), value '1' if class has >0 genes on negative strand
my $tot_ngenomes  = 0;  # total number of genomes
my @cls_A         = (); # all classes, in order
my %cdsavglen_HA  = (); # key: class name, value: [0..$i..$ncds-1] array of summed, then average, 
                        # cds lengths for this class, gene $i
my $genome_avglen = 0;  # average length of all genomes
my $tot_na        = 0;  # total number of genomes with >= 1 anomalies
my $tot_act       = 0;  # total number of anomalies

my $N4 = undef; # threshold for anomaly 4: any genome in a class with < this number of genomes is an anomaly 4

my $frc_haspos = undef; # fraction of all genomes with >= 1 gene, that have at least one gene on positive strand
my $frc_hasneg = undef; # fraction of all genomes with >= 1 gene, that have at least one gene on negative strand

# Do 2 passes over the file. On the first pass collect stats (CDS lengths, etc.)
# And on the second pass detect anomalies based on those stats.
# 
for(my $p = 1; $p <= 2; $p++) { 
  if($p == 2) { # convert our total lengths into averages
    my $n_haspos = 0; # number of genomes with >= 1 gene, that have at least one gene on positive strand
    my $n_hasneg = 0; # number of genomes with >= 1 gene, that have at least one gene on negative strand
    my $n_hasany = 0; # number of genomes with >= 1 gene
    foreach my $cls (@cls_A) { 
      my $ncds = scalar(@{$cdsavglen_HA{$cls}});
      for(my $i = 0; $i < $ncds; $i++) { 
        $cdsavglen_HA{$cls}[$i] /= $cls_ngenome_H{$cls};
        # printf("cls: $cls; gene $i; %.1f\n", $cdsavglen_HA{$cls}[$i]);
      }
      if($ncds > 0) { 
        $n_hasany += $cls_ngenome_H{$cls};
        if($cls_haspos_H{$cls}) { $n_haspos += $cls_ngenome_H{$cls}; }
        if($cls_hasneg_H{$cls}) { $n_hasneg += $cls_ngenome_H{$cls}; }
      }
    }
    if($n_hasany == 0) { die "ERROR no genomes have any annotated genes"; }
    $frc_haspos = $n_haspos / $n_hasany;
    $frc_hasneg = $n_hasneg / $n_hasany;
    
    # set thresholds, now that we know number of genomes
    $N4 = (1. - $F4) * $tot_ngenomes;

    # get average genome length (relevant for a7)
    $genome_avglen /= $tot_ngenomes;
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
        my @el_A  = split(/\s+/, $line);
        my @expA = ("#accn", "#cds", "#pos", "#neg", "#both", "#unkn", "strand-string", "cls", "tot-len", "g1");
        my $nexp = 10;
        if(scalar(@el_A) < $nexp) { die "ERROR fewer than expected columns in genome section"; }
        for($i = 1; $i <= $nexp; $i++) { 
          if($el_A[($i-1)] ne $expA[($i-1)]) { die "ERROR genome section header line, token $i is unexpected: %s ne %s", $el_A[($i-1)], $expA[($i-1)]; } 
        }
      }    
      
      # based on where we are, look for anomalies
      if($in_genome_section) { 
        # we already checked that the column headers should be 
        #accn       #cds   #pos   #neg  #both  #unkn  strand-string  cls  tot-len       g1     g2     g3     g4     g5     g6     g7     g8
        #NC_003977      7      7      0      0      0  +++++++          1     3215     2532   1203    846    681    465    639    552
        my @el_A  = split(/\s+/, $line);
        my $nel  = scalar(@el_A);
        my ($accn, $ncds, $npos, $nneg, $nboth, $nunkn, $strand_str, $cls, $totlen) = 
            ($el_A[0], $el_A[1], $el_A[2], $el_A[3], $el_A[4], $el_A[5], $el_A[6], $el_A[7], $el_A[8]);
        my @cdslen_A   = @el_A;
        my $ncds2      = $nel - $npregene;
        if($ncds == 0) { # genomes with 0 genes have no 'strand-string' token, fix that:
          $strand_str = "";
          $cls        = $el_A[6];
          $totlen     = $el_A[7];
          $ncds2++;
          @cdslen_A = (); # no CDS lengths are in @el_A
        }
        else { # remove non-CDS length values from @cdslen_A
          splice(@cdslen_A, 0, 9); 
        }

        # if we're in the first pass, gather statistics
        if($p == 1) { 
          if(! exists $cls_ngenome_H{$cls}) { 
            push(@cls_A, $cls);
            $cls_ngenome_H{$cls} = 0;
            $cls_haspos_H{$cls}  = ($strand_str =~ m/\+/) ?  1 : 0;
            $cls_hasneg_H{$cls}  = ($strand_str =~ m/\-/) ?  1 : 0;
          }
          $cls_ngenome_H{$cls}++;
          $tot_ngenomes++;
          if($ncds ne $ncds2) { die "ERROR discrepancy between stated and inferred number of CDS: $ncds ne $ncds2"; }
          # store lengths
          if(! exists $cdsavglen_HA{$cls}) { 
            @{$cdsavglen_HA{$cls}} = ();
            for($i = 0; $i < $ncds; $i++) { $cdsavglen_HA{$cls}[$i] = 0; }
          }
          for($i = 0; $i < $ncds; $i++) { 
            $cdsavglen_HA{$cls}[$i] += $el_A[($npregene+$i)];
          }
          $genome_avglen += $totlen;
        }
        
        # if we're in the second pass, look for anomalies
        else { 
          my %have_a_H = ();
          my $idx_str;
          my $offset;
          $have_a_H{"1"} = check_a1($ncds);
          $have_a_H{"2"} = check_a2($nboth);
          $have_a_H{"3"} = check_a3($nunkn);
          $have_a_H{"4"} = check_a4($N4, $cls_ngenome_H{$cls});
          $have_a_H{"5"} = 0;
          $have_a_H{"6"} = 0;
          if($cls_ngenome_H{$cls} > 0) { 
            $have_a_H{"5"} = check_a5ora6($F5and6, $frc_haspos, $cls_haspos_H{$cls});
            $have_a_H{"6"} = check_a5ora6($F5and6, $frc_hasneg, $cls_hasneg_H{$cls});
          }
          $have_a_H{"7"} = check_a7($F7, $totlen, $genome_avglen);
          $have_a_H{"8"} = 0;
          $have_a_H{"9"} = 0;
          if($ncds > 0) { 
            ($have_a_H{"8"}, $idx_str) = check_a8($F8, \@cdslen_A, \@{$cdsavglen_HA{$cls}});
            ($have_a_H{"9"}, $offset)  = check_a9(\@cdslen_A, \@{$cdsavglen_HA{$cls}});
          }

          my $na = 0;
          for($i = 1; $i <= $na_types; $i++) { 
            if($have_a_H{"$i"}) { 
              printf("%-10s  anomaly-#%d  %s", $accn, $i, "has " . $desc_H{$i}); 
              # add anomaly specific info:
              if($i == 4) { 
                printf(" (%d gene(s), order: $strand_str)", length($strand_str)); 
              }
              if($i == 7) { 
                printf(" (%d nucleotides %s)", $totlen, 
                       (($totlen < $genome_avglen) ? 
                        sprintf("< %.1f nt (%.3f * %.1f)", ((1. - $F7) * $genome_avglen), (1. - $F7), $genome_avglen) : 
                        sprintf("< %.1f nt (%.3f * %.1f)", ((1. + $F7) * $genome_avglen), (1. + $F7), $genome_avglen)));
              }
              if($i == 8) { 
                printf(" (class $cls:$strand_str, anomalous CDS lengths: $idx_str)");
              }
              if($i == 9) { 
                printf(" (class $cls:$strand_str, offset: $offset)");
              }
              printf("\n");
              $act_A[$i]++; 
              $tot_act++;
              $na++;
            }
          }
          $na_A[$na]++;
          if($na > 0) { 
            $tot_na++; 
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
printf("#\n");
printf("#anomaly-#   count  description\n");
printf("#---------  ------  -----------\n");
for($i = 1; $i <= $na_types; $i++) { 
  printf("%-10d  %6d  %s\n", $i, $act_A[$i], "genomes have " . $desc_H{"$i"});
}
printf("%-10s  %6d  %s\n", "total", $tot_na, "-");
printf("#\n");
printf("##-of-anomalies   count\n");
printf("#--------------   -----\n");
for($i = 0; $i <= $na_types; $i++) { 
  printf("%-15d  %6d\n", $i, $na_A[$i]);
}
printf("%-15s  %6d\n", "total", $tot_ngenomes);

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
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

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
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

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
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($nunkn) = (@_);
  
  return ($nunkn > 0) ? '1' : '0';
}

# Subroutine: check_a4()       
# Synopsis:   Check for anomaly 4: genome is in a class that 
#             includes less than F1 fraction of all genomes.
# Args:       $N4:                F4 * total number of genomes 
#             $ngenomes_in_class: number of genomes in this genome's class
#
# Returns:    '1' if $ngenomes_in_class >= $N1, else '0'

sub check_a4 {
  my $sub_name = "check_a4()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($N4, $ngenomes_in_class) = (@_);

  return ($ngenomes_in_class < $N4) ? '1' : '0';
}

# Subroutine: check_a5ora6()       
# Synopsis:   Check for anomaly 5 or 6: genome has >= 1 gene on a strand
#             on which > $F5 fraction of genomes have 0 genes.
# 
# Args:       $F5:      fractional threshold for anomalies 5 and 6    
#             $frc_has: fraction of genomes with genes with >= 1 on strand of interest (positive or negative)
#             $cls_has: '1' if this genome is in a class with >= 1 gene on strand of interest (positive or negative)
#
# Returns:    '1' if (($frc_has < (1. - $F5)) && $cls_has = 1), else '0'

sub check_a5ora6 {
  my $sub_name = "check_a5ora6()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($F5, $frc_has, $cls_has) = (@_);
  
  return (($frc_has < (1. - $F5)) && ($cls_has == 1)) ? '1' : '0';
}

# Subroutine: check_a7()       
# Synopsis:   Check for anomaly 7: genome length deviates from mean by more than
#                                  F7 fraction. 
# Args:       $F7:       F7 fraction
#             $cur_len:  length of current genome
#             $avg_len: mean length of all genomes
#
# Returns:    '1' if ($cur_len < ((1. - $F7) * $avg_len)) or
#                    ($cur_len < ((1. + $F7) * $avg_len))
#             else '0'

sub check_a7 {
  my $sub_name = "check_a7()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($F7, $cur_len, $avg_len) = (@_);

  if   ($cur_len < ((1. - $F7) * $avg_len)) { return 1; }
  elsif($cur_len > ((1. + $F7) * $avg_len)) { return 1; }
  else                                       { return 0; }
}

# Subroutine: check_a8()       
# Synopsis:   Check for anomaly 8: any gene length deviates from mean by more than
#                                  F8 fraction. 
# Args:       $F8:         F8 fraction
#             $cur_len_AR: ref to array of CDS lengths for current genome
#             $avg_len_AR: average length of all genes in all genomes of this class
#
# Returns:    Two values: 
#             First value: 
#               '1' if ($cur_len_AR->[$i] < ((1. - $F8) * $avg_len_AR->[$i])) or
#               '1' if ($cur_len_AR->[$i] > ((1. + $F8) * $avg_len_AR->[$i]))
#                else '0'
#             Second value:
#             "" if first returned value is '0'
#             string of ($i+1) values for which condition for first value holds
#             if more than one, each is separated by a ','.

sub check_a8 {
  my $sub_name = "check_a8()";
  my $nargs_exp = 3;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($F8, $cur_len_AR, $avg_len_AR) = (@_);

  my $n = scalar(@{$cur_len_AR}); 
  if($n != (scalar(@{$avg_len_AR}))) { die "ERROR in check_a8, number of genes in the two arrays differs"; }

  my $ret_val1 = 0;
  my $ret_val2 = "";
  for(my $i = 0; $i < $n; $i++) { 
    if($cur_len_AR->[$i] < ((1. - $F8) * $avg_len_AR->[$i])) { $ret_val1 = 1; $ret_val2 .= ($i+1) . ","; }
  }

  # remove final ',' if it exists
  $ret_val2 =~ s/\,$//;
  
  return ($ret_val1, $ret_val2);
}

# Subroutine: check_a9()       
# Synopsis:   Check for anomaly 9: shifting gene order minimizes average CDS length deviation
# Args:       $cur_len_AR: ref to array of CDS lengths for current genome
#             $avg_len_AR: average length of all genes in all genomes of this class
#
# Returns:    Two values: 
#             First value: 
#               '1' if a shift makes gene lengths closer to average
#               '0' if not
#             Second value:
#             0 if first returned value is '0'
#             $i: 1 <= $i <= scalar(@cur_len_AR); number of genes to shift by to get closest match

sub check_a9 {
  my $sub_name = "check_a9()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cur_len_AR, $avg_len_AR) = (@_);

  my $n = scalar(@{$cur_len_AR}); 
  if($n != (scalar(@{$avg_len_AR}))) { die "ERROR in check_a9, number of genes in the two arrays differs"; }

  my @copy_A = @{$cur_len_AR};

  my ($i, $j); # counters

  # set min diff as 10 times sum of all averages
  my $min_idx = -1;
  my $min_diff = 0; 
  for(my $i = 0; $i < $n; $i++) {
    $min_diff += $avg_len_AR->[$i];
  }
  $min_diff *= 10;

  for(my $i = 0; $i < $n; $i++) { 
    my $diff = 0.;
    for(my $j = 0; $j < $n; $j++) {
      $diff += abs($copy_A[$j] - $avg_len_AR->[$j]); 
    }
    if($diff < $min_diff) { $min_diff = $diff; $min_idx = $i; }

    # shift array by 1 to prepare for next comparison
    my $tmp = pop(@copy_A); # take off last element
    unshift(@copy_A, $tmp); # put it at the beginning
  }
  if($min_idx == -1) { die "ERROR in check_a9() unable to find a minimum index"; }

  my $ret_val = ($min_idx == 0) ? 0 : 1;

  return ($ret_val, $min_idx);
}
