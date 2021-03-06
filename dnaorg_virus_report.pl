#!/usr/bin/env perl
# EPN, Thu Jun 4 14:45:12 2015
#
# This script uses BioEasel's SqFile module.
#
use strict;
use warnings;
use Getopt::Long;
use Time::HiRes qw(gettimeofday);
use Bio::Easel::SqFile;

# examine an output file from dnaorg_compare_genomes and find and report 'anomalies'

my $executable         = $0;
my $do_purge           = 0;
my $do_purge_anom_only = 0;
# threshold related variables
my $F4             = undef;
my $F5and6         = undef;
my $F7             = undef;
my $F8a            = undef;
my $F8b            = undef;
my $default_F4     = 0.99;
my $default_F5and6 = 0.99;
my $default_F7     = 2;
my $default_F8a    = 2;
my $default_F8b    = 0.5;

my $outfile_zero   = undef;
my $outfile_anom   = undef;
my $outfile_plist  = undef;

my $esl_test_cds_script = "/panfs/pan1/dnaorg/programs/esl-test-cds-translate-vs-fetch.pl";


my $usage  = "\ndnaorg_virus_report.pl\n";
$usage .= "\t<output from dnaorg_compare_genomes.pl>\n";
$usage .= "\n";
$usage .= " BASIC OPTIONS:\n";
$usage .= "  -purge    : remove all gene seqs from a genome with >= 1 anomaly or >= 1 CDS test failure from their CDS files created by dnaorg_compare_genomes.pl\n";
$usage .= "  -panom    : only remove all gene seqs from a genome with >= 1 anomaly, don't do the CDS tests\n";
$usage .= "\n";
$usage .= " OPTIONS CONTROLLING THRESHOLDS:\n";
$usage .= "  -F4     <f>: set threshold for anomaly 4 to <f> (gene count/strand order differs from >= <f> fraction of genomes)                      [df: $default_F4]\n";
$usage .= "  -F5and6 <f>: set threshold for anomaly 5 and 6 to <f> (class has >= 1 gene on a strand which <f> fraction of genomes do not)           [df: $default_F5and6]\n";
$usage .= "  -F7     <d>: set threshold for anomaly 7 to <f> (genome length deviates by more than <d> standard errors from mean)                    [df: $default_F7]\n";
$usage .= "  -F8a    <d>: set threshold a for anomaly 8 to <f> (> F8b fraction of CDS lengths deviate by more than <f> std errors from class mean)  [df: $default_F8a]\n";
$usage .= "  -F8b    <f>: set threshold b for anomaly 8 to <f> (> <f> fraction of CDS lengths deviate by more than F8a std errors from class mean)  [df: $default_F8b]\n";
$usage .= "\n";
$usage .= " OPTIONS FOR CREATING ADDITIONAL OUTPUT FILES:\n";
$usage .= "  -ozero <s>: save a list of all accessions with 0   anomalies to file <s>\n";
$usage .= "  -oanom <s>: save a list of all accessions with >=1 anomalies to file <s>\n";
$usage .= "  -plist <s>: with -purge, create a file named <s> with a list of the newly created fasta files with anomalous seqs purged.\n";

$usage .= "\n";

my $i; # a counter
&GetOptions( "purge"    => \$do_purge,
             "panom"    => \$do_purge_anom_only,
             "F4=s"     => \$F4, 
             "F5and6=s" => \$F5and6, 
             "F7=s"     => \$F7, 
             "F8a=s"    => \$F8a,
             "F8b=s"    => \$F8b,
             "ozero=s"  => \$outfile_zero,
             "oanom=s"  => \$outfile_anom, 
             "plist=s"  => \$outfile_plist);

if(scalar(@ARGV) != 1) { die $usage; }
my ($infile) = (@ARGV);

# store options used, so we can output them 
my $opts_used_short = "";
my $opts_used_long  = "";
if($do_purge) { 
  $opts_used_short .= "-purge ";
  $opts_used_long  .= "# option:  purging gene seqs from anomalous/test failing seqs from sequence files [-purge]\n";
}
if($do_purge_anom_only) { 
  $opts_used_short .= "-panom ";
  $opts_used_long  .= "# option:  only purge sequences that have anomalies, don't do CDS tests [-panom]\n";
}
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
if(defined $F8a) { 
  $opts_used_short .= "-F8a $F8a";
  $opts_used_long  .= "# option:  setting F8a anomaly 8 threshold to $F8a [-F8a]\n";
}
if(defined $F8a) { 
  $opts_used_short .= "-F8b $F8b";
  $opts_used_long  .= "# option:  setting F8a anomaly 8 threshold to $F8b [-F8b]\n";
}
if(defined $outfile_zero) { 
  $opts_used_short .= "-ozero $outfile_zero";
  $opts_used_long  .= "# option:  saving list of accessions with 0 anomalies to file $outfile_zero [-ozero]\n";
}
if(defined $outfile_anom) { 
  $opts_used_short .= "-oanom $outfile_anom";
  $opts_used_long  .= "# option:  saving list of accessions with >=1 anomalies to file $outfile_anom [-oanom]\n";
}
if(defined $outfile_plist) { 
  $opts_used_short .= "-pfile $outfile_plist";
  $opts_used_long  .= "# option:  saving list of new fasta files created to file $outfile_plist [-plist]\n";
}


# check for incompatible option values/combinations:
if(defined $F4 && ($F4 < 0. || $F4 > 1.0)) { 
  die "ERROR with -F4 <f>, <f> must be a number between 0 and 1."; 
}
if(defined $F5and6 && ($F5and6 < 0. || $F5and6 > 1.0)) { 
  die "ERROR with -F5and6 <f>, <f> must be a number between 0 and 1."; 
}
if(defined $F7 && ($F7 < 0.)) { 
  die "ERROR with -F7 <f>, <f> must be a positive number."; 
}
if(defined $F8a && ($F8a < 0. || $F8a > 1.0)) { 
  die "ERROR with -F8a <f>, <f> must be a positive number."; 
}
if(defined $F8b && ($F8b < 0. || $F8b > 1.0)) { 
  die "ERROR with -F8b <f>, <f> must be a number between 0 and 1."; 
}
if(defined $outfile_plist && (! $do_purge)) { 
  die "ERROR with -plist only makes sense in combination with -purge"; 
}  
if($do_purge_anom_only && (! $do_purge)) { 
  die "ERROR with -panom only makes sense in combination with -purge"; 
}  

# set defaults for variables not set by the user via options
if(! defined $F4)     { $F4     = $default_F4;     }
if(! defined $F5and6) { $F5and6 = $default_F5and6; }
if(! defined $F7)     { $F7     = $default_F7;     }
if(! defined $F8a)    { $F8a    = $default_F8a;    }
if(! defined $F8b)    { $F8b    = $default_F8b;    }

# define anomaly descriptions:
my %desc_H = ();
$desc_H{"1"} = "0 CDS annotations.";
$desc_H{"2"} = "at least one gene with coding sequence on both strands";
$desc_H{"3"} = "at least one gene on an unknown strand";
$desc_H{"4"} = sprintf("a CDS order/strand string that differs from >= %.3f fraction of genomes", $F4); 
$desc_H{"5"} = sprintf("a CDS on the positive strand, when >= %.3f fraction of all genomes do not", $F5and6);
$desc_H{"6"} = sprintf("a CDS on the negative strand, when >= %.3f fraction of all genomes do not", $F5and6);
$desc_H{"7"} = sprintf("total length that deviates by more than %d standard errors from the mean", $F7);
$desc_H{"8"} = sprintf("> %.3f fraction of CDS lengths that deviate by more than %d standard errors from the mean for that class", $F8b, $F8a);
#$desc_H{"9"} = sprintf("shifting gene order descreases summed CDS length deviation from mean");

my $na_types = 8;
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

printf("#\n"); 
printf("###############################################################\n"); 
printf("# Summary of class definition and membership from compare file:\n"); 
printf("#\n"); 


if(defined $outfile_zero) { 
  open(OUTZERO, ">" . $outfile_zero) || die "ERROR unable to open $outfile_zero for writing"; 
}
if(defined $outfile_zero) { 
  open(OUTANOM, ">" . $outfile_anom) || die "ERROR unable to open $outfile_anom for writing"; 
}
if(defined $outfile_plist) { 
  open(OUTPFILE, ">" . $outfile_plist) || die "ERROR unable to open $outfile_plist for writing"; 
}

##############################################
# Process the dnaorg_compare_genomes.pl output
##############################################
# variables that let us know where we are in the file
my $in_genome_section  = 0;
my $in_summary_section = 0;
my $in_fetch_section   = 0;

my %cls_ngenome_H     = (); # key: class name (usually integer), value number of genomes in class
my %cls_haspos_H      = (); # key: class name (usually integer), value '1' if class has >0 genes on positive strand
my %cls_hasneg_H      = (); # key: class name (usually integer), value '1' if class has >0 genes on negative strand
my %cls_sstr_H        = (); # key: class name (usually integer), value strand string for that class
my %cls_ncds_H        = (); # key: class name (usually integer), number of CDS annotations for that class
my $tot_ngenomes      = 0;  # total number of genomes
my @cls_A             = (); # all classes, in order
my %cds_len_HAA       = (); # key: class name $cls, value: [0..$i..$cls_cds_H{$cls}-1] array of $cls_ngenome_H{$cls}
                            # lengths for class $cls, cds $i
my %cds_len_mean_HA   = (); # key: class name, value: [0..$i..$cls_cds_H{$cls}-1] mean cds length
                            # for this class, cds $i
my %cds_len_stderr_HA = (); # key: class name, value: [0..$i..$cls_cds_H{$cls}-1] standard error in cds length
                            # for this class, cds $i
my @genome_len_A      = (); # [0..$i..$tot_ngenomes-1] length of genome $i
my $genome_len_mean   = 0;  # mean length of all genomes
my $genome_len_stderr = 0;  # standard error in length of all genomes
my $tot_na            = 0;  # total number of genomes with >= 1 anomalies
my $tot_act           = 0;  # total number of anomalies

my $N4 = undef; # threshold for anomaly 4: any genome in a class with < this number of genomes is an anomaly 4

my $frc_haspos = undef; # fraction of all genomes with >= 1 gene, that have at least one gene on positive strand
my $frc_hasneg = undef; # fraction of all genomes with >= 1 gene, that have at least one gene on negative strand

my @out_A               = (); # the output lines that list the anomalies
my @out_extra_info_A    = (); # the extra info substrings of the out_A lines, used along with extra_info_ct_H() for determining singleton's
my @out_anomaly_A       = (); # the anomaly index of the out_A lines
my @out_class_A         = (); # the class index of the out_A lines
my %extra_info_ct_H     = (); # key: an extra info string, value: number of times we've seen that string.

my $a5_is_possible = 0; # set to '1' if anomaly a5 can be observed
my $a6_is_possible = 0; # set to '1' if anomaly a6 can be observed

# variables related to -purge
my %to_purge_H = (); # key: accession; value: 1 to purge this sequence
my @purged1_file_A   = (); # array of file names of round 1 purged fasta files created
my @purged1_output_A = (); # array of strings to output describing round 1 purged fasta files created or not created (because all seqs were purged) 
my @purged2_file_A   = (); # array of file names of round 2 purged fasta files created
my @purged2_output_A = (); # array of strings to output describing round 2 purged fasta files created or not created (because all seqs were purged) 


# Do 2 passes over the file. On the first pass collect stats (CDS lengths, etc.)
# And on the second pass detect anomalies based on those stats.
# 
for(my $p = 1; $p <= 2; $p++) { 
  if($p == 2) { # convert our total lengths into means
    my $n_haspos = 0; # number of genomes with >= 1 gene, that have at least one gene on positive strand
    my $n_hasneg = 0; # number of genomes with >= 1 gene, that have at least one gene on negative strand
    my $n_hasany = 0; # number of genomes with >= 1 gene
    foreach my $cls (@cls_A) { 
      @{$cds_len_mean_HA{$cls}} = ();
      my $ncds = scalar(@{$cds_len_HAA{$cls}});
      for(my $i = 0; $i < $ncds; $i++) { 
        $cds_len_mean_HA{$cls}[$i]   = mean(\@{$cds_len_HAA{$cls}[$i]});
        $cds_len_stderr_HA{$cls}[$i] = standard_error(\@{$cds_len_HAA{$cls}[$i]});
        #printf("cls: $cls; gene $i; %.1f\n", $cds_len_mean_HA{$cls}[$i]);
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

    # get mean and std err for genome length (relevant for a7)
    $genome_len_mean   = mean(\@genome_len_A);
    $genome_len_stderr = standard_error(\@genome_len_A);
    # printf("genome_len_mean:   $genome_len_mean\n");
    # printf("genome_len_stderr: $genome_len_stderr\n");

    # determine if anomalies a5 and a6 are possible:
    # a5 is only possible if >= F5and6 fraction genomes lack any gene on positive strand
    $a5_is_possible = ($frc_haspos < (1. - $F5and6)) ? 1 : 0;
    # a6 is only possible if >= F5and6 fraction genomes lack any gene on negative strand
    $a6_is_possible = ($frc_hasneg < (1. - $F5and6)) ? 1 : 0;
  }
  
  open(IN, $infile) || die "ERROR unable to open $infile on pass $p";
  while(my $line = <IN>) { 
    if($line =~ m/\w/) { 
      # potentially update variables which tell us where we are
      if((! $in_genome_section) && (! $in_summary_section)) { 
        if($line !~ m/^\#/) { 
          $in_genome_section  = 1;
          $in_summary_section = 0;
          $in_fetch_section   = 0;
        }
      }
      elsif($in_genome_section) { 
        if($line =~ /^\# Number-of-classes\:\s+(\d+)/) { 
          $in_genome_section  = 0;
          $in_summary_section = 1;
          $in_fetch_section   = 0;
        }
      }
      elsif($in_summary_section) { 
        if($line =~ m/^# Fetching/) { 
          $in_genome_section  = 0;
          $in_summary_section = 0;
          $in_fetch_section   = 1;
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
        my @cur_cds_len_A   = @el_A;
        my $ncds2      = $nel - $npregene;
        if($ncds == 0) { # genomes with 0 genes have no 'strand-string' token, fix that:
          $strand_str = "";
          $cls        = $el_A[6];
          $totlen     = $el_A[7];
          $ncds2++;
          @cur_cds_len_A = (); # no CDS lengths are in @el_A
        }
        else { # remove non-CDS length values from @cur_cds_len_A
          splice(@cur_cds_len_A, 0, 9); 
        }

        # if we're in the first pass, gather statistics
        if($p == 1) { 
          if(! exists $cls_ngenome_H{$cls}) { 
            push(@cls_A, $cls);
            $cls_ngenome_H{$cls} = 0;
            $cls_sstr_H{$cls}    = $strand_str;
            $cls_ncds_H{$cls}    = length($strand_str);
            $cls_haspos_H{$cls}  = ($strand_str =~ m/\+/) ?  1 : 0;
            $cls_hasneg_H{$cls}  = ($strand_str =~ m/\-/) ?  1 : 0;
          }
          $cls_ngenome_H{$cls}++;
          $genome_len_A[$tot_ngenomes] = $totlen;
          $tot_ngenomes++;
          if($ncds ne $ncds2) { die "ERROR discrepancy between stated and inferred number of CDS: $ncds ne $ncds2"; }
          # store lengths
          if(! exists $cds_len_HAA{$cls}) { 
            # initialize both array dimensions of this hash of 2D arrays
            @{$cds_len_HAA{$cls}} = ();
            for($i = 0; $i < $ncds; $i++) { @{$cds_len_HAA{$cls}[$i]} = (); }
          }
          for($i = 0; $i < $ncds; $i++) { 
            push(@{$cds_len_HAA{$cls}[$i]}, $el_A[($npregene+$i)]);
          }
        }
        
        # if we're in the second pass, look for anomalies
        else { 
          my %have_a_H = ();
          my $ndev = 0;
          my $dev_str = "";
          # my $offset;
          $have_a_H{"1"} = check_a1($ncds);
          $have_a_H{"2"} = check_a2($nboth);
          $have_a_H{"3"} = check_a3($nunkn);
          $have_a_H{"4"} = check_a4($N4, $cls_ngenome_H{$cls});
          $have_a_H{"5"} = 0;
          $have_a_H{"6"} = 0;
          if($cls_ngenome_H{$cls} > 0) { 
            if($a5_is_possible) { 
              $have_a_H{"5"} = check_a5ora6($F5and6, $frc_haspos, $cls_haspos_H{$cls});
            }
            if($a6_is_possible) { 
              $have_a_H{"6"} = check_a5ora6($F5and6, $frc_hasneg, $cls_hasneg_H{$cls});
            }
          }
          $have_a_H{"7"} = check_a7($F7, $totlen, $genome_len_mean, $genome_len_stderr);
          $have_a_H{"8"} = 0;
          $have_a_H{"9"} = 0;
          if($ncds > 0) { 
            ($have_a_H{"8"}, $ndev, $dev_str) = check_a8($F8a, $F8b, \@cur_cds_len_A, \@{$cds_len_mean_HA{$cls}}, \@{$cds_len_stderr_HA{$cls}});
            # ($have_a_H{"9"}, $other_cls, $cds_offset)  = check_for_cds_shifts_a9(\@cur_cds_len_A, \%cds_len_mean_HA, \%cds_sstr_H, $cls);
          }
          
          my $na = 0;
          for($i = 1; $i <= $na_types; $i++) { 
            if($have_a_H{"$i"}) { 
              my $outline = sprintf("%-10s  anomaly-#%d  %9s  %s", $accn, $i, "", "has " . $desc_H{$i}); 
              # add anomaly specific info:
              my $extra_info = "";
              if($i == 4 || $i == 5 || $i == 6) { 
                $extra_info = sprintf(" (class $cls:%d gene(s), order: $strand_str)", $cls_ncds_H{$cls});
              }
              if($i == 7) { 
                $extra_info = sprintf(" (%d nucleotides %s)", $totlen, 
                                    (($totlen < $genome_len_mean) ? 
                                     sprintf("< %.1f nt (%.1f - (%d * %.3f))", ($genome_len_mean - ($F7 * $genome_len_stderr)), $genome_len_mean, $F7, $genome_len_stderr) :
                                     sprintf("> %.1f nt (%.1f + (%d * %.3f))", ($genome_len_mean + ($F7 * $genome_len_stderr)), $genome_len_mean, $F7, $genome_len_stderr)));
              }
              if($i == 8) { 
                $extra_info = " (class $cls:$strand_str, $ndev anomalous CDS lengths: $dev_str)";
              }
              if(! exists $extra_info_ct_H{$extra_info}) { 
                $extra_info_ct_H{$extra_info} = 1; 
              }
              else { 
                $extra_info_ct_H{$extra_info}++;
              }
              $outline .= $extra_info . "\n";
              $act_A[$i]++; 
              $tot_act++;
              $na++;
              push(@out_A, $outline);
              push(@out_extra_info_A, $extra_info); # this is necessary only so we can determine which have unique extra_info substrings in their $outline,
                                                    # these are 'SINGLETON's.
              push(@out_class_A, $cls);
              push(@out_anomaly_A, $i);
            }
          }
          if($na > 0) { # add a blank line in between each genome
            push(@out_A, "\n");
            push(@out_extra_info_A, "");
            push(@out_anomaly_A, -1);
            push(@out_class_A, -1);
            if($do_purge) { # store this accession, we'll need it to purge anomalous seqs later
              $to_purge_H{$accn} = 1;
            }
            if(defined $outfile_anom) { 
              print OUTANOM "$accn\n"; 
            }
          }              
          else { # no anomalies for this accession
            if(defined $outfile_zero) { 
              print OUTZERO "$accn\n"; 
            }
          }

          $na_A[$na]++;
          if($na > 0) { 
            $tot_na++; 
          }
        }
      }
      elsif($in_summary_section) { 
        # echo this part 
        if($p == 1) { 
          print $line; 
        }
      }
      ########################################
      # Possibly purge sequences from the fasta files
      elsif($in_fetch_section && $do_purge) { 
        if($p == 2) { # only do this on the second pass
          my $fa_file = $line;
          chomp $fa_file;
          my $have_fa_file = 0;
          # Fetching 315 full genome sequences for class  1 ... done. [Hepatitis-C_r9.NC_004102/Hepatitis-C_r9.NC_004102.c1.fg.fa]
          # Fetching 425 CDS sequences for class  1 gene  1 ... done. [FMDV_r26.NC_004004/FMDV_r26.NC_004004.c1.g1.fa]
          if($fa_file =~ s/^\# Fetching.+CDS sequences.+done.\s+\[//) { 
            $fa_file =~ s/\]//;
            $have_fa_file = 1;
          }
          elsif($fa_file =~ s/^\# Fetching.+full genome sequences for.+done.\s+\[//) { 
            $fa_file =~ s/\]//;
            $have_fa_file = 1;
          }
          if($have_fa_file) { 
            my $output_fa_file    = $fa_file; 
            my $output_plist_file = $fa_file; 
            if($output_fa_file !~ m/\.fa$/) { die "ERROR $fa_file does not end in .fa"; }
            $output_fa_file    =~ s/\.fa$/.purged1.fa/;
            $output_plist_file =~ s/\.fa$/.purged1.list/;
            my ($nkept, $npurged) = purgeSeqs($fa_file, $output_fa_file, $output_plist_file, \%to_purge_H, \@purged1_file_A, \@purged1_output_A);
            if($nkept > 0) { 
              if(defined $outfile_plist) { 
                print OUTPFILE "$output_fa_file\n"; 
              }
            }
          }
        }
      }
    }
  }
  close(IN);
}

#############################
# If we're purging sequences, perform round 2 of the purging:
#   Go back over the round 1 purged fasta files and run $esl_test_cds_script on them.
#   Then parse its output and purge any sequences that failed >= 1 test to create
#   the round 2 purged fasta files.
#   
#   If -panom was used, we skip this round 2 purging.
#
if($do_purge && (! $do_purge_anom_only)) { 
  %to_purge_H = (); # reinitialize this hash, we'll fill it with names of sequences to remove in round 2
  foreach my $fa_file (@purged1_file_A) { 
    if($fa_file !~ m/\.fg\.purged1\./) { # skip the full genomes, we can't do CDS tests on those
      my $cds_test_output = $fa_file;
      $cds_test_output =~ s/\.fa$/\.cds-test/;
      my $cmd = "perl $esl_test_cds_script -incompare -subset $fa_file > $cds_test_output";
      runCommand($cmd, 0);
      
      # now parse the output
      open(IN, $cds_test_output) || die "ERROR unable to open freshly created file $cds_test_output for reading";
      while(my $line = <IN>) { 
        chomp $line;
        ##protein-accession  nt-accession  mincoord  maxcoord  tr-len  nexon  str  incomplete?  start   num-Ns  num-oth     T1   T2   T3   T4   T5   T6   T7   T8   T9  pass/fail
        #CAD62370.1          AJ539138          1092      8090    6999      1    +           no      1        1        1      0    0    0    0    0    0    0    0    0    pass     
        if($line !~ m/^\#/ && $line =~ m/\s+fail/) { 
          if($line =~ /\S+\s+(\S+)/) { 
            $to_purge_H{$1} = 1;
          }
          else { 
            die "ERROR unable to parse esl_test_cds_translate_vs_fetch.pl output line $line"; 
          }
        }
      }
    } # end of 'if($fa_file !~ m/\.fg\.purged1\./)'
  } # end of 'foreach my $fa_file (@purged1_file_A)'
  # now we have a full list of sequences to purge in %to_purge_H
  # go back through each file an remove any sequences in %to_purge_H
  foreach my $fa_file (@purged1_file_A) { 
    my $output_fa_file = $fa_file;
    if($output_fa_file =~ s/\.purged1\./.purged2./) { 
      my $output_plist_file = $output_fa_file;
      $output_plist_file =~ s/\.fa/\.list/;
      my ($nkept, $npurged) = purgeSeqs($fa_file, $output_fa_file, $output_plist_file, \%to_purge_H, \@purged2_file_A, \@purged2_output_A);
    }
    else { die "ERROR round 2 purging $fa_file file name unparseable"; }
  }
}


# go back and identify singletons and modify output lines for singletons
my @act_singleton_A = ();
my $tot_nsingleton = 0;
for($i = 1; $i <= $na_types; $i++) { 
  $act_singleton_A[$i] = 0;
}
for(my $l = 0; $l < scalar(@out_A); $l++) { 
  my $line       = $out_A[$l];
  my $extra_info = $out_extra_info_A[$l];
  my $anomaly    = $out_anomaly_A[$l];
  my $class      = $out_class_A[$l];
  my $is_singleton = 0; # changed to '1' below if we determine this is a singleton
  if($anomaly == 7) { # special case, can't use extra_info_ct to determine if it's a singleton
    if($cls_ngenome_H{$class} == 1) { 
      $is_singleton = 1;
    }
  }
  elsif($extra_info =~ m/\w/ && $extra_info_ct_H{$extra_info} == 1) { 
    $is_singleton = 1;
  }
  if($is_singleton) { 
    $act_singleton_A[$anomaly]++;
    $tot_nsingleton++;
    $line =~ s/             has /  SINGLETON  has /;
  }
  $out_A[$l] = $line;
}


# print summary table that lists counts of each anomaly
printf("#\n");
printf("###############################################################\n"); 
printf("# Description and number of occurences of each anomaly:\n");
printf("#\n");
printf("#   SINGLETONs are numbered in the '(sngl)' column and\n");
printf("#   explained at end of this output.\n");
printf("#\n");
printf("#               count\n");
printf("#           -------------\n");
printf("#anomaly-#   total (sngl)  description\n");
printf("#---------  -------------  -----------\n");
for($i = 1; $i <= $na_types; $i++) { 
  # one exception here
  if(($i == 5 && (! $a5_is_possible)) || 
     ($i == 6 && (! $a6_is_possible))) { 
    printf("%-10d  %6s %6s  %s\n", $i, "N/A", "", "genomes have " . $desc_H{"$i"});
  }
  else { 
    printf("%-10d  %6d %6s  %s\n", $i, $act_A[$i], sprintf("(%d)", $act_singleton_A[$i]), "genomes have " . $desc_H{"$i"});
  }
}
printf("%-10s  %6d %6s  anomalies exist in %d genomes of %d (%6.4f)\n", "total", $tot_act, "", $tot_na, $tot_ngenomes, $tot_na / $tot_ngenomes); 
printf("#\n");
printf("%-10s  %6d %6s  %s\n", "none", $tot_ngenomes - $tot_na, "", sprintf("genomes have zero anomalies (%6.4f)%s", ($tot_ngenomes - $tot_na) / $tot_ngenomes, (defined $outfile_zero) ? " [saved to file: $outfile_zero]" : ""));
printf("#\n");


# print summary table that lists counts of genomes that have N anomalies
#printf("##-of-anomalies   count  fraction\n");
#printf("#--------------   -----  --------\n");
#for($i = 0; $i <= $na_types; $i++) { 
#  printf("%-15d  %6d    %6.4f\n", $i, $na_A[$i], $na_A[$i] / $tot_ngenomes);
#}
#printf("%-15s  %6d    %6.4f\n", "total", $tot_ngenomes, 1.0);

# print individual anomalies
# (we already added SINGLETON annotation above)
printf("#\n");
printf("###############################################################\n"); 
printf("# Individual anomalies:\n");
printf("#\n");
foreach my $line (@out_A) { 
  print $line;
}

# print tail, explanation of SINGLETONs:

printf("###############################################################\n"); 
print("#\n");
print("# Explanation of SINGLETON (sngl) annotation:\n");
print("# A4 SINGLETONs are the only genomes in their class with anomaly 4\n");
print("# A5 SINGLETONs are the only genomes in their class with anomaly 5\n");
print("# A6 SINGLETONs are the only genomes in their class with anomaly 6\n");
print("# A7 SINGLETONs are the only genomes in their class with anomaly 7\n");
print("# A8 SINGLETONs are the only genomes in their class with anomaly 8 and their particular set of CDS indices that deviate in length from the mean\n");
print("#\n");
print("# Explanation of class definition:\n");
print("# All genomes with the same \"strand-string\" are placed in the same class.\n");
print("# A genome's strand-string is a string specifying the strand and order of all CDS annotation for that genome.\n");
print("# For example, a genome with 3 CDS, from 1..99, 300..398, and 600..698 all on the positive strand would have a strand string of \"+++\".\n");
print("# The same genome but with CDS #2 on the negative strand would have a strand string of \"+-+\".\n");
print("#\n");
print("###############################################################\n"); 

if(scalar(@purged1_output_A) > 0) { 
  printf("#\n");
  printf("# Newly created fasta files after purging sequences with >= 1 anomaly.\n");
  foreach my $line (@purged1_output_A) { 
    print $line; 
  }
  printf("#\n");
  print("###############################################################\n"); 
}
if(scalar(@purged2_output_A) > 0) { 
  printf("#\n");
  printf("# Newly created fasta files after purging sequences that had 0 anomalies but failed >= 1 CDS tests.\n");
  foreach my $line (@purged2_output_A) { 
    print $line; 
  }
  printf("#\n");
  print("###############################################################\n"); 
}

my $n = 0;
if(defined $outfile_zero) { 
  close(OUTZERO); 
  if($n == 0) { print("#\n"); }
  $n++; 
  printf("# List of accessions with zero anomalies saved to $outfile_zero.\n");
}
if(defined $outfile_anom) { 
  close(OUTANOM); 
  if($n == 0) { print("#\n"); }
  $n++; 
  printf("# List of accessions with >= 1 anomalies saved to $outfile_anom.\n");
}
if(defined $outfile_plist) { 
  close(OUTPFILE); 
  if($n == 0) { print("#\n"); }
  $n++; 
  printf("# List of newly created fasta files with anomalous sequences purged saved to $outfile_plist\n"); 
}
if($n > 0) { print("#\n###############################################################\n"); }


#############
# SUBROUTINES
#############
# Subroutine: mean()
# Synopsis:   Given an array of values, return the mean.
# Args:       $AR: ref to an array of the values
#
# Returns:    mean: sum(x) / n

sub mean { 
  my $sub_name = "mean()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($AR) = (@_);

  my $n = scalar(@{$AR});
  my $mean = 0.;
  for(my $i = 0; $i < $n; $i++) { 
    $mean += $AR->[$i];
  }
  $mean /= $n;

  return $mean;
}
# Subroutine: estimated_variance()
# Synopsis:   Given an array of values, return the estimated variance.
# Args:       $AR: ref to an array of the values
#
# Returns:    estimated variance: sum{(x-u)^2} / n

sub estimated_variance {
  my $sub_name = "estimated_variance()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($AR) = (@_);

  my $mean = mean($AR);

  my $sum = 0.;
  my $n = scalar(@{$AR});
  for(my $i = 0; $i < $n; $i++) { 
    $sum += ($AR->[$i] - $mean) * ($AR->[$i] - $mean);
  }
  
  return $sum / $n;
}

# Subroutine: standard_error()
# Synopsis:   Given an array of values, return the standard error.
# Args:       $AR: ref to an array of the values
#
# Returns:    standard error: sqrt of estimated variance: sqrt(sum{(x-u)^2} / n)

sub standard_error { 
  my $sub_name = "standard_error()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($AR) = (@_);

  return sqrt(estimated_variance($AR));
}

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
#             F7 standard errors.
# Args:       $F7:        F7 fraction
#             $cur_len:   length of current genome
#             $mean_len:  mean length of all genomes
#             $stderr:    standard error in genome length distribution
#
# Returns:    '1' if ($cur_len < ((1. - $F7) * $mean_len)) or
#                    ($cur_len < ((1. + $F7) * $mean_len))
#             else '0'

sub check_a7 {
  my $sub_name = "check_a7()";
  my $nargs_exp = 4;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($F7, $cur_len, $mean_len, $stderr) = (@_);

  if   ($cur_len < ($mean_len - ($F7 * $stderr))) { return 1; }
  elsif($cur_len > ($mean_len + ($F7 * $stderr))) { return 1; }
  else                                      { return 0; }
}

# Subroutine: check_a8()       
# Synopsis:   Check for anomaly 8: more than F8b fraction of CDS' deviate from mean
#             by more than F8a standard errors from the mean length.
# Args:       $F8a:         F8a standard error threshold
#             $F8b:         F8b fraction threshold
#             $cur_len_AR:  ref to array of CDS lengths for current genome
#             $mean_len_AR: ref to array of mean length of CDS lengths
#             $stderr_AR:   ref to array of std err in CDS length distributions
#
# Returns:    Three return values: 
#             First return value: '1' if more than F8b fraction of CDS' deviate 
#             from mean by more than F8a standard errors from the mean length.
#             else '0'
#             Second return value:
#             "" if first returned value is '0'
#             string of ($i+1) values for each CDS $i for which length deviates
#             from mean by more than F8b standard errors.

sub check_a8 {
  my $sub_name = "check_a8()";
  my $nargs_exp = 5;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($F8a, $F8b, $cur_len_AR, $mean_len_AR, $stderr_AR) = (@_);

  my $n = scalar(@{$cur_len_AR}); 
  if($n != (scalar(@{$mean_len_AR}))) { die "ERROR in check_a8, number of genes in the two arrays differs"; }

  my $ndeviates    = 0;
  my $deviates_str = "";

  for(my $i = 0; $i < $n; $i++) { 
    #printf("mean_len_AR->[$i]: $mean_len_AR->[$i]; stderr: $stderr_AR->[$i]; cur: $cur_len_AR->[$i]\n");
    if(($cur_len_AR->[$i] < ($mean_len_AR->[$i] - ($F8a * $stderr_AR->[$i]))) ||  # too low
       ($cur_len_AR->[$i] > ($mean_len_AR->[$i] + ($F8a * $stderr_AR->[$i])))) {  # too high
      $ndeviates++;
      $deviates_str .= ($i+1) . ","; 
    }
  }

  # remove final ',' if it exists
  $deviates_str =~ s/\,$//;

  my $ret_val1 = ($ndeviates > ($F8b * $n)) ? 1 : 0;
  my $ret_val2 = $ndeviates;
  my $ret_val3 = ($ndeviates > ($F8b * $n)) ? $deviates_str : "";
  
  return ($ret_val1, $ret_val2, $ret_val3);
}

# Subroutine: check_a9()       
# Synopsis:   Check for anomaly 9: shifting gene order minimizes mean CDS length deviation
# Args:       $cur_len_AR: ref to array of CDS lengths for current genome
#             $mean_len_AR: mean length of all genes in all genomes of this class
#
# Returns:    Two values: 
#             First value: 
#               '1' if a shift makes gene lengths closer to mean
#               '0' if not
#             Second value:
#             0 if first returned value is '0'
#             $i: 1 <= $i <= scalar(@cur_len_AR); number of genes to shift by to get closest match

sub check_a9 {
  my $sub_name = "check_a9()";
  my $nargs_exp = 2;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }

  my ($cur_len_AR, $mean_len_AR) = (@_);

  my $n = scalar(@{$cur_len_AR}); 
  if($n != (scalar(@{$mean_len_AR}))) { die "ERROR in check_a9, number of genes in the two arrays differs"; }

  my @copy_A = @{$cur_len_AR};

  my ($i, $j); # counters

  # set min diff as 10 times sum of all means
  my $min_idx = -1;
  my $min_diff = 0; 
  for(my $i = 0; $i < $n; $i++) {
    $min_diff += $mean_len_AR->[$i];
  }
  $min_diff *= 10;

  for(my $i = 0; $i < $n; $i++) { 
    my $diff = 0.;
    for(my $j = 0; $j < $n; $j++) {
      $diff += abs($copy_A[$j] - $mean_len_AR->[$j]); 
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

# Subroutine: runCommand()
# Args:       $cmd:            command to run, with a "system" command;
#             $be_verbose:     '1' to output command to stdout before we run it, '0' not to
#
# Returns:    amount of time the command took, in seconds
# Dies:       if $cmd fails

sub runCommand {
  my $sub_name = "runCommand()";
  my $nargs_exp = 2;

  my ($cmd, $be_verbose) = @_;

  if($be_verbose) { 
    print ("Running cmd: $cmd\n"); 
  }

  my ($seconds, $microseconds) = gettimeofday();
  my $start_time = ($seconds + ($microseconds / 1000000.));
  system($cmd);
  ($seconds, $microseconds) = gettimeofday();
  my $stop_time = ($seconds + ($microseconds / 1000000.));

  if($? != 0) { die "ERROR command failed:\n$cmd\n"; }

  return ($stop_time - $start_time);
}

# Subroutine: purgeSeqs()
# Synopsis:   Given an input sequence file and a hash of sequence names that may or may not
#             be in the input file, create a new sequence file with all sequences from 
#             the input file except those that have names that exist as keys in the hash.
#             If $out_file_AR is defined push the output file (if it has >= 1 seqs)
#             to @{$out_file_AR} and remove it if it has 0 seqs. 
#             If $out_text_AR is defined push an informative message about the file to that
#             array.
#
# Args:       $input_fa_file:     fasta file to read, we'll output some of these seqs to a new file
#             $output_fa_file:    fasta file to create, with a subset of the seqs in $input_fa_file
#             $output_plist_file: text file to create with a list of sequences removed
#             $to_purge_HR:       ref to hash of sequences to remove (purge) from $input_fa_file
#                                 key: <seqname>, value: '1' (or any defined value).
#             $out_file_AR:       ref to array of output file names to add to if we create a non-empty file here
#             $out_text_AR:       ref to array of text describing the output file we created
#                                 (or deleted if if has 0 seqs).
# Returns:    2 values:
#               $nkept:   number of sequences from $input_fa_file written to $output_fa_file
#               $npurged: number of sequences removed from $input_fa_file
# Dies:       if we have a problem reading the file $input_fa_file

sub purgeSeqs {
  my $sub_name = "purgeSeqs()";
  my $nargs_exp = 3;
  
  my ($input_fa_file, $output_fa_file, $output_plist_file, $to_purge_HR, $out_file_AR, $out_text_AR) = @_;

  my $nkept   = 0;
  my $npurged = 0;

  open(OUTFA, ">" . $output_fa_file) || die "ERROR unable to open $output_fa_file for writing"; 
  open(OUTPL, ">" . $output_plist_file) || die "ERROR unable to open $output_plist_file for writing"; 
  my $sqfile = Bio::Easel::SqFile->new({ fileLocation => $input_fa_file });
  my $nseq = $sqfile->nseq_ssi;
  # check each sequence to see if it should be purged, if so, purge it
  for(my $i = 1; $i <= $nseq; $i++) { 
    my $seq = $sqfile->fetch_consecutive_seqs(1, "", 80, undef); # get the next sequence in the file
    if($seq =~ /^\>(\S+)\s+/) { 
      my $name = $1;
      # strip the version from the name, which *may* be in "accession.version" format, e.g.: KF693782.1
      my $name_accn_only = $name;
      stripVersion(\$name_accn_only);
      if(! exists $to_purge_H{$name_accn_only}) { 
        print OUTFA $seq; # if it is not in the purge list, output it to OUTFA
        $nkept++;
      }
      else { 
        print OUTPL $name . "\n";
        $npurged++; 
      }
    }
    else { 
      die "ERROR unable to determine name for sequence $i in file $input_fa_file";
    }
  }
  $sqfile->close_sqfile;
  if(-e $input_fa_file . ".ssi") { unlink $input_fa_file . ".ssi"; } # clean up the index file
  close(OUTFA);
  close(OUTPL);


  if($nkept == 0) { 
    if($npurged == 0) { die "ERROR didn't read any sequences in file $input_fa_file"; }
    if(defined $out_text_AR) { 
      push(@{$out_text_AR}, sprintf("# All %4d sequence(s) in $input_fa_file were purged [no output fasta file created]\n", $npurged));  
    }
    if(-e $output_fa_file) { 
      if(-s $output_fa_file) { die "ERROR didn't write any sequences to $output_fa_file, but it is not empty"; }
      unlink $output_fa_file;
    }
  }
  else { # we did output >= 1 sequences to $output_fa_file
    if(defined $out_text_AR) { 
      push(@{$out_text_AR}, sprintf("# Created file $output_fa_file [%4d sequences kept; %4d sequences purged%s]\n", $nkept, $npurged, (defined $output_plist_file) ? ", listed in $output_plist_file" : ""));
    }
    if(defined $out_file_AR) { 
      push(@{$out_file_AR}, $output_fa_file);
    }
  }
  return ($nkept, $npurged);
}

# Subroutine: stripVersion()
# Purpose:    Given a ref to an accession.version string, remove the version.
# Args:       $accver_R: ref to accession version string
# Returns:    Nothing, $$accver_R has version removed
sub stripVersion {
  my $sub_name  = "stripVersion()";
  my $nargs_exp = 1;
  if(scalar(@_) != $nargs_exp) { die "ERROR $sub_name entered with wrong number of input args"; }
  
  my ($accver_R) = (@_);

  $$accver_R =~ s/\.[0-9]*$//; # strip version

  return;
}
