$edir = "/panfs/pan1/dnaorg/programs";

$ctr = 0;
while($line = <>) { 
  chomp $line;
  $accn = $line;
  $accn =~ s/\.ntlist//;
  $jobname = $accn;
  $errfile = $accn . ".err";
  printf("qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n \"perl $edir/dnaorg_virus_report.pl -ozero $accn/$accn.pass-list -oanom $accn/$accn.anom-list -purge -plist $accn/$accn.purge-list $accn/$accn.compare > $accn.report\"\n");
}
