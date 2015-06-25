$edir = "/panfs/pan1/dnaorg/programs";

$ctr = 0;
while($line = <>) { 
#NC_003977 Hepatitis-B_r1
  chomp $line;
  ($accn, $name) = split(/\s+/, $line);
  $name_accn = $name . "." . $accn;
  $jobname = $name_accn;
  $errfile = $name_accn . ".report.err";
  printf("qsub -N $jobname -b y -v SGE_FACILITIES -P unified -S /bin/bash -cwd -V -j n -o /dev/null -e $errfile -m n \"perl $edir/dnaorg_virus_report.pl -ozero $name_accn/$name_accn.pass-list -oanom $name_accn/$name_accn.amom-list -purge -plist $name_accn/$name_accn.purge-list $name_accn/$name_accn.compare > $name_accn.report\"\n");
}
