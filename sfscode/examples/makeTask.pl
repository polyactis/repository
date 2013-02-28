#!/usr/bin/perl -w
use strict;

my $NSIMS = 1000;   #TOTAL NUMBER OF SIMULATIONS YOU WANT TO RUN
my $NJOBS = 10;     #TOTAL NUMBER OF JOBS (PROCESSORS) YOU WANT TO USE
my $FOLD = 0;       #START IN $FOLD+1
my $THETA = 0.0015; #SET NON-VARYING PARAMETERS

# SET PARAMETERS YOU WOULD LIKE TO VARY IN AN ARRAY AS FOLLOWS:
my @RHO = (0, 0.001);

my %SEEDS = ();
#NOW LOOP OVER VALUES OF ALL PARAMETERS
foreach my $r (@RHO){
  #YOU CAN INCLUDE MULTIPLE foreach LOOPS IF VARYING MORE THAN 1 PARAMETER
  $FOLD++;
  unless(chdir $FOLD){
    mkdir $FOLD;
    chdir $FOLD;
  }
  print "$FOLD\n";

  #DEFAULT OUTPUT FILE IS "tasklist"
  open(OUT,">tasklist") or die "cannot write to tasklist\n";

  # UPDATE DIRECTORY TO INDICATE THE LOCATION OF SFS_CODE AND UPDATE COMMAND LINE AS NECESSARY FOR YOUR SET OF SIMULATIONS
  my $cmd = "time /DIRECTORY/sfs_code 1 ";
  $cmd .= int($NSIMS/$NJOBS+0.5);
  $cmd .= " -t $THETA -r $r ";
  
  for(my $i=0; $i<$NJOBS; $i++){
    my $seed = int(rand(40000000));
    if(exists($SEEDS{$seed})){ #ENSURE UNIQUE SEEDS FOR EACH TASK
      redo;
    }
    print OUT "$cmd -s $seed\n";
  }
  
  close(OUT);
  chdir "..";
}
