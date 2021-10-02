#!/usr/bin/perl

use strict;


my $TableCSV           = "Tables_S2.tsv";
my $TableBED           = "Tables_S2.bed";
my $TableBEDSRT        = "Tables_S2.srt.bed";
my $TableBEDSRTDIST    = "Tables_S2.srt.dist.bed";
my $TableBEDSRTDISTUNQ = "Tables_S2.srt.dist.uniq.bed";
my $TablePAIRS         = "Tables_S2.pairs.txt";



open(BED, ">$TableBED") || die "Can't open $TableBED for writing: $!\n";

# Tsga14,chr6,30653457,30693749,,P,,,1,3,mat,50-60,Mest Cluster

my %type;
open(CSV, $TableCSV) || die "Can't open $TableCSV for reading: $!\n";
while( my $line = <CSV>)  
  {
  	chomp($line);
    my @entry = split(/\t/, $line);

    if(length($entry[0]) > 0)
      { 
        #print $line, "\n";

        my $score = $entry[11];

        if($score =~ m/[0-9]+-[0-9]+/)
          { 
             my @scr = split(/-/, $score);
             $score = ($scr[1]+$scr[0]) / 2;
             print BED "$entry[1]\t$entry[2]\t$entry[3]\t$entry[0]\t$score\t.\n";
             $type{$entry[0]} = $entry[9];
          }
      }
  }
close CSV;


close BED;

system("bedtools sort -i ${TableBED} > ${TableBEDSRT}");
system("bedtools closest -k 1 -N -D a -a ${TableBEDSRT} -b ${TableBEDSRT} > ${TableBEDSRTDIST}");
#system("bedtools closest -k 100 -N -d -a ${TableBEDSRT} -b ${TableBEDSRT} > ${TableBEDSRTDIST}");


#my %UniqueList;
# chr1    36530534        36535739        Ankrd23 55      .       chr1    63143595        63176822        Ndufs1  55      .       26607857

open(BEDDISTUNQ, ">$TableBEDSRTDISTUNQ") || die "Can't open $TableBEDSRTDISTUNQ for writing: $!\n";
open(BEDDIST, "$TableBEDSRTDIST") || die "Can't open $TableBEDSRTDIST for reading: $!\n";
my %seen;
while( my $line = <BEDDIST>)  
  {
    chomp($line);
    my @entry = split(/\t/, $line);
    my $outline="";

    my $id = "";
    if( $entry[1] < $entry[7] )
        { 
          $outline      = $line;
          $id           = $entry[3] . "_" .$entry[9];
          $seen{ $id } += 1;
        }
    else{ 
          $outline      = "$entry[6]\t$entry[7]\t$entry[8]\t$entry[9]\t$entry[10]\t$entry[11]\t$entry[0]\t$entry[1]\t$entry[2]\t$entry[3]\t$entry[4]\t$entry[5]\t$entry[12]";
          $id           = $entry[9] . "_" .$entry[3];
          $seen{ $id } += 1;
        }

    if($seen{ $id } <= 1)
      { 
    	print BEDDISTUNQ $outline, "\n"; 
      }
  }
close BEDDIST;
close BEDDISTUNQ;


open(PAIRS, ">$TablePAIRS") || die "Can't open $TablePAIRS for writing: $!\n";

open(BEDDISTUNQ, "$TableBEDSRTDISTUNQ") || die "Can't open $TableBEDSRTDISTUNQ for reading: $!\n";
while( my $line = <BEDDISTUNQ>)  
  {
    chomp($line);
    my @entry = split(/\t/, $line);
    print PAIRS "$entry[3]\t$entry[4]\t", $type{$entry[3]}, "\t$entry[9]\t$entry[10]\t", $type{$entry[9]}, "\t$entry[12]\n";
  }
close BEDDISTUNQ;
close PAIRS;





print "_END_OF_SCRIPT_\n";