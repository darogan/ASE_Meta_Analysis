#!/usr/bin/env perl
use warnings;
use strict;
use Getopt::Long;
use FindBin qw($Bin);
use lib "$FindBin::Bin/../source";
use CF::Constants;
use CF::Helpers;

##########################################################################
# Copyright 2014, Philip Ewels (phil.ewels@babraham.ac.uk)               #
# Copyright 2020, Russell HAmilton (rsh46@cam.ac.uk)                     #
#                                                                        #
# This file is part of Cluster Flow.                                     #
#                                                                        #
# Cluster Flow is free software: you can redistribute it and/or modify   #
# it under the terms of the GNU General Public License as published by   #
# the Free Software Foundation, either version 3 of the License, or      #
# (at your option) any later version.                                    #
#                                                                        #
# Cluster Flow is distributed in the hope that it will be useful,        #
# but WITHOUT ANY WARRANTY; without even the implied warranty of         #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
# GNU General Public License for more details.                           #
#                                                                        #
# You should have received a copy of the GNU General Public License      #
# along with Cluster Flow.  If not, see <http://www.gnu.org/licenses/>.  #
##########################################################################

# Module requirements
my %requirements = (
	'cores' 	=> '2',
	'memory' 	=> '8G',
	'modules' 	=> ['SNPsplit', 'samtools'],
	'references'=> 'SNPsplit',
	'time' 		=> sub {
		my $cf = $_[0];
		my $num_files = $cf->{'num_starting_merged_aligned_files'};
		$num_files = ($num_files > 0) ? $num_files : 1;
                return CF::Helpers::minutes_to_timestamp ($num_files * 6 * 60);
	}
);

# Help text
my $helptext = "".("-"x15)."\n SNPsplit Module\n".("-"x15)."\n
SNPsplit: A tool to determine allele-specific alignments from high-throughput sequencing experiments that have been aligned to N-masked genomes \n\n";

# Setup
my %cf = CF::Helpers::module_start(\%requirements, $helptext);


print CF::Helpers::human_readable_to_bytes($cf{'memory'}), "\n";


# MODULE
# Check that we have a genome defined
if(!defined($cf{'refs'}{'SNPsplit'})){
	die "\n\n###CF Error: No SNPsplit ref path found in run file $cf{run_fn} for job $cf{job_id}. Exiting..";
} else {
	warn "\nAligning against SNPsplit path: $cf{refs}{SNPsplit}\n\n";
}

# Set up optional parameters
my $SE_or_PE       = (defined($cf{'params'}{'paired'})) ? '--paired' : '';
my $sorted         = (defined($cf{'params'}{'sorted'})) ? '--no_sort' : '';
my $sorted_paired  = (defined($cf{'params'}{'sorted,paired'})) ? '--no_sort --paired' : '';


open (RUN,'>>',$cf{'run_fn'}) or die "###CF Error: Can't write to $cf{run_fn}: $!";

# Print version information about the module.
warn "---------- SNPsplit version information ----------\n";
warn `SNPsplit --version`;
warn "\n------- End of SNPsplit version information ------\n";

# we want e.g. SNPsplit --paired --no_sort --snp_file <SNPFILE> <BAM>

foreach my $file (@{$cf{'prev_job_files'}}){
                my $timestart = time;
                my $filetype; 
                my $command = "";
                my $prefix = $file;
            
                if ($prefix =~ s/\.([sb]am$)//){
                        $filetype = $1;
                        warn "\nGuessing file $file is a $filetype file\n";
                } else {
                        warn "\n Cant determine file-type for $file. Not the expected bam... \n";
                }

                my $output_fn = $prefix.".SNPsplit_report.txt";

                $command .= "SNPsplit $SE_or_PE $sorted $sorted_paired --snp_file $cf{refs}{SNPsplit} $file";

                warn "\n###CFCMD $command\n\n";


                if(!system ($command)){
                        # SNPsplit worked - print out resulting filenames
                        my $duration =  CF::Helpers::parse_seconds(time - $timestart);
                        warn "###CF SNPsplit successfully exited, took $duration..\n";
                        if(-e $output_fn){
                                print RUN "$cf{'job_id'}\t$output_fn\n";
                        } else {
                                warn "\n###CF Error! SNPsplit output file $output_fn not found..\n";
                        }
                } else {
                        warn "\n###CF Error! SNPsplit failed, exited in an error state for input file '$file': $? $!\n\n";
                }
}





close (RUN);
