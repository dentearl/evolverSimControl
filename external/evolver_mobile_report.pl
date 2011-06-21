#!/usr/bin/perl -w
# Copyright (C) 2008-2011 by
# George Asimenos, Robert C. Edgar, Serafim Batzoglou and Arend Sidow.
# 
# All rights reserved. Reproduced and distributed here with permission.
# 
##############################
use strict;
use warnings;
use Getopt::Long;

my $MOBILES_LOG; # location of the mobiles.log logfile, from the handle_mobiles.pl step
GetOptions('mobilesLog=s' => \$MOBILES_LOG);
my %argHash=(
             mobilesLog => $MOBILES_LOG);
for my $key (keys %argHash){
    if ($argHash{$key} eq ""){
        &usage("Please specify missing argument, --$key.\n");
    }
}
sub usage{
    my $message = $_[0];
    print $message;
    print "USAGE: $0 <intra step evo.[chr].log files> --mobilesLog <logFile>\n\n";
    exit(2);
}

my %stats;
my $state = 0;
my %rpg;
sub read_mobiles_log
{
    my $logFile = $_[0];
	local *FH;
	open (FH, '<', $logFile) || die("Cannot open logfile: $logFile");

	while(<FH>)
	{
		chomp;
		if (/^  - RPG #(\d+): (\d+) copies/)
		{
			$rpg{'RPG.' . $1} = $2;
		}
	}
	
	close (FH);
}

sub dotted
{
	my $x = "$a,$b";
	if ($x =~ /^([^.]+)\.(\d+),([^.]+)\.(\d+)/)
	{
		my $z = $1 cmp $3;
		return ($z || ($2 <=> $4));
	}
	else
	{
		return $a cmp $b;
	}
}

read_mobiles_log($MOBILES_LOG);

while (<>)
{
	chomp;
	if (/^Mobile elements:$/)
	{
		$state = 1;
	}
	elsif (/of genome/)
	{
		$state = 0;
	}
	elsif ($state == 1)
	{
		s/^\s+//;
		s/\s+$//;
		next unless (/^\d+/);
		my ($INDEX, $NAME, $RATE, $RATE_OR_RATEPERBASE, $PCT, $AVGDEL, $STDDEV, $DELETES, $ACCEPTED, $ACCBASES, $SEQUENCE) = split(/\s+/);
		
		$stats{$NAME}->{'success'} += $ACCEPTED;
		$stats{$NAME}->{'bases'} += $ACCBASES;
	}
}


print "RPG Report:\n\n";
printf "%9s %10s %9s %9s %9s\n", "", "Desired", "Actual", "", "";
printf "%9s %10s %9s %9s %9s\n", "Name", "Copies", "Copies", "Diff", "Bases";
printf "%9s=%10s=%9s=%9s=%9s\n", "=" x 9, "=" x 10, "=" x 9, "=" x 9, "=" x 9;

my $t_req = 0;
my $t_success = 0;
my $t_diff = 0;
my $t_bases = 0;

foreach my $rpg (sort dotted keys %rpg)
{
	defined($stats{$rpg}) || die("Invalid rpg $rpg");
	$t_req += $rpg{$rpg};
	$t_success += $stats{$rpg}->{'success'};
	$t_diff += $stats{$rpg}->{'success'} - $rpg{$rpg};
	$t_bases += $stats{$rpg}->{'bases'};
	printf "%9s %10u %9u %9d %9d\n", $rpg, $rpg{$rpg}, $stats{$rpg}->{'success'}, $stats{$rpg}->{'success'} - $rpg{$rpg}, $stats{$rpg}->{'bases'};
}

printf "%9s=%10s=%9s=%9s=%9s\n", "=" x 9, "=" x 10, "=" x 9, "=" x 9, "=" x 9;
printf "%9s %10u %9u %9d %9d\n", "Total", $t_req, $t_success, $t_diff, $t_bases;


print "\n\nME Report:\n\n";
printf "%19s %9s %11s\n", "Name", "Copies", "Bases";
printf "%19s=%9s=%11s\n", "=" x 19, "=" x 9, "=" x 11;

my $m_success = 0;
my $m_bases = 0;

foreach my $rpg (sort dotted keys %stats)
{
	if (!($rpg =~ /^RPG/))
	{
		$m_success += $stats{$rpg}->{'success'};
		$m_bases += $stats{$rpg}->{'bases'};
		printf "%19s %9u %11d\n", $rpg, $stats{$rpg}->{'success'}, $stats{$rpg}->{'bases'};
	}
}

printf "%19s=%9s=%11s\n", "=" x 19, "=" x 9, "=" x 11;
printf "%19s %9u %11d\n", "Total", $m_success, $m_bases;
