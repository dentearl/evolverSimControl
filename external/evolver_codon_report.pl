#!/usr/bin/perl -w
# Copyright (C) 2008-2011 by
# George Asimenos, Robert C. Edgar, Serafim Batzoglou and Arend Sidow.
# 
# All rights reserved. Reproduced and distributed here with permission.
# 
##############################
use warnings;
use strict;

my %tr = (
	'AAA' => 'K',
	'AAC' => 'N',
	'AAG' => 'K',
	'AAT' => 'N',
	'ACA' => 'T',
	'ACC' => 'T',
	'ACG' => 'T',
	'ACT' => 'T',
	'AGA' => 'R',
	'AGC' => 'S',
	'AGG' => 'R',
	'AGT' => 'S',
	'ATA' => 'I',
	'ATC' => 'I',
	'ATG' => 'M',
	'ATT' => 'I',
	'CAA' => 'Q',
	'CAC' => 'H',
	'CAG' => 'Q',
	'CAT' => 'H',
	'CCA' => 'P',
	'CCC' => 'P',
	'CCG' => 'P',
	'CCT' => 'P',
	'CGA' => 'R',
	'CGC' => 'R',
	'CGG' => 'R',
	'CGT' => 'R',
	'CTA' => 'L',
	'CTC' => 'L',
	'CTG' => 'L',
	'CTT' => 'L',
	'GAA' => 'E',
	'GAC' => 'D',
	'GAG' => 'E',
	'GAT' => 'D',
	'GCA' => 'A',
	'GCC' => 'A',
	'GCG' => 'A',
	'GCT' => 'A',
	'GGA' => 'G',
	'GGC' => 'G',
	'GGG' => 'G',
	'GGT' => 'G',
	'GTA' => 'V',
	'GTC' => 'V',
	'GTG' => 'V',
	'GTT' => 'V',
	'TAA' => '*',
	'TAC' => 'Y',
	'TAG' => '*',
	'TAT' => 'Y',
	'TCA' => 'S',
	'TCC' => 'S',
	'TCG' => 'S',
	'TCT' => 'S',
	'TGA' => '*',
	'TGC' => 'C',
	'TGG' => 'W',
	'TGT' => 'C',
	'TTA' => 'L',
	'TTC' => 'F',
	'TTG' => 'L',
	'TTT' => 'F'
);

my %ff = (
	'ACA' => 'T',
	'ACC' => 'T',
	'ACG' => 'T',
	'ACT' => 'T',
	'CCA' => 'P',
	'CCC' => 'P',
	'CCG' => 'P',
	'CCT' => 'P',
	'CGA' => 'R',
	'CGC' => 'R',
	'CGG' => 'R',
	'CGT' => 'R',
	'CTA' => 'L',
	'CTC' => 'L',
	'CTG' => 'L',
	'CTT' => 'L',
	'GCA' => 'A',
	'GCC' => 'A',
	'GCG' => 'A',
	'GCT' => 'A',
	'GGA' => 'G',
	'GGC' => 'G',
	'GGG' => 'G',
	'GGT' => 'G',
	'GTA' => 'V',
	'GTC' => 'V',
	'GTG' => 'V',
	'GTT' => 'V',
	'TCA' => 'S',
	'TCC' => 'S',
	'TCG' => 'S',
	'TCT' => 'S'
);

my $genome1 = "Ancestral";
my $genome2 = "Evolved";

$genome1 = $ARGV[0] if(defined($ARGV[0]));
$genome2 = $ARGV[1] if(defined($ARGV[1]));

sub uniq
{
	my %saw;
	@saw{@_} = ();
	return (sort keys %saw);
}

my @cod = sort keys %tr;
my @aa = uniq(values %tr);

my %cod1;
my %aa1;
my %cod2;
my %aa2;
my %cod1cod2;
my %aa1aa2;
my $e_tot = 0;
my $e_syn = 0;
my $e_non = 0;
my $ff_tot = 0;
my $ff_mut = 0;

foreach my $c (@cod)
{
	$cod1{$c} = 0;
	$cod2{$c} = 0;
	foreach my $c_ (@cod)
	{
		$cod1cod2{$c.$c_} = 0;
	}
}

foreach my $a (@aa)
{
	$aa1{$a} = 0;
	$aa2{$a} = 0;
	foreach my $a_ (@aa)
	{
		$aa1aa2{$a.$a_} = 0;
	}
}

while(<STDIN>)
{
	chomp;
	my @line = split(/\s+/);
	die("Invalid line: $_") unless (@line == 2);
	my $c1 = $line[0];
	my $c2 = $line[1];
	die("Invalid codons: $_") unless (defined($tr{$c1}) && defined($tr{$c2}));

	my $a1 = $tr{$c1};
	my $a2 = $tr{$c2};

	$cod1{$c1}++;
	$cod2{$c2}++;
	$aa1{$a1}++;
	$aa2{$a2}++;
	$cod1cod2{$c1.$c2}++;
	$aa1aa2{$a1.$a2}++;
}

foreach my $c1 (@cod)
{
	my $a1 = $tr{$c1};
	foreach my $c2 (@cod)
	{
		my $n = $cod1cod2{$c1.$c2};
		next unless ($n > 0);

		my $a2 = $tr{$c2};

		$e_tot += $n;
		if ($c1 ne $c2)
		{
			if ($a1 eq $a2) { $e_syn += $n } else { $e_non += $n; }
		}

		if (defined($ff{$c1}) && defined($ff{$c2}) && ($ff{$c1} eq $ff{$c2}))
		{
			$ff_tot += $n;
			$ff_mut += $n if ($c1 ne $c2);
		}
	}
}

sub mydiv
{
	my ($fmt, $x, $y) = @_;
	return sprintf($fmt, ($x / (1.0 * $y))) if ($y != 0);

	if ($x == 0)
	{
		return "n/a";
	}
	else
	{
		return "inf";
	}
}

sub showline
{
	my ($label, $anc, $anc_tot, $ev, $ev_tot) = @_;

	my $anc_pc = $anc * 100.0 / $anc_tot;
	my $ev_pc = $ev * 100.0 / $ev_tot;
	my $dcount = $ev - $anc;
	my $dpct = $ev_pc - $anc_pc;
	my $chg = "";

	if ($anc == 0)
	{
		if ($ev == 0) { $chg = "+0.0%" } else { $chg = "inf" }
	}
	else
	{
		$chg = sprintf("%+.1f%%", ($ev - $anc) * 100.0 / $anc);
	}
	printf "%5s  %10d  %5.1f%%  %10d  %5.1f%%  %10d  %+6.1f%%  %7s\n", $label, $anc, $anc_pc, $ev, $ev_pc, $dcount, $dpct, $chg;
}

if ($e_tot == 0)
{
	print "No aligned codons; nothing to report.\n";
	exit;
}

print "Codon Alignment Report\n";
print "======================\n";
print "\n";
print "Genome 1: $genome1\n";
print "Genome 2: $genome2\n";
print "\n";
print "\n";
print "Evolutionary statistics:\n";
print "------------------------\n";
print "\n";
printf "%10d  %5.1f%%  Total aligned codons\n", $e_tot, 100.0;
printf "%10d  %5.1f%%  Identical codons\n", $e_tot - $e_syn - $e_non, ($e_tot - $e_syn - $e_non) * 100.0 / $e_tot;
printf "%10d  %5.1f%%  Synonymous changes\n", $e_syn, $e_syn * 100.0 / $e_tot;
printf "%10d  %5.1f%%  Non-Synonymous changes\n", $e_non, $e_non * 100.0 / $e_tot;
print "\n";
printf "%10s          Ks/Ka approx.\n", mydiv("%.2f", $e_syn, $e_non);
print "\n";
printf "%10d          4D sites\n", $ff_tot;
printf "%10d          4D differencies\n", $ff_mut;
printf "%10s          4D mutation rate\n", mydiv("%.4f", $ff_mut, $ff_tot);
print "\n";
print "\n";
print "Frequencies of aligned codons:\n";
print "------------------------------\n";
print "\n";
print "          Genome1    Gen1.    Genome2    Gen2.                             \n";
print "Codon       Count     Pct       Count     Pct      dCount     dPct    % Chg\n";
print "=====  ==========  ======  ==========  ======  ==========  =======  =======\n";
foreach my $c (@cod)
{
	my $label = $tr{$c};
	$label .= (defined($ff{$c}) ? "-" : " ");
	$label .= $c;
	showline($label, $cod1{$c}, $e_tot, $cod2{$c}, $e_tot);
}
showline("Total", $e_tot, $e_tot, $e_tot, $e_tot);
print "\n";
print "\n";
print "Frequencies of aligned amino acids:\n";
print "-----------------------------------\n";
print "\n";
print "          Genome1    Gen1.    Genome2    Gen2.                             \n";
print "   AA       Count     Pct       Count     Pct      dCount     dPct    % Chg\n";
print "=====  ==========  ======  ==========  ======  ==========  =======  =======\n";
foreach my $a (@aa)
{
	showline($a, $aa1{$a}, $e_tot, $aa2{$a}, $e_tot);
}
showline("Total", $e_tot, $e_tot, $e_tot, $e_tot);

print "\n";
print "\n";
print "Codon correspondence frequencies:\n"; 
print "---------------------------------\n"; 
print "\n";
print 'G1\G2: ', join("    ", @cod), "\n";

foreach my $c1 (@cod)
{
	printf "%3s", $c1;
	foreach my $c2 (@cod)
	{
		printf "%7d", $cod1cod2{$c1.$c2};
	}
	print "\n";
}

print "\n";
print "\n";
print "Amino acid correspondence frequencies:\n"; 
print "--------------------------------------\n"; 
print "\n";
print 'G1\G2:   ', join("      ", @aa), "\n";

foreach my $a1 (@aa)
{
	printf "%3s", $a1;
	foreach my $a2 (@aa)
	{
		printf "%7d", $aa1aa2{$a1.$a2};
	}
	print "\n";
}

