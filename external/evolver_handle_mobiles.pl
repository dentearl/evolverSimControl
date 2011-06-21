#!/usr/bin/perl -w
# Copyright (C) 2008-2011 by
# George Asimenos, Robert C. Edgar, Serafim Batzoglou and Arend Sidow.
# 
# All rights reserved. Reproduced and distributed here with permission.
# 
##############################

use strict;
use warnings;
use IO::Handle;
use Getopt::Long;

my $EVO_BIN        = "";
my $CVT_BIN        = "";
my $SCR_DIR        = "";
my $GENOME_NAME    = ""; # this is the string stored in the .rev file,
                         # not the path to the cycle
my $PARENT_DIR     = "";
my $STEPSIZE       = "";
my $MES_CFG        = "";
my $ME_FA          = "";
my $ME_GFF         = "";
my $LTR_FA         = "";
my $MODEL_TXT      = "";

GetOptions('evo=s'       => \$EVO_BIN
           ,'cvt=s'      => \$CVT_BIN
           ,'py=s'       => \$SCR_DIR
           ,'parentDir=s'=> \$PARENT_DIR
           ,'genome=s'   => \$GENOME_NAME
           ,'stepSize=f' => \$STEPSIZE
           ,'mescfg=s'   => \$MES_CFG
           ,'mefa=s'     => \$ME_FA
           ,'megff=s'    => \$ME_GFF
           ,'ltr=s'      => \$LTR_FA
           ,'model=s'    => \$MODEL_TXT);
my %argHash=(
             evo      => $EVO_BIN,
             cvt      => $CVT_BIN,
             py       => $SCR_DIR,
             parentDir=> $PARENT_DIR,
             genome   => $GENOME_NAME,
             stepSize => $STEPSIZE,
             mescfg   => $MES_CFG,
             mefa     => $ME_FA,
             megff    => $ME_GFF,
             ltr      => $LTR_FA,
             model    => $MODEL_TXT);
for my $key (keys %argHash){
    if ($argHash{$key} eq ""){
        &usage("Please specify missing argument, --$key.\n");
    }
}

my %fate;
########################################

sub usage{
    my $message = $_[0];
    print $message;
    print "USAGE: $0 --evo <Path_to_evo> -- cvt <Path_to_cvt> --py <Path_to_py>\n",
        "\t--genome <ingenome> --stepSize <branch length> --mescfg <model file .cfg>\n",
        "\t--mefa <mobile fasta> --megff <moblie annotiation .gff> --ltr <ltr fasta> --model <model .txt>\n(human-readable log to stdout; machine-readable to stderr)\n\n";
    exit(2);
}

sub fate
{
	my $class = shift;
	my $title = shift;
	my $label = shift;
	push(@{$fate{"$class.$title"}}, $label);
}

sub cmd
{
	my $s = shift;
	my $r = 0xffff & system($s);
	($r == 0) || die("ERROR: $0: cmd='$s', status=" . ($r>>8));
}

sub fix_geneix
{
	my $in = shift;
	my $out = shift;
	cmd("python $SCR_DIR/evolver_gff_fixgeneix.py $in >$out 2>/dev/null");
}

sub read_gff
{
	my $fname = shift;
	my $gff = shift;

	local *FH;
	open(FH, '<', $fname) || die("ERROR: $0: Cannot open $fname");

	while(<FH>)
	{
		chomp;
		my @s = split(/\t/);
		my $t = shift(@s);
		push(@{$gff->{$t}}, [ @s ]);
	}
	
	close(FH) || die("ERROR: $0: Cannot close $fname");
}

sub write_gff
{
	my $fname = shift;
	my $gff = shift;

	local *FH;
	open(FH, '>', $fname) || die("ERROR: $0: Cannot open $fname");

	foreach my $g (sort(keys(%$gff)))
	{
		foreach my $s (@{$gff->{$g}})
		{
			print(FH $g, "\t", join("\t", @$s), "\n");
		}
	}
	
	close(FH) || die("ERROR: $0 Cannot close $fname");
}

sub read_conf
{
	my $fname = shift;
	my $conf = shift;

	local *FH;
	open(FH, '<', $fname) || die("ERROR: $0: Cannot open $fname");

	while(<FH>)
	{
		s/^\s+//;
		s/\s+$//;
		/^#/ && next;
		/^$/ && next;

		my @s = split(/\s+/);
		my $d = lc($s[0]);
		my $l = scalar(@s);
		if (($d eq 'rate') && ($l == 12))
		{
			$conf->{'rate'}->{$s[1]} = { 'relins' => $s[2], 'avgdel' => $s[3], 'stddev' => $s[4], 'pct' => $s[5], 'minlen' => $s[6], 'maxlen' => $s[7], 'maxcount' => $s[8], 'duprate' => $s[9], 'delrate' => $s[10], 'lifetime' => $s[11] };
		}
		elsif (($d eq 'totalinsertrate') && ($l == 2))
		{
			$conf->{'totins'} = $s[1];
		}
		elsif (($d eq 'ltrclass') && ($l == 4))
		{
			$conf->{'ltr'}->{$s[1]} = { 'minlen' => $s[2], 'maxlen' => $s[3] };
		}
		elsif (($d eq 'branchlengthfactor') && ($l == 2))
		{
			$conf->{'blf'} = $s[1];
		}
		elsif (($d eq 'rpgheader') && ($l == 4))
		{
			$conf->{'rpg_avgdel'} = $s[1];
			$conf->{'rpg_stddev'} = $s[2];
			$conf->{'rpg_pct'} = $s[3];
		}
		elsif (($d eq 'polyatail') && ($l == 2))
		{
			$conf->{'rpg_polya'} = $s[1];
		}
		elsif (($d eq 'countpertick') && ($l == 2))
		{
			$conf->{'rpg_count'} = $s[1];
		}
		elsif (($d eq 'maxrpgsize') && ($l == 2))
		{
			$conf->{'rpg_maxsize'} = $s[1];
		}
		else
		{
			print("WARNING: Ignoring line: '$_' in config file $fname\n");
		}
	}
	
	close(FH) || die("ERROR $0: Cannot close $fname");

	#
	# Validation
	#
	defined($conf->{'totins'}) || die("ERROR: $0: TotalInsertRate missing from configuration file");
	defined($conf->{'blf'}) || die("ERROR: $0: BranchLengthFactor missing from configuration file");
	defined($conf->{'rpg_avgdel'}) || die("ERROR: $0: RPGHeader missing from configuration file");
	defined($conf->{'rpg_polya'}) || die("ERROR: $0: PolyATail missing from configuration file");
	defined($conf->{'rpg_count'}) || die("ERROR: $0: CountPerTick missing from configuration file");
	defined($conf->{'rpg_maxsize'}) || die("ERROR: $0: MaxRPGSize missing from configuration file");
	#
	# LTR classes need to exist
	#
	foreach (keys(%{$conf->{'ltr'}}))
	{
		defined($conf->{'rate'}->{$_}) || die("ERROR: $0: Invalid LTR class '$_'");
	}

	print("* Loaded config file: ", scalar(keys(%{$conf->{'rate'}})), " classes (", scalar(keys(%{$conf->{'ltr'}})), " LTR-like)\n");
}

sub write_pretty
{
	local *F = shift;
	my $s = shift;

	my $l = length($s);
	while ($l > 60)
	{
		print(F substr($s, 0, 60), "\n");
		$s = substr($s, 60);
		$l -= 60;
	}

	($l > 0) && print(F "$s\n");
}

sub read_fasta_flat
{
	my $fname = shift;
	my $seq = shift;
	
	$$seq = '';
	
	local *FH;
	open(FH, '<', $fname) || die("ERROR: $0: Cannot open $fname");

	my $title = '';

	while(<FH>)
	{
		s/^\s+//;
		s/\s+$//;
		/^$/ && next;
		
		if (/>(.*)$/)
		{
			($1 eq '') && die("ERROR: $0: Invalid header '$_' in file $fname");
			$title = $1;
		}
		else
		{
			($title eq '') && die("ERROR: $0: Sequence with no header in file $fname");
			$$seq .= $_;
		}
	}
	
	close(FH) || die("ERROR: $0: Cannot close $fname");
}

sub read_fasta
{
	my $fname = shift;
	my $seq = shift;
	
	%$seq = ();
	
	local *FH;
	open(FH, '<', $fname) || die("ERROR: $0: Cannot open $fname");

	my $title = '';

	while(<FH>)
	{
		s/^\s+//;
		s/\s+$//;
		/^$/ && next;
		
		if (/>(.*)$/)
		{
			($1 eq '') && die("ERROR: $0: Invalid header '$_' in file $fname");
			defined($seq->{$1}) && die("ERROR: $0: Duplicate header '$_' in file $fname");
			$title = $1;
		}
		else
		{
			($title eq '') && die("ERROR: $0: Sequence with no header in file $fname");
			$seq->{$title} .= $_;
		}
	}
	
	close(FH) || die("ERROR: $0: Cannot close $fname");
}

sub nicesort
{
	if ($a =~ /^([^.]+)\.(\d+)/)
	{
		my ($a_title, $a_id) = ($1, $2);
		if ($b =~ /^([^.]+)\.(\d+)/)
		{
			my ($b_title, $b_id) = ($1, $2);
			if ($a_title eq $b_title)
			{
				return $a_id <=> $b_id;
			}
			else
			{
				return $a_title cmp $b_title;
			}
		}
		else
		{
			return $a cmp $b;
		}
	}
	else
	{
		return $a cmp $b;
	}
}

sub write_fasta
{
	my $fname = shift;
	my $seq = shift;

	local *FH;
	open(FH, '>', $fname) || die("ERROR: $0: Cannot open $fname");

	foreach (sort nicesort (keys(%$seq)))
	{
		my $s = $seq->{$_};
		print(FH ">$_\n");
		write_pretty(*FH, $s);
	}

	close(FH) || die("ERROR: $0: Cannot close $fname");
}

sub decompose
{
	my $inh = shift;
	my $outh = shift;
	
	%$outh = ();
	
	foreach (keys(%$inh))
	{
		if (/(.+)\.([^.:]+):ACTIVE$/)
		{
			$outh->{$1}->{$2}->{'active'} = 1;
			$outh->{$1}->{$2}->{'seq'} = $inh->{$_};
		}
		elsif (/(.+)\.([^.]+)$/)
		{
			$outh->{$1}->{$2}->{'active'} = 0;
			$outh->{$1}->{$2}->{'seq'} = $inh->{$_};
		}
		else
		{
			die("ERROR: $0: Invalid mobile element header $_");
		}
	}
}

sub make_header
{
	my $mes = shift;
	my $class = shift;
	my $title = shift;
	if ($mes->{$class}->{$title}->{'active'})
	{
		return "$class.$title:ACTIVE";
	}
	else
	{
		return "$class.$title";
	}
}

sub recompose
{
	my $inh = shift;
	my $outh = shift;

	%$outh = ();
	
	foreach my $class (keys(%$inh))
	{
		foreach my $title (keys(%{$inh->{$class}}))
		{
			my $t = make_header($inh, $class, $title);
			$outh->{$t} = $inh->{$class}->{$title}->{'seq'};
		}
	}
}

sub print_mes
{
	my $h = shift;
	foreach my $class (sort(keys(%$h)))
	{
		print("    $class: ", join(', ', sort {$a <=> $b} (keys(%{$h->{$class}}))), "\n");
	}
}

sub validate_mes
{
	my $me = shift;
	my $cfg = shift;
	foreach my $class (keys(%$me))
	{
		defined($cfg->{'rate'}->{$class}) || die("ERROR: $0: Don't know anything about ME class '$class'");
	}
}

sub validate_ltrs
{
	my $me = shift;
	my $ltr = shift;
	my $cfg = shift;
	foreach my $class (keys(%$ltr))
	{
		defined($cfg->{'ltr'}->{$class}) || die("ERROR: $0: Class '$class' appears in LTR file but isn't defined as LTR in config");
		foreach my $title (keys(%{$ltr->{$class}}))
		{
			defined($me->{$class}->{$title}) || die("ERROR: $0: Element $class.$title appears in LTR file but not in ME file");
		}
	}
	foreach my $class (keys(%{$cfg->{'ltr'}}))
	{
		defined($me->{$class}) && (defined($ltr->{$class}) || die("ERROR: $0: Class '$class' configured as LTR but does not appear in LTR file"));
	}
}

sub validate_gff
{
	my $gff = shift;
	my $seqs = shift;
	
	foreach (keys(%$gff))
	{
		defined($seqs->{$_}) || die("ERROR: $0: Invalid element '$_' in GFF file");
	}
}

sub birth_death
{
	my $cfg = shift;
	
	my $mes = shift;
	my $ltrs = shift;
	my $gff = shift;

	my $out_mes = shift;
	my $out_ltrs = shift;
	my $out_gff = shift;

	my $branch = shift;

	%$out_mes = ();
	%$out_ltrs = ();
	%$out_gff = ();

	print("* Performing birth/death process...\n");
	foreach my $class (keys(%$mes))
	{
		my $maxtitle = 0;
		foreach my $title (keys(%{$mes->{$class}}))
		{
			if ($maxtitle < $title)
			{
				$maxtitle = $title;
			}
		}
		foreach my $title (keys(%{$mes->{$class}}))
		{
			my $roll = rand();
			my $prob = ($cfg->{'rate'}->{$class}->{'duprate'} + $cfg->{'rate'}->{$class}->{'delrate'}) * $branch;
			if ($roll < $prob)
			{
				my $prob = ($cfg->{'rate'}->{$class}->{'delrate'} * $branch) / $prob;
				$roll = rand();
				if ($roll < $prob)
				{
					if ($mes->{$class}->{$title}->{'active'} == 0)
					{
						fate($class, $title, "bddel");
						print("  - Deleting element $class.$title\n");
						next;
					}
					print("  - WARNING: Active element $class.$title was going to be deleted; will duplicate instead.\n");
				}
				++$maxtitle;
				fate($class, $title, "dupto:$class.$maxtitle");
				print("  - Duplicating element $class.$title -> $class.$maxtitle\n");
				fate($class, $maxtitle, "dupfrom:$class.$title");
				$out_mes->{$class}->{$maxtitle}->{'seq'} = $mes->{$class}->{$title}->{'seq'};
				$out_mes->{$class}->{$maxtitle}->{'active'} = 0;
				my $oldhead = make_header($mes, $class, $title);
				my $newhead = make_header($out_mes, $class, $maxtitle);
				$out_gff->{$newhead} = [ @{$gff->{$oldhead}} ] if (defined($gff->{$oldhead}));
				if (defined($ltrs->{$class}))
				{
					$out_ltrs->{$class}->{$maxtitle}->{'seq'} = $ltrs->{$class}->{$title}->{'seq'};
					$out_ltrs->{$class}->{$maxtitle}->{'active'} = $ltrs->{$class}->{$title}->{'active'};
				}
			}
			my $t = make_header($mes, $class, $title);
			$out_mes->{$class}->{$title}->{'seq'} = $mes->{$class}->{$title}->{'seq'};
			$out_mes->{$class}->{$title}->{'active'} = $mes->{$class}->{$title}->{'active'};
			$out_gff->{$t} = [ @{$gff->{$t}} ] if (defined($gff->{$t}));
			if (defined($ltrs->{$class}))
			{
				$out_ltrs->{$class}->{$title}->{'seq'} = $ltrs->{$class}->{$title}->{'seq'};
				$out_ltrs->{$class}->{$title}->{'active'} = $ltrs->{$class}->{$title}->{'active'};
			}
		}
		# The following should now never happen as the active element is never deleted
		# but it's still there because in the beginning there is no active element
		if (scalar(keys(%{$out_mes->{$class}})) == 0)
		{
			my $title = '';
			my $r = int(rand(scalar(keys(%{$mes->{$class}}))));
			do
			{
				$title = each(%{$mes->{$class}});
			} while ($r--);
			
			print("  - Class $class ran out of items; resurrecting $class.$title\n");
			fate($class, $title, "resurrected");
			$out_mes->{$class}->{$title}->{'seq'} = $mes->{$class}->{$title}->{'seq'};
			$out_mes->{$class}->{$title}->{'active'} = $mes->{$class}->{$title}->{'active'};
			my $t = make_header($mes, $class, $title);
			$out_gff->{$t} = [ @{$gff->{$t}} ] if (defined($gff->{$t}));
			if (defined($ltrs->{$class}))
			{
				$out_ltrs->{$class}->{$title}->{'seq'} = $ltrs->{$class}->{$title}->{'seq'};
				$out_ltrs->{$class}->{$title}->{'active'} = $ltrs->{$class}->{$title}->{'active'};
			}
		}
	}
}

sub remove_excess
{
	my $cfg = shift;
	my $mes = shift;
	my $ltrs = shift;
	my $gff = shift;

	foreach my $class (sort(keys(%$mes)))
	{
		if (scalar(keys(%{$mes->{$class}})) > $cfg->{'rate'}->{$class}->{'maxcount'})
		{
			my $extra = scalar(keys(%{$mes->{$class}})) - $cfg->{'rate'}->{$class}->{'maxcount'};
			while ($extra--)
			{
				my $title = '';
				my $r = 0;
				do
				{
					$r = int(rand(scalar(keys(%{$mes->{$class}}))));
					do
					{
						$title = each(%{$mes->{$class}});
					} while ($r--);
				} while ($mes->{$class}->{$title}->{'active'});
				my $header = make_header($mes, $class, $title);
				print("  - Too many elements in class $class; deleting $header\n");
				fate($class, $title, "del2many");
				delete($mes->{$class}->{$title});
				delete($gff->{$header});
				if (defined($ltrs->{$class}))
				{
					delete($ltrs->{$class}->{$title});
				}
			}
		}
	}
}

sub evolve_mes
{
	my $cfg = shift;
	my $br = shift;
	my $mes = shift;
	my $me_gff = shift;
	my $me_fa = shift;
	my $ltrs = shift;
	my $ltr_fa = shift;
	my $evolved_me_gff = shift;
	my $evolved_me_fa = shift;
	my $evolved_ltr_fa = shift;
	my $model = shift;

	$br *= $cfg->{'blf'};
	
	#
	# Evolve MEs
	#
	cmd("$CVT_BIN -fromfasta $me_fa -torev $me_fa.tmp_mobiles.rev -genome mobiles >& /dev/null");
	unlink($evolved_me_fa);
	unlink($evolved_me_gff);
	foreach my $class (sort(keys(%$mes)))
	{
		foreach my $title (sort {$a <=> $b} (keys(%{$mes->{$class}})))
		{
			my $t = make_header($mes, $class, $title);
			my $seed = int(rand(1000000));
			print("* Evolving ME $class.$title...\n");
			cmd("$EVO_BIN -inseq $me_fa.tmp_mobiles.rev -branchlength $br -model $model -inannots $me_gff -outseq $me_fa.tmp_mobiles.out.rev -outannots $me_fa.tmp_mobiles.gff -chrname $t -seed $seed -log $me_fa.evo.log >& /dev/null");
			cmd("$CVT_BIN -fromrev $me_fa.tmp_mobiles.out.rev -tofasta $evolved_me_fa.tmp >& /dev/null");
			cmd("cat $evolved_me_fa.tmp >> $evolved_me_fa");
			cmd("cat $me_fa.tmp_mobiles.gff >> $evolved_me_gff");
			unlink("$evolved_me_fa.tmp");
			unlink("$me_fa.tmp_mobiles.gff");
			unlink("$me_fa.tmp_mobiles.out.rev");
			unlink("$me_fa.evo.log");
		}
	}
	unlink("$me_fa.tmp_mobiles.rev");

	#
	# Evolve LTRs
	#
	cmd("$CVT_BIN -fromfasta $ltr_fa -torev $ltr_fa.tmp_mobiles.rev -genome mobiles >& /dev/null");
	unlink($evolved_ltr_fa);
	foreach my $class (sort(keys(%$ltrs)))
	{
		foreach my $title (sort(keys(%{$ltrs->{$class}})))
		{
			my $t = make_header($ltrs, $class, $title);
			my $seed = int(rand(1000000));
			print("* Evolving LTR $class.$title...\n");
			cmd("$EVO_BIN -inseq $ltr_fa.tmp_mobiles.rev -branchlength $br -model $model -outseq $ltr_fa.tmp_mobiles.out.rev -chrname $t -seed $seed -log $ltr_fa.evo.log >& /dev/null");
			cmd("$CVT_BIN -fromrev $ltr_fa.tmp_mobiles.out.rev -tofasta $evolved_ltr_fa.tmp >& /dev/null");
			cmd("cat $evolved_ltr_fa.tmp >> $evolved_ltr_fa");
			unlink("$evolved_ltr_fa.tmp");
			unlink("$ltr_fa.tmp_mobiles.out.rev");
			unlink("$ltr_fa.evo.log");
		}
	}
	unlink("$ltr_fa.tmp_mobiles.rev");
}

sub read_genomic_gff
{
	my $gff_file = shift;
	my $gff = shift;
	print("* Reading genome gff entries...");

	local *FH;
	open(FH, '<', $gff_file) || die("ERROR: $0: Cannot open $gff_file");

	ENTRY: while(<FH>)
	{
		chomp;
		my @line = split(/\t/);
		next unless (($line[2] eq 'CDS') || ($line[2] eq 'UTR'));
		my @fields = split(/\s*;\s*/, $line[8]);
		foreach my $f (@fields)
		{
			next unless (length($f) > 0);
			my @attr = split(/\s+/, $f);
			if ($attr[0] eq 'gene_index')
			{
				if ($attr[1] =~ /^\d+$/)
				{
					push(@{$gff->{$line[0] . ':' . $attr[1]}}, [ (@line) ]);
					next ENTRY;
				}
				else
				{
					die("ERROR: $0: Invalid gene_index " . $attr[1]);
				}
			}
		}
		die("ERROR: $0: Missing gene_index");
	}
	
	close(FH);
	print(" ", scalar(keys(%$gff)), " entries loaded.\n");
}

sub filter_length
{
	my $maxlen = shift;
	my $gff = shift;

	print("* Deleting long genes...");
	my $num_deleted = 0;
	my $count = scalar(keys(%$gff)); # Reset each() pointer
	while (my ($gid, $recs) = each(%$gff))
	{
		my $l = 0;
		foreach my $r (@$recs)
		{
			$l += ${$r}[4] - ${$r}[3] + 1;
		}
		if ($l > $maxlen)
		{
			++$num_deleted;
			delete($gff->{$gid});
		}
	}
	print(" ", $num_deleted, " deleted, ", scalar(keys(%$gff)), " remaining.\n");
}

sub output_random_gene
{
	my $fname = shift;
	my $gff = shift;
	(scalar(keys(%$gff)) == 0) && die("ERROR: $0: No more genes");

	my $gid;
	my $g;
	my $r = int(rand(scalar(keys(%$gff))));
	do
	{
		($gid, $g) = each(%$gff);
	} while ($r--);
	
	my @lines = @$g;

	(scalar(@lines) == 0) && die("ERROR: $0: Invalid gene");

	die("ERROR: $0: Missing strand info") unless ((${$lines[0]}[6] eq '+') || (${$lines[0]}[6] eq '-'));
	@lines = (${$lines[0]}[6] eq '+') ? (sort {${$a}[3] <=> ${$b}[3]} @lines) : (sort {${$b}[3] <=> ${$a}[3]} @lines);
	
	local *FH;
	open(FH, '>', $fname) || die("ERROR: $0: Cannot open $fname");
	foreach my $line (@lines)
	{
		print(FH join("\t", @$line), "\n");
	}
	close(FH);

	delete($gff->{$gid});
}

sub activate
{
	my $cfg = shift;
	my $branch = shift;
	my $mes = shift;
	my $ltrs = shift;
	my $gff = shift;
	my $outseqs = shift;

	my $rsum = 0;
	
	foreach my $class (keys(%$mes))
	{
		$rsum += $cfg->{'rate'}->{$class}->{'relins'};
	}

	my $globrate = $cfg->{'totins'};
	print("* Active element selection:\n");

	foreach my $class (sort(keys(%$mes)))
	{
		my $activetitle = '';
		foreach my $title (keys(%{$mes->{$class}}))
		{
			if ($mes->{$class}->{$title}->{'active'})
			{
				$activetitle = $title;
				last;
			}
		}

		if ($activetitle ne '')
		{
			my $p = $branch / $cfg->{'rate'}->{$class}->{'lifetime'};
			my $roll = rand();
			if ($roll < $p)
			{
				if (scalar(keys(%{$mes->{$class}})) == 1)
				{
					print("  - It was time for active element $class.$activetitle to die, but it will be kept as there are no other elements in that class.\n");
					fate($class, $activetitle, "remainactivedie");
				}
				else
				{
					print("  - Deleting $class.$activetitle\n");
					fate($class, $activetitle, "deactivate");
					my $header = make_header($mes, $class, $activetitle);
					delete($gff->{$header});
					delete($mes->{$class}->{$activetitle});
					if (defined($ltrs->{$class}))
					{
						delete($ltrs->{$class}->{$activetitle});
					}
					$activetitle = '';
				}
			}
			else
			{
				print("  - Keeping $class.$activetitle\n");
				fate($class, $activetitle, "remainactive");
			}
		}

		my $title = $activetitle;
		if ($title eq '')
		{
			my $r = int(rand(scalar(keys(%{$mes->{$class}}))));
			do
			{
				$title = each(%{$mes->{$class}});
			} while ($r--);
			print("  - Activating $class.$title\n");
			fate($class, $title, "activate");
		}

		my $old_header = make_header($mes, $class, $title);
		$mes->{$class}->{$title}->{'active'} = 1;
		my $new_header = make_header($mes, $class, $title);
		if ($old_header ne $new_header)
		{
			if (defined($gff->{$old_header}))
			{
				$gff->{$new_header} = $gff->{$old_header};
				delete($gff->{$old_header});
			}
		}
		my $seq = $mes->{$class}->{$title}->{'seq'};
		if (defined($ltrs->{$class}))
		{
			$seq = $ltrs->{$class}->{$title}->{'seq'} . $seq . $ltrs->{$class}->{$title}->{'seq'};
		}
		my $rate = $cfg->{'rate'}->{$class}->{'relins'};
		my $finalrate = $rate * $globrate / $rsum;
		my $avgdel = int($cfg->{'rate'}->{$class}->{'avgdel'} * length($seq));
		my $stddev = int($cfg->{'rate'}->{$class}->{'stddev'} * length($seq));
		my $pct = $cfg->{'rate'}->{$class}->{'pct'};

		my $ratestr = sprintf("%.5f", $finalrate);
		my $fa_header = "$class.$title; rateperbase \"$ratestr\"; avgdel \"$avgdel\"; stddev \"$stddev\"; pct \"$pct\";";
		$outseqs->{$fa_header} = $seq;
	}
}

sub rpgpick
{
	my $max = shift;

	if ($max > 100)
	{
		$max = 100;
	}

	my $sum = 0;
	for (my $i = 1; $i <= $max; $i++)
	{
		$sum += 1.0 / ($i * $i);
	}

	my $roll = rand();

	$roll = $roll * $sum;

	my $cdf = 0;
	for (my $i = 1; $i <= $max; $i++)
	{
		$cdf += 1.0 / ($i * $i);
		if ($cdf >= $roll)
		{
			return $i;
		}
	}
	
	die("ERROR: $0: Numerical instability!");
}

sub getgenomesize
{
    # here, $genome must be the *name* of the genome (i.e., the string that describes
    # it, e.g., 'hg18'), not the path to the child cycle..
	my $fname = shift;
	my $genome = shift;
    
    if (! -e $fname){
        die("ERROR: $0: ERROR: $0: Unable to locate $fname.")
    }
	local *FH;
	open(FH, '<', $fname) || die("ERROR: $0: Cannot open $fname");
    my $gSize;
    my @s;
	while(<FH>)
	{
		chomp;
		@s = split(/\t/);
        $gSize = scalar(@s);
		($gSize != 2) && die("ERROR: $0: Invalid genome size from $fname ($gSize != 2)");

		if ($s[1] eq $genome)
		{
			return $s[0];
		}
	}
	close(FH) || die("ERROR: $0: Cannot close $fname");
	die("ERROR: $0: Invalid genome size (size:[@s], genome we seek:[$genome] file: $fname). Verify that genome names are the same.");
}

sub activate_rpg
{
	my $cfg = shift;
	my $branch = shift;
	my $ingff = shift;
	my $inseq = shift;
	my $ingen = shift;
	my $seqs = shift;

	my %gff;
	
	my $rpgs_left = int($cfg->{'rpg_count'} * $branch);
	($rpgs_left == 0) && return;

	print("* Distributing a total of $rpgs_left RPG copies...\n");

	read_genomic_gff($ingff, \%gff);
	filter_length($cfg->{'rpg_maxsize'}, \%gff);
	
	cmd("$CVT_BIN -showgenomesizes $inseq 2>/dev/null >rpg_tmp_genomesizes.txt");
	my $genomesize = getgenomesize('rpg_tmp_genomesizes.txt', $ingen);
	unlink('rpg_tmp_genomesizes.txt');

	my $rpg_index = 1;
	while ($rpgs_left > 0)
	{
		my $copy_count = rpgpick($rpgs_left);
		$rpgs_left -= $copy_count;
		output_random_gene('rpg_tmp.gff', \%gff);
		cmd("$CVT_BIN -xgffseqs $inseq -gff rpg_tmp.gff -out rpg_tmp.fa -genome $ingen 2>/dev/null");

		my $s;
		read_fasta_flat('rpg_tmp.fa', \$s);
		$s .= 'A' x $cfg->{'rpg_polya'};

		my $l = length($s);
		
		unlink('rpg_tmp.gff');
		unlink('rpg_tmp.fa');
		
		my $rate = $copy_count / ($genomesize * $branch);
		my $avgdel = int($cfg->{'rpg_avgdel'} * $l);
		my $stddev = int($cfg->{'rpg_stddev'} * $l);
		my $pct = $cfg->{'rpg_pct'};
		print("  - RPG #$rpg_index: $copy_count copies (rate=$rate)\n");

		my $header = "RPG.$rpg_index; rate \"$rate\"; avgdel \"$avgdel\"; stddev \"$stddev\"; pct \"$pct\";";
		
		$seqs->{$header} = $s;
		++$rpg_index;
	}
	close(FH);
}

sub randomdna
{
	my $length = shift;

	my @DNA = ('A', 'C', 'G', 'T');
	my $s = "";

	for (my $i = 0; $i < $length; ++$i)
	{
		$s .= $DNA[int(rand(4))];
	}

	return $s;
}

sub boundcheck
{
	my $cfg = shift;
	my $mes = shift;
	my $ltrs = shift;
	my $gff = shift;

	foreach my $class (keys(%$mes))
	{
		my $med = int(($cfg->{'rate'}->{$class}->{'minlen'} + $cfg->{'rate'}->{$class}->{'maxlen'}) / 2.0);
		foreach my $title (keys(%{$mes->{$class}}))
		{
			my $del = 0;
			my $header = make_header($mes, $class, $title);
			my $l = length($mes->{$class}->{$title}->{'seq'});
			
			if ($l < $cfg->{'rate'}->{$class}->{'minlen'})
			{
				my $n = $med - $l;
				$mes->{$class}->{$title}->{'seq'} .= randomdna($n);
				print("  - Appending $n nucleotides to ME $class.$title\n");
				fate($class, $title, "append:$n");
			}
			elsif ($l > $cfg->{'rate'}->{$class}->{'maxlen'})
			{
				if (defined($gff->{$header}))
				{
					print("  - Deleting ME $class.$title as its length ($l) is out of range (", $cfg->{'rate'}->{$class}->{'minlen'}, "-", $cfg->{'rate'}->{$class}->{'maxlen'}, ") and it has GFF entries.\n");
					fate($class, $title, "deleted:rangegff");
					if ($mes->{$class}->{$title}->{'active'})
					{
						print("  - WARNING: An active element just got deleted.\n");
					}
					$del = 1;
					delete($mes->{$class}->{$title});
					delete($gff->{$header});
					if (defined($ltrs->{$class}))
					{
						delete($ltrs->{$class}->{$title});
					}
					next;
				}
				else
				{
					my $n = $l - $med;
					$mes->{$class}->{$title}->{'seq'} = substr($mes->{$class}->{$title}->{'seq'}, $n);
					print("  - Deleting $n nucleotides from ME $class.$title\n");
					fate($class, $title, "remove:$n");
				}
			}

			if (defined($cfg->{'ltr'}->{$class}))
			{
				my $ll = length($ltrs->{$class}->{$title}->{'seq'});
				my $lmed = int(($cfg->{'ltr'}->{$class}->{'minlen'} + $cfg->{'ltr'}->{$class}->{'maxlen'}) / 2.0);
				if ($ll < $cfg->{'ltr'}->{$class}->{'minlen'})
				{
					my $n = $lmed - $ll;
					$ltrs->{$class}->{$title}->{'seq'} .= randomdna($n);
					print("  - Appending $n nucleotides to LTR $class.$title\n");
				}
				elsif ($ll > $cfg->{'ltr'}->{$class}->{'maxlen'})
				{
					my $n = $ll - $lmed;
					$ltrs->{$class}->{$title}->{'seq'} = substr($ltrs->{$class}->{$title}->{'seq'}, $n);
					print("  - Deleting $n nucleotides from LTR $class.$title\n");
				}
			}
		}
		if (scalar(keys(%{$mes->{$class}})) == 0)
		{
			print("* Deleting class $class because all its elements got out of range...\n");
			delete($mes->{$class});
			delete($ltrs->{$class});
		}
	}
}

sub validate_branch
{
	my $cfg = shift;
	my $branch = shift;
	foreach my $class (keys(%{$cfg->{'rate'}}))
	{
		my $duprate = $cfg->{'rate'}->{$class}->{'duprate'};
		my $delrate = $cfg->{'rate'}->{$class}->{'delrate'};
		if ((($duprate + $delrate) * $branch) > 1)
		{
			print("WARNING: Branch too long! The calculated birth+death probability for class ", $class, " exceeds 1.0 -- it will be adjusted!\n");
		}
		if ($cfg->{'rate'}->{$class}->{'lifetime'} < $branch)
		{
			print("WARNING: Branch too long! The lifetime for class ", $class, " is shorter (", $cfg->{'rate'}->{$class}->{'lifetime'}, ").\n");
		}
	}
}

sub main
{
	STDOUT->autoflush(1);
	STDERR->autoflush(1);
	if(@_){ 
        usage("Wrong Number of Arguments!\n");
    }
#	my $genome = shift;
#	my $branch = shift;
# 	$evo = shift;
# 	$cvt = shift;
# 	$py = shift;
#   print "evo: $evo\ncvt: $cvt\npy: $py\n";
    my $branch = $STEPSIZE;

    ##############################
	#my $cfgfile = 'ME.cfg';
    #my $in_fa_file = 'ME_input.fa';
	#my $in_ltr_file = 'ME_input_ltrs.fa';
	#my $in_gff_file = 'ME_input.gff';
	#my $model_file = 'ME.model.txt';
    my $cfgfile     = $MES_CFG;
    my $in_fa_file  = $ME_FA;
    my $in_ltr_file = $LTR_FA;
    my $in_gff_file = $ME_GFF;
    my $model_file  = $MODEL_TXT;
    
	my $tmp_fa_file          = 'ME_intermediate.fa';
	my $tmp_ltr_file         = 'ME_intermediate_ltrs.fa';
	my $tmp_unfixed_gff_file = 'ME_intermediate.unfixed.gff';
	my $tmp_gff_file         = 'ME_intermediate.gff';
    
	my $evo_fa_file          = 'ME_evolved.fa';
	my $evo_ltr_file         = 'ME_evolved_ltrs.fa';
	my $evo_unfixed_gff_file = 'ME_evolved.unfixed.gff';
	my $evo_gff_file         = 'ME_evolved.gff';

	my $out_fa_file    = 'ME_output.fa';
	my $out_ltr_file   = 'ME_output_ltrs.fa';
	my $out_gff_file   = 'ME_output.gff';
	my $active_fa_file = 'mes.fa';

    ####
    # these infiles correspond to the parent genome.
	my $ingff = "$PARENT_DIR/annots.gff";
	my $inseq = "$PARENT_DIR/seq.rev";

	my %conf;
	my %seqs;
	my %mes;
	my %ltrseqs;
	my %ltrs;
	my %gff;

	read_conf($cfgfile, \%conf);
	validate_branch(\%conf, $branch);
	read_fasta($in_fa_file, \%seqs);
	decompose(\%seqs, \%mes);
	validate_mes(\%mes, \%conf);
	read_fasta($in_ltr_file, \%ltrseqs);
	decompose(\%ltrseqs, \%ltrs);
	validate_mes(\%ltrs, \%conf);
	validate_ltrs(\%mes, \%ltrs, \%conf);
	read_gff($in_gff_file, \%gff);
	validate_gff(\%gff, \%seqs);

	print("* Loaded the following elements from the input:\n");
	print_mes(\%mes);

	foreach my $class (keys(%mes))
	{
		foreach my $title (keys(%{$mes{$class}}))
		{
			fate($class, $title, "init");
		}
	}

	my %outmes;
	my %outseqs;
	my %outltrs;
	my %outltrseqs;
	my %outgff;

	birth_death(\%conf, \%mes, \%ltrs, \%gff, \%outmes, \%outltrs, \%outgff, $branch);
	
	recompose(\%outmes, \%outseqs);
	write_fasta($tmp_fa_file, \%outseqs);
	recompose(\%outltrs, \%outltrseqs);
	write_fasta($tmp_ltr_file, \%outltrseqs);
	write_gff($tmp_unfixed_gff_file, \%outgff);
	fix_geneix($tmp_unfixed_gff_file, $tmp_gff_file);

	print("* Mobile elements after the birth-death process:\n");
	print_mes(\%outmes);

	evolve_mes(\%conf, $branch, \%outmes, $tmp_gff_file, $tmp_fa_file, \%outltrs, $tmp_ltr_file, $evo_unfixed_gff_file, $evo_fa_file, $evo_ltr_file, $model_file);
	fix_geneix($evo_unfixed_gff_file, $evo_gff_file);

	my %finalmes;
	my %finalseqs;
	my %finalltrs;
	my %finalltrseqs;
	my %finalgff;

	read_fasta($evo_fa_file, \%finalseqs);
	decompose(\%finalseqs, \%finalmes);
	validate_mes(\%finalmes, \%conf);
	read_fasta($evo_ltr_file, \%finalltrseqs);
	decompose(\%finalltrseqs, \%finalltrs);
	validate_mes(\%finalltrs, \%conf);
	validate_ltrs(\%finalmes, \%finalltrs, \%conf);
	read_gff($evo_gff_file, \%finalgff);
	validate_gff(\%finalgff, \%finalseqs);

	boundcheck(\%conf, \%finalmes, \%finalltrs, \%finalgff);
	remove_excess(\%conf, \%finalmes, \%finalltrs, \%finalgff);
	%finalseqs = ();
	%finalltrseqs = ();
	
	my %activeseqs;
	activate(\%conf, $branch, \%finalmes, \%finalltrs, \%finalgff, \%activeseqs);

	recompose(\%finalmes, \%finalseqs);
	write_fasta($out_fa_file, \%finalseqs);
	recompose(\%finalltrs, \%finalltrseqs);
	write_fasta($out_ltr_file, \%finalltrseqs);
	write_gff($out_gff_file, \%finalgff);

	activate_rpg(\%conf, $branch, $ingff, $inseq, $GENOME_NAME, \%activeseqs);
	write_fasta($active_fa_file, \%activeseqs);

	foreach my $key (keys(%fate))
	{
		print(STDERR $key, " ", join(',', @{$fate{$key}}), "\n");
	}
}

main(@ARGV);
