#!/usr/bin/perl -w
use strict;
use Cwd qw(cwd);

$|=1;

my $corr_seq_gaps   = 0;
my $make_reg_labels = 0;
   $make_reg_labels = 0 if $corr_seq_gaps == 0;

my $OAS_txt = shift(@ARGV);
   $OAS_txt =~ s/\\/\//g;

my ($OAS_path, $OAS_name);
if($OAS_txt =~ /\//){   ($OAS_path) = $OAS_txt =~ /(.*)\/.*/   } else {   ($OAS_path) = cwd;   }
if($OAS_txt =~ /\//){   ($OAS_name) = $OAS_txt =~ /.*\/(.*)\.txt/i   } else {   ($OAS_name) = $OAS_txt =~ /(.*)\.txt/i   }
print "$OAS_name in $OAS_path\n\n";


open(OAS, $OAS_txt) || die $!;

my $MUT_txt = $OAS_txt;
	$MUT_txt =~ s/\.txt/-MUT\.txt/i;

if($OAS_txt eq $MUT_txt){ 
	print "!! $OAS_txt is not a txt file !!\n"; exit;
} else {
	open( MUT, ">$MUT_txt") || die "Cannot write $MUT_txt!\n";
					#print MUT join("\t", ("OASid", "v_call", "d_call", "j_call", "region", "pos", "aa_gln", "aa_seq"))."\n";
	print MUT join("\t", ("OASid", "region", "pos", "aa_gln", "aa_seq", "trypepid"))."\n";
}

# NEW 28. Oct 2022: Write regions to REG

my $REG_txt = $OAS_txt;
	$REG_txt =~ s/\.txt/-REG\.txt/i;

if($OAS_txt eq $REG_txt){ 
	print "!! $OAS_txt is not a txt file !!\n"; exit;
} else {
	open( REG, ">$REG_txt") || die "Cannot write $REG_txt!\n";
					#print REG join("\t", ("OASid", "v_call", "d_call", "j_call", "region", "pos", "aa_gln", "aa_REG"))."\n";
	print REG join("\t", ("OASid", "v_call", "d_call", "j_call", "FWR1_aa", "CDR1_aa", "FWR2_aa", "CDR2_aa", "FWR3_aa", "CDR3_aa", "FWR4_aa"))."\n";
}

# NEW 28. Oct 2022: Write peptide~region alignment + MUT to PEP

my $PEP_txt = $OAS_txt;
	$PEP_txt =~ s/\.txt/-PEP\.txt/i;

if($OAS_txt eq $PEP_txt){ 
	print "!! $OAS_txt is not a txt file !!\n"; exit;
} else {
	open( PEP, ">$PEP_txt") || die "Cannot write $PEP_txt!\n";
					#print PEP join("\t", ("OASid", "v_call", "d_call", "j_call", "PEPion", "pos", "aa_gln", "aa_PEP"))."\n";
	print PEP join("\t", ("OASPEPid", "region", "reg_limits","overlap", "oltype", "olmut"))."\n";
}

my $OAS_fasta = $OAS_txt;
	$OAS_fasta =~ s/\.txt/\.fasta/i;

if($OAS_txt eq $OAS_fasta){ 
	print "!! $OAS_txt is not a txt file !!\n"; exit;
} else {
	open(FASTA, ">$OAS_fasta") || die "Cannot write $OAS_fasta!\n";
}

my $OAS_trypep = $OAS_txt;
	$OAS_trypep =~ s/\.txt/\.peptides\.txt/i;

if($OAS_txt eq $OAS_trypep){ 
	print "!! $OAS_trypep is not a txt file !!\n"; exit;
} else {
	open(ANALYSIS, ">$OAS_trypep") || die "Cannot write $OAS_trypep!\n";
	print ANALYSIS join(
	    "\t",
	    (   "pepid",        "pre4",  "pepseq", "post4", "pep_len", "curr_startpos", "curr_endpos", "pre_dist",
	        "post_dist", "n_Met", "n_Cys",  "n_Trp", "n_Gln", "n_Asn", "n_nt_Gln", "penalty", "aacomp", "MCWQN"
	    )
	) . "\n";

}

my	$OAS_trystat  = $OAS_path."/OAS-digest-stats.txt";

if($OAS_txt eq $OAS_trystat){ 
	print "!! $OAS_trystat is not a txt file !!\n"; exit;
} else {
	unless(-e $OAS_trystat){
		open(TRYSTAT, ">$OAS_trystat") || die "Cannot write $OAS_trystat!\n";
		print TRYSTAT join("\t", (   "OASds", "n_IGseq", "n_pep", "n_unique", "n_aacomp", "n_MCWQN" ) ). "\n";
		close(TRYSTAT);
	}

}


my $header = <OAS>; chomp($header);
my @H = split(/\t/, $header);

# checking columns
print "> Checking column nanes of $OAS_txt:\n";
my %col_names = (	0=>"OASid",
						1=>"locus",
						2=>"stop_codon",
						3=>"vj_in_frame",
						4=>"v_frameshift",
						5=>"productive",
						6=>"rev_comp",
						7=>"complete_vdj",
						8=>"v_call",
						9=>"d_call",
						10=>"j_call",
						11=>"sequence_alignment_aa",
						12=>"germline_alignment_aa",
						13=>"v_sequence_alignment_aa",
						14=>"v_germline_alignment_aa",
						15=>"j_sequence_alignment_aa",
						16=>"j_germline_alignment_aa",
						17=>"fwr1_aa",
						18=>"cdr1_aa",
						19=>"fwr2_aa",
						20=>"cdr2_aa",
						21=>"fwr3_aa",
						22=>"fwr4_aa",
						23=>"cdr3_aa",
						24=>"junction_aa",
						25=>"Redundancy");

foreach (0..$#H){
	
	if($H[$_] eq $col_names{$_}){
		#print " - [$_] ".$H[$_]." ~ ".$col_names{$_}." -> ok.\n";
	} else {
		print " - [$_] ".$H[$_]." ~ ".$col_names{$_}." -> ";
		print "check failed!\n";
		exit;
	}

}

#exit;
my $lcnt = 0;
my $report_step = 1;

my %len_mismatch = ();
my $gap_first_aa = 0;
my $entry_cnt = 0;
my $pep_count = 0;

my %key_aa20 = ();
my %key_MCWQN = ();
my %key_pepseq =();

while(<OAS>){
	$entry_cnt++;

	chomp;
	my @L = split(/\t/, $_);
	my @SEQ = split(//, $L[11]);
	my @GLN = split(//, $L[12]);
	my @MUT = ();

	my @POS_region = ();

	my @regions        = ("FWR1", "CDR1", "FWR2", "CDR2", "FWR3", "CDR3", "FWR4");
	my @peptides       = ($L[17], $L[18], $L[19], $L[20], $L[21], $L[23], $L[22] );
	my @region_lengths = (length($L[17]), length($L[18]), length($L[19]), length($L[20]), length($L[21]), length($L[23]), length($L[22]) );
	my @region_starts  = ();
	my @region_stops   = ();
	my @region_labels  = ();
	my $last           = 1;

	print REG join("\t", ($L[0], @L[8..10], @region_lengths) )."\n";

	my %aacomp_cache = ();

	# reading regions

	my @entry_description = ();

	foreach my $i (0..$#region_lengths){
		my $curr_start = $last;
		my $curr_stop  = $last + $region_lengths[$i] - 1;
		my $check_seq  = join("", @SEQ[($curr_start-1)..($curr_stop-1)]);

		unless($check_seq eq $peptides[$i]){
			if($corr_seq_gaps == 1){

				if($curr_stop <= $#SEQ){
					my $check_seq1 = join("", @SEQ[($curr_start-1+1)..($curr_stop-1+1)]);
					if($check_seq1 = $peptides[$i]){
						$curr_stop++;
						$region_lengths[$i]++;
						$peptides[$i] = join("", @SEQ[($curr_start-1)..($curr_stop-1)]);
						$gap_first_aa++;
						print "G"

					} else {
						print "Not fixed by assuming gap. :-(\n\n"; exit;
					}
				} else {
					print "!! CANNOT FIX THE PROBLEM; EXITING ...\n";
					exit;
				}
			} else {
				print "!! SeqLength does not correspond to sum of region lengths !!\n";
				print "!! Use read-OAS.pl with seq-check paramater.              !!\n";
				exit;
			}

		}

		push(@region_starts, $curr_start);
		push(@region_stops , $curr_stop );
		if($make_reg_labels == 1){
			if($region_lengths[$i] == 0){
				push(@region_labels, "");
			} elsif ($region_lengths[$i] == 1) {
				push(@region_labels, "|");
			} elsif ($region_lengths[$i] == 2) {
				push(@region_labels, "<>");
			} elsif ($region_lengths[$i] < 6){
				push(@region_labels, "<"."~"x($region_lengths[$i]-2).">");
			} else {
				my $label_n_front = int(($region_lengths[$i]-6)/2);
				my $label_n_tail  = $region_lengths[$i] - 6 - $label_n_front;
				push(@region_labels, "<"."~"x($label_n_front).$regions[$i]."~"x($label_n_tail).">");
			}
		}

		foreach (1..$region_lengths[$i]) {push(@POS_region, $regions[$i] )}

		push(@entry_description , $region_starts[$#region_starts]."-".$regions[$i]."-".$region_stops[$#region_stops]);
		$last = $curr_stop+1;
	}

	# Generatin peptide list/table
	my @tryclv = (1, &try_cleavage_sites4($L[11])) ;
	push(@tryclv, length($L[11]) ) if $tryclv[$#tryclv] != length($L[11]);

	my $misclv = 0;

	my $curr_startpos = 1;
	my %pep_to_pos = ();

	my %pepreg = ();
	my @pepregids = ();
	foreach my $itryclv (1..$#tryclv){
		$pep_count++;
		my $curr_endpos = $tryclv[$itryclv];
		my $pepid = join(".", ($L[0], $curr_startpos, $curr_endpos));

		my $pepseq = substr( $L[11], $curr_startpos-1, $curr_endpos - $curr_startpos + 1);
		my $pep_len = $curr_endpos - $curr_startpos + 1;


		# NEW VERSION 4:
		my $pre_len = $curr_startpos - 1;
		$pre_len = 4 if $pre_len > 4;

		my $post_len = length($L[11]) - $curr_endpos;
		$post_len = 4 if $post_len > 4;

		my $pre4 = substr( $L[11], $curr_startpos - $pre_len - 1, $pre_len );
		if ( $pre_len < 4 ) {
		  $pre4 = "-" x ( 4 - $pre_len ) . $pre4;
		}
		my $post4 = substr( $L[11], $curr_endpos, $post_len );
		if ( $post_len < 4 ) {
		  $post4 .= "-" x ( 4 - $post_len );
		}

		my $pre_dist = -1;
		$pre_dist = $curr_startpos - $tryclv[ $itryclv - 2 ] if $itryclv > 1;

		my $post_dist = -1;
		$post_dist = $tryclv[$itryclv + 1] - $curr_endpos if $itryclv < $#tryclv;

		my $n_Met    = ( $pepseq =~ tr/M// );
		my $n_Cys    = ( $pepseq =~ tr/C// );
		my $n_Trp    = ( $pepseq =~ tr/W// );
		my $n_Gln    = ( $pepseq =~ tr/Q// );
		my $n_Asn    = ( $pepseq =~ tr/N// );
		my $n_nt_Gln = 0;
		$n_nt_Gln = 1 if $pepseq =~ /^Q/;
		my $penalty =  $n_Met + $n_Cys + $n_Trp + $n_Gln + $n_Asn + $n_nt_Gln;

		
		my $curr_key_MCWQN = join("-", ($n_Met, $n_Cys, $n_Trp, $n_Gln, $n_Asn));
		$key_MCWQN{$curr_key_MCWQN}++;

		my ($aacomp);

		if(defined $aacomp_cache{$pepseq} ){
			$aacomp = $aacomp_cache{$pepseq};
		} else {
			$aacomp = &aacomp($pepseq);
		}

		$key_aa20{$aacomp}++;
		$key_pepseq{$pepseq}++;

		print ANALYSIS join(
		  "\t",
		  (   $pepid,        $pre4,  $pepseq, $post4, $pep_len, $curr_startpos, $curr_endpos, $pre_dist,
		      $post_dist, $n_Met, $n_Cys,  $n_Trp, $n_Gln, $n_Asn, $n_nt_Gln, $penalty, $aacomp, $curr_key_MCWQN
		  )
		) . "\n";


		foreach my $seqpos ($curr_startpos..$curr_endpos){
			#print "   - $seqpos ~ pep = ".join(".", ($L[0], $curr_startpos, $curr_endpos))."\n";
			$pep_to_pos{$seqpos} = $pepid;
		}

		foreach my $ireg (0..$#region_starts){
			my $c_reg_start = $region_starts[$ireg];
			my $c_reg_stop  = $region_stops[$ireg];
			my $reg_limits = join("-", ($c_reg_start, $c_reg_stop));
			my $pepregid = join("~", ($pepid, $regions[$ireg]));

			if ($curr_startpos <= $c_reg_start && $curr_endpos >= $c_reg_stop){
				push(@pepregids, $pepregid);
				#print join("\t", ($pepid, $regions[$ireg], $reg_limits,$c_reg_stop - $c_reg_start + 1), "FR")."\n";
				$pepreg{$pepregid} = join("\t", ($pepid, $regions[$ireg], $reg_limits,$c_reg_stop - $c_reg_start + 1), "FR");
			} elsif ($curr_startpos >= $c_reg_start && $curr_startpos <= $c_reg_stop && $curr_endpos >= $c_reg_start && $curr_endpos <= $c_reg_stop){
				push(@pepregids, $pepregid);
				#print join("\t", ($pepid, $regions[$ireg], $reg_limits,$curr_endpos - $curr_startpos + 1), "FP")."\n";
				$pepreg{$pepregid} = join("\t", ($pepid, $regions[$ireg], $reg_limits,$curr_endpos - $curr_startpos + 1), "FP");
			} elsif ($curr_startpos >= $c_reg_start && $curr_startpos <= $c_reg_stop){
				push(@pepregids, $pepregid);
				#print join("\t", ($pepid, $regions[$ireg], $reg_limits, $c_reg_stop - $curr_startpos + 1), "R")."\n";
				$pepreg{$pepregid} = join("\t", ($pepid, $regions[$ireg], $reg_limits, $c_reg_stop - $curr_startpos + 1), "R");
			} elsif ($curr_endpos >= $c_reg_start && $curr_endpos <= $c_reg_stop){
				push(@pepregids, $pepregid);
				#print join("\t", ($pepid, $regions[$ireg], $reg_limits, $curr_endpos - $c_reg_start + 1), "L")."\n";
				$pepreg{$pepregid} = join("\t", ($pepid, $regions[$ireg], $reg_limits, $curr_endpos - $c_reg_start + 1), "L");
			}
		}

		$curr_startpos = $curr_endpos + 1;
	}

	my %pepreg_mut = ();
	if($#POS_region != $#SEQ){
		print "> Incomplete region assignment:\n";
		print join("", @region_labels)."\n" if $make_reg_labels == 1;
		print $L[12]."\n";
		print $L[11]."\n";

	} elsif ($#SEQ == $#GLN){
		foreach my $ipos (0..$#SEQ){
			my $pos = $ipos + 1;
			if($SEQ[$ipos] ne $GLN[$ipos]){
				push(@MUT, $SEQ[$ipos]);
				print MUT join("\t", ($L[0], $POS_region[$ipos], $pos, $GLN[$ipos], $SEQ[$ipos], $pep_to_pos{$pos}))."\n";

				my $pepregid = $pep_to_pos{$pos}."~".$POS_region[$ipos];
				$pepreg_mut{$pepregid}++;
			} else {
#				print ".";
			}
		}
	} else {
		print "D";
		my $diff = $#GLN-$#SEQ;
		$len_mismatch{$diff}++;
	}

	foreach my $pepregid (@pepregids){
		#print "O [$pepregid]".$pepreg{$pepregid}."\t";
		my $n_olmut = 0;
		   $n_olmut = $pepreg_mut{$pepregid} if defined $pepreg_mut{$pepregid};

		print PEP join("\t", ($pepreg{$pepregid}, $n_olmut) )."\n";
	}

	my $accnr = $L[0];
	   $accnr =~ s/\:\:/-/g;
	my $VDR = join("_", @L[8..10]);
		$VDR =~ s/__/_/g;

	print FASTA ">IG|$accnr|$VDR ".join(":", @entry_description)."\n";
	print FASTA $L[11]."\n";

   if($entry_cnt >= $report_step){
      print ".";   
      if($report_step < 2**13) {
         $report_step *= 2;   
      } else {
         $report_step += 2**13;   
      }
      
   }


}

print ".. done. [$entry_cnt entries; len mismatches: ";
for my $diff (sort {$a <=> $b} keys(%len_mismatch)){
	print "diff=$diff->n=".$len_mismatch{$diff}."; ";
}

print " excluded; AA1 gaps = $gap_first_aa].\n";
close FASTA;
close MUT;

print " = Keys: ".scalar(keys(%key_aa20))." AA20 keys and ".scalar(keys(%key_MCWQN))." MCWQN keys for ".scalar(keys(%key_pepseq))." pep seqs.\n";
open(TRYSTAT, ">>$OAS_trystat") || die "Cannot write $OAS_trystat!\n";
print TRYSTAT join("\t", (   $OAS_name, $entry_cnt, $pep_count, scalar(keys(%key_pepseq)), scalar(keys(%key_aa20)), scalar(keys(%key_MCWQN)) ) ). "\n";
close(TRYSTAT);
print " = stats written to $OAS_trystat.\n";
print ". done\n";

#  ___ _   _ ___ ___
# / __| | | | _ ) __|
# \__ \ |_| | _ \__ \
# |___/\___/|___/___/


sub try_cleavage_sites4 {

    my @seq = split(//, $_[0]);
    my @clv_pos = ();
    push(@seq, ".");
    foreach my $i (1..($#seq-1)){
      push(@clv_pos, $i+1) if(($seq[$i] eq "K" || $seq[$i] eq "R") &&   $seq[$i+1] ne "P" );
    }
    return sort { $a <=> $b } @clv_pos;
}


sub aacomp {
     my ($aacomp);
     my %CNT = ("A" => 0, "C" => 0, "D" => 0, "E" => 0, "F" => 0, "G" => 0, "H" => 0, "I" => 0, "K" => 0, "L" => 0, "M" => 0, "N" => 0, "P" => 0, "Q" => 0, "R" => 0, "S" => 0, "T" => 0, "V" => 0, "W" => 0, "Y" => 0);
     my @PEP = split(//, $_[0]);
     foreach my $aa (@PEP){  $CNT{"$aa"}++  }
	  my @return;
	  foreach my $aa (sort{$a cmp $b} keys(%CNT)){
	  	push(@return, $CNT{$aa});
	  }
     return join("-", @return);
  }
