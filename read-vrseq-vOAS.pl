#!/usr/bin/perl -w
use strict;
$|=1;
my $parse_ANARCI = 0;
my $run_flg = 1;
my $skip_flg = 0;
my $seq_chk = 0;
my $clean_up = 0;

use Cwd qw(cwd);
my $currdir = cwd;


my @OAS_index_locs = ("dataset.index");

my ($OASID_index_file);
foreach my $loc ( @OAS_index_locs ){
   $OASID_index_file = $loc if -e $loc;
}
 

my %OASID_index = ();
if(-e $OASID_index_file){
   print "< reading OASID index (filename -  OASid): `$OASID_index_file`\n";
   open(OASIDX, $OASID_index_file) || die "Cannot read $OASID_index_file!\n";
   <OASIDX>;
   while(<OASIDX>){
      chomp;
      my ($f, $id, $oasid) = split(/\t/, $_);
      $OASID_index{$f} = $oasid;
   }
   close OASIDX;
} else {
   print "!! OAS index `OASID.index` (filename <-> OASid) expected in working directory cannot be read!!\n";
   $OASID_index{"NA"} = ".";
   #exit;
}


my $z7 = "c:\\Progra~1\\7-Zip\\7z.exe e #OASgz -o#OASpath";

my @entry_columnns = ("Run", "Link", "Author", "Species", "Longitudinal", "Age", "Disease", "Vaccine", "BType", "BSource", "Subject", "Chain", "Unique sequences", "Isotype", "Total sequences");

my @colex = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 22, 24, 30, 32, 34, 36, 38, 40, 42, 44, 46, 49, 94); # , 95, 96 not on list -> ANARCI entries

my ($OAS_dataset_log_file);
my @OAS_dataset_log_line = ();;

my $sample = 0;
my %SAMPLE = ();
my %POOL = ();
my $cnt_red_gr1 = 0;

foreach (0..$#ARGV){
   $ARGV[$_] =~  s/\\/\//g;
   if($ARGV[$_] =~ /skip/i ){
      $run_flg = 0;
      print " + PARA: skip-parsing active: no reading and parsing of data table; just reading of file information.\n";
   } elsif($ARGV[$_] =~ /anarci/i ){
      $parse_ANARCI = 1;
      print " + PARA: ANARCI parsing ON.\n";
      print "!! ANARCI PARSING IS NOT TESTED YET !!\n";
      <STDIN>;
   } elsif($ARGV[$_] =~ /seq/i && $ARGV[$_] =~ /check/i){
      print " + PARA: check and exclude sequences = ON.\n";
      $seq_chk = 1;
   } elsif($ARGV[$_] =~ /sample=/i){
      ($sample) = $ARGV[$_] =~ /sample=(.*)/i;
      print " + sampling set of $sample entries (red>1) <- ".$ARGV[$_]."\n";
      exit if $sample < 1;
   } elsif ($ARGV[$_] =~ /clean/i && $ARGV[$_] =~ /up/i){
      print " + PARA: cleanup = ON.\n";
      $clean_up = 1;
   }
}

my %meta = ();
my ($OAS_file, $out_file);
my $input_is_archive = 0;
my ($OASgz);
my @cleanup_files = ();

my $entry_index = "";

my ($curr_OASNAME, $curr_OASID, $OAStxt, $OAStxtN);

if(defined $ARGV[0] && -e $ARGV[0]){
   $OAS_file = shift(@ARGV);

   unless($OAS_file =~ /\//){
      $OAS_file = $currdir."/".$OAS_file;
   }


   print " = Input file `$OAS_file`\n";
   my ($OAS_path) =  $OAS_file =~ /(.*)\/.*/;

   if($OAS_file =~ /\.csv\.gz$/i){
      print " = Input file is packed archive.\n";
      $OASgz = $OAS_file;
      $OAS_file =~ s/\.csv\.gz$/\.csv/i;

      print " = OAS csv File: $OAS_file";
      print "(exists)" if -e $OAS_file;
      print ".\n";

      $out_file = $OAS_file;
      $out_file =~ s/\.csv$/\.##tbl\.txt/i;
      $OAStxt = $out_file;
      $OAStxt =~ s/##tbl/OAS/;
      
      print " = OAS TXT File: $OAStxt";
      print "(exists)" if -e $OAStxt;
      print ".\n";

      $OAStxtN = $OAStxt;
      $OAStxtN =~ s/\.txt/.N$sample\.txt/i;
      if($OAStxtN eq $OAStxt){
         print "ERROR creating valid OAStxtN for subsampling of $sample entries.\n";
         exit
      }

      if($sample > 0){
         print " = OAS TXT SubSampling File: $OAStxtN";
         print "(exists)" if -e $OAStxtN;
         print ".\n";
      }
   }

   my ($OAS_dataset) =  $OAS_path =~ /.*\/(.*)/;
      ($curr_OASNAME) =  $OAS_file =~ /.*\/(.*)\./;

   print " = OAS dataset (from dir): $OAS_dataset\n";
   print " = OAS name (from filename): $curr_OASNAME\n";

   # READING EXISTING ENTRIES FIRST

   $entry_index = $OAS_path."/_OAS-entry-index.txt";

   if(defined $OASID_index{$curr_OASNAME}){
      $curr_OASID = $OASID_index{$curr_OASNAME};
      print "< $curr_OASID from $curr_OASNAME (lookup in $OASID_index_file)\n";
   } else{
      $curr_OASID = $curr_OASNAME;
   }

   $OAS_dataset_log_file = $OAS_path."/".$OAS_dataset.".read-OAS.log";
   print " > logging to $OAS_dataset_log_file.\n";
   unless(-e $OAS_dataset_log_file){
      print "> writing new log file template (`$OAS_dataset_log_file`).\n";
      open(OASLOG, ">$OAS_dataset_log_file") || die "Cannot write to $OAS_dataset_log_file!\n";
      print OASLOG join("\t", ("dataset", "OASid", "name", "file", "n_entries", "n_seqchk_failed", "n_red_gr1"))."\n";
      close OASLOG;
   }
   push(@OAS_dataset_log_line, ($OAS_dataset, $curr_OASID, $curr_OASNAME, $OAS_file));

   print " + ENTRY INDEX: $entry_index\n";
   if(-e $entry_index){
      print "< reading existing OAS ENTRY index (`$entry_index`).\n";
      open(IN, $entry_index) || die "Cannot read entry file $entry_index!\n";
      <IN>;

      while(<IN>){
         chomp;
         my($id, undef) = split(/\t/, $_);
         if($id eq $curr_OASID){
            $skip_flg = 1;
            print " + $id ALREADY IN entry_index !!\n";
         }
         
      }
      close IN;

   } else {
      print "> writing new ENTRY INDEX template (`$entry_index`).\n";
      open(OUT, ">".$entry_index) || die "Cannot write entry file $entry_index!\n";
      print OUT join("\t", ("OASid", "OASname", @entry_columnns))."\n";
      close OUT;
   }

   # unpacking, if required

   if(-e $OAStxt && $skip_flg == 1){
      print "> OAS text file $OAStxt already exists and entry indexed. Skipping unpacking.\n";
   } elsif (-e $OAS_file){
      print "> Unpacked file $OAS_file already exists. Unpacking skipped.\n";
   } else {
      print "# Unpacking $OASgz ...";
      
      my $cmd_7z = $z7;
         $cmd_7z =~ s/#OASgz/$OASgz/;
         $cmd_7z =~ s/#OASpath/$OAS_path/;
      my @log_7z = qx($cmd_7z);
      if(-e $OAS_file){
         print ". done.\n";
         push(@cleanup_files, $OAS_file);
      } else {
         print ". FAILED!\n";
         print join("\n", @log_7z)."\n";
         exit;
      }

   }
} else {
   print "!! Usage: read-OAS.pl oas-file.csv\n\n";
   if(defined $ARGV[0]){
      print " --> cannot access ".$ARGV[0]."\n";
   } else {
      print " --> no input file defined.\n";
   }
   
   exit;
}

if(-e $OAStxt && $skip_flg == 1 && ($sample > 0 && -e $OAStxtN )){
   print ". finished (delete $OAStxt is reprocessing is required).\n";
   exit;
}

open(OAS, $OAS_file) || die "Cannot read or access $OAS_file (OAS file).\n";

my $header = <OAS>;
   $header =~ s/[\{\}]//g;
   $header =~ s/et al\.\,/et al\./g;
   $header =~ s/\"//g;
chomp($header);

my @H = split(/\, /, $header);

print "# Reading OAS ENTRY INFORMATION (from `$OAS_file`) ... ";

my @entry_v = ();
my %entry = ();
foreach (0..$#H){
   $H[$_] =~ s/\"//g;
   my ($p, $v) = $H[$_] =~ /(.*)\: (.*)/;
   $v =~ s/^\s+//;
   $v =~ s/\s+$//;
   $entry{$p} = $v;
}

for my $e (0..$#entry_columnns){
   if(defined $entry{$entry_columnns[$e]}){
      print ".";
   } else {
      print " [$e] ".$entry_columnns[$e]." -> ";
      print "ENTRY NOT FOUND.\n";
      print "(in ".join(", ", keys(%entry)).")\n";
      exit;
   }
   push(@entry_v, $entry{$entry_columnns[$e]});
}
print " done.\n";

if($skip_flg == 0){
   print "> adding ENTRY INFORMATION to index (`$entry_index`)\n";
   open(OUT, ">>".$entry_index) || die "Cannot write entry file $entry_index!\n";
   print OUT join("\t", ($curr_OASID, $curr_OASNAME, @entry_v))."\n";
   close OUT;
}

print "# Reading OAS DATA $curr_OASID (from `$OAS_file`) ... ";

$header = <OAS>;
$header =~ s/[\{\}]//g;
chomp($header);
@H = split(/\,/, $header);

my $lcnt = 0; 
my $ecnt = 0; 
my $report_step = 1;
my $OAS_entry_cnt = 0;

my $curr_out = $out_file;
   $curr_out =~ s/##tbl/OAS/;

if(-e $OAStxt && ($sample > 0 && -e $OAStxtN )){
      $run_flg = 0;
      print "!! Skipping processing because OAStxt and subsampling files already exist!\n";
}

my $cnt_seq_check_failed = 0;

if($run_flg == 1 && $parse_ANARCI == 0){
   if($seq_chk == 0){
      open(OUT, ">".$curr_out) || die "Cannot write entry file $curr_out!\n";
      print OUT join("\t", ("OASid", @H[@colex]))."\n";

      while (<OAS>){
         $OAS_entry_cnt++;
         my $curr_ID =  $curr_OASID."::".sprintf("%06s", $OAS_entry_cnt);
         chomp;
         $lcnt++;
         my($line1t94, $ANARCI, $ANARCHIstat) = $_ =~ /(.*)\,\"\{(.*)\}\}\"\,(.*)/;
         my @L = split(/\,/, $line1t94);
         print OUT join("\t", ($curr_ID, @L[@colex]))."\n";

         if($L[94] > 1){
            $cnt_red_gr1++;
            $SAMPLE{$OAS_entry_cnt} = $L[94];
         } else {
            $POOL{$OAS_entry_cnt} = $L[94];
         }

      }
   } else {
      print "> Writing to OAS.txt ($curr_out)\n";
      open(OUT, ">".$curr_out) || die "Cannot write entry file $curr_out!\n";
      print OUT join("\t", ("OASid", @H[@colex]))."\n";

      my $curr_excluded = $curr_out;
         $curr_excluded =~ s/\.txt$/\.fail\.txt/i;
      
      print "> Preparing OAS.fail.txt ($curr_excluded)\n";
      open(FAIL, ">".$curr_excluded) || die "Cannot write entry file $curr_out!\n";
      print FAIL join("\t", ("OASid", @H[@colex]))."\n";


      while (<OAS>){
         $OAS_entry_cnt++;
         my $curr_ID =  $curr_OASID."::".sprintf("%06s", $OAS_entry_cnt);
         chomp;
         $lcnt++;
         my($line1t94, $ANARCI, $ANARCHIstat) = $_ =~ /(.*)\,\"\{(.*)\}\}\"\,(.*)/;

         my @L = split(/\,/, $line1t94);

         my $len_aligned = length($L[13]);
         my $len_regions = length($L[34]) + length($L[36]) + length($L[38]) + length($L[40]) + length($L[42]) + length($L[44]) + length($L[46]);

         if($len_aligned == $len_regions){

            print OUT join("\t", ($curr_ID, @L[@colex]))."\n";
            $ecnt++;
            if($L[94] > 1){
               $cnt_red_gr1++;
               $SAMPLE{$OAS_entry_cnt} = $L[94];
            } else {
               $POOL{$OAS_entry_cnt} = $L[94];
            }

         } else {
           $cnt_seq_check_failed++;
           print "!";
           print FAIL join("\t", ($curr_ID, @L[@colex]))."\n";
         }

      }

      close OUT;
      close FAIL;
      if($cnt_seq_check_failed == 0){
         push(@cleanup_files, $curr_excluded);
      }
   }
} elsif ($run_flg == 1 && $parse_ANARCI == 0) {

} else {
   print "...FILE ALREADY EXISTS-> SKIPPING...";
}


close OAS;

print "done ($lcnt).\n";

if($sample > 0){
   if($sample > $ecnt){
      print "!! Data file contains less entries ($ecnt) than sample size specificed (N = $sample) !!\n";
      $sample = $ecnt;
      print "!! sample reduced to number of entries (N = $sample) !!\n";
   }

   print "# Sampling $sample lines (red > 1 are inlcuded always; N = $cnt_red_gr1)\n";

   my @LINES = ();
   if($sample == $cnt_red_gr1){
      @LINES = sort {$a <=> $b} keys(%SAMPLE);
      print " - sampling all entries red > 1!\n";
      <STDIN>;
   } elsif ($sample < $cnt_red_gr1 ){
      @LINES = sort {$SAMPLE{$b} <=> $SAMPLE{$a}} keys(%SAMPLE);
      @LINES = @LINES[0..($sample-1)];
      print " - sampling $sample entries with highest redundancy.\n";

      <STDIN>;
   } elsif ($sample > $cnt_red_gr1){
      @LINES = sort {$a <=> $b} keys(%SAMPLE);
      my @rand_sample = sort {rand > 0.5} keys(%POOL);
      print "[Y] ".scalar(@rand_sample)."\n";
      my $n_rnd_samples = $#rand_sample+1;
      push(@LINES, @rand_sample);
      my $n_lines = $#LINES;
      print "[Y] ".$n_lines."\n";
      
      if($sample > $ecnt){
         print "!! Data file contains less entries ($lcnt) than sample size specificed (N = $sample) !!\n";
         <STDIN>;   
         @LINES = @LINES[0..($lcnt-1)];
      } else {
         @LINES = @LINES[0..($sample-1)];
      }

      
      print " - sampling $cnt_red_gr1 entries with redundancy > 1 and ".$n_rnd_samples." random samples.\n";

   }

my @undef_lines = ();
foreach (0..$#LINES){
   push(@undef_lines, $_) if !defined($LINES[$_]);
}

if(scalar(@undef_lines) > 0){
   print " !!> undefined lines: ".join(" ", @undef_lines)." (".scalar(@undef_lines).")\n";
   <STDIN>

}

   @LINES = sort {$a <=> $b} @LINES;
   open(TXT, $curr_out) || die "Cannot write entry file $curr_out!\n";
   open(SAMPLE, ">".$OAStxtN) || die "Cannot write sample file $OAStxtN!\n";
   
   my $H = <TXT>;
   print SAMPLE $H;
   my $OAS_entry_line = 0;
   my $scnt = 0;
   while(my $L = <TXT>){

      $OAS_entry_line++;
      if(scalar(@LINES) > 0 && $OAS_entry_line == $LINES[0]){
         shift(@LINES);
         print SAMPLE $L;
         $scnt++;
      }
   }

   close TXT;
   close SAMPLE;
   print " = $scnt lines written to $OAStxtN.\n";
}


push(@OAS_dataset_log_line, ($OAS_entry_cnt, $cnt_seq_check_failed));

if(-e $OAS_dataset_log_file){
   open(OASLOG, ">>$OAS_dataset_log_file") || die "Cannot write to $OAS_dataset_log_file!\n";
   push(@OAS_dataset_log_line, $cnt_red_gr1);
   print OASLOG join("\t", @OAS_dataset_log_line )."\n";
   close OASLOG;
   } else {
   print "Failed to write log:\n";
   print join("\t", @OAS_dataset_log_line)."\n";
}

if($clean_up == 1){
   foreach my $file_to_delete(@cleanup_files){
      print "! DELETING $file_to_delete\n";
      unlink($file_to_delete);
   }
} else {
   print "~ Cleanup of skipped (of ".join(", ", @cleanup_files).")\n";
}