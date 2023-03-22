
###############################################################################
 #
 #  This file is part of canu, a software program that assembles whole-genome
 #  sequencing reads into contigs.
 #
 #  This software is based on:
 #    'Celera Assembler' (http://wgs-assembler.sourceforge.net)
 #    the 'kmer package' (http://kmer.sourceforge.net)
 #  both originally distributed by Applera Corporation under the GNU General
 #  Public License, version 2.
 #
 #  Canu branched from Celera Assembler at its revision 4587.
 #  Canu branched from the kmer project at its revision 1994.
 #
 #  This file is derived from:
 #
 #    src/pipelines/ca3g/Defaults.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-FEB-27 to 2015-SEP-21
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-OCT-21
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Defaults;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(getCommandLineOptions addCommandLineOption addCommandLineError writeLog caExit caFailure getNumberOfCPUs getPhysicalMemorySize getAllowedResources diskSpace printOptions printHelp setParametersFromFile setParametersFromCommandLine checkParameters getGlobal setGlobal setGlobalIfUndef showErrorRates setErrorRate setDefaults);

use strict;
use Carp qw(cluck);
use Sys::Hostname;
use Text::Wrap;

use Filesys::Df;  #  for diskSpace()

my %global;
my %synops;
my %synnam;

my $cLineOpts = "";
my $specLog   = "";



#  Return the second argument, unless the first argument is found in
#  %global, in which case return that.
#
sub getGlobal ($) {
    my $var = shift @_;

    $var =~ tr/A-Z/a-z/;

    caFailure("parameter '$var' is not known", undef) if (!exists($global{$var}));

    return($global{$var});
}


sub setGlobalSpecialization ($@) {
    my $val = shift @_;

    foreach my $var (@_) {
        $global{$var} = $val;
    }

    return(1);
}


sub setGlobal ($$) {
    my $VAR = shift @_;
    my $var = $VAR;
    my $val = shift @_;
    my $set = 0;

    $var =~ tr/A-Z/a-z/;
    $val = undef  if ($val eq "");  #  Set to undefined, the default for many of the options.

    #  Map 'true'/'false' et al. to 0/1.

    $val = 0  if (($val =~ m/^false$/i) || ($val =~ m/^f$/i));
    $val = 1  if (($val =~ m/^true$/i)  || ($val =~ m/^t$/i));

    #  Translate from generic to specialized var

    foreach my $alg ("ovl", "mecat2asmpw") {
        foreach my $opt ("gridoptions") {
            $set += setGlobalSpecialization($val, ("${opt}cor${alg}", "${opt}obt${alg}", "${opt}utg${alg}"))  if ($var eq "${opt}${alg}");
        }

        foreach my $opt ("memory", "threads", "concurrency") {
            $set += setGlobalSpecialization($val, ("cor${alg}${opt}", "obt${alg}${opt}", "utg${alg}${opt}"))  if ($var eq "${alg}${opt}");
        }
    }

    foreach my $opt ("overlapper") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    #  e.g., corOvlHashBlockLength
    foreach my $opt ("ovlerrorrate", "ovlhashblocklength", "ovlrefblocksize", "ovlrefblocklength", "ovlhashbits", "ovlhashload", "ovlmersize", "ovlmerthreshold", "ovlmerdistinct", "ovlmertotal", "ovlfrequentmers") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    #  e.g., cormecat2asmpwBlockSize
    foreach my $opt ("mecat2asmpwblocksize", "mecat2asmpwmersize", "mecat2asmpwrealign", "mecat2asmpwsensitivity") {
        $set += setGlobalSpecialization($val, ("cor${opt}", "obt${opt}", "utg${opt}"))  if ($var eq "${opt}");
    }

    return  if ($set > 0);

    if ($var eq "errorrate") {
        setErrorRate($val);
        return;
    }

    #  If we got a parameter we don't understand, we should be parsing command line options or
    #  reading spec files, and we can let the usual error handling handle it.

    addCommandLineError("ERROR:  Paramter '$VAR' is not known.\n")   if (!exists($global{$var}));

    $global{$var} = $val;
}



sub setGlobalIfUndef ($$) {
    my $var = shift @_;
    my $val = shift @_;

    $var =~ tr/A-Z/a-z/;
    $val = undef  if ($val eq "");  #  Set to undefined, the default for many of the options.

    return  if (defined($global{$var}));

    $global{$var} = $val;
}



sub getCommandLineOptions () {
    return($cLineOpts);
}



sub addCommandLineOption ($) {
    if ($cLineOpts =~ m/\s$/) {
        $cLineOpts .= "$_[0]";
    } else {
        $cLineOpts .= " $_[0]";
    }
}



sub addCommandLineError($) {
    $global{'errors'} .= shift @_;
}



sub writeLog ($) {
    my $wrk = shift @_;

    my $time = time();
    my $host = hostname();
    my $pid  = $$;

    open(F, "> $wrk/canu-logs/${time}_${host}_${pid}_canu");
    print F $specLog;
    close(F);
}



#  Use caExit() for transient errors, like not opening files, processes that die, etc.
sub caExit ($$) {
    my  $msg   = shift @_;
    my  $log   = shift @_;

    print STDERR "================================================================================\n";
    print STDERR "Don't panic, but a mostly harmless error occurred and canu failed.\n";
    print STDERR "\n";

    #  Really should pass in $wrk
    if (defined($log)) {
        my  $df = diskSpace($log);

        print STDERR "Disk space available:  $df GB\n";
        print STDERR "\n";
    }

    if (-e $log) {
        print STDERR "Last 50 lines of the relevant log file ($log):\n";
        print STDERR "\n";
        system("tail -n 50 $log");
        print STDERR "\n";
    }

    print STDERR "canu failed with '$msg'.\n";
    print STDERR "\n";

    exit(1);
}


#  Use caFailure() for errors that definitely will require code changes to fix.
sub caFailure ($$) {
    my  $msg   = shift @_;
    my  $log   = shift @_;

    print STDERR "================================================================================\n";
    print STDERR "Please panic.  canu failed, and it shouldn't have.\n";
    print STDERR "\n";
    print STDERR "Stack trace:\n";
    print STDERR "\n";
    cluck;
    print STDERR "\n";

    if (-e $log) {
        print STDERR "Last few lines of the relevant log file ($log):\n";
        print STDERR "\n";
        system("tail -n 50 $log");
    }

    print STDERR "\n";
    print STDERR "canu failed with '$msg'.\n";

    exit(1);
}


#
#  Host management - these really belong in 'Execution.pm' (or 'Utilities.pm') but can't go there
#  (Execution.pm) and be used here too.
#

sub getNumberOfCPUs () {
    my $os   = $^O;
    my $ncpu = 1;

    #  See http://stackoverflow.com/questions/6481005/obtain-the-number-of-cpus-cores-in-linux

    if ($os eq "freebsd") {
        $ncpu = int(`/sbin/sysctl -n hw.ncpu`);
    }

    if ($os eq "darwin") {
        $ncpu = int(`/usr/bin/getconf _NPROCESSORS_ONLN`);
    }

    if ($os eq "linux") {
        $ncpu = int(`getconf _NPROCESSORS_ONLN`);
    }

    return($ncpu);
}


sub getPhysicalMemorySize () {
    my $os     = $^O;
    my $memory = 1;

    if ($os eq "freebsd") {
        $memory = `/sbin/sysctl -n hw.physmem` / 1024 / 1024 / 1024;
    }

    if ($os eq "darwin") {
        $memory = `/usr/sbin/sysctl -n hw.memsize` / 1024 / 1024 / 1024;
    }

    if ($os eq "linux") {
        open(F, "< /proc/meminfo");        #  Way to go, Linux!  Make it easy on us!
        while (<F>) {
            if (m/MemTotal:\s+(\d+)/) {
                $memory = $1 / 1024 / 1024;
            }
        }
        close(F);
    }

    return(int($memory + 0.5));  #  Poor man's rounding
}



sub dirname ($) {
    my $d = shift @_;

    return($d)  if (-d $d);

    my @d = split '/', $d;
    pop @d;

    $d = join('/', @d);

    return($d);
}



sub diskSpace ($) {
    my $wrk   = dirname($_[0]);
    my $df    = df($wrk, 1024);

    my $total = int(10 * $df->{blocks} / 1048576) / 10;
    my $used  = int(10 * $df->{used}   / 1048576) / 10;
    my $free  = int(10 * $df->{bfree}  / 1048576) / 10;
    my $avail = int(10 * $df->{bavail} / 1048576) / 10;

    #print STDERR "Disk space: total $total GB, used $used GB, free $free GB, available $avail GB\n";

    return (wantarray) ? ($total, $used, $free, $avail) : $avail;
}


sub printOptions () {
    my $pretty = 0;

    foreach my $k (sort values %synnam) {
        my $o = $k;
        my $u = $synops{$k};

        next   if (length($u) == 0);

        if ($pretty == 0) {
            $o = substr("$k                                    ", 0, 40);

        } else {
            $Text::Wrap::columns = 100;

            $o = "$o\n";
            $u = wrap("    ", "    ", $u) . "\n";
        }

        print "$o$u\n";
    }
}


sub printHelp () {

    return   if (!exists($global{'errors'}));

    print "\n";
    print "usage: canu [-correct | -trim | -assemble] \\\n";
    print "            [-s <assembly-specifications-file>] \\\n";
    print "             -p <assembly-prefix> \\\n";
    print "             -d <assembly-directory> \\\n";
    print "             genomeSize=<number>[g|m|k] \\\n";
    print "             errorRate=0.X \\\n";
    print "            [other-options] \\\n";
    print "            [-pacbio-raw | -pacbio-corrected | -nanopore-raw | -nanopore-corrected] *fastq\n";
    print "\n";
    print "  By default, all three stages (correct, trim, assemble) are computed.\n";
    print "  To compute only a single stage, use:\n";
    print "    -correct  - generate corrected reads\n";
    print "    -trim     - generate trimmed reads\n";
    print "    -assemble - generate an assembly\n";
    print "\n";
    print "  The assembly is computed in the (created) -d <assembly-directory>, with most\n";
    print "  files named using the -p <assembly-prefix>.\n";
    print "\n";
    print "  The genome size is your best guess of the genome size of what is being assembled.\n";
    print "  It is used mostly to compute coverage in reads.  Fractional values are allowed: '4.7m'\n";
    print "  is the same as '4700k' and '4700000'\n";
    print "\n";
    print "  The errorRate is not used correctly (we're working on it).  Set it to 0.06 and\n";
    print "  use the various utg*ErrorRate options.\n";
    print "\n";
    print "  A full list of options can be printed with '-options'.  All options\n";
    print "  can be supplied in an optional sepc file.\n";
    print "\n";
    print "  Reads can be either FASTA or FASTQ format, uncompressed, or compressed\n";
    print "  with gz, bz2 or xz.  Reads are specified by the technology they were\n";
    print "  generated with:\n";
    print "    -pacbio-raw         <files>\n";
    print "    -pacbio-corrected   <files>\n";
    print "    -nanopore-raw       <files>\n";
    print "    -nanopore-corrected <files>\n";
    print "\n";
    print "Complete documentation at http://canu.readthedocs.org/en/latest/\n";
    print "\n";
    print "$global{'errors'}";
    print "\n";

    exit(1);
}



sub makeAbsolute ($) {
    my $var = shift @_;
    my $val = getGlobal($var);
    if (defined($val) && ($val !~ m!^/!)) {
        $val = "$ENV{'PWD'}/$val";
        setGlobal($var, $val);
        $val =~ s/\\\"/\"/g;
        $val =~ s/\"/\\\"/g;
        $val =~ s/\\\$/\$/g;
        $val =~ s/\$/\\\$/g;

        addCommandLineOption("\"$var=$val\"");
    }
}



sub fixCase ($) {
    my $var = shift @_;
    my $val = getGlobal($var);

    if (defined($val)) {
        $val =~ tr/A-Z/a-z/;
        setGlobal($var, $val);
    }
}



sub setParametersFromFile ($@) {
    my $specFile  = shift @_;
    my @fragFiles = @_;

    #  Client should be ensuring that the file exists before calling this function.
    die "specFile '$specFile' not found.\n"  if (! -e "$specFile");

    $specLog .= "\n";
    $specLog .= "###\n";
    $specLog .= "###  Reading options from '$specFile'\n";
    $specLog .= "###\n";
    $specLog .= "\n";

    open(F, "< $specFile") or caExit("can't open '$specFile' for reading: $!", undef);

    while (<F>) {
        $specLog .= $_;

        s/^\s+//;
        s/\s+$//;

        next if (m/^\s*\#/);
        next if (m/^\s*$/);

        if (-e $_) {
            my $xx = $_;
            $xx = "$ENV{'PWD'}/$xx" if ($xx !~ m!^/!);
            if (-e $xx) {
                push @fragFiles, $xx;
            } else {
                addCommandLineError("ERROR:  File not found '$_' after appending absolute path.\n");
            }
        } elsif (m/\s*(\w*)\s*=([^#]*)#*.*$/) {
            my ($var, $val) = ($1, $2);
            $var =~ s/^\s+//; $var =~ s/\s+$//;
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            undef $val if ($val eq "undef");
            setGlobal($var, $val);
        } else {
            addCommandLineError("ERROR:  File not found or unknown specFile option line '$_'.\n");
        }
    }
    close(F);

    return(@fragFiles);
}



sub setParametersFromCommandLine(@) {
    my @specOpts = @_;

    if (scalar(@specOpts) > 0) {
        $specLog .= "\n";
        $specLog .= "###\n";
        $specLog .= "###  Reading options from the command line.\n";
        $specLog .= "###\n";
        $specLog .= "\n";
    }

    foreach my $s (@specOpts) {
        $specLog .= "$s\n";

        if ($s =~ m/\s*(\w*)\s*=(.*)/) {
            my ($var, $val) = ($1, $2);
            $var =~ s/^\s+//; $var =~ s/\s+$//;
            $val =~ s/^\s+//; $val =~ s/\s+$//;
            setGlobal($var, $val);
        } else {
            addCommandLineError("ERROR:  Misformed command line option '$s'.\n");
        }
    }
}



sub setExecDefaults ($$) {
    my $tag         = shift @_;
    my $name        = shift @_;

    $global{"gridOptions${tag}"}   = undef;
    $synops{"gridOptions${tag}"}   = "Grid engine options applied to $name jobs";

    $global{"${tag}Memory"}        = undef;
    $synops{"${tag}Memory"}        = "Amount of memory, in gigabytes, to use for $name jobs";

    $global{"${tag}Threads"}       = undef;
    $synops{"${tag}Threads"}       = "Number of threads to use for $name jobs";

    $global{"${tag}Concurrency"}   = undef;
    $synops{"${tag}Concurrency"}   = "If grid not enabled, number of $name jobs to run at the same time; default is n_proc / n_threads";
}



sub showErrorRates ($) {
    my $prefix = shift @_;

    print STDERR "${prefix}\n";
    print STDERR "${prefix}genomeSize          -- ", getGlobal("genomeSize"), "\n";
    print STDERR "${prefix}errorRate           -- ", getGlobal("errorRate"), "\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}corOvlErrorRate     -- ", getGlobal("corOvlErrorRate"), "\n";
    print STDERR "${prefix}obtOvlErrorRate     -- ", getGlobal("obtOvlErrorRate"), "\n";
    print STDERR "${prefix}utgOvlErrorRate     -- ", getGlobal("utgOvlErrorRate"), "\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}obtErrorRate        -- ", getGlobal("obtErrorRate"), "\n";
    print STDERR "${prefix}\n";
    print STDERR "${prefix}utgGraphErrorRate   -- ", getGlobal("utgGraphErrorRate"), "\n";
    print STDERR "${prefix}utgBubbleErrorRate  -- ", getGlobal("utgBubbleErrorRate"), "\n";
    print STDERR "${prefix}utgMergeErrorRate   -- ", getGlobal("utgMergeErrorRate"), "\n";
    print STDERR "${prefix}utgRepeatErrorRate  -- ", getGlobal("utgRepeatErrorRate"), "\n";
    print STDERR "${prefix}\n";
    #print STDERR "${prefix}corErrorRate        -- ", getGlobal("corErrorRate"), "\n";
    print STDERR "${prefix}cnsErrorRate        -- ", getGlobal("cnsErrorRate"), "\n";
}



#  Defaults are set for yeast:
#    trimming   errorRate = 0.009  obtOvlErrorRate = 0.06  obtErrorRate = 0.035
#    assembly   errorRate = 0.009  utgOvlErrorRate = 0.06  bogart 0.035
#
sub setErrorRate ($@) {
    my $er      = shift @_;
    my $verbose = shift @_;

    print STDERR "-- Set errorRate to $er (verbose='$verbose')\n"  if (defined($verbose));

    #  Can NOT call setGlobal() for this, because it calls setErrorRate()!.
    $global{"errorrate"} = $er;
    setGlobal("corOvlErrorRate",    $er * 3);  #  Not used, except for realigning
    setGlobal("obtOvlErrorRate",    $er * 3);  #  Generally must be smaller than utgGraphErrorRate
    setGlobal("utgOvlErrorRate",    $er * 3);

    setGlobal("obtErrorRate",       $er * 3);

    setGlobal("utgGraphErrorRate",  $er * 2);
    setGlobal("utgBubbleErrorRate", $er * 2 + 0.5 * $er);  #  Not tested!
    setGlobal("utgMergeErrorRate",  $er * 2 - 0.5 * $er);
    setGlobal("utgRepeatErrorRate", $er * 2);

    #  Removed, is usually set in CorrectReads, can be set from command line directly.
    #setGlobal("corErrorRate",       $er * 10);  #  Erorr rate used for raw sequence alignment/consensus
    setGlobal("cnsErrorRate",       $er * 2.5);

    showErrorRates("--  ")  if (defined($verbose));
}



sub setOverlapDefaults ($$$) {
    my $tag     = shift @_;  #  If 'cor', some parameters are loosened for raw pacbio reads
    my $name    = shift @_;
    my $default = shift @_;  #  Sets ${tag}Overlapper

    #  Which overlapper to use.

    $global{"${tag}Overlapper"}               = $default;
    $synops{"${tag}Overlapper"}               = "Which overlap algorithm to use for $name";

    #  OverlapInCore parameters.

    $global{"${tag}OvlHashBlockLength"}       = undef;
    $synops{"${tag}OvlHashBlockLength"}       = "Amount of sequence (bp) to load into the overlap hash table";

    $global{"${tag}OvlRefBlockSize"}          = undef;
    $synops{"${tag}OvlRefBlockSize"}          = "Number of reads to search against the hash table per batch";

    $global{"${tag}OvlRefBlockLength"}        = 0;
    $synops{"${tag}OvlRefBlockLength"}        = "Amount of sequence (bp) to search against the hash table per batch";

    $global{"${tag}OvlHashBits"}              = ($tag eq "cor") ? 18 : 22;
    $synops{"${tag}OvlHashBits"}              = "Width of the kmer hash.  Width 22=1gb, 23=2gb, 24=4gb, 25=8gb.  Plus 10b per ${tag}OvlHashBlockLength";

    $global{"${tag}OvlHashLoad"}              = 0.75;
    $synops{"${tag}OvlHashLoad"}              = "Maximum hash table load.  If set too high, table lookups are inefficent; if too low, search overhead dominates run time; default 0.75";

    $global{"${tag}OvlMerSize"}               = ($tag eq "cor") ? 19 : 22;
    $synops{"${tag}OvlMerSize"}               = "K-mer size for seeds in overlaps";

    $global{"${tag}OvlMerThreshold"}          = "auto";
    $synops{"${tag}OvlMerThreshold"}          = "K-mer frequency threshold; mers more frequent than this count are ignored; default 'auto'";

    $global{"${tag}OvlMerDistinct"}           = undef;
    $synops{"${tag}OvlMerDistinct"}           = "K-mer frequency threshold; the least frequent fraction of distinct mers can seed overlaps";

    $global{"${tag}OvlMerTotal"}              = undef;
    $synops{"${tag}OvlMerTotal"}              = "K-mer frequency threshold; the least frequent fraction of all mers can seed overlaps";

    $global{"${tag}OvlFrequentMers"}          = undef;
    $synops{"${tag}OvlFrequentMers"}          = "Do not seed overlaps with these kmers (fasta format)";

    #  mecat2asmpw parameters.

    $global{"${tag}mecat2asmpwVersion"}              = "1.5b1";
    $synops{"${tag}mecat2asmpwVersion"}              = "Version of the mecat2asmpw jar file to use";

    $global{"${tag}mecat2asmpwFilterThreshold"}      = "0.000005";
    $synops{"${tag}mecat2asmpwFilterThreshold"}      = "Value between 0 and 1. kmers which comprise more than this percentage of the input are downweighted";

    $global{"${tag}mecat2asmpwBlockSize"}            = 1500;
    $synops{"${tag}mecat2asmpwBlockSize"}            = "Number of reads per 1GB; memory * blockSize = the size of  block loaded into memory per job";

    $global{"${tag}mecat2asmpwMerSize"}              = ($tag eq "cor") ? 16 : 22;
    $synops{"${tag}mecat2asmpwMerSize"}              = "K-mer size for seeds in mecat2asmpw";

    $global{"${tag}mecat2asmpwReAlign"}              = undef;
    $synops{"${tag}mecat2asmpwReAlign"}              = "Compute actual alignments from mecat2asmpw overlaps; 'raw' from mecat2asmpw output, 'final' from overlap store; uses either obtErrorRate or ovlErrorRate, depending on which overlaps are computed";

    $global{"${tag}mecat2asmpwSensitivity"}          = "normal";
    $synops{"${tag}mecat2asmpwSensitivity"}          = "Coarse sensitivity level: 'normal' or 'high'; default 'normal'";
}



sub setDefaults () {

    #####  General Configuration Options (aka miscellany)

    $global{"canuIteration"}               = 1;  #  See documentation in Execution.pm
    $global{"canuIterationMax"}            = 2;

    $global{"showNext"}                    = undef;
    $synops{"showNext"}                    = "Don't run any commands, just report what would run";

    $global{"pathMap"}                     = undef;
    $synops{"pathMap"}                     = "File with a hostname to binary directory map";

    $global{"shell"}                       = "/bin/sh";
    $synops{"shell"}                       = "Command interpreter to use; sh-compatible (e.g., bash), NOT C-shell (csh or tcsh); default '/bin/sh'";

    #$global{"java"}                        = "java";
    #$synops{"java"}                        = "Java interpreter to use; at least version 1.8; default 'java'";

    #####  Cleanup options

    $global{"saveOverlaps"}                = 0;
    $synops{"saveOverlaps"}                = "Save intermediate overlap files, almost never a good idea";

    $global{"saveReadCorrections"}         = 0;
    $synops{"saveReadCorrections"}         = "Save intermediate read correction files, almost never a good idea";

    $global{"saveMerCounts"}               = 0;
    $synops{"saveMerCounts"}               = "Save full mer counting results, sometimes useful";

    #####  Error Rates

    $global{"errorRate"}                   = undef;
    $synops{"errorRate"}                   = "The expected error rate in the input reads";

    $global{"corOvlErrorRate"}             = undef;
    $synops{"corOvlErrorRate"}             = "Overlaps above this error rate are not computed";

    $global{"obtOvlErrorRate"}             = undef;
    $synops{"obtOvlErrorRate"}             = "Overlaps at or below this error rate are used to trim reads";

    $global{"utgOvlErrorRate"}             = undef;
    $synops{"utgOvlErrorRate"}             = "Overlaps at or below this error rate are used to trim reads";

    #$global{"utgErrorRate"}                = undef;
    #$synops{"utgErrorRate"}                = "Overlaps at or below this error rate are used to construct unitigs (BOG and UTG)";

    $global{"utgGraphErrorRate"}           = undef;
    $synops{"utgGraphErrorRate"}           = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgBubbleErrorRate"}          = undef;
    $synops{"utgBubbleErrorRate"}          = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgMergeErrorRate"}           = undef;
    $synops{"utgMergeErrorRate"}           = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"utgRepeatErrorRate"}          = undef;
    $synops{"utgRepeatErrorRate"}          = "Overlaps at or below this error rate are used to construct unitigs (BOGART)";

    $global{"corErrorRate"}                = undef;
    $synops{"corErrorRate"}                = "Only use raw alignments below this error rate to construct corrected reads";

    $global{"cnsErrorRate"}                = undef;
    $synops{"cnsErrorRate"}                = "Consensus expects alignments at about this error rate";

    #####  Minimums and maximums

    $global{"minReadLength"}               = 500;
    $synops{"minReadLength"}               = "Reads shorter than this length are not loaded into the assembler; default 1000";

    $global{"minOverlapLength"}            = 500;
    $synops{"minOverlapLength"}            = "Overlaps shorter than this length are not computed; default 500";

    $global{"minMemory"}                   = undef;
    $synops{"minMemory"}                   = "Minimum amount of memory needed to compute the assembly (do not set unless prompted!)";
    $global{"minThreads"}                  = undef;
    $synops{"minThreads"}                  = "Minimum number of compute threads suggested to compute the assembly";

    $global{"maxMemory"}                   = undef;
    $synops{"maxMemory"}                   = "Maximum memory to use by any component of the assembler";
    $global{"maxThreads"}                  = undef;
    $synops{"maxThreads"}                  = "Maximum number of compute threads to use by any component of the assembler";

    #####  Stopping conditions

    $global{"stopBefore"}                  = undef;
    $synops{"stopBefore"}                  = "Tell canu when to halt execution";

    $global{"stopAfter"}                   = undef;
    $synops{"stopAfter"}                   = "Tell canu when to halt execution";

    #####  Grid Engine configuration, internal parameters.  These are filled out in canu.pl, right after this function returns.

    $global{"availableHosts"}                       = undef;  #  Internal list of cpus-memory-nodes describing the grid

    $global{"gridEngine"}                           = undef;
    $global{"gridEngineSubmitCommand"}              = undef;
    $global{"gridEngineHoldOption"}                 = undef;
    $global{"gridEngineHoldOptionNoArray"}          = undef;
    $global{"gridEngineSyncOption"}                 = undef;
    $global{"gridEngineNameOption"}                 = undef;
    $global{"gridEngineArrayOption"}                = undef;
    $global{"gridEngineArrayName"}                  = undef;
    $global{"gridEngineOutputOption"}               = undef;
    $global{"gridEnginePropagateCommand"}           = undef;
    $global{"gridEngineThreadsOption"}              = undef;
    $global{"gridEngineMemoryOption"}               = undef;
    $global{"gridEngineMemoryUnits"}                = undef;
    $global{"gridEngineNameToJobIDCommand"}         = undef;
    $global{"gridEngineNameToJobIDCommandNoArray"}  = undef;
    $global{"gridEngineTaskID"}                     = undef;
    $global{"gridEngineArraySubmitID"}              = undef;
    $global{"gridEngineJobID"}                      = undef;

    #####  Grid Engine Pipeline

    $global{"useGrid"}                     = 1;
    $synops{"useGrid"}                     = "If 'true', enable grid-based execution; if 'false', run all jobs on the local machine; if 'remote', create jobs for grid execution but do not submit; default 'true'";

    #####  Grid Engine configuration, for each step of the pipeline

    $global{"gridOptions"}                 = undef;
    $synops{"gridOptions"}                 = "Grid engine options applied to all jobs";

    $global{"gridOptionsExecutive"}        = undef;
    $synops{"gridOptionsExecutive"}        = "Grid engine options applied to the canu executive script";

    $global{"gridOptionsJobName"}          = undef;
    $synops{"gridOptionsJobName"}          = "Grid jobs job-name suffix";

    #####  Grid Engine configuration and parameters, for each step of the pipeline (memory, threads)

    setExecDefaults("meryl",  "mer counting");

    setExecDefaults("cor",    "read correction");

    setExecDefaults("corovl",  "overlaps for correction");
    setExecDefaults("obtovl",  "overlaps for trimming");
    setExecDefaults("utgovl",  "overlaps for unitig construction");

    setExecDefaults("cormecat2asmpw", "mecat2asmpw overlaps for correction");
    setExecDefaults("obtmecat2asmpw", "mecat2asmpw overlaps for trimming");
    setExecDefaults("utgmecat2asmpw", "mecat2asmpw overlaps for unitig construction");

    setExecDefaults("ovb",    "overlap store bucketizing");
    setExecDefaults("ovs",    "overlap store sorting");

    setExecDefaults("red",    "read error detection");
    setExecDefaults("oea",    "overlap error adjustment");

    setExecDefaults("bat",    "unitig construction");
    setExecDefaults("cns",    "unitig consensus");

    #####  Overlapper

    setOverlapDefaults("cor", "correction",             "mecat2asmpw");  #  Overlaps computed for correction
    setOverlapDefaults("obt", "overlap based trimming", "ovl");   #  Overlaps computed for trimming
    setOverlapDefaults("utg", "unitig construction",    "ovl");   #  Overlaps computed for unitigging

    ##### Overlap Store

    $global{"ovlStoreMemory"}              = undef;
    $synops{"ovlStoreMemory"}              = "How much memory, in gigabytes, to use when constructing overlap stores";

    $global{"ovlStoreThreads"}             = 1;
    $synops{"ovlStoreThreads"}             = "Unused, only one thread supported";

    $global{"ovlStoreConcurrency"}         = undef;
    $synops{"ovlStoreConcurrency"}         = "Unused, only one process supported";

    $global{"ovlStoreMethod"}              = "sequential";
    $synops{"ovlStoreMethod"}              = "Use the 'sequential' or 'parallel' algorithm for constructing an overlap store; default 'sequential'";

    $global{"ovlStoreSlices"}              = undef;
    $synops{"ovlStoreSlices"}              = "How many pieces to split the sorting into, for the parallel store build";

    #####  Mers

    $global{"merylMemory"}                 = undef;
    $synops{"merylMemory"}                 = "Amount of memory, in gigabytes, to use for mer counting";

    $global{"merylThreads"}                = undef;
    $synops{"merylThreads"}                = "Number of threads to use for mer counting";

    $global{"merylConcurrency"}            = undef;
    $synops{"merylConcurrency"}            = "Unused, there is only one process";

    #####  Overlap Based Trimming

    $global{"obtErrorRate"}                = undef;
    $synops{"obtErrorRate"}                = "Stringency of overlaps to use for trimming";

    $global{"trimReadsOverlap"}            = 1;
    $synops{"trimReadsOverlap"}            = "Minimum overlap between evidence to make contiguous trim; default '1'";

    $global{"trimReadsCoverage"}           = 1;
    $synops{"trimReadsCoverage"}           = "Minimum depth of evidence to retain bases; default '1'";

    #$global{"splitReads..."}               = 1;
    #$synops{"splitReads..."}               = "";

    #####  Fragment/Overlap Error Correction

    $global{"enableOEA"}                   = 1;
    $synops{"enableOEA"}                   = "Do overlap error adjustment - comprises two steps: read error detection (RED) and overlap error adjustment (OEA); default 'true'";

    $global{"redBatchSize"}                = undef;
    $synops{"redBatchSize"}                = "Number of reads per fragment error detection batch";

    $global{"redBatchLength"}              = undef;
    $synops{"redBatchLength"}              = "Number of bases per fragment error detection batch";

    $global{"oeaBatchSize"}                = undef;
    $synops{"oeaBatchSize"}                = "Number of reads per overlap error correction batch";

    $global{"oeaBatchLength"}              = undef;
    $synops{"oeaBatchLength"}              = "Number of bases per overlap error correction batch";

    #####  Unitigger & BOG & bogart Options

    $global{"unitigger"}                   = "bogart";
    $synops{"unitigger"}                   = "Which unitig algorithm to use; only 'bogart' supported; default 'bogart'";

    $global{"genomeSize"}                  = undef;
    $synops{"genomeSize"}                  = "An estimate of the size of the genome";

    $global{"batOptions"}                  = undef;
    $synops{"batOptions"}                  = "Advanced options to bogart";

    $global{"batMemory"}                   = undef;
    $synops{"batMemory"}                   = "Approximate maximum memory usage, in gigabytes, default is the maxMemory limit";

    $global{"batThreads"}                  = undef;
    $synops{"batThreads"}                  = "Number of threads to use; default is the maxThreads limit";

    $global{"batConcurrency"}              = undef;
    $synops{"batConcurrency"}              = "Unused, only one process supported";

    #####  Unitig Filtering Options (also set in bogart/bogart.C)

    $global{"contigFilter"}                = "2 1000 0.75 0.75 2";
    $synops{"contigFilter"}                = "Parameters to filter out 'unassembled' unitigs: minReads; minLength; singleReadSpan; lowCovSpan, lowCovDepth";

    #####  Consensus Options

    $global{"cnsPartitions"}               = undef;
    $synops{"cnsPartitions"}               = "Partition consensus into N jobs";

    $global{"cnsPartitionMin"}             = undef;
    $synops{"cnsPartitionMin"}             = "Don't make a consensus partition with fewer than N reads";

    $global{"cnsMaxCoverage"}              = 0;
    $synops{"cnsMaxCoverage"}              = "Limit unitig consensus to at most this coverage; default '0' = unlimited";

    $global{"cnsConsensus"}                = "pbdagcon";
    $synops{"cnsConsensus"}                = "Which consensus algorithm to use; 'pbdagcon' (fast, reliable); 'utgcns' (multialignment output); 'quick' (single read mosaic); default 'pbdagcon'";

    #####  Correction Options

    $global{"corPartitions"}               = undef;
    $synops{"corPartitions"}               = "Partition read correction into N jobs";

    $global{"corPartitionMin"}             = undef;
    $synops{"corPartitionMin"}             = "Don't make a read correction partition with fewer than N reads";

    $global{"corMinEvidenceLength"}        = undef;
    $synops{"corMinEvidenceLength"}        = "Limit read correction to only overlaps longer than this; default: unlimited";

    $global{"corMaxEvidenceErate"}         = undef;
    $synops{"corMaxEvidenceErate"}         = "Limit read correction to only overlaps at or below this fraction error; default: unlimited";

    $global{"corMaxEvidenceCoverageGlobal"}= "1.0x";
    $synops{"corMaxEvidenceCoverageGlobal"}= "Limit reads used for correction to supporting at most this coverage; default: '1.0x' = 1.0 * estimated coverage";

    $global{"corMaxEvidenceCoverageLocal"} = "2.0x";
    $synops{"corMaxEvidenceCoverageLocal"} = "Limit reads being corrected to at most this much evidence coverage; default: '2.0x' = 2.0 * estimated coverage";

    $global{"corOutCoverage"}              = 40;
    $synops{"corOutCoverage"}              = "Only correct the longest reads up to this coverage; default 40";

    $global{"corMinCoverage"}              = 4;
    $synops{"corMinCoverage"}              = "Minimum number of bases supporting each corrected base, if less than this sequences are split; default 4";

    $global{"corFilter"}                   = "expensive";
    $synops{"corFilter"}                   = "Method to filter short reads from correction; 'quick' or 'expensive'; default 'expensive'";

    $global{"corConsensus"}                = "falconpipe";
    $synops{"corConsensus"}                = "Which consensus algorithm to use; only 'falcon' and 'falconpipe' are supported; default 'falconpipe'";

    $global{"falconSense"}                 = undef;
    $synops{"falconSense"}                 = "Path to fc_consensus.py or falcon_sense.bin";

    #  Convert all the keys to lowercase, and remember the case-sensitive version

    foreach my $k (keys %global) {
        (my $l = $k) =~ tr/A-Z/a-z/;

        if (! exists($synnam{$l})) {
            $synnam{$l} = $k;

            if (!exists($global{$l})) {
                $global{$l} = $global{$k};
                delete $global{$k};
            }

            #print "$k -> $l\n";
        }
    }

    #  If this is set, it breaks the consensus.sh and overlap.sh scripts.  Good grief!  Why
    #  are you running this in a task array!?

    if (exists($ENV{getGlobal("gridEngineTaskID")})) {
        undef $ENV{getGlobal("gridEngineTaskID")};
        print STDERR "ENV: ", getGlobal("gridEngineTaskID"), " needs to be unset, done.\n";
    }

    #  Finally, set the global default error rate

    setErrorRate(0.01);
}



sub checkParameters () {

    #
    #  Fiddle with filenames to make them absolute paths.
    #

    makeAbsolute("pathMap");

    makeAbsolute("corOvlFrequentMers");
    makeAbsolute("obtOvlFrequentMers");
    makeAbsolute("utgOvlFrequentMers");

    #
    #  Adjust case on some of them
    #

    fixCase("corOverlapper");
    fixCase("obtOverlapper");
    fixCase("utgOverlapper");

    fixCase("corConsensus");
    fixCase("cnsConsensus");

    fixCase("corFilter");

    fixCase("unitigger");
    fixCase("stopBefore");
    fixCase("stopAfter");

    #
    #  Check for inconsistent parameters
    #

    if (getGlobal("minReadLength") < getGlobal("minOverlapLength")) {
        my $mr = getGlobal("minReadLength");
        my $mo = getGlobal("minOverlapLength");

        addCommandLineError("ERROR:  minReadLength=$mr must be at least minOverlapLength=$mo.\n");

        #  Or we can just reset one or the other....
        #print STDERR "-- WARNING: minReadLength reset from $mr to $mo (limited by minOverlapLength)\n";
        #setGlobal("minOverlapLength", $mo);
    }

    #
    #  Check for invalid usage
    #

    foreach my $tag ("cor", "obt", "utg") {
        if ((getGlobal("${tag}Overlapper") ne "mecat2asmpw") &&
            (getGlobal("${tag}Overlapper") ne "ovl")) {
            addCommandLineError("ERROR:  Invalid '${tag}Overlapper' specified (" . getGlobal("${tag}Overlapper") . "); must be 'mecat2asmpw' or 'ovl'\n");
        }
    }

    if ((getGlobal("unitigger") ne "unitigger") &&
        (getGlobal("unitigger") ne "bogart")) {
        addCommandLineError("ERROR:  Invalid 'unitigger' specified (" . getGlobal("unitigger") . "); must be 'unitigger' or 'bogart'\n");
    }

    if ((getGlobal("corConsensus") ne "utgcns") &&
        (getGlobal("corConsensus") ne "falcon") &&
        (getGlobal("corConsensus") ne "falconpipe")) {
        addCommandLineError("ERROR:  Invalid 'corConsensus' specified (" . getGlobal("corConsensus") . "); must be 'utgcns' or 'falcon' or 'falconpipe'\n");
    }

    if ((getGlobal("cnsConsensus") ne "quick") &&
        (getGlobal("cnsConsensus") ne "pbdagcon") &&
        (getGlobal("cnsConsensus") ne "utgcns")) {
        addCommandLineError("ERROR:  Invalid 'cnsConsensus' specified (" . getGlobal("cnsConsensus") . "); must be 'quick', 'pbdagcon', or 'utgcns'\n");
    }


    if ((!defined("lowCoverageAllowed") &&  defined("lowCoverageDepth")) ||
        ( defined("lowCoverageAllowed") && !defined("lowCoverageDepth"))) {
        addCommandLineError("ERROR:  Invalid 'lowCoverageAllowed' and 'lowCoverageDepth' specified; both must be set\n");
    }

    #if ((getGlobal("cleanup") ne "none") &&
    #    (getGlobal("cleanup") ne "light") &&
    #    (getGlobal("cleanup") ne "heavy") &&
    #    (getGlobal("cleanup") ne "aggressive")) {
    #    addCommandLineError("ERROR:  Invalid cleaup specified (" . getGlobal("cleanup") . "); must be 'none', 'light', 'heavy' or 'aggressive'\n");
    #}

    if ((getGlobal("corFilter") ne "quick") &&
        (getGlobal("corFilter") ne "expensive")) {
        addCommandLineError("ERROR:  Invalid 'corFilter' specified (" . getGlobal("corFilter") . "); must be 'quick' or 'expensive'\n");
    }


    if ((getGlobal("useGrid") ne "0") &&
        (getGlobal("useGrid") ne "1") &&
        (getGlobal("useGrid") ne "remote")) {
        addCommandLineError("ERROR:  Invalid 'useGrid' specified (" . getGlobal("useGrid") . "); must be 'true', 'false' or 'remote'\n");
    }

    if (defined(getGlobal("stopBefore"))) {
        my $ok = 0;
        my $st = getGlobal("stopBefore");
        $st =~ tr/A-Z/a-z/;

        my $failureString = "ERROR:  Invalid stopBefore specified (" . getGlobal("stopBefore") . "); must be one of:\n";

        my @stopBefore = ("gatekeeper",
                          "meryl",
                          "trimReads",
                          "splitReads",
                          "unitig",
                          "consensusConfigure");

        foreach my $sb (@stopBefore) {
            $failureString .= "ERROR:      '$sb'\n";
            $sb =~ tr/A-Z/a-z/;
            if ($st eq $sb) {
                $ok++;
                setGlobal('stopBefore', $st);
            }
        }

        addCommandLineError($failureString)   if ($ok == 0);
    }

    if (defined(getGlobal("stopAfter"))) {
        my $ok = 0;
        my $st = getGlobal("stopAfter");
        $st =~ tr/A-Z/a-z/;

        my $failureString = "ERROR:  Invalid stopAfter specified (" . getGlobal("stopAfter") . "); must be one of:\n";

        my @stopAfter = ("gatekeeper",
                         "meryl",
                         "mecat2asmpwConfigure",
                         "overlapStoreConfigure",
                         "overlapStore",
                         "unitig",
                         "consensusConfigure",
                         "consensusCheck",
                         "consensusLoad",
                         "consensusFilter");

        foreach my $sa (@stopAfter) {
            $failureString .= "ERROR:      '$sa'\n";
            $sa =~ tr/A-Z/a-z/;
            if ($st eq $sa) {
                $ok++;
                setGlobal('stopAfter', $st);
            }
        }

        addCommandLineError($failureString)   if ($ok == 0);
    }

    addCommandLineError("ERROR:  Required parameter 'errorRate' is not set\n")    if (! defined(getGlobal("errorRate")));
    addCommandLineError("ERROR:  Required parameter 'genomeSize' is not set\n")   if (! defined(getGlobal("genomeSize")));

    #
    #  Java?  Need JRE 1.8.
    #
=pod
    if ((getGlobal("corOverlapper") eq "mecat2asmpw") ||
        (getGlobal("obtOverlapper") eq "mecat2asmpw") ||
        (getGlobal("utgOverlapper") eq "mecat2asmpw")) {
        my $java       = getGlobal("java");
        my $versionStr = "unknown";
        my $version    = 0;

        open(F, "$java -showversion 2>&1 |");
        while (<F>) {
            #  First word is either "java" or "openjdk" or ...
            if (m/^.*\s+version\s+\"(\d+.\d+)(.*)\"$/) {
                $versionStr = "$1$2";
                $version    =  $1;
            }
        }
        close(F);

        print STDERR "-- Detected Java(TM) Runtime Environment '$versionStr' (from '$java').\n";

        addCommandLineError("ERROR:  mecat2asmpw overlapper requires java version at least 1.8.0; you have $versionStr\n")   if ($version < 1.8);
    }
=cut
    #
    #  Falcon?  Need to find it.
    #

    if ((getGlobal("corConsensus") eq "falcon") ||
        (getGlobal("corConsensus") eq "falconpipe")) {
        my $falcon = getGlobal("falconSense");

        addCommandLineError("ERROR:  Didn't find falcon program with option falconSense='$falcon'")   if ((defined($falcon)) && (! -e $falcon));
    }

    #
    #  Set default error rates based on the per-read error rate.
    #

    setGlobalIfUndef("corOvlErrorRate",      3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("obtOvlErrorRate",      3.0 * getGlobal("errorRate"));
    setGlobalIfUndef("utgOvlErrorRate",      3.0 * getGlobal("errorRate"));

    setGlobalIfUndef("ovlErrorRate",         2.5 * getGlobal("errorRate"));
    setGlobalIfUndef("utgGraphErrorRate",    2.0 * getGlobal("errorRate"));
    setGlobalIfUndef("utgBubbleErrorRate",   2.0 * getGlobal("errorRate") + 0.5 * getGlobal("errorRate"));
    setGlobalIfUndef("utgMergeErrorRate",    2.0 * getGlobal("errorRate") - 0.5 * getGlobal("errorRate"));
    setGlobalIfUndef("utgRepeatErrorRate",   2.0 * getGlobal("errorRate"));

    setGlobalIfUndef("corsErrorRate",        10.0 * getGlobal("errorRate"));
    setGlobalIfUndef("cnsErrorRate",         2.5 * getGlobal("errorRate"));
}


1;
