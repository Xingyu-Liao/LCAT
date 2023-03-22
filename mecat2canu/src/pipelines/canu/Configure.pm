
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

package canu::Configure;

require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(displayMemoryValue displayGenomeSize configureAssembler);

use strict;
use Carp qw(cluck);
use Sys::Hostname;

use canu::Defaults;


#  This is called to expand parameter ranges for memory and thread parameters.
#  Examples of valid ranges:
#
#  no units  - 1-4:2                      - assumes 'g' in the adjust (if memory)
#  one unit  - 1g-4:2 1-4g:2 1-4:2g       - all the others are set to 'g'
#  two units - 1g-4g:2 1g-4:2g 1-4g:2g    - bgn/end are the same, stp uses end
#  all three - 1g-4g:2g                   - use as is
#
#  Quirks:  1g-2000m will increment every 1m.
#           1g-2000m:1g only adds 1g.
#           1g-2048m:1g adds 1 and 2g.

sub expandRange ($$) {
    my $var = shift @_;
    my $val = shift @_;

    my @v = split ',', $val;
    my @r;

    foreach my $v (@v) {
        my $bgn;  my $bgnu;
        my $end;  my $endu;
        my $stp;  my $stpu;

        #  Decode the range.

        if      ($v =~ m/^(\d+\.{0,1}\d*)([kKmMgGtT]{0,1})$/) {
            $bgn = $1;  $bgnu = $2;
            $end = $1;  $endu = $2;
            $stp =  1;  $stpu = $2;
        } elsif ($v =~ m/^(\d+\.{0,1}\d*)([kKmMgGtT]{0,1})-(\d+\.{0,1}\d*)([kKmMgGtT]{0,1})$/) {
            $bgn = $1;  $bgnu = $2;
            $end = $3;  $endu = $4;
            $stp =  1;  $stpu = $4;
        } elsif ($v =~ m/^(\d+\.{0,1}\d*)([kKmMgGtT]{0,1})-(\d+\.{0,1}\d*)([kKmMgGtT]{0,1}):(\d+\.{0,1}\d*)([kKmMgGtT]{0,1})$/) {
            $bgn = $1;  $bgnu = $2;
            $end = $3;  $endu = $4;
            $stp = $5;  $stpu = $6;
        } else {
            caExit("can't parse '$var' entry '$v'", undef);
        }

        #  Undef things that are null.  The code that follows this was written assuming undef.

        $bgnu = undef   if ($bgnu eq "");
        $endu = undef   if ($endu eq "");
        $stpu = undef   if ($stpu eq "");

        #  Process the range

        my $def = defined($bgnu) + defined($endu) + defined($stpu);

        #  If no units, this could be a memory or a thread setting.  Don't use units.
        if      ($def == 0) {
        }

        #  If only one unit specified, set the others to the same.
        elsif ($def == 1) {
            if    (defined($bgnu))  { $endu = $stpu = $bgnu;  }
            elsif (defined($endu))  { $bgnu = $stpu = $endu;  }
            elsif (defined($stpu))  { $bgnu = $endu = $stpu;  }
        }

        #  If two units specified, set the unset as:
        #    bgn or end unset - set based on the other range
        #    stp unset        - set on end if stp<end otherwise bgn
        elsif ($def == 2) {

            if ((!defined($bgnu) && ($endu ne $stpu)) ||
                (!defined($endu) && ($bgnu ne $stpu)) ||
                (!defined($stpu) && ($bgnu ne $endu))) {
                print STDERR "--\n";
                print STDERR "-- WARNING: incomplete and inconsistent units on '$var=$val'.\n";
            }
                
            $bgnu = $endu  if (!defined($bgnu));
            $endu = $bgnu  if (!defined($endu));
            $stpu = $endu  if (!defined($stpu) && ($stp <= $end));
            $stpu = $bgnu  if (!defined($stpu) && ($stp >  $end));
        }

        #  Nothing to do if all three are set!
        elsif ($def == 3) {
        }

        my $b = adjustMemoryValue("$bgn$bgnu");
        my $e = adjustMemoryValue("$end$endu");
        my $s = adjustMemoryValue("$stp$stpu");

        for (my $ii=$b; $ii<=$e; $ii += $s) {
            push @r, $ii;
        }
    }

    #print "$var = ";
    #foreach my $r (@r) {
    #    print "$r ";
    #}
    #print "\n";

    return(@r);
}


#  Side effect!  This will RESET the $global{} parameters to the computed value.  This lets
#  the rest of canu - in particular, the part that runs the jobs - use the correct value.  Without
#  resetting, I'd be making code changes all over the place to support the values returned.

sub getAllowedResources ($$$$) {
    my $tag  = shift @_;  #  Variant, e.g., "cor", "utg"
    my $alg  = shift @_;  #  Algorithm, e.g., "mecat2asmpw", "ovl"
    my $err  = shift @_;  #  Report of things we can't run.
    my $all  = shift @_;  #  Report of things we can run.

    #  If no grid, or grid not enabled, everything falls under 'lcoal'.

    my $class = ((getGlobal("useGrid") ne "0") && (defined(getGlobal("gridEngine")))) ? "grid" : "local";

    #  Figure out limits.

    my $maxMemory    = getGlobal("maxMemory");
    my $maxThreads   = getGlobal("maxThreads");

    my $taskMemory   = getGlobal("${tag}${alg}Memory");   #  Algorithm limit, "utgovlMemory", etc.
    my $taskThreads  = getGlobal("${tag}${alg}Threads");  #

    #  The task limits MUST be defined.

    caExit("${tag}${alg}Memory is not defined", undef)   if (!defined($taskMemory));
    caExit("${tag}${alg}Threads is not defined", undef)  if (!defined($taskThreads));

    #  If the maximum limits aren't set, default to 'unlimited' (for the grid; we'll effectively filter
    #  by the number of jobs we can fit on the hosts) or to the current hardware limits.

    $maxMemory  = (($class eq "grid") ? 1024 * 1024 : getPhysicalMemorySize())  if (!defined($maxMemory));    #  1 PB memory!
    $maxThreads = (($class eq "grid") ? 1024        : getNumberOfCPUs())        if (!defined($maxThreads));   #  1 k  cores!

    #  Build a list of the available hardware configurations we can run on.  If grid, we get this
    #  from the list of previously discovered hosts.  If local, it's just this machine.

    my @gridCor;  #  Number of cores
    my @gridMem;  #  GB's of memory
    my @gridNum;  #  Number of nodes

    if ($class eq "grid") {
        my @grid = split '\0', getGlobal("availableHosts");

        if (scalar(@grid) == 0) {
            caExit("invalid useGrid (" . getGlobal("useGrid") . ") and gridEngine (" . getGlobal("gridEngine") . "); found no execution hosts - is grid available from this host?", undef);
        }

        foreach my $g (@grid) {
            my ($cpu, $mem, $num) = split '-', $g;

            push @gridCor, $cpu;
            push @gridMem, $mem;
            push @gridNum, $num;
        }
    } else {
        push @gridCor, $maxThreads;
        push @gridMem, $maxMemory;
        push @gridNum, 1;
    }

    #  The task usually has multiple choices, and we have a little optimization problem to solve.  For each
    #  pair of memory/threads, compute three things:
    #    a) how many processes we can get running
    #    b) how many cores we can get running
    #    c) how much memory we can consume
    #  We then (typically) want to maximize the number of cores we can get running.
    #  Other options would be number of cores * amount of memory.

    my @taskMemory  = expandRange("${tag}${alg}Memory",  $taskMemory);
    my @taskThreads = expandRange("${tag}${alg}Threads", $taskThreads);

    #  Filter out task settings that can't be run based on the gridMemory/gridThreads or masterMemory/masterThreads setting.
    #  (actually, this just reports those that would be filtered; the actual filtering is inline in the algorithm)

    my $ignoreM;
    my $ignoreT;

    foreach my $m (@taskMemory) {
        $m = adjustMemoryValue($m);
    }

    foreach my $m (@taskMemory) {
        next  if ($m <= $maxMemory);
        $ignoreM .= ","  if (defined($ignoreM));
        $ignoreM .= "${m}g";
    }
    foreach my $t (@taskThreads) {
        next  if ($t <= $maxThreads);
        $ignoreT .= ","  if (defined($ignoreT));
        $ignoreT .= "$t";
    }

    #  Too verbose with long value lists
    #
    #if      (defined($ignoreM) && defined($ignoreT)) {
    #    $err .= "-- Can't use ${tag}${alg}Memory=$ignoreM and ${tag}${alg}Threads=$ignoreT because of maxMemory=${maxMemory}g and maxThreads=$maxThreads limits.\n";
    #
    #} elsif (defined($ignoreM)) {
    #    $err .= "-- Can't use ${tag}${alg}Memory=$ignoreM because of maxMemory=${maxMemory}g limit.\n";
    #
    #} elsif (defined($ignoreT)) {
    #    $err .= "-- Can't use ${tag}${alg}Threads=$ignoreT because of maxThreads=$maxThreads limit.\n";
    #}

    #  Find task memory/thread settings that will maximize the number of cores running.  This used
    #  to also compute best as 'cores * memory' but that is better handled by ordering the task
    #  settings parameters.  The example below will pick the largest (last) configuration that
    #  maximizes core utilization:
    #
    #    taskThreads = 4,8,32,64
    #    taskMemory  = 16g,32g,64g

    my ($bestCores,  $bestCoresM,  $bestCoresT)  = (0, undef, undef);

    foreach my $m (@taskMemory) {
        foreach my $t (@taskThreads) {

            next  if ($m > $maxMemory);   #  Bail if either of the suggest settings are
            next  if ($t > $maxThreads);  #  larger than the maximum allowed.

            my $processes = 0;
            my $cores     = 0;
            my $memory    = 0;

            #  For a job using $m GB memory and $t threads, we can compute how many processes will
            #  fit on each node in our set of available machines.  The smaller of the two is then
            #  the number of processes we can run on this node.

            for (my $ii=0; $ii<scalar(@gridCor); $ii++) {
                my $np_cpu = $gridNum[$ii] * int($gridCor[$ii] / $t);  #  Each process uses $t cores, node has $gridCor[$ii] cores available.
                my $np_mem = $gridNum[$ii] * int($gridMem[$ii] / $m);  #  Same idea.

                my $np = ($np_cpu < $np_mem) ? $np_cpu : $np_mem;

                $processes += $np;
                $cores     += $np * $t;
                $memory    += $np * $m;
            }

            #  Save the best one seen so far.

            if ($bestCores <= $cores) {
                $bestCores  = $cores;
                $bestCoresM = $m;
                $bestCoresT = $t;
            }
        }
    }

    if (!defined($bestCoresM)) {
        print STDERR "--\n";
        print STDERR "-- Task $tag$alg can't run on any available machines.\n";
        print STDERR "-- It is requesting ", getGlobal("${tag}${alg}Memory"), " GB memory and ", getGlobal("${tag}${alg}Threads"), " threads.\n";
        print STDERR "-- See above for hardware limits.\n";
        print STDERR "--\n";

        caExit("task $tag$alg failed to find a configuration to run on", undef);
    }

    $taskMemory  = $bestCoresM;
    $taskThreads = $bestCoresT;

    #  Check for stupidity.

    caExit("invalid taskThread=$taskMemory; maxMemory=$maxMemory", undef)     if ($taskMemory > $maxMemory);
    caExit("invalid taskThread=$taskThreads; maxThreads=$maxThreads", undef)  if ($taskThreads > $maxThreads);

    #  Reset the global values for later use.

    setGlobal("${tag}${alg}Memory",  $taskMemory);
    setGlobal("${tag}${alg}Threads", $taskThreads);

    #  Finally, reset the concurrency (if we're running locally) so we don't swamp our poor workstation.

    my $concurrent = undef;

    if ($class eq "local") {
        my $nc = int($maxThreads / $taskThreads);

        if (($taskThreads * getGlobal("${tag}${alg}Concurrency") > $maxThreads)) {
            $err .= "-- Reset concurrency from ", getGlobal("${tag}${alg}Concurrency"), " to $nc.\n";
            setGlobal("${tag}${alg}Concurrency", $nc);
        }

        if (!defined(getGlobal("${tag}${alg}Concurrency"))) {
            setGlobal("${tag}${alg}Concurrency", $nc);
        }

        $concurrent = getGlobal("${tag}${alg}Concurrency");
    }

    #  And report.

    my $nam;

    if    ($alg eq "bat")      {  $nam = "bogart (unitigger)"; }
    elsif ($alg eq "cns")      {  $nam = "utgcns (consensus"; }
    elsif ($alg eq "cor")      {  $nam = "falcon_sense (read correction)"; }
    elsif ($alg eq "meryl")    {  $nam = "meryl (k-mer counting)"; }
    elsif ($alg eq "oea")      {  $nam = "overlap error adjustment"; }
    elsif ($alg eq "ovb")      {  $nam = "overlap store parallel bucketizer"; }
    elsif ($alg eq "ovlStore") {  $nam = "overlap store sequential building"; }
    elsif ($alg eq "ovs")      {  $nam = "overlap store parallel sorting"; }
    elsif ($alg eq "red")      {  $nam = "read error detection (overlap error adjustment)"; }
    elsif ($alg eq "mecat2asmpw")     {  $nam = "mecat2asmpw (overlapper)"; }
    elsif ($alg eq "ovl")      {  $nam = "overlapper"; }
    else {
        caFailure("unknown task '$alg' in getAllowedResources().", undef);
    }

    $all .= "-- Allowed to";
    $all .= " run " . substr("   $concurrent", -3) . " job" . (($concurrent == 1) ? " " : "s") . " concurrently,"  if (defined($concurrent));
    $all .= " run under grid control,"                                                                             if (!defined($concurrent));
    $all .= " and use up to " . substr("   $taskThreads", -3) . " compute thread" . (($taskThreads == 1) ? " " : "s");
    $all .= " and " . substr("   $taskMemory", -4) . " GB memory for stage '$nam'.\n";

    return($err, $all);
}




#  Converts number with units to gigabytes.  If no units, gigabytes is assumed.
sub adjustMemoryValue ($) {
    my $val = shift @_;

    return(undef)                     if (!defined($val));

    return($1)                        if ($val =~ m/^(\d+\.{0,1}\d*)$/);
    return($1 / 1024 / 1024)          if ($val =~ m/^(\d+\.{0,1}\d*)[kK]$/);
    return($1 / 1024)                 if ($val =~ m/^(\d+\.{0,1}\d*)[mM]$/);
    return($1)                        if ($val =~ m/^(\d+\.{0,1}\d*)[gG]$/);
    return($1 * 1024)                 if ($val =~ m/^(\d+\.{0,1}\d*)[tT]$/);
    return($1 * 1024 * 1024)          if ($val =~ m/^(\d+\.{0,1}\d*)[pP]$/);

    die "Invalid memory value '$val'\n";
}


#  Converts gigabytes to number with units.
sub displayMemoryValue ($) {
    my $val = shift @_;

    return(($val * 1024 * 1024)        . "k")   if ($val < adjustMemoryValue("1m"));
    return(($val * 1024)               . "m")   if ($val < adjustMemoryValue("1g"));
    return(($val)                      . "g")   if ($val < adjustMemoryValue("1t"));
    return(($val / 1024)               . "t");
}


#  Converts number with units to bases.
sub adjustGenomeSize ($) {
    my $val = shift @_;

    return(undef)               if (!defined($val));

    return($1)                  if ($val =~ m/^(\d+\.{0,1}\d*)$/i);
    return($1 * 1000)           if ($val =~ m/^(\d+\.{0,1}\d*)[kK]$/i);
    return($1 * 1000000)        if ($val =~ m/^(\d+\.{0,1}\d*)[mM]$/i);
    return($1 * 1000000000)     if ($val =~ m/^(\d+\.{0,1}\d*)[gG]$/i);
    return($1 * 1000000000000)  if ($val =~ m/^(\d+\.{0,1}\d*)[tT]$/i);

    die "Invalid genome size '$val'\n";
}


#  Converts bases to number with units.
sub displayGenomeSize ($) {
    my $val = shift @_;

    return(($val))                        if ($val < adjustGenomeSize("1k"));
    return(($val / 1000)          . "k")  if ($val < adjustGenomeSize("1m"));
    return(($val / 1000000)       . "m")  if ($val < adjustGenomeSize("1g"));
    return(($val / 1000000000)    . "g")  if ($val < adjustGenomeSize("1t"));
    return(($val / 1000000000000) . "t");
}





#
#  If minMemory or minThreads isn't defined, pick a reasonable pair based on genome size.
#

sub configureAssembler () {

    #  Parse units on things the user possibly set.

    setGlobal("genomeSize", adjustGenomeSize(getGlobal("genomeSize")));

    setGlobal("minMemory",  adjustMemoryValue(getGlobal("minMemory")));
    setGlobal("maxMemory",  adjustMemoryValue(getGlobal("maxMemory")));


    #  For overlapper and mecat2asmpw, allow larger maximums for larger genomes.  More memory won't help
    #  smaller genomes, and the smaller minimums won't hurt larger genomes (which are probably being
    #  run on larger machines anyway, so the minimums won't be used).

    #  For uncorrected overlapper, both memory and thread count is reduced.  Memory because it is
    #  very CPU bound, and thread count because it can be quite unbalanced.

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("corOvlMemory", "2-6");     setGlobalIfUndef("corOvlThreads", "1");
        setGlobalIfUndef("obtOvlMemory", "4-8");     setGlobalIfUndef("obtOvlThreads", "1-8");
        setGlobalIfUndef("utgOvlMemory", "4-8");     setGlobalIfUndef("utgOvlThreads", "1-8");

        setGlobalIfUndef("cormecat2asmpwMemory", "4-6");   setGlobalIfUndef("cormecat2asmpwThreads", "1-16");
        setGlobalIfUndef("obtmecat2asmpwMemory", "4-6");   setGlobalIfUndef("obtmecat2asmpwThreads", "1-16");
        setGlobalIfUndef("utgmecat2asmpwMemory", "4-6");   setGlobalIfUndef("utgmecat2asmpwThreads", "1-16");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("500m")) {
        setGlobalIfUndef("corOvlMemory", "2-6");     setGlobalIfUndef("corOvlThreads", "1");
        setGlobalIfUndef("obtOvlMemory", "4-8");     setGlobalIfUndef("obtOvlThreads", "1-8");
        setGlobalIfUndef("utgOvlMemory", "4-8");     setGlobalIfUndef("utgOvlThreads", "1-8");

        setGlobalIfUndef("cormecat2asmpwMemory", "8-13");   setGlobalIfUndef("cormecat2asmpwThreads", "1-16");
        setGlobalIfUndef("obtmecat2asmpwMemory", "8-13");   setGlobalIfUndef("obtmecat2asmpwThreads", "1-16");
        setGlobalIfUndef("utgmecat2asmpwMemory", "8-13");   setGlobalIfUndef("utgmecat2asmpwThreads", "1-16");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("2g")) {
        setGlobalIfUndef("corOvlMemory", "2-8");     setGlobalIfUndef("corOvlThreads", "1");
        setGlobalIfUndef("obtOvlMemory", "4-12");    setGlobalIfUndef("obtOvlThreads", "1-8");
        setGlobalIfUndef("utgOvlMemory", "4-12");    setGlobalIfUndef("utgOvlThreads", "1-8");

        setGlobalIfUndef("cormecat2asmpwMemory", "16-32");   setGlobalIfUndef("cormecat2asmpwThreads", "4-16");
        setGlobalIfUndef("obtmecat2asmpwMemory", "16-32");   setGlobalIfUndef("obtmecat2asmpwThreads", "4-16");
        setGlobalIfUndef("utgmecat2asmpwMemory", "16-32");   setGlobalIfUndef("utgmecat2asmpwThreads", "4-16");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("5g")) {
        setGlobalIfUndef("corOvlMemory", "2-8");     setGlobalIfUndef("corOvlThreads", "1");
        setGlobalIfUndef("obtOvlMemory", "4-16");    setGlobalIfUndef("obtOvlThreads", "1-8");
        setGlobalIfUndef("utgOvlMemory", "4-16");    setGlobalIfUndef("utgOvlThreads", "1-8");

        setGlobalIfUndef("cormecat2asmpwMemory", "16-48");  setGlobalIfUndef("cormecat2asmpwThreads", "4-16");
        setGlobalIfUndef("obtmecat2asmpwMemory", "16-48");  setGlobalIfUndef("obtmecat2asmpwThreads", "4-16");
        setGlobalIfUndef("utgmecat2asmpwMemory", "16-48");  setGlobalIfUndef("utgmecat2asmpwThreads", "4-16");

    } else {
        setGlobalIfUndef("corOvlMemory", "2-8");     setGlobalIfUndef("corOvlThreads", "1");
        setGlobalIfUndef("obtOvlMemory", "4-16");    setGlobalIfUndef("obtOvlThreads", "1-8");
        setGlobalIfUndef("utgOvlMemory", "4-16");    setGlobalIfUndef("utgOvlThreads", "1-8");

        setGlobalIfUndef("cormecat2asmpwMemory", "32-64");  setGlobalIfUndef("cormecat2asmpwThreads", "4-16");
        setGlobalIfUndef("obtmecat2asmpwMemory", "32-64");  setGlobalIfUndef("obtmecat2asmpwThreads", "4-16");
        setGlobalIfUndef("utgmecat2asmpwMemory", "32-64");  setGlobalIfUndef("utgmecat2asmpwThreads", "4-16");
    }

    #  Overlapper block sizes probably don't need to be modified based on genome size.

    setGlobalIfUndef("corOvlHashBlockLength",   2500000);   setGlobalIfUndef("corOvlRefBlockSize",   20000);   setGlobalIfUndef("corOvlRefBlockLength", 0);
    setGlobalIfUndef("obtOvlHashBlockLength", 100000000);   setGlobalIfUndef("obtOvlRefBlockSize", 2000000);   setGlobalIfUndef("obtOvlRefBlockLength", 0);
    setGlobalIfUndef("utgOvlHashBlockLength", 100000000);   setGlobalIfUndef("utgOvlRefBlockSize", 2000000);   setGlobalIfUndef("utgOvlRefBlockLength", 0);

    #  The sequential overlap store build is mostly memory agnostic.  With lots of overlaps, smaller
    #  memory sizes can run out of open file handles.  This really should be using the number of
    #  overlaps found.

    if      (getGlobal("genomeSize") < adjustGenomeSize("100m")) {
        setGlobalIfUndef("ovlStoreMemory", 4);
    } else {
        setGlobalIfUndef("ovlStoreMemory", 8);
    }

    #  The parallel store build really doesn't change much.  Also should be based on the number of overlaps found.

    if      (getGlobal("genomeSize") < adjustGenomeSize("100m")) {
        setGlobalIfUndef("ovbMemory",   "2-4");     setGlobalIfUndef("ovbThreads",   "1");
        setGlobalIfUndef("ovsMemory",   "4-16");    setGlobalIfUndef("ovsThreads",   "1");
        setGlobalIfUndef("ovlStoreSlices", 16);

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("1g")) {
        setGlobalIfUndef("ovbMemory",   "2-4");     setGlobalIfUndef("ovbThreads",   "1");
        setGlobalIfUndef("ovsMemory",   "8-24");    setGlobalIfUndef("ovsThreads",   "1");
        setGlobalIfUndef("ovlStoreSlices", 64);

    } else {
        setGlobalIfUndef("ovbMemory",   "2-4");     setGlobalIfUndef("ovbThreads",   "1");
        setGlobalIfUndef("ovsMemory",   "8-32");    setGlobalIfUndef("ovsThreads",   "1");
        setGlobalIfUndef("ovlStoreSlices", 128);
    }

    #  Correction and consensus are likewise somewhat invariant.

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("cnsMemory",     "8-32");     setGlobalIfUndef("cnsThreads",      "1-4");
        setGlobalIfUndef("corMemory",     "6-16");     setGlobalIfUndef("corThreads",      "1-4");
        setGlobalIfUndef("cnsPartitions", "8");        setGlobalIfUndef("cnsPartitionMin", "15000");
        setGlobalIfUndef("corPartitions", "128");       setGlobalIfUndef("corPartitionMin", "5000");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("1g")) {
        setGlobalIfUndef("cnsMemory",     "16-48");    setGlobalIfUndef("cnsThreads",      "2-8");
        setGlobalIfUndef("corMemory",     "10-16");    setGlobalIfUndef("corThreads",      "2-8");
        setGlobalIfUndef("cnsPartitions", "64");       setGlobalIfUndef("cnsPartitionMin", "20000");
        setGlobalIfUndef("corPartitions", "256");      setGlobalIfUndef("corPartitionMin", "15000");

    } else {
        setGlobalIfUndef("cnsMemory",     "16-64");    setGlobalIfUndef("cnsThreads",      "2-8");
        setGlobalIfUndef("corMemory",     "10-16");    setGlobalIfUndef("corThreads",      "2-8");
        setGlobalIfUndef("cnsPartitions", "256");      setGlobalIfUndef("cnsPartitionMin", "25000");
        setGlobalIfUndef("corPartitions", "512");      setGlobalIfUndef("corPartitionMin", "25000");
    }

    #  Meryl too, basically just small or big.  This should really be using the number of bases
    #  reported from gatekeeper.

    if      (getGlobal("genomeSize") < adjustGenomeSize("100m")) {
        setGlobalIfUndef("merylMemory", "4-8");     setGlobalIfUndef("merylThreads", "1-4");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("1g")) {
        setGlobalIfUndef("merylMemory", "16-64");    setGlobalIfUndef("merylThreads", "1-16");

    } else {
        setGlobalIfUndef("merylMemory", "64-256");   setGlobalIfUndef("merylThreads", "1-32");
    }

    #  Overlap error adjustment
    #
    #  Configuration is primarily done though memory size.  If that blows up for some reason,
    #  the actual number of reads (batchSize) or bases (batchLength) can be restricted.  I expect
    #  those to be used only from the command line, so they're left unset here.
    #
    #setGlobalIfUndef("redBatchSize", "");    setGlobalIfUndef("redBatchLength", "");
    #setGlobalIfUndef("oeaBatchSize", "");    setGlobalIfUndef("oeaBatchLength", "");

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("redMemory",   "2-8");    setGlobalIfUndef("redThreads",   "1-4");
        setGlobalIfUndef("oeaMemory",   "2");      setGlobalIfUndef("oeaThreads",   "1");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("500m")) {
        setGlobalIfUndef("redMemory",   "4-12");    setGlobalIfUndef("redThreads",   "1-6");
        setGlobalIfUndef("oeaMemory",   "2");       setGlobalIfUndef("oeaThreads",   "1");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("2g")) {
        setGlobalIfUndef("redMemory",   "4-16");    setGlobalIfUndef("redThreads",   "1-8");
        setGlobalIfUndef("oeaMemory",   "2");       setGlobalIfUndef("oeaThreads",   "1");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("5g")) {
        setGlobalIfUndef("redMemory",   "4-32");    setGlobalIfUndef("redThreads",   "1-8");
        setGlobalIfUndef("oeaMemory",   "2");       setGlobalIfUndef("oeaThreads",   "1");

    } else {
        setGlobalIfUndef("redMemory",   "4-32");    setGlobalIfUndef("redThreads",   "1-8");
        setGlobalIfUndef("oeaMemory",   "2");       setGlobalIfUndef("oeaThreads",   "1");
    }

    #  And bogart.

    if      (getGlobal("genomeSize") < adjustGenomeSize("40m")) {
        setGlobalIfUndef("batMemory",   "2-16");        setGlobalIfUndef("batThreads",   "1-4");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("500m")) {
        setGlobalIfUndef("batMemory",   "8-64");        setGlobalIfUndef("batThreads",   "2-8");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("2g")) {
        setGlobalIfUndef("batMemory",   "32-256");      setGlobalIfUndef("batThreads",   "4-16");

    } elsif (getGlobal("genomeSize") < adjustGenomeSize("5g")) {
        setGlobalIfUndef("batMemory",   "128-512");     setGlobalIfUndef("batThreads",   "8-32");

    } else {
        setGlobalIfUndef("batMemory",   "256-1024");    setGlobalIfUndef("batThreads",   "16-64");
    }

    #  Finally, use all that setup to pick actual values for each component.

    my $err;
    my $all;

    ($err, $all) = getAllowedResources("",    "bat",      $err, $all);
    ($err, $all) = getAllowedResources("cor", "mecat2asmpw",     $err, $all);
    ($err, $all) = getAllowedResources("obt", "mecat2asmpw",     $err, $all);
    ($err, $all) = getAllowedResources("utg", "mecat2asmpw",     $err, $all);
    ($err, $all) = getAllowedResources("",    "red",      $err, $all);
    ($err, $all) = getAllowedResources("",    "oea",      $err, $all);
    ($err, $all) = getAllowedResources("",    "cns",      $err, $all);
    ($err, $all) = getAllowedResources("",    "ovlStore", $err, $all);
    ($err, $all) = getAllowedResources("",    "ovb",      $err, $all);
    ($err, $all) = getAllowedResources("",    "ovs",      $err, $all);
    ($err, $all) = getAllowedResources("cor", "ovl",      $err, $all);
    ($err, $all) = getAllowedResources("obt", "ovl",      $err, $all);
    ($err, $all) = getAllowedResources("utg", "ovl",      $err, $all);
    ($err, $all) = getAllowedResources("",    "meryl",    $err, $all);
    ($err, $all) = getAllowedResources("",    "cor",      $err, $all);

    print STDERR "--\n" if (defined($err));
    print STDERR $err   if (defined($err));
    print STDERR "--\n";
    print STDERR $all;
}


