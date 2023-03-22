
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
 #    src/pipelines/ca3g/Overlapmecat2asmpw.pm
 #
 #  Modifications by:
 #
 #    Brian P. Walenz from 2015-MAR-27 to 2015-SEP-21
 #      are Copyright 2015 Battelle National Biodefense Institute, and
 #      are subject to the BSD 3-Clause License
 #
 #    Brian P. Walenz beginning on 2015-NOV-03
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #    Sergey Koren beginning on 2015-NOV-20
 #      are a 'United States Government Work', and
 #      are released in the public domain
 #
 #  File 'README.licenses' in the root directory of this distribution contains
 #  full conditions and disclaimers for each license.
 ##

package canu::Overlapmecat2asmpw;


require Exporter;

@ISA    = qw(Exporter);
@EXPORT = qw(mecat2asmpwConfigure mecat2asmpwPrecomputeCheck mecat2asmpwCheck);


use strict;

use File::Path qw(make_path remove_tree);

use canu::Defaults;
use canu::Execution;
use canu::Gatekeeper;
use canu::HTML;

#  Map long reads to long reads with mecat2asmpw.

#  Problems:
#  - mecat2asmpw .dat output can't be verified.  if the job is killed, a partial .dat is output.  It needs to write to
#    .dat.WORKING, then let the script rename it when the job finishes successfully.  There is no output file name option.



sub mecat2asmpwConfigure ($$$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $bin     = getBinDirectory();

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    my $path = "$wrk/1-overlapper";

    caFailure("invalid type '$typ'", undef)  if (($typ ne "partial") && ($typ ne "normal"));

    goto allDone   if (skipStage($WRK, $asm, "$tag-mecat2asmpwConfigure") == 1);
    goto allDone   if (-e "$path/precompute.sh") && (-e "$path/mecat2asmpw.sh");
    goto allDone   if (-e "$path/ovljob.files");
    goto allDone   if (-e "$wrk/$asm.ovlStore");

    print STDERR "--\n";
    print STDERR "-- OVERLAPPER (mecat2asmpw) (correction)\n"  if ($tag eq "cor");
    print STDERR "-- OVERLAPPER (mecat2asmpw) (trimming)\n"    if ($tag eq "obt");
    print STDERR "-- OVERLAPPER (mecat2asmpw) (assembly)\n"    if ($tag eq "utg");
    print STDERR "--\n";

    make_path("$path") if (! -d "$path");

    #  mecat2asmpw parameters - filterThreshold needs to be a string, else it is printed as 5e-06.

    my $numHashes       = (getGlobal("${tag}mecat2asmpwSensitivity") eq "normal") ? "512"                                    : "768";
    my $minNumMatches   = (getGlobal("${tag}mecat2asmpwSensitivity") eq "normal") ?   "3"                                    :   "2";
    my $threshold       = (getGlobal("${tag}mecat2asmpwSensitivity") eq "normal") ?   "0.04"                                 :   "0.04";
    my $filterThreshold = (getGlobal("${tag}mecat2asmpwSensitivity") eq "normal") ?   getGlobal("${tag}mecat2asmpwFilterThreshold") :   getGlobal("${tag}mecat2asmpwFilterThreshold"); #  Also set in Meryl.pm

    #  Constants.

    my $merSize       = getGlobal("${tag}mecat2asmpwMerSize");

    my $numReads      = getNumberOfReadsInStore($wrk, $asm);
    my $memorySize    = getGlobal("${tag}mecat2asmpwMemory");
    my $blockPerGb    = getGlobal("${tag}mecat2asmpwBlockSize");
    if ($numHashes > 768) { 
       $blockPerGb = int($blockPerGb / 2);
    }

    my $blockSize = 200000;

    print STDERR "-- Given $memorySize GB, can fit $blockSize reads per block.\n";

    #  Divide the reads into blocks of ovlHashBlockSize.  Each one of these blocks is used as the
    #  table in mecat2asmpw.  Several of these blocks are used as the queries.

    my @blocks;    #  Range of reads to extract for this block
    my @blockBgn;  #  First read in the block
    my @blockLen;  #  Number of reads in the block

    my @hashes;    #  One for each job, the block that is the hash table
    my @skipSelf;  #  One for each job, jobs that would search block N against hash block N need to be handled special
    my @convert;   #  One for each job, flags to the mecat2asmpw-ovl conversion program

    push @blocks,   "no zeroth block, makes the loop where this is used easier";
    push @blockBgn, "no zeroth block";
    push @blockLen, "no zeroth block";
    push @hashes,   "no zeroth job";
    push @skipSelf, "no zeroth job";
    push @convert,  "no zeroth job";

    #for (my $bgn=1; $bgn < $numReads; $bgn += $blockSize) {
    #   my $end = $bgn + $blockSize - 1;
    #   $end = $numReads  if ($end > $numReads);

    #    print STDERR "BLOCK ", scalar(@blocks), " reads from $bgn through $end\n";

    #    push @blocks, "-b $bgn -e $end";

    #   push @blockBgn, $bgn;
    #    push @blockLen, $end - $bgn + 1;
    #}

	open(L2, "< $wrk/filelist") or caFailure("can't open '$wrk/filelist' for reading: $!", undef);
    #for (my $bgn=1; $bgn < $numReads; $bgn += $blockSize) {
	my($range, @tmp, $bgn, $end);
	while( <L2> ){
		#my $end = $bgn + $blockSize - 1;
        #$end = $numReads  if ($end > $numReads);
		chomp;
		@tmp = split(/\t/,$_);
        $bgn = $tmp[0];
		$end = $tmp[1];
		
        print STDERR "BLOCK ", scalar(@blocks), " reads from $bgn through $end\n";
        push @blocks, "-b $bgn -e $end";
        push @blockBgn, $bgn;
        push @blockLen, $end - $bgn + 1;
		
    }
    close(L2);
	
	
	
	
	
    #  Each mecat2asmpw job will process one block against a set of other blocks.  We'll pick, arbitrarily,
    #  to use num_blocks/4 for that size, unless it is too small.

    my $numBlocks = scalar(@blocks);
    my $qryStride = ($numBlocks < 16) ? (2) : int($numBlocks / 4);

    print STDERR "-- For $numBlocks blocks, set stride to $qryStride blocks.\n";
    print STDERR "-- Logging partitioning to '$path/partitioning.log'.\n";

    open(L, "> $path/partitioning.log") or caExit("can't open '$path/partitioning.log' for writing: $!\n", undef);

    #  Make queries.  Each hask block needs to search against all blocks less than or equal to it.
    #  Each job will search at most $qryStride blocks at once.  So, queries could be:
    #  1:  1 vs 1,2,3  (with self-allowed, and quert block 1 implicitly included)
    #  2:  1 vs 4,5    (with no-self allowed)
    #  3:  2 vs 2,3,4
    #  4:  2 vs 5
    #  5:  3 vs 3,4,5
    #  6:  4 vs 4,5
    #  7:  5 vs 5

    #make_path("$path/queries");

    #for (my $bid=1; $bid < $numBlocks; $bid++) {

        #  Note that we never do qbgn = bid; the self-self overlap is special cased.

        #for (my $qbgn = $bid; $qbgn < $numBlocks; $qbgn += $qryStride) {

            #  mecat2asmpw's read labeling is dumb.  It dumps all reads into one array, labeling the hash
            #  reads 0 to N, and the query reads N+1 to N+1+M.  This makes it impossible to tell if
            #  we are comparing a read against itself (e.g., when comparing hash block 1 against
            #  query blocks 1,2,3,4 -- compared to hash block 1 against query blocks 5,6,7,8).
            #
            #  So, there is a giant special case to compare the hash block against itself enabled by
            #  default.  This needs to be disabled for the second example above.
            #
            #  This makes for ANOTHER special case, that of the last block.  There are no query
            #  blocks, just the hash block, and we need to omit the flag...do we?

            #my $andSelf = "";

            #if ($bid == $qbgn) {
                #  The hash block bid is in the range and we need to compute hash-to-hash
                #  overlaps, and exclude the block from the queries.
                #push @skipSelf, "";
                #$qbgn++;
                #$andSelf = " (and self)";
            #} else {
                #  The hash block bid isn't in the range of query blocks, don't allow hash-to-hash
                #  overlaps.
                #push @skipSelf, "--no-self";
            #}

            #my $qend = $qbgn + $qryStride - 1;                 #  Block bid searches reads in dat files from
            #$qend = $numBlocks-1   if ($qend >= $numBlocks);   #  qbgn to qend (inclusive).


            #my $job = substr("000000" . scalar(@hashes), -6);  #  Unique ID for this compute

            #  Make a place to save queries.  If this is the last-block-special-case, make a directory,
            #  but don't link in any files.  Without the directory, we'd need even more special case
            #  code down in mecat2asmpw.sh to exclude the -q option for this last block.

            #make_path("$path/queries/$job");

            #if ($qbgn < $numBlocks) {
                #print L "Job ", scalar(@hashes), " computes block $bid vs blocks $qbgn-$qend$andSelf,\n";

                #for (my $qid=$qbgn; $qid <= $qend; $qid++) {
                    #my $qry = substr("000000" . $qid, -6);             #  Name for the query block

                    #symlink("$path/blocks/$qry.dat", "$path/queries/$job/$qry.dat");
                #}

            #} else {
                #print L "Job ", scalar(@hashes), " computes block $bid vs itself.\n";
                #$qbgn = $bid;  #  Otherwise, the @convert -q value is bogus
            #}

            #  This is easy, the ID of the hash.

            #push @hashes, substr("000000" . $bid, -6);  #  One new job for block bid with qend-qbgn query files in it

            #  Annoyingly, if we're against 'self', then the conversion needs to know that the query IDs
            #  aren't offset by the number of hash reads.

            #if ($andSelf eq "") {
                #push @convert, "-h $blockBgn[$bid] $blockLen[$bid] -q $blockBgn[$qbgn]";
            #} else {
                #push @convert, "-h $blockBgn[$bid] 0 -q $blockBgn[$bid]";
            #}
        #}
    #}

    #close(L);

    #  The ignore file is created in Meryl.pm


    #  The seed length is the shortest read such that all reads longer than this sum to 50x genome size.
    #  genomeSize must be set (canu should be failing early if it isn't).

    my $seedLength = 500;

    {
        my @readLengths;

        open(F, "< $wrk/$asm.gkpStore/reads.txt") or caExit("can't open '$wrk/$asm.gkpStore/reads.txt' for reading: $!", undef);
        while (<F>) {
            my @v = split '\s+', $_;
            push @readLengths, $v[2];
        }
        close(F);

        @readLengths = sort { $b <=> $a } @readLengths;

        my $readLengthSum = 0;
        my $targetSum     = getGlobal("corOutCoverage") * getGlobal("genomeSize");

        foreach my $l (@readLengths) {
            $readLengthSum += $l;

            if ($readLengthSum > $targetSum) {
                $seedLength = $l;
                last;
            }
        }

        print STDERR "-- Computed seed length $seedLength from desired output coverage ", getGlobal("corOutCoverage"), " and genome size ", getGlobal("genomeSize"), "\n";
    }

    #  Create a script to generate precomputed blocks, including extracting the reads from gkpStore.

    #getAllowedResources($tag, "mecat2asmpw");

   # my $javaPath = getGlobal("java");

    open(F, "> $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for writing: $!", undef);
	##########
	##### add:
	##########
	make_path("$path/blocks") if (! -d "$path/blocks");
	open(FFF, "> $path/blocks/ovlprep") or caFailure("can't open '$path/blocks/ovlprep' for writing: $!", undef);
    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "jobid=\$" . getGlobal("gridEngineTaskID") . "\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need " . getGlobal("gridEngineTaskID") . " set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    for (my $ii=1; $ii < scalar(@blocks); $ii++) {
        print F "if [ \$jobid -eq $ii ] ; then\n";
        print F "  rge=\"$blocks[$ii]\"\n";
        print F "  job=\"", substr("000000" . $ii, -6), "\"\n";
        print F "fi\n";
        print F "\n";
		##########
		##### add:
		##########
		print FFF " -allreads -allbases $blocks[$ii]\n";
    }
    print F "\n";
    print F "if [ x\$job = x ] ; then\n";
    print F "  echo Job partitioning error.  jobid \$jobid is invalid.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ ! -d $path/blocks ]; then\n";
    print F "  mkdir $path/blocks\n";
    print F "fi\n";
    print F "\n";
    #print F "if [ -e $path/blocks/\$job.dat ]; then\n";
    #print F "  echo Job previously completed successfully.\n";
    #print F "  exit\n";
    #print F "fi\n";
    #print F "\n";
    print F "#  If the fasta exists, our job failed, and we should try again.\n";
    print F "if [ -e \"$path/blocks/\$job.fasta\" ] ; then\n";
    print F "  rm -f $path/blocks/\$job.dat\n";
    print F "fi\n";
    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    print F "\$bin/gatekeeperDumpFASTQ \\\n";
    print F "  -G $wrk/$asm.gkpStore \\\n";
    print F "  \$rge \\\n";
    print F "  -nolibname \\\n";
    print F "  -fasta \\\n";
    print F "  -o $path/blocks/\$job \\\n";
    print F "|| \\\n";
    print F "mv -f $path/blocks/\$job.fasta $path/blocks/\$job.fasta.FAILED\n";
    print F "\n";
    print F "\n";
    print F "if [ ! -e \"$path/blocks/\$job.fasta\" ] ; then\n";
    print F "  echo Failed to extract fasta.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
	print F "\$bin/fastaconvert -D$path/blocks/\$job.fasta -O$path/blocks/\$job.fa\n";
	print F "rm $path/blocks/\$job.fasta\n";
	print F "mv $path/blocks/\$job.fa $path/blocks/\$job.fasta\n";
    print F "echo Starting mecat2asmpw precompute.\n";
    print F "\n";
    print F "#  So mecat2asmpw writes its output in the correct spot.\n";
    #print F "cd $path/blocks\n";
    #print F "\n";
    #print F "$javaPath -server -Xmx", getGlobal("${tag}mecat2asmpwMemory"), "g \\\n";
    #print F "  -jar \$bin/mecat2asmpw-" . getGlobal("${tag}mecat2asmpwVersion") . ".jar \\\n";
    #print F "  --weighted -k $merSize \\\n";
    #print F "  --num-hashes $numHashes \\\n";
    #print F "  --num-min-matches $minNumMatches \\\n";
    #print F "  --threshold $threshold \\\n";
    #print F "  --filter-threshold $filterThreshold \\\n";
    #print F "  --min-store-length " . ($seedLength-1) . " \\\n";
    #print F "  --num-threads ", getGlobal("${tag}mecat2asmpwThreads"), " \\\n";
    #print F "  -f $wrk/0-mercounts/$asm.ms$merSize.frequentMers.ignore \\\n"   if (-e "$wrk/0-mercounts/$asm.ms$merSize.frequentMers.ignore");
    #print F "  -p $path/blocks/\$job.fasta \\\n";
    #print F "  -q $path/blocks \\\n";
    #print F "|| \\\n";
    #print F "mv -f $path/blocks/\$job.dat $path/blocks/\$job.dat.FAILED\n";
    #print F "\n";
    #print F "if [ ! -e \"$path/blocks/\$job.dat\" ] ; then\n";
    #print F "  echo mecat2asmpw failed.\n";
    #print F "  exit 1\n";
    #print F "fi\n";
    #print F "\n";
    #print F "#  Clean up, remove the fasta input\n";
    #print F "rm -f $path/blocks/\$job.fasta\n";
    #print F "\n";
    print F "exit 0\n";
    close(F);close(FFF);

    #  Create a script to run mecat2asmpw.
    open(F, "> $path/mecat2asmpw.sh") or caFailure("can't open '$path/mecat2asmpw.sh' for writing: $!", undef);

    print F "#!" . getGlobal("shell") . "\n";
    print F "\n";
    print F "jobid=\$" . getGlobal("gridEngineTaskID") . "\n";
    print F "if [ x\$jobid = x -o x\$jobid = xundefined -o x\$jobid = x0 ]; then\n";
    print F "  jobid=\$1\n";
    print F "fi\n";
    print F "if [ x\$jobid = x ]; then\n";
    print F "  echo Error: I need " . getGlobal("gridEngineTaskID") . " set, or a job index on the command line.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    for (my $ii=1; $ii < scalar(@blocks); $ii++) {
        print F "if [ \$jobid -eq $ii ] ; then\n";
        print F "  blk=\"$blocks[$ii]\"\n";
        print F "  slf=\"$skipSelf[$ii]\"\n";
        print F "  cvt=\"-h 1 0 -q 1\"\n";
        print F "  qry=\"", substr("000000" . $ii, -6), "\"\n";
        print F "fi\n";
        print F "\n";
    }
	
	##########
	##### add:
	##########
	my $cvt="-h 0 0 -q 0";
	my $tmp="\$jobid";
	my $qry=substr("000000".$tmp, -6);
	#print F "cvt=\"-h 0 0 -q 0\"";
	#print F "\n";
	#print F "qry=\"$qry\"";
	
    print F "\n";
    print F "if [ x\$qry = x ]; then\n";
    print F "  echo Error: Job index out of range.\n";
    print F "  exit 1\n";
    print F "fi\n";
    print F "\n";
    print F "if [ -e $path/results/\$qry.ovb.gz ]; then\n";
    print F "  echo Job previously completed successfully.\n";
    print F "  exit\n";
    print F "fi\n";
    print F "\n";
    print F "if [ ! -d $path/results ]; then\n";
    print F "  mkdir $path/results\n";
    print F "fi\n";
    print F "\n";
    print F "echo Running block \$blk in query \$qry\n";

    print F "\n";
    print F getBinDirectoryShellCode();
    print F "\n";
    #print F "if [ ! -e \"$path/results/\$qry.mecat2asmpw\" ] ; then\n";
    #print F "  $javaPath -server -Xmx", getGlobal("${tag}mecat2asmpwMemory"), "g \\\n";
    #print F "    -jar \$bin/mecat2asmpw-" . getGlobal("${tag}mecat2asmpwVersion"). ".jar \\\n";
    #print F "    --weighted -k $merSize \\\n";
    #print F "    --num-hashes $numHashes \\\n";
    #print F "    --num-min-matches $minNumMatches \\\n";
    #print F "    --threshold $threshold \\\n";
    #print F "    --filter-threshold $filterThreshold \\\n";
    #print F "    --min-store-length " . ($seedLength-1) . " \\\n";
    #print F "    --num-threads ", getGlobal("${tag}mecat2asmpwThreads"), " \\\n";
    #print F "    -f $wrk/0-mercounts/$asm.ms$merSize.frequentMers.ignore \\\n"   if (-e "$wrk/0-mercounts/$asm.ms$merSize.frequentMers.ignore");
    #print F "    -s $path/blocks/\$blk.dat \$slf \\\n";
    #print F "    -q $path/queries/\$qry \\\n";
    #print F "  > $path/results/\$qry.mecat2asmpw.WORKING \\\n";
    #print F "  && \\\n";
    #print F "  mv -f $path/results/\$qry.mecat2asmpw.WORKING $path/results/\$qry.mecat2asmpw\n";
    #print F "fi\n";
	
	##########
	##### add:
	##########
    	
       if ($tag eq "obt"){
	print F "\$bin/mecat2trimpw -P$wrk/1-overlapper/blocks -T".getGlobal("${tag}mecat2asmpwThreads")." -S\$jobid -E".(scalar(@blocks)-1)."\n";
	}else{
	print F "\$bin/mecat2asmpw -P$wrk/1-overlapper/blocks -T".getGlobal("${tag}mecat2asmpwThreads")." -S\$jobid -E".(scalar(@blocks)-1)."\n";
	}
	
        #print F "\$bin/mecat2asmpw -P$wrk/1-overlapper/blocks -T".getGlobal("${tag}mecat2asmpwThreads")." -S\$jobid -E".(scalar(@blocks)-1)."\n";
	print F "cat $wrk/1-overlapper/blocks/\$\{jobid\}_*.r >$wrk/1-overlapper/results/\$qry.mecat2asmpw"."\n";
	print F "rm -f $wrk/1-overlapper/blocks/\$\{jobid\}_*.r";
    print F "\n";
    print F "if [   -e \"$path/results/\$qry.mecat2asmpw\" -a \\\n";
    print F "     ! -e \"$path/results/\$qry.ovb.gz\" ] ; then\n";
    print F "  \$bin/mecat2asmpwConvert \\\n";
    print F "    \$cvt \\\n";
    print F "    -o $path/results/\$qry.mecat2asmpw.ovb.gz \\\n";
    print F "    $path/results/\$qry.mecat2asmpw\n";
    print F "fi\n";
    print F "\n";

    if (getGlobal("${tag}mecat2asmpwReAlign") eq "raw") {
        print F "if [ -e \"$path/results/\$qry.mecat2asmpw.ovb.gz\" ] ; then\n";
        print F "  \$bin/overlapPair \\\n";
        print F "    -G $wrk/$asm.gkpStore \\\n";
        print F "    -O $path/results/\$qry.mecat2asmpw.ovb.gz \\\n";
        print F "    -o $path/results/\$qry.ovb.gz \\\n";
        print F "    -partial \\\n"  if ($typ eq "partial");
        print F "    -erate ", getGlobal("obtOvlErrorRate"), " \\\n"  if ($typ eq "partial");
        print F "    -erate ", getGlobal("utgOvlErrorRate"), " \\\n"  if ($typ eq "normal");
        print F "    -memory " . getGlobal("${tag}mecat2asmpwMemory") . " \\\n";
        print F "    -t " . getGlobal("${tag}mecat2asmpwThreads") . " \n";
        print F "fi\n";
    } else {
        print F "mv -f \"$path/results/\$qry.mecat2asmpw.ovb.gz\" \"$path/results/\$qry.ovb.gz\"\n";
    }

    print F "\n";
    print F "\n";
    #print F "rm -rf $path/queries/\$qry\n";
    print F "\n";
    print F "exit 0\n";

    close(F);

    if (-e "$path/precompute.sh") {
        my $numJobs = 0;
        open(F, "< $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for reading: $!", undef);
        while (<F>) {
            $numJobs++   if (m/^\s+job=/);
        }
        close(F);

        print STDERR "-- Configured $numJobs mecat2asmpw precompute jobs.\n";
    }

    if (-e "$path/mecat2asmpw.sh") {
        my $numJobs = 0;
        open(F, "< $path/mecat2asmpw.sh") or caFailure("can't open '$path/mecat2asmpw.sh' for reading: $!", undef);
        while (<F>) {
            $numJobs++  if (m/^\s+qry=\"(\d+)\"$/);
        }
        close(F);

        print STDERR "-- Configured $numJobs mecat2asmpw overlap jobs.\n";
    }

  finishStage:
    emitStage($WRK, $asm, "$tag-mecat2asmpwConfigure");
    buildHTML($WRK, $asm, $tag);
    stopAfter("mecat2asmpwConfigure");

  allDone:
}


sub mecat2asmpwPrecomputeCheck ($$$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $attempt = getGlobal("canuIteration");

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    my $path    = "$wrk/1-overlapper";

    goto allDone   if (skipStage($WRK, $asm, "$tag-mecat2asmpwPrecomputeCheck", $attempt) == 1);
    goto allDone   if (-e "$path/precompute.files");
    goto allDone   if (-e "$wrk/$asm.ovlStore");

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/precompute.sh") or caFailure("can't open '$path/precompute.sh' for reading: $!", undef);
    while (<F>) {
        if (m/^\s+job=\"(\d+)\"$/) {
            if (-e "$path/blocks/$1.fasta") {
                push @successJobs, "$path/blocks/$1.fasta\n";
            } else {
                $failureMessage .= "--   job $path/blocks/$1.fasta FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }
    close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " mecat2asmpw precompute jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to precompute mecat2asmpw indices.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        print STDERR "-- mecat2asmpw precompute attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

        emitStage($WRK, $asm, "$tag-mecat2asmpwPrecomputeCheck", $attempt);
        buildHTML($WRK, $asm, $tag);

        submitOrRunParallelJob($WRK, $asm, "${tag}mecat2asmpw", $path, "precompute", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- All ", scalar(@successJobs), " mecat2asmpw precompute jobs finished successfully.\n";

    open(L, "> $path/precompute.files") or caExit("failed to open '$path/precompute.files'", undef);
    print L @successJobs;
    close(L);

    setGlobal("canuIteration", 0);
    emitStage($WRK, $asm, "$tag-mecat2asmpwPrecomputeCheck");
    buildHTML($WRK, $asm, $tag);
    stopAfter("mecat2asmpwPrecompute");

  allDone:
}

sub mecat2asmpwCheck ($$$$) {
    my $WRK     = shift @_;  #  Root work directory (the -d option to canu)
    my $wrk     = $WRK;      #  Local work directory
    my $asm     = shift @_;
    my $tag     = shift @_;
    my $typ     = shift @_;
    my $attempt = getGlobal("canuIteration");

    $wrk = "$wrk/correction"  if ($tag eq "cor");
    $wrk = "$wrk/trimming"    if ($tag eq "obt");
    $wrk = "$wrk/unitigging"  if ($tag eq "utg");

    my $path    = "$wrk/1-overlapper";

    goto allDone   if (skipStage($WRK, $asm, "$tag-mecat2asmpwCheck", $attempt) == 1);
    goto allDone   if (-e "$path/mecat2asmpw.files");
    goto allDone   if (-e "$wrk/$asm.ovlStore");

    #  Figure out if all the tasks finished correctly.

    my $currentJobID   = 1;
    my @mecat2asmpwJobs;
    my @successJobs;
    my @failedJobs;
    my $failureMessage = "";

    open(F, "< $path/mecat2asmpw.sh") or caExit("failed to open '$path/mecat2asmpw.sh'", undef);
    while (<F>) {
        if (m/^\s+qry=\"(\d+)\"$/) {
            if      (-e "$path/results/$1.ovb.gz") {
                push @mecat2asmpwJobs,    "$path/results/$1.mecat2asmpw\n";
                push @successJobs, "$path/results/$1.ovb.gz\n";

            } elsif (-e "$path/results/$1.ovb") {
                push @mecat2asmpwJobs,    "$path/results/$1.mecat2asmpw\n";
                push @successJobs, "$path/results/$1.ovb\n";

            } elsif (-e "$path/results/$1.ovb.bz2") {
                push @mecat2asmpwJobs,    "$path/results/$1.mecat2asmpw\n";
                push @successJobs, "$path/results/$1.ovb.bz2\n";

            } elsif (-e "$path/results/$1.ovb.xz") {
                push @mecat2asmpwJobs,    "$path/results/$1.mecat2asmpw\n";
                push @successJobs, "$path/results/$1.ovb.xz\n";

            } else {
                $failureMessage .= "--   job $path/results/$1.ovb FAILED.\n";
                push @failedJobs, $currentJobID;
            }

            $currentJobID++;
        }
    }
    close(F);

    #  Also find the queries symlinks so we can remove those.  And the query directories, because
    #  the last directory can be empty, and so we'd never see it at all if only finding files.

    #open(F, "find $path/queries -print |");
    #while (<F>) {
    #    push @mecat2asmpwJobs, $_;
    #}
    #close(F);

    #  Failed jobs, retry.

    if (scalar(@failedJobs) > 0) {

        #  If not the first attempt, report the jobs that failed, and that we're recomputing.

        if ($attempt > 1) {
            print STDERR "--\n";
            print STDERR "-- ", scalar(@failedJobs), " mecat2asmpw jobs failed:\n";
            print STDERR $failureMessage;
            print STDERR "--\n";
        }

        #  If too many attempts, give up.

        if ($attempt > getGlobal("canuIterationMax")) {
            caExit("failed to compute mecat2asmpw overlaps.  Made " . ($attempt-1) . " attempts, jobs still failed", undef);
        }

        #  Otherwise, run some jobs.

        print STDERR "-- mecat2asmpw attempt $attempt begins with ", scalar(@successJobs), " finished, and ", scalar(@failedJobs), " to compute.\n";

        emitStage($WRK, $asm, "$tag-mecat2asmpwCheck", $attempt);
        buildHTML($WRK, $asm, $tag);
        submitOrRunParallelJob($WRK, $asm, "${tag}mecat2asmpw", $path, "mecat2asmpw", @failedJobs);
        return;
    }

  finishStage:
    print STDERR "-- Found ", scalar(@successJobs), " mecat2asmpw overlap output files.\n"; 

    open(L, "> $path/mecat2asmpw.files") or caExit("failed to open '$path/mecat2asmpw.files'", undef);
    print L @mecat2asmpwJobs;
    close(L);

    open(L, "> $path/ovljob.files") or caExit("failed to open '$path/ovljob.files'", undef);
    print L @successJobs;
    close(L);

    setGlobal("canuIteration", 0);
    emitStage($WRK, $asm, "$tag-mecat2asmpwCheck");
    buildHTML($WRK, $asm, $tag);
    stopAfter("mecat2asmpwOverlap");

  allDone:
}

