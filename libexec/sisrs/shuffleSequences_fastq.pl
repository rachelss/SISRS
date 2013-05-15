#!/usr/bin/perl
#script is part of the Velvet de novo assembly package by Daniel Zerbino
#https://github.com/dzerbino/velvet/blob/master/contrib/shuffleSequences_fasta/shuffleSequences_fastq.pl

$filenameA = $ARGV[0];
$filenameB = $ARGV[1];
$filenameOut = $ARGV[2];

open $FILEA, "< $filenameA";
open $FILEB, "< $filenameB";

open $OUTFILE, "> $filenameOut";

while(<$FILEA>) {
print $OUTFILE $_;
$_ = <$FILEA>;
print $OUTFILE $_;
$_ = <$FILEA>;
print $OUTFILE $_;
$_ = <$FILEA>;
print $OUTFILE $_;

$_ = <$FILEB>;
print $OUTFILE $_;
$_ = <$FILEB>;
print $OUTFILE $_;
$_ = <$FILEB>;
print $OUTFILE $_;
$_ = <$FILEB>;
print $OUTFILE $_;
}