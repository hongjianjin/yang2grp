#!/usr/bin/perl -w
use strict;

my $spec = shift(@ARGV);
my $chromSizeDir="chromSize";
open(IN, "<".$chromSizeDir.$spec.".sizes") || die "failed open\n";
my %size;
while (my $line = <IN>) {
  chomp $line;
  my @a = split(/\t/, $line);
  $a[0] =~ s/chr//;
  my $t = "chr$a[0]";
  if ($a[0] eq 'MT') {$t = "chrM";}
  $size{$t} = $a[1];
  $size{$a[0]} = $a[1];
}
close (IN);

while (my $line = <>) {
  chomp $line;
  my @a = split(/\s+/, $line);
  if ($a[0] =~ /NT/) {next;}
  if ($a[0] =~ /GL/) {next;}
  if ($a[0] =~ /KB/) {next;}
  if ($a[0] =~ /KI/) {next;}
  if ($a[0] =~ /JH/) {next;}
  $a[0] =~ s/MT/M/;
  if ($a[4] < 1) {
  next;
  }
  if ($a[1] < 0) {$a[1] = 0;}
  if ($a[2] > $size{$a[0]}) {$a[2] = $size{$a[0]};}
  if ($line =~ /^chr/) {
    print join("\t", @a), "\n";
  }
  else {
    print "chr", join("\t", @a), "\n";
  }
}
close (IN);
