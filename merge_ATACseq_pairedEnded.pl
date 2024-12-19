#!/usr/bin/perl
use strict;

if ($#ARGV < 1) {die "$0 dir sample\n";}
my $dir = shift(@ARGV);
my $sample = shift(@ARGV);
my (%start, %end, %qual);
open(IN, "<$dir/$sample.bed") || die "failed open $dir/$sample.bed\n";
open(OUT, ">$dir/$sample-frag.bed") || die "failed open $dir/$sample-frag.bed for writing\n";
open(OUT1, ">$dir/$sample-frag-free.bed") || die "failed open $dir/$sample-frag-free.bed for writing\n";
open(OUT2, ">$dir/$sample-frag-mono.bed") || die "failed open $dir/$sample-frag-mono.bed for writing\n";
open(OUT3, ">$dir/$sample-frag-di.bed") || die "failed open $dir/$sample-frag-di.bed for writing\n";
open(OUT4, ">$dir/$sample-frag-tri.bed") || die "failed open $dir/$sample-frag-tri.bed for writing\n";
my $chrFLAG=0;
while (my $line = <IN>) {
  #1	3000604	3000698	D52608P1:194:C5F5JACXX:3:2113:8031:39563/2	60	+
  chomp $line;
  my @a = split(/\t/, $line);
  my @b = split(/\//, $a[3]);
  if($_ = ~/^chr/ && $chrFLAG==0){ #added by HJ, on 6/11/2018 
	$chrFLAG=1;
  }
  if (!defined($start{$a[0]}{$b[0]}) || $start{$a[0]}{$b[0]} > $a[1]) {$start{$a[0]}{$b[0]} = $a[1];}
  if (!defined($end{$a[0]}{$b[0]}) || $end{$a[0]}{$b[0]} < $a[2]) {$end{$a[0]}{$b[0]} = $a[2];}
  $qual{$a[0]}{$b[0]} = $a[4];
}
close (IN);

#nature paper
#free: <100
#mono: 180-247
#di: 315-473
#tri: 558-615
my @chr = (1..22, 'X', 'Y');
if ($chrFLAG ==1) {
  @chr=map "chr".$_,@chr; #added by HJ, on 6/11/2018
}
foreach my $chr (@chr) {
  if (!%{$start{$chr}}) {next;}
  #$chr = "chr$chr";
  my @rec = ();  # 615 - 1000
  my @rec1 = (); # <100
  my @rec2 = (); # 150-247
  my @rec3 = (); # 316-473
  my @rec4 = (); # 558-615
  foreach my $rd (keys %{$start{$chr}}) {
	my $size = $end{$chr}{$rd} - $start{$chr}{$rd};
	if ($size > 1000) {next;}
	push @rec, "$start{$chr}{$rd}\t$chr\t$start{$chr}{$rd}\t$end{$chr}{$rd}\t$rd\t$qual{$chr}{$rd}\t+";
	if ($size < 100) {
	  push @rec1, "$start{$chr}{$rd}\t$chr\t$start{$chr}{$rd}\t$end{$chr}{$rd}\t$rd\t$qual{$chr}{$rd}\t+";
	}
	if ($size >= 150 && $size <= 250) {
	  push @rec2, "$start{$chr}{$rd}\t$chr\t$start{$chr}{$rd}\t$end{$chr}{$rd}\t$rd\t$qual{$chr}{$rd}\t+";
	}
	if ($size >= 316 && $size <= 473) {
	  push @rec3, "$start{$chr}{$rd}\t$chr\t$start{$chr}{$rd}\t$end{$chr}{$rd}\t$rd\t$qual{$chr}{$rd}\t+";
	}
	if ($size >= 558 && $size <= 615) {
	  push @rec4, "$start{$chr}{$rd}\t$chr\t$start{$chr}{$rd}\t$end{$chr}{$rd}\t$rd\t$qual{$chr}{$rd}\t+";
	}
  }
  my $str=""; #added by HJ, on 6/11/2018
  my @srec = sort {$a<=>$b;} @rec;
  foreach my $rec (@srec) {
	my @a = split(/\t/, $rec, 2);
	$str= ( $chrFLAG ==1 ? $a[1] : chr$a[1] ); #added by HJ, on 6/11/2018
	print OUT "$str\n";
  }
  @srec = sort {$a<=>$b;} @rec1;
  foreach my $rec (@srec) {
	my @a = split(/\t/, $rec, 2);
	$str= ( $chrFLAG ==1 ? $a[1] : chr$a[1] );
	print OUT1 "$str\n";
  }
  @srec = sort {$a<=>$b;} @rec2;
  foreach my $rec (@srec) {
	my @a = split(/\t/, $rec, 2);
	$str= ( $chrFLAG ==1 ? $a[1] : chr$a[1] );
	print OUT2 "$str\n";
  }
  @srec = sort {$a<=>$b;} @rec3;
  foreach my $rec (@srec) {
	my @a = split(/\t/, $rec, 2);
	$str= ( $chrFLAG ==1 ? $a[1] : chr$a[1] );
	print OUT3 "$str\n";
  }
  @srec = sort {$a<=>$b;} @rec4;
  foreach my $rec (@srec) {
	my @a = split(/\t/, $rec, 2);
	$str= ( $chrFLAG ==1 ? $a[1] : chr$a[1] );
	print OUT4 "$str\n";
  }
}
close (OUT);
close (OUT1);
close (OUT2);
close (OUT3);
close (OUT4);
########################
#update 6/11/2018 to take care of reads with'chr' in its chromosome notation