#!/usr/bin/perl
use strict;
#use Data::Dumper;
while (<>) {
  chomp;
  my @t=split(/\t/);
  my ($chr, $strand, $exons)=split(/:/, $t[1]);
  my @ex=map { [ split('-',$_) ] } (split(',', $exons));
  my $len=0;
  map { $len+=($$_[1]-$$_[0]+1) } @ex;
  print ">$t[0] $strand $chr:$ex[0]->[0]-$ex[-1]->[1] $len ",
   scalar(@ex), " $t[3] $t[4]\n";
  my @genes=split(/,/, $t[3]);
  my %og; #overlapping genes
  map { $og{$_}={} } @genes;
  my @txs=split(/,/,$t[2]);
  foreach my $tx (@txs) {
    my ($c,$t,$g)=split(/\|/,$tx);
    if (my $gh=$og{$g}) { #overlapping gene hash
      my $tl=$$gh{$c};
      if ($tl) {
        $$gh{$c}=$tl.','.$t;
      } else {
        $$gh{$c}=$t;
      }
    }
  }
  foreach my $g (@genes) {
    print sprintf("%26s  ", $g);
    my $ch=$og{$g}; #hash ref
    my @c=keys(%$ch);
    for (my $i=0;$i<@c;$i++) {
      print sprintf("%26s  ", ' ') if $i;
      print $c[$i]," : ",${$ch}{$c[$i]}."\n";
    }
  }
}

