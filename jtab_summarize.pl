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
  my $excount = scalar(@ex);
  my $ovlcodes='';
  #add introns? $intron_lens
  my $exon_lens = join(',', (map { $$_[1]-$$_[0]+1 } @ex) );
  #my @genes=split(/,/, $t[3]);
  my @txs=split(/,/,$t[2]);
  foreach my $tx (@txs) {
    my ($c,$t,$g)=split(/\|/,$tx);
    $ovlcodes.=$c unless ($c=~m/[iyxs]/ || index($ovlcodes,$c)>=0);
  }
  $ovlcodes='.' unless length($ovlcodes);
  print join("\t", ">$t[0]",$strand, "$chr:$ex[0]->[0]-$ex[-1]->[1]",$len,
    $excount,$t[3], $ovlcodes, $t[4], $exon_lens)."\n";

}

