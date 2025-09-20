#!/usr/bin/perl 


use strict;


my %remove=();
open(INF, "tagged_transgenes.txt");
while(<INF>){
	chomp $_;
	($_=~/^Gene/) && next;
	my @a=split("\t",$_);
	$remove{$a[1]}=1;
print $a[1]."\n";
}
close(INF);


my @files=("human/hs_bedfiles/exon.bed","human/hs_bedfiles/intron.bed","mouse/mm_bedfiles/exon.bed","mouse/mm_bedfiles/intron.bed");
foreach my $file(@files){

	my $output=$file;
	$output =~ s/\.bed/_filtered\.bed/;

	open(INF, $file) || die "Mal:$!\n";
	open(OUTF, ">$output") || die "Mal:$!\n";
	while(<INF>){
		chomp $_;
		my @a=split("\t",$_);
		my $gene=$a[3];
		$gene =~ s/\..*//g;
		($remove{$gene}) && next;
		print OUTF $_."\n";
	}
	close(INF);
}


=cut

