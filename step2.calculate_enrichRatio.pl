#!/usr/bin/perl
die "perl peak.bed sample.reads.bed control.reads.bed\n" if(@ARGV != 3);
my $peak_bed=shift;
my $sample_reads_bed=shift;
my $control_reads_bed=shift;

`bedtools intersect -a $peak_bed -b $sample_reads_bed -c > tmp.peak_sample.count`;
`bedtools intersect -a $peak_bed -b $control_reads_bed -c > tmp.peak_control.count`;
my %peak_sample_count=read_count("tmp.peak_sample.count");
my %peak_control_count=read_count("tmp.peak_control.count");

open(PB,$peak_bed) || die;
while(my $line=<PB>){
	chomp $line;
	my @sub=split/\s+/,$line;
	if($peak_sample_count{$sub[3]} <= 0){
		die;
	}
	print $line,"\t",$peak_sample_count{$sub[3]},"\t",$peak_control_count{$sub[3]},"\n";
}

sub read_count{
	my $file=shift;
	my %hash;
	open(IN,$file) || die;
	while(my $line=<IN>){
		chomp $line;
		my @sub=split/\s+/,$line;
		$hash{$sub[3]}=$sub[-1];
	}
	close IN;
	return %hash;
}
