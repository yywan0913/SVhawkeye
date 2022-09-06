#!/usr/bin/perl -w
use strict;

die "Usage:\n$0 cpu_umber shell_file\n" if (@ARGV!=2);
my $cpu=$ARGV[0];
my $shell_file=$ARGV[1];
my @cmd;
open IN,$shell_file;
while (<IN>){
	chomp;
	s/&$//g;
    next if (/^#/);
	next if(! /\S+/);
	push @cmd,$_;
}
close IN;

for (my $i=0; $i<=$#cmd; $i++) {
	my $cmd=$cmd[$i];
#	print "$i\n";
	defined (my $id=fork ()) || die "Can not fork:$!";
	if ( $id !=0 ) {
		wait if($i+1 >= $cpu);
	}else{
		exec $cmd;
		exit 0;
	}
    sleep 0.5;
}
while (wait != -1) { sleep 1; }
0;
