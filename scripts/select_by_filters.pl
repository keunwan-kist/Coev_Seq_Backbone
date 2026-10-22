use strict;
use Data::Dumper;

if(@ARGV<1){
	die "[].pl [pdb]\n";
}

my $pdb = $ARGV[0];

#cav 225.994
#deg 10.6623
#geom 1
#rmsd 1.07917
#score_filter -1149.28


my @filts = ("cav","deg","geom","score_filter");

my %filters;
my $on=0;
open(IN,$pdb) or die "can't open file:$pdb\n";
while(<IN>){
	chomp;
	next if(/^$/);
	if($on == 1){
		my ($term,$val)=split;
		
		
		$filters{$term}=$val;
	}
	
	if(/\#END_POSE_ENERGIES_TABLE/){
		$on=1;
	}


}
close IN;


#print Dumper \%filters;
my @scores;

for my $f (@filts){
	#if($f eq "geom" and $filters{$f} == 0){
	#	exit;
	#}
	#next if($f eq "geom");
	push @scores, $filters{$f};
}
print "$pdb @scores\n";
