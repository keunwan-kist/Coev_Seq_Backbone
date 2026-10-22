use strict;

if(@ARGV<1){
	die "[].pl [summary.AM1_1557g2_rn_0001_coev_Test_out.1.0.pos_all_connected]\n";
}

my $summary = $ARGV[0];
my %S;
my $tot=0;
my $ratio = 0.015;
open(IN,$summary) or die "can't open file:$_\n";
while(<IN>){
	chomp;
	# AM1_1557g2_rn_0001_coev_Test_out.1.0.pos_all_connected.0.AM1_1557g2_rn_0001_0001.pdb -553.953 419.610000
	my ($pdb,$energy,$coev_score) = split;
	$S{$pdb}{energy}=$energy;
	$S{$pdb}{coev_score}=$coev_score;

	# AM1_1557g2_rn_0001_coev_Test_out.1.0.pos_all_connected.232.AM1_1557g2_rn_0001_0001.pdb
	$pdb=~/(.+all_connected)/;
	$S{$pdb}{path}="CoevOpt_str/$1/$pdb";
	$tot++;
}
close IN;


my $cnt=0;
for my $p (sort {$S{$b}{coev_score}<=>$S{$a}{coev_score}} keys %S){

	print "$S{$p}{path} $S{$p}{energy} $S{$p}{coev_score}\n";
	$cnt++;

	last if($cnt > $tot*$ratio);

}

#print "\n\n";
$cnt=0;
for my $p (sort {$S{$a}{energy}<=>$S{$b}{energy}} keys %S){

	#print "$S{$p}{path} $S{$p}{energy} $S{$p}{coev_score}\n";
	$cnt++;

	last if($cnt > $tot*$ratio);

}
