use strict;

if(@ARGV<2){
	die "[].pl [pmrf infer result] [score_energy.txt]\n";
}


my $result = $ARGV[0];
my $energy = $ARGV[1];

my %Coev;

my @lines;
open(IN,$result) or die "can't open file:$result\n";
while(<IN>){
	chomp;
	push @lines,$_;
}
close IN;

#   0    IQDIFRATTQDVRELLKCDRVVVYRFNPDWSGEFVAEATAPNILPLSELEVPTVWEDTYLQENQGGKYRDNATTVVADIYQESYRDCHLEILEWFQARAYLVVPIFQGETLWGLLAAYQNSAPRQWHKVELYLLKQLANQLGVALQQAELLKQVR 419.610000

for(my $i=8;$i<@lines;$i++){
	$lines[$i]=~s/^\s+//;
	my @val=split(/\s+/,$lines[$i]);
	#next unless(@val < 6 or $val[1]=~/^[-\w+]/);
	$Coev{$val[0]}{score}=$val[2];
	$Coev{$val[0]}{seq}=$val[1];
	#print "$val[-1] $val[1]\n";
}



use Data::Dumper;

open(IN,$energy) or die "can't open file:$energy\n";
while(<IN>){
	chomp;
	#AM1_1557g2_rn_0001_coev_Test_out.1.0.pos_all_connected.0.AM1_1557g2_rn_0001_0001.pdb 121.127 10.3419 1 -553.953

	my ($id,$cav,$deg,$geom,$ene)=split;
	next if($geom==0);
	$id=~/all_connected\.(\d+)\./;
	my $coev_score = $Coev{$1}{score};

	print "$id $ene $coev_score\n";

}
close IN;




