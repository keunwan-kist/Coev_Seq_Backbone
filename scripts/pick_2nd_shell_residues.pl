use strict;
use Data::Dumper;

if(@ARGV<2){
	die "[].pl [coev_net_cytoscape_with_str/slr1393.3.0.sif] [binding_site_pos.txt] [depth : default 1]\n";
}

my $sif = $ARGV[0];
my $pos_info = $ARGV[1];
my $depth = $ARGV[2] || 1;

my %P;
open(IN,$pos_info) or die "can't open file:$_\n";
while(<IN>){
	chomp;
	next if(/^$/);
	$P{$_}=1;
}
close IN;


my %Map;
my %AA;
open(IN,$sif) or die "can't open file:$sif\n";
while(<IN>){
	chomp;
	next if(/^$/);

#	L1 aa V146
	my ($a1,$int,$a2)=split;	
	
	my $aa1=$a1;
	my $aa2=$a2;
	
	$a1=~s/^\w//g;
	$a2=~s/^\w//g;
	
	$AA{$a1}=$aa1;
	$AA{$a2}=$aa2;

	$Map{$aa1}{$aa2}=1;
	$Map{$aa2}{$aa1}=1;	# undirected network 
}
close IN;


my %shell;
my %start_shell;
for my $res (sort {$a<=>$b} keys %P){		# number index
	my $aa = $AA{$res};
	next if(not defined $aa);
	$shell{$aa}=1;	
	$start_shell{$aa}=1;
}

for(my $i=0;$i<$depth;$i++){	# second shell residues

	for my $aa (keys %shell){
		for my $int_res (keys %{$Map{$aa}}){
			#print "$int_res\n";
			$shell{$int_res}=1;
		}
	}
}


for my $res (sort { $a <=> $b } keys %shell){
	next if(defined $start_shell{$res});
	my $i = $res;
	$i=~s/^\w//;
	print "$i\t$res\n";
}


