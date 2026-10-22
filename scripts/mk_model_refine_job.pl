use strict;
use File::Basename;

if(@ARGV<1){
	die "[].pl [resfile list txt] [rosetta xml: default fast_Design_New.xml] [input pdb: default AM1_1557g2_rn_0001.pdb]\n";
}
my $list = $ARGV[0];
my $xml  = $ARGV[1] || "fast_Design_New.xml";		# RosettaScripts protocol
my $pdb  = $ARGV[2] || "AM1_1557g2_rn_0001.pdb";	# input structure
my @L;
open(IN,$list) or die "can't open file:$_\n";
while(<IN>){
	chomp;
	push @L,$_;	
}
close IN;

my $cur=`pwd`;chomp $cur;

# keep absolute paths as given, resolve relative ones against the current dir
my $xml_abs = ($xml =~ m{^/}) ? $xml : "$cur/$xml";
my $pdb_abs = ($pdb =~ m{^/}) ? $pdb : "$cur/$pdb";

for my $resfile (@L){
	my ($name,$path,$suffix)=fileparse($resfile,".resfile");
		
	#`mkdir -p $cur/$path/$name`;
	#`rmdir $cur/$path/$name`;

	# CoevNetwork_opt/CoevOpt_str/slr1393_nostr_coev_Test_out.2.5.mut_single.all_connected.89.resfile 
	
	#`rm $cur/$path/*slr1393_0001.pdb`;
	#`rm $cur/$path/*score.sc`;

	print "cd $cur/$path;rosetta_scripts.linuxgccrelease -parser:protocol $xml_abs -s $pdb_abs -ex1 -ex2 -use_input_sc -linmem_ig 10 -resfile $cur/$resfile -nstruct 1 -out:prefix $name. -mute all\n";
	

}
# 12 x 5 x 30 = 1800
# 1  - 1120 sec 
# 96 - 107520 sec 
