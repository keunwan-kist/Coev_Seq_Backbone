use strict;

use File::Basename;
use Data::Dumper;
use Storable;

if(@ARGV<4){
	die "[].pl [pdb] [pair output from SMRF] [Zscore cutoff: default 1.5] [binding site list file] [structure|no_structure]\n";
}

my $pdb = $ARGV[0];
my $pair = $ARGV[1];

#########################################
my $DIST_CUT = 5.0;
my $ZSCORE_CUT = $ARGV[2] || 1.5 ;
##########################################

my %BI;
if(defined $ARGV[3]){
	open(IN,$ARGV[3]) or die "can't open file:$ARGV[3]\n";
	while(<IN>){
		chomp;
		next if(/^$/);
		$BI{$_}=1;
	}
	close IN;
}

my $str_opt;
if(defined $ARGV[4]){
	$str_opt =  $ARGV[4];	# consider only close residues or not? 
}


my ($name,$path,$suffix) = fileparse($pdb,".pdb");



my %Res;
my %PDB;
my %HET;
open(IN,$pdb) or die "can't open file:$pdb\n";
while(<IN>){
	chomp;
	if(/^ATOM/){
		my $h = substr($_,13,1);
		next if($h eq "H");
		my $ele = substr($_,12,4); $ele=~s/\s+//g;
		next if($ele eq "C");
		next if($ele eq "O");
		next if($ele eq "N");
		
		my $resi = substr($_,22,4);$resi=~s/\s+//g;
		my $x = substr($_,30,8);$x=~s/\s+//g;
		my $y = substr($_,38,8);$y=~s/\s+//g;
		my $z = substr($_,46,8);$z=~s/\s+//g;
		
		$Res{$resi}{$ele}{x}=$x;
		$Res{$resi}{$ele}{y}=$y;
		$Res{$resi}{$ele}{z}=$z;
		
		my $aa3 = substr($_,17,3);
		my $aa1 = &one_from_three($aa3);
		$PDB{$resi}{aa3}=$aa3;
		$PDB{$resi}{aa1}=$aa1;
	
	}elsif(/^HETATM/){
		my $h = substr($_,13,1);
		next if($h eq "H");
		my $ele = substr($_,12,4); $ele=~s/\s+//g;
		my $resi = substr($_,22,4);$resi=~s/\s+//g;
		my $x = substr($_,30,8);$x=~s/\s+//g;
		my $y = substr($_,38,8);$y=~s/\s+//g;
		my $z = substr($_,46,8);$z=~s/\s+//g;
		
		my $aa3 = substr($_,17,3);
		
		$HET{$aa3}{$ele}{x}=$x;
		$HET{$aa3}{$ele}{y}=$y;
		$HET{$aa3}{$ele}{z}=$z;
	}
}
close IN;

my %RI;
my %Pocket;
my @res= sort {$a<=>$b} keys %Res;
for(my $i=0;$i<@res;$i++){
	my @ele_i = keys %{$Res{$res[$i]}};
	for(my $j=$i+1;$j<@res;$j++){
		my @ele_j = keys %{$Res{$res[$j]}};
		
		my $int_flag=0;
		for my $ei (@ele_i){
			my $xi=$Res{$res[$i]}{$ei}{x};
			my $yi=$Res{$res[$i]}{$ei}{y};
			my $zi=$Res{$res[$i]}{$ei}{z};
			for my $ej (@ele_j){
				my $xj=$Res{$res[$j]}{$ej}{x};
				my $yj=$Res{$res[$j]}{$ej}{y};
				my $zj=$Res{$res[$j]}{$ej}{z};
				
				my $d = sqrt( ($xi-$xj)**2 + ($yi-$yj)**2 + ($zi-$zj)**2  );
				if($d < $DIST_CUT){
					$int_flag=1;last;
				}
			}
			last if($int_flag);
		}
		$RI{$res[$i]}{$res[$j]}=1 if($int_flag);	
	}

	for my $aa3 (keys %HET){
		my $int_flag=0;
		for my $het_ele (keys %{$HET{$aa3}}){
			my $xh=$HET{$aa3}{$het_ele}{x};
			my $yh=$HET{$aa3}{$het_ele}{y};
			my $zh=$HET{$aa3}{$het_ele}{z};
			for my $ei (@ele_i){
				my $xi=$Res{$res[$i]}{$ei}{x};
				my $yi=$Res{$res[$i]}{$ei}{y};
				my $zi=$Res{$res[$i]}{$ei}{z};

				my $d = sqrt( ($xi-$xh)**2 + ($yi-$yh)**2 + ($zi-$zh)**2  );
				if($d < $DIST_CUT){
					$Pocket{$res[$i]}=$d;
					$int_flag=1;last;
				}			
			}
			last if($int_flag);
		} # het ele end	
	}	# aa3 end
}

######################## output files ###########################
my $sif = "$name.$ZSCORE_CUT.sif";
my $edge_att = "$name.$ZSCORE_CUT.edges.attrs";
my $edge_att2 = "$name.$ZSCORE_CUT.edges2.attrs";
my $node_att = "$name.nodes.attrs";

open(O1,">$sif") or die "can't open file:$sif\n";
open(O2,">$edge_att") or die "can't open file:$edge_att\n";
open(O3,">$node_att") or die "can't open file:$node_att\n";
open(O4,">$edge_att2") or die "can't open file:$edge_att2\n";
##################################################################

print O3 "PocketPos\n";
if(keys %BI > 1){
	%Pocket=();		# if file is provided, use the info as a pocket def
	for my $i (sort {$a<=>$b} keys %BI){
		$Pocket{$i}=1;
		#print "$i\n";
		print O3 "$PDB{$i}{aa1}$i = BindingSite\n";
	}
}else{
	for my $i (sort {$a<=>$b} keys %Pocket){
		#print "$i\n";
		print O3 "$PDB{$i}{aa1}$i = BindingSite\n";
	}
}

print O2 "CoevZscore\n";
print O4 "LongRangeInt\n";
open(IN,$pair) or die "can't open file:$pair\n";
<IN>;<IN>;<IN>;
while(<IN>){
	chomp;
#	Pos1 Pos2      Score    Z-score
#	---- ---- ---------- ----------
#	   2    3   0.127695   2.985262

	my $line=$_;
	$line=~s/^\s+//;
	my ($pos1,$pos2,$score,$Zscore)=split;

	next if($Zscore < $ZSCORE_CUT);

	my $pdb_pos1=$pos1;
	my $pdb_pos2=$pos2;
	
	my $aa_pos1=$PDB{$pdb_pos1}{aa1};
	my $aa_pos2=$PDB{$pdb_pos2}{aa1};

	# only for structually-close positions..  
	if(defined $str_opt and $str_opt eq "structure"){	
		next unless(defined	$RI{$pdb_pos1}{$pdb_pos2} or defined $RI{$pdb_pos2}{$pdb_pos1});	
	}else{
		#do nothing.. able to consider long-range interaction???
	}
	
	#print "$aa_pos1$pdb_pos1 $aa_pos2$pdb_pos2 $score $Zscore\n";
	
	print O1 "$aa_pos1$pdb_pos1 aa $aa_pos2$pdb_pos2\n";
	print O2 "$aa_pos1$pdb_pos1 (aa) $aa_pos2$pdb_pos2 = $Zscore\n";
	
	if(abs($pdb_pos1-$pdb_pos2) > 7){	# 2 turns of helix, long range interaction 
		print O4 "$aa_pos1$pdb_pos1 (aa) $aa_pos2$pdb_pos2 = LongRange\n";
	}
}
close IN;


close O1;
close O2;
close O3;


sub one_from_three {
    my $three = shift;
    my $one = ".";
    $one="A" if ($three eq "ALA");
    $one="C" if ($three eq "CYS");
    $one="D" if ($three eq "ASP");
    $one="E" if ($three eq "GLU");
    $one="F" if ($three eq "PHE");
    $one="G" if ($three eq "GLY");
    $one="H" if ($three eq "HIS");
    $one="I" if ($three eq "ILE");
    $one="K" if ($three eq "LYS");
    $one="L" if ($three eq "LEU");
    $one="M" if ($three eq "MET");
    $one="N" if ($three eq "ASN");
    $one="P" if ($three eq "PRO");
    $one="Q" if ($three eq "GLN");
    $one="R" if ($three eq "ARG");
    $one="S" if ($three eq "SER");
    $one="T" if ($three eq "THR");
    $one="V" if ($three eq "VAL");
    $one="W" if ($three eq "TRP");
    $one="Y" if ($three eq "TYR");
    return $one;
}

