use strict;

use File::Basename;
use Data::Dumper;

if(@ARGV<1){
	die "[].pl [pdb]\n";
}

my $pdb = $ARGV[0];
my $msa = $ARGV[1];

my $DIST_CUT = 7.0;

my ($name,$path,$suffix) = fileparse($pdb,".pdb");


my %Res;
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
			
	}
}
close IN;

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
		print "$res[$i] $res[$j]\n" if($int_flag);	
	}
}



