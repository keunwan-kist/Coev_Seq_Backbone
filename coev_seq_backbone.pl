use strict;

use Storable;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

my $usage = "Usage: $0 -a3m <a3m fn: first ref seq> -pdb <pdb> -chain <A> -dist <dist cut: 5.0> -pair <SMRF .pair file> -zscore <Zscore cutoff: 0.5> -bsite <binding site file> -resf_outfn <filename> -resf_by_hotspot <0> -ncf_len <5> -ncf_var_std <1.5> -aa_ratio <0.1> -depth <10> -coev_pair_seq_dist <5> -exclude_bsite <1>\n";

sub print_usage {
    print <<"USAGE";
Usage: perl script.pl -a3m <file> -pdb <file> -pair <file> -bsite <file>

Options:
-a3m         A3m file 
-pdb         PDB file (required)
-pair        PMRF coev pair file (required)
-bsite       Binding site file (required)
-resf_outfn  RESFILE out file
-h           Show help message

Example:
  perl coev_seq_backbone.pl -a3m psiblast_nr90.a3m -pdb crtw_alphafold.pdb -pair hhblits_alignment_files/psiblast_crtw_str_nr90.pair -resf_by_hotspot 1 -zscore 0.5 -bsite bs_4A.txt

USAGE
}

# ---- Special Arguments ----
my %args = (
	zscore  => 0.5,
	ncf_len => 5,
	aa_ratio => 0.1, 
	depth => 10, 
	ncf_var_std => 1.5,
	ncf_len => 5,
  dist => 7,
	exclude_bsite => 1,
  resf_by_hotspot => 0,
	coev_pair_seq_dist => 5,
	chain => "A",
);

GetOptions("a3m:s"=>\$args{a3m},
           "pdb:s"=>\$args{pdb},
           "chain:s"=>\$args{chain},
           "dist:s"=>\$args{dist},
           "pair:s"=>\$args{pair},
           "zscore:s"=>\$args{zscore},
           "bsite:s"=>\$args{bsite},
           "ncf_len:i"=>\$args{ncf_len},
           "aa_ratio:s"=>\$args{aa_ratio},
           "resf_outfn:s"=>\$args{resf_outfn},
           "depth:s"=>\$args{depth},
           "resf_by_hotspot:i"=>\$args{resf_by_hotspot},
           "coev_pair_seq_dist:i"=>\$args{coev_pair_seq_dist},
           "exclude_bsite:i"=>\$args{exclude_bsite},
           "ncf_var_std:s"=>\$args{ncf_var_std},
					 'h:s'   => \$args{h},
					 'help:s'   => \$args{help},
) or die "Incorrect usage!\n";

my $a3m_ali = $args{a3m} or (print_usage() and die $usage);
my $pair = $args{pair} or (print_usage() and die $usage);
my $ZSCORE_CUT  = $args{zscore};
my $pdb = $args{pdb};
my $chain = $args{chain};
my $DIST_CUT = $args{dist};
my $bsite = $args{bsite};
my $NCF_CUT = $args{ncf_len};
my $coev_pair_seq_dist = $args{coev_pair_seq_dist};			# sequence distance to be a coev pair, default 5 (only applied for pymol printing) 
my $AA_RATIO_CUT = $args{aa_ratio};	# aa ratio cut, evo profile for coev pos
my $NCF_VAR_STD = $args{ncf_var_std};			# variation for remodelled loop length, allowed until original_length x 1.5 (default) fold  
my $DEPTH = $args{depth};			# def of coev depth from hotspots 
my ($ali_fn,$ali_path,$ali_suf)=fileparse($a3m_ali,".a3m");
my ($pdb_fn,$pdb_path,$pdb_suf)=fileparse($pdb,".pdb");
my $resf_outfn = $args{resf_outfn} || "$ali_fn.$pdb_fn.cb.resf";
my $resf_by_hotspot = $args{resf_by_hotspot};	# is resf defined by hotspot-connected backbone? or 0 - all coev backbone 
my $exclude_bsite = $args{exclude_bsite};	# disabling designs for binding sites? or 0 - designing binding sites  
my ($name,$path,$suffix) = fileparse($pdb,".pdb");

if(defined $args{h} or defined $args{help}) {
	print_usage();
	exit;
}
print "Arguments received:\n";
foreach my $key (keys %args) {
	print "$key => $args{$key}\n";
}
print "resf_outfn -> $resf_outfn\n";
print "=====================================\n\n";

my %BI;
if(defined $bsite){
	open(IN,$bsite) or die "can't open file:$bsite\n";
	while(<IN>){
		chomp;
		next if(/^$/);
		$BI{$_}=1;
	}
	close IN;
}

my (%Res,%PDB,%HET,%RI,%Pocket);

if(defined $pdb){		# if pdb option is given 

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

			$HET{"$aa3 $resi"}{$ele}{x}=$x;
			$HET{"$aa3 $resi"}{$ele}{y}=$y;
			$HET{"$aa3 $resi"}{$ele}{z}=$z;
		}
	}
	close IN;

	## pdb loading finish here 
	## now calculating residue-residue contacts and Pocket info from pdb 

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
}	# PDB option END


die "No binding site info, please use -bsite or put HETATM in -pdb option!\n" if(keys %BI < 1 and keys %Pocket < 1);
# now we have a Pocket or BI to define binding site 
for my $ri (keys %Pocket){
	$BI{$ri}=1;		# combine Pocket + BI 
}

## Reading coev pair info 
my %CoevPair;
my %CoevPos;

open(IN,$pair) or die "can't open file:$pair\n";
<IN>;<IN>;
while(<IN>){
	chomp;
	my $line=$_;
	$line=~s/^\s+//;
	my ($pos1,$pos2,$score,$Zscore)=split;

	next if($Zscore < $ZSCORE_CUT);

	my $pdb_pos1=$pos1;		# index pdb must start with 1 
	my $pdb_pos2=$pos2;
	
	my $res_pos1=$PDB{$pdb_pos1}{aa1};
	my $res_pos2=$PDB{$pdb_pos2}{aa1};

	# only for structually-close positions..  
	if(defined $pdb){	
		next unless(defined	$RI{$pdb_pos1}{$pdb_pos2} or defined $RI{$pdb_pos2}{$pdb_pos1});	
	}else{
		#do nothing.. to consider long-range interaction?
	}
	
	my $seq_dist = abs($pdb_pos1-$pdb_pos2);

	# $CoevPair{"$res_pos1 $pdb_pos1"}{"$res_pos2 $pdb_pos2"}{seq_dist} = $seq_dist;  
	# $CoevPair{"$res_pos1 $pdb_pos1"}{"$res_pos2 $pdb_pos2"}{coev_zscore} = $Zscore;  

	$CoevPair{$pdb_pos1}{$pdb_pos2}{seq_dist} = $seq_dist;  
	$CoevPair{$pdb_pos1}{$pdb_pos2}{coev_zscore} = $Zscore;  
	$CoevPair{$pdb_pos2}{$pdb_pos1}{seq_dist} = $seq_dist;		# same 
	$CoevPair{$pdb_pos2}{$pdb_pos1}{coev_zscore} = $Zscore;		# same

	$CoevPos{$pdb_pos1}=$res_pos1;
	$CoevPos{$pdb_pos2}=$res_pos2;
}
close IN;


# defining coev-backbone connected to the binding hot spots.. 
#my %shell;
#for my $pdb_idx (keys %BI){   # initializing shell with binding sites 
#	$shell{$pdb_idx}{$pdb_idx} = "START";
#}
#for(my $d=1;$d<=$DEPTH;$d++){  # second shell residues
#	for my $p1 (keys %shell){
#    for my $p2 (keys %{$CoevPair{$p1}}){		# p2 is added to the shell with depth info 
#      $shell{$p2}{$p1}="$d $CoevPair{$p1}{$p2}{coev_zscore}" if(not defined $shell{$p2}{$p1});
#    }
#  }
#}

my %shell;
my @current_boundary;  # 현재 depth에서 주변을 탐색할 residue 목록

# 1. 초기 binding site 세팅
for my $pdb_idx (keys %BI) {
	$shell{$pdb_idx} = "START"; 
	push @current_boundary, $pdb_idx;
}

# 2. BFS search
for (my $d = 1; $d <= $DEPTH; $d++) {
	my @next_boundary; # 다음 depth를 위한 임시 배열

	for my $p1 (@current_boundary) {
		# CoevPair에 p1 데이터가 안전하게 존재하는지 확인
		next unless exists $CoevPair{$p1} && ref $CoevPair{$p1} eq 'HASH';

		for my $p2 (keys %{$CoevPair{$p1}}) {
			# p2가 이미 binding site이거나, 이미 shell에 등록된 노드라면 건너뜀 (최단거리 유지)
			next if exists $BI{$p2} || exists $shell{$p2};
			$shell{$p2} = "$d $p1:$p2 $CoevPair{$p1}{$p2}{coev_zscore}";
			push @next_boundary, $p2;
		}
	}
	@current_boundary = @next_boundary;
}

#print Dumper \%shell;exit;

# reading alignment info 
my %Ali;
my $ali_cnt=0;
my $id;
open(IN,$a3m_ali) or die "can't open file:$_\n";
while(<IN>){
	chomp;
	if(/^>(.+)/){
		$id = $1;
		
		my $seq = <IN>;	# seq must be in one line!
		chomp $seq;
		
		if($ali_cnt==0){
			$Ali{ref}{id}=$id;
			$Ali{ref}{seq}=$seq;
		}else{
			$Ali{"s$ali_cnt"}{id}=$id;
			$Ali{"s$ali_cnt"}{seq}=$seq;
		}
		$ali_cnt++;
	}
}
close IN;

print STDERR "$ali_cnt sequences are read from $a3m_ali!\n"; 

my $ref_seq = $Ali{ref}{seq};
my $ref_seq_len = length $ref_seq;

my @coev_bb_index;
for(my $i=0;$i<$ref_seq_len;$i++){
	my $idx = $i+1;

	if(defined $BI{$idx}){		# binding sites are always a coev-bb 
		push @coev_bb_index,2;
	}elsif(defined $CoevPair{$idx}){
		push @coev_bb_index,1;
	}else{
		push @coev_bb_index,0;
	}
}
print "Coev BB Map by Zscore $ZSCORE_CUT\n@coev_bb_index\n";

my @coev_bb_index_hotspot;
for(my $i=0;$i<$ref_seq_len;$i++){
	my $idx = $i+1;
	
	if(exists $BI{$idx}){		# binding sites are always a coev-bb 
		push @coev_bb_index_hotspot,2;
	}elsif(exists $shell{$idx}){
		push @coev_bb_index_hotspot,1;
	}else{
		push @coev_bb_index_hotspot,0;
	}
}
print "Hotspot-connected Coev BB Map by Zscore $ZSCORE_CUT\n@coev_bb_index_hotspot\n";


## here, checking evo profile for coev-backbone 
my $Ali_n = 0;
my %Coev_aa_profile;
for my $ss (keys %Ali){

	my $ali_seq = $Ali{$ss}{seq};
	my $ali_seq_id = $Ali{$ss}{id};

	$ali_seq =~ s/[a-z]//g;		# delete all insertions 
	my $ali_seq_len = length $ali_seq;

	die "$ss $ali_seq_id length is not the same to ref_seq!\nA $ali_seq\nR $ref_seq\n" if($ali_seq_len != $ref_seq_len);

	for my $resi (keys %CoevPos){
		my $aa = $CoevPos{$resi};
		my $ali_aa = substr($ali_seq,$resi-1,1);
		push @{$Coev_aa_profile{$resi}{aas}},$ali_aa;
	}
	$Ali_n++;
}

my (%AA, %AA20);
my @aa_keys = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","-","X");
my @values = ("0","0") x 11;
@AA{@aa_keys} = @values;
my @aa_keys20 = ("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y");
my @values20 = ("0","0") x 10;
@AA20{@aa_keys20} = @values20;

for my $i (keys %Coev_aa_profile){

	for my $aa (@{$Coev_aa_profile{$i}{aas}}){
		next if(not defined $AA{uc $aa});
		$Coev_aa_profile{$i}{ratio}{$aa}++;
	}

	for my $aa (@aa_keys){
		if(defined $Coev_aa_profile{$i}{ratio}{$aa}){
			$Coev_aa_profile{$i}{ratio}{$aa} /= $Ali_n;
		}else{
			$Coev_aa_profile{$i}{ratio}{$aa} = 0;
		}
	}
}

##############################################
## now writing a resfile for coev backbone
################################################
open(OUT,">$resf_outfn") or die "can't open file:$resf_outfn\n";
print OUT "NATRO\nstart\n";	# 2 lines have to be fixed


my @resf_aa_pos;
if($resf_by_hotspot == 1){
	print STDERR "\nResfile is defined only for hotspot-connected coev backbone!\n";
	@resf_aa_pos = @coev_bb_index_hotspot;		# only for hotspot-connected backbone
}else{
	print STDERR "\nResfile is defined only for all coev backbone!\n";
	@resf_aa_pos = @coev_bb_index;			# for all coev backbone
}

my @final_design_aa_pos;
for(my $i=0;$i<@resf_aa_pos;$i++){
	
	next if($resf_aa_pos[$i] == 0);		# next if not coev backbone
	next if($exclude_bsite==1 and $resf_aa_pos[$i] == 2);		# bsite is not allowed to design.. 

	my $resi = $i+1;
	my $aa = $CoevPos{$resi};
	my %cur_aas;

	for my $aa (keys %{$Coev_aa_profile{$resi}{ratio}}){
		next if(not defined $AA20{$aa});		# next if not stardard aa, e.g. "-", "x"
		$cur_aas{$aa}=1 if($Coev_aa_profile{$resi}{ratio}{$aa} > $AA_RATIO_CUT);
	}

	my $ori_aa = substr($ref_seq,$i,1);
	$cur_aas{$ori_aa}=1 if(defined $AA20{$ori_aa}); # put original aa in ref seq to resf
	
	my $cur_aas_n = keys %cur_aas;
	
	if($cur_aas_n <= 1){	# if only ref aa
		print OUT "$resi $chain NATAA\n";
	}else{
		my $cur_aas_str = join("",keys %cur_aas);
		print OUT "$resi $chain PIKAA $cur_aas_str\n";
	}
	push @final_design_aa_pos,$resi;
}
close OUT;


my $final_design_aa_pos_pymol_sele = join("+",@final_design_aa_pos);

print STDERR "FILE: $ali_fn.$pdb_fn.cb.resf resfile is generated!\n";
print "select Design_Sites, resi $final_design_aa_pos_pymol_sele\n";

my @binding_site_ary = sort {$a<=>$b} keys %BI;
my $BI_pos_pymol_sele = join("+",@binding_site_ary);
print "select BI_Sites, name ca and resi $BI_pos_pymol_sele\n\n";

# coevolving among the selected sites
my @coevpair_ary;
my %coevpair_ary_site;
my @coev_to_bs = (@final_design_aa_pos,@binding_site_ary);
for(my $i=0;$i<@coev_to_bs;$i++){
	my $resi = $coev_to_bs[$i];
	for(my $j=$i+1;$j<@coev_to_bs;$j++){
		my $resj = $coev_to_bs[$j];
		next if(abs($resj-$resi) < $coev_pair_seq_dist);	# only for distant region
		if(defined $CoevPair{$resi}{$resj}){		# coevoling..
			push @coevpair_ary,"($resi,$resj)";
			$coevpair_ary_site{$resi}=1;
			$coevpair_ary_site{$resj}=1;
		}
	}
}
my $coevpair_ary_str = join(",",@coevpair_ary);
print "\ncoevolving_pairs = [$coevpair_ary_str]\n";
print 'run C:\D\Programs\pymol_rc\pymol_scripts_python3-version\coevpair_line.py'."\n";
print "[example in pymol]: draw_coevo_lines_cgo( obj_name=\"chCom\", chain_id=\"$chain\", pairs=coevolving_pairs, line_color=cmd.get_color_tuple(\"cyan\"), line_radius=0.5, cgo_object_name=\"coevo_lines\" )\n";
my @coevpair_ary_sites = sort {$a<=>$b} keys %coevpair_ary_site;
print "select coevoing_pairs_site,resi ",join("+",@coevpair_ary_sites),"\n\n";

### printing all coev pairs
my @all_coevpair_ary;
my @all_coevpair_ary_site;
for my $i (keys %CoevPair){
	for my $j (keys %{$CoevPair{$i}}){
		next if(abs($j-$i) < $coev_pair_seq_dist);	# seq distance
		push @all_coevpair_ary,"($i,$j)";
	}
}
my $all_coevpair_ary_str = join(",",@all_coevpair_ary);
print "\nall_coevolving_pairs = [$all_coevpair_ary_str]\n";
print "draw_coevo_lines_cgo( obj_name=\"chCom\", chain_id=\"$chain\", pairs=all_coevolving_pairs, line_color=cmd.get_color_tuple(\"white\"), line_radius=0.5, cgo_object_name=\"all_coevo_lines\" )\n\n";


## okay, lets check non-coevolving fragments (NCFs)
## !! NCFs has to be defined from all coev backbone, not hotspot-connected one? 

my %ncf_ary;
my $ncf_cnt=0;
my @ncf_ary_pymol_sele;

for(my $i=0;$i<@coev_bb_index;$i++){

	my $idx_i = $i+1;
	
	if($coev_bb_index[$i] == 0){	# if not a coev-bb, it becomes NCF 
		my $ncf_len=1;
		for(my $j=$i+1;$j<@coev_bb_index;$j++){
			my $idx_j = $j+1;
			if($coev_bb_index[$j] == 0){
				$ncf_len++;
			}else{
				if($ncf_len >= $NCF_CUT){
					$ncf_cnt++;
					$ncf_ary{$ncf_cnt} = "$idx_i $j";			# real index starting from 1, j+1 == 1, so use j  
					print "NCF$ncf_cnt $idx_i-$j $ncf_len ($coev_bb_index[$i],$coev_bb_index[$j])\n";
					push @ncf_ary_pymol_sele, "resi $idx_i-$j";
				}
				$i=$j;	# updating i index
				last;
			}
		}
	}
}

my $ncf_ary_pymol_sele_str = join(" or ",@ncf_ary_pymol_sele);
print "select NCF_regions, $ncf_ary_pymol_sele_str\n";

# lets write a contigmap for RFdiffusion 

my @contig_ary;
my $start_idx=1;
my $contig_n=0;
for my $i (sort {$a<=>$b} keys %ncf_ary){
	my ($s,$e)=split(/\s+/,$ncf_ary{$i});
	my $len = $e-$s+1;
	my $s_1=$s-1;
	my $contig = "$chain$start_idx-$s_1";
	push @contig_ary,$contig;
	push @contig_ary,"$len-$len";			# temp, consider NCF_VAR_STD, ...  

	$start_idx = $e+1;
	$contig_n++;

	if($contig_n == keys %ncf_ary){		# for last one 
		if($e==$ref_seq_len){		
			# do nothing if e is the last seq
		}else{
			my $last_contig = "$chain$start_idx-$ref_seq_len";
			push @contig_ary,$last_contig;
		}
	}
}

my $full_contig=join("/",@contig_ary);

print "\nDefault contig directly from $pdb\n";
print "contigmap.contigs=[$full_contig]\n";


my @ref_seq_aa_ary = split(//,$ref_seq);

my %ncf_var_info;

for my $i (sort {$a<=>$b} keys %ncf_ary){
	
	my ($s,$e)=split(/\s+/,$ncf_ary{$i});
	my $ref_frg_seq = substr($ref_seq,$s-1,$e-$s+1);
	my $ref_frg_seq_len = length $ref_frg_seq;
	my $min_len = 100000;	
	my $max_len = 0;
	
	$ncf_var_info{$i}{original_length}=$ref_frg_seq_len;

	for my $ss (keys %Ali){

		my $ali_seq = $Ali{$ss}{seq};
		my $ali_seq_len = length $ali_seq;
		my $ali_seq_id = $Ali{$ss}{id};
	
		my %ref_idx_map;
		my @gapped_ref_seq;
		my $matched_idx_in_ali_seq=0;

		for(my $j=0;$j<$ali_seq_len;$j++){
			my $aa = substr($ali_seq,$j,1);
			if($aa =~/[a-z]/){		# insertion in aligned seq
				push @gapped_ref_seq,".";
			}else{
				push @gapped_ref_seq,$ref_seq_aa_ary[$matched_idx_in_ali_seq];		# insertion
					$ref_idx_map{$matched_idx_in_ali_seq+1}=$j+1;		# index mapping, this is real index starting from 1 
					$matched_idx_in_ali_seq++;
			}
		}

		my $gapped_ref_seqs = join("",@gapped_ref_seq);	# now ref has a "." so same length to ali_seq
		# print "$ref_seq\n$ali_seq\n$gapped_ref_seqs\n";

		my $new_s = $ref_idx_map{$s};
		my $new_e = $ref_idx_map{$e};
		my $ncf_frg_seq = substr($ali_seq,$new_s-1,$new_e-$new_s+1);		# s, e is real index starting from 1

		# checking stem region is not a gap 
		my $stem_s = substr($ali_seq,$new_s-1-1,1);
		my $stem_e = substr($ali_seq,$new_e-1+1,1);
	
		next if($stem_s eq "-" or $stem_e eq "-");	# no gap in the stem region.. 

		my $ori_ncf_frg_seq = $ncf_frg_seq;
		$ncf_frg_seq=~s/[\-\.]//g;	# deleting gaps..
		my $ncf_frg_seq_len = length $ncf_frg_seq;

		push @{$ncf_var_info{$i}{all}{$ncf_frg_seq_len}},"$ncf_frg_seq:$ori_ncf_frg_seq:$ss:$ali_seq_id";
		$min_len = $ncf_frg_seq_len if($min_len > $ncf_frg_seq_len);
		$max_len = $ncf_frg_seq_len if($max_len < $ncf_frg_seq_len);
	}

	$ncf_var_info{$i}{range}="$min_len-$max_len";
	#print "min $min_len max $max_len original $ref_frg_seq_len\n";

	my $min_len_flex = int($ref_frg_seq_len/$NCF_VAR_STD); 
	$min_len_flex = $min_len if($min_len > $min_len_flex);		# native length is preferred 
	
	my $max_len_flex = int($ref_frg_seq_len*$NCF_VAR_STD); 
	$max_len_flex = $max_len if($max_len < $max_len_flex);		# native length is preferred 

	$ncf_var_info{$i}{flex_range}="$min_len_flex-$max_len_flex";

}

# print Dumper \%ncf_var_info;				# how to interprete lower letter in a3m?

my @contig_ary2;
my $start_idx2=1;
my $contig_n2=0;
for my $i (sort {$a<=>$b} keys %ncf_ary){
	my ($s,$e)=split(/\s+/,$ncf_ary{$i});
	my $len = $e-$s+1;
	my $s_1=$s-1;
	my $contig = "$chain$start_idx2-$s_1";
	push @contig_ary2,$contig;
	push @contig_ary2,$ncf_var_info{$i}{flex_range};			# consider NCF_VAR_STD, ...  

	$start_idx2 = $e+1;
	$contig_n2++;

	if($contig_n2 == keys %ncf_ary){		# for last one 
		if($e==$ref_seq_len){		
			# do nothing if e is the last seq
		}else{
			my $last_contig = "$chain$start_idx2-$ref_seq_len";
			push @contig_ary2,$last_contig;
		}
	}
}

my $full_contig2=join("/",@contig_ary2);

print "\nFlexible contigs directly from $a3m_ali\n";
print "contigmap.contigs=[$full_contig2]\n";



# to do : contig with binder 
# to do : use evolutionary info to design sequences with proteinMPNN
# to do : connecting rf_diffusion output pdb (flexible contig) with resf file
# to do : alphafold model checking protocol..for proteinMPNN design 


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

