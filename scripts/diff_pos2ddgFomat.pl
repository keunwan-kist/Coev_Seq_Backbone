use strict;
use Data::Dumper;
use File::Basename;

=pod
----------------------------------
PMRF psiblast_slrpsiblast2_no_str.mrf
SEQ LQNIFRATSDEVRHLLSCDRVVVYRFNPDWSGEFIHESVAQMWEPLKDLQNNFPLWQDTYLQENEGGKYRNHESTAVGDVETAGYRDCHLDILRRFEIRAYLTVPVFVGEQLWGLLGAYQNGAPRHWQAREIHLLHQIANQLGVAVYQAQLLARFQ
POPULATION 500
MutationRate 0.05
SurviveRate 0.2
Iterations 5000
----------------------------------
   0    LQNIFRATSDEVRHLLSCDRVVVYRFNPDWSGEFVHESVAQMWEPLKDLQNNFPLWQDTYLQENEGGKYRNHESTAVGDVETAGYRDCHLDILRRFEIRAYLTVPVFVGEQLWGLLGAYQNGAPRHWQAREIHLLHQIANQLGVAVYQAQLLARFQ    513.380000
   1    LQNIFRATSDEVRHLLSCDRVVVYRFNPDWSGEFVHESVAQMWEPLKDLQNNFPLWQDTYLQENEGGKYRNHESTAVGDVETAGYRDCHLDILRRFEIRAYLTVPVFVGEQLWGLLGAYQNGAPRHWQAREIHLLHQIANQLGVAVYQAQLLARFQ    513.380000

=cut

if(@ARGV<2){
	die "[].pl [slr1393_coev_slrpsiblast2_0.05.rand56.txt] [slr1393.fasta]\n";
}

my $result = $ARGV[0];
my ($prefix,$path,$suffix) = fileparse($result,".txt");

my $fasta = $ARGV[1];
my $query_seq;

open(IN,$fasta) or die "can't open file:$_\n";
while(<IN>){
	chomp;
	if(/^\>/){
	}else{
		$query_seq = $_;
	}
}
close IN;



my %Coev;

open(IN,$result) or die "can't open file:$_\n";
while(<IN>){
	chomp;
	my $line = $_; 
	$line=~s/^\s+//;
	if($line=~/^\d+/){
		my ($rank,$seq,$coev_score)=split(/\s+/,$line);
		#$Coev{$rank}{seq}=$seq;
		#$Coev{$rank}{score}=$coev_score;
		$Coev{$seq}=$coev_score;
	}
}
close IN;

my @query_seqs = split(//,$query_seq);

my $cnt=0;
my $cur=`pwd`;chomp $cur;
for my $seq (sort {$a<=>$b} keys %Coev){
	my $score = $Coev{$seq};
	print "$query_seq\n";
	print "$seq\n";
	my @seqs=split(//,$seq);
	`mkdir -p $cur/CoevOpt_str`;
	`mkdir -p $cur/CoevOpt_str/$prefix/`;
	open(OUT,">$cur/CoevOpt_str/$prefix/$prefix.$cnt.resfile") or die "can't open file:$_\n";
	
	print OUT "NATRO\n";           
	print OUT "EX 1 EX 2\n";       
	print OUT "USE_INPUT_SC\n";    
	print OUT "start\n";

	for(my $i=0;$i<@query_seqs;$i++){
		if($query_seqs[$i] ne $seqs[$i]){
			my $idx = $i+1;
			print OUT "$idx A PIKAA $seqs[$i]\n";
			print "$idx $query_seqs[$i] !=  $seqs[$i]\n";
		}
	}
	close OUT;
	$cnt++;
}

print STDERR "$prefix $cnt UNIQ sequences\n";

