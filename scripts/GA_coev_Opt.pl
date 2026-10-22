use EA::GA_COEV_ar;
use Getopt::Long;
use Data::Dumper;

my $usage = "Usage: $0 -pmrf <PMRF model file> -ref_seq <ref sequence file> -N <population size :default 100> -mut_rate <mutation rate : default 0.05> -survive <survive rate : default 0.2> -iter <iteration number : default 100> -allowed_pos <mutatable residue index file> -fixed_pos <critical binding residue index file>\n";


# ---- Special Arguments ----
my %args=();
GetOptions("pmrf:s"=>\$args{pmrf},
     			 "ref_seq:s"=>\$args{ref_seq},
     			 "N:i"=>\$args{N},
     			 "mut_rate:i"=>\$args{mut_rate},
     			 "survive:i"=>\$args{survive},
     			 "iter:i"=>\$args{iter},
     			 "fixed_pos:s"=>\$args{fixed_pos},
     			 "allowed_pos:s"=>\$args{allowed_pos},
) or die "Incorrect usage!\n";

my $PMRF_Model = $args{pmrf} or die $usage;
my $seq_fn = $args{ref_seq} or die $usage;
my $fixed_site_fn = $args{fixed_pos} or die $usage;
my $allowed_site_fn = $args{allowed_pos} or die $usage;
my $pop_size = $args{N} || 100;
my $mut_rate = $args{mut_rate} || 0.05;
my $survive_rate = $args{survive} || 0.2;
my $iter = $args{iter} || 100;


######################################################################
# read fasta file 
open(IN,$seq_fn) or die "can't open file:$seq_fn\n";
my $target='';
while(<IN>){
	chomp;
	next if(/^$/);
	next if(/^>/);
	$target .= $_;
}
close IN;

my %restrict_pos;
open(IN,$fixed_site_fn) or die "can't open file:$fixed_site_fn\n"; 
while(<IN>){
	chomp;
	my @val=split;
	$restrict_pos{$val[0]} = 1;		# hash key - residue number 
}
close IN;

my %allowed_pos;
open(IN,$allowed_site_fn) or die "can't open file:$allowed_site_fn\n"; 
while(<IN>){
	chomp;
	my @val=split;
	$allowed_pos{$val[0]} = 1;		# hash key - residue number 
}
close IN;

######################################################################

# create an array of allowed tokens
@token = ();
for ('A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y') {
	push @token, $_;
}

# New GA object that sets the alphabet and the fitness function
$p = EA::GA_COEV_ar->new(\&fitness_function,\@token,$target,\%allowed_pos,\%restrict_pos);

#$p->set_mut_rate($mut_rate);
#$p->set_survive_rate($survive_rate);

# initialise the population
$p->init_pop($pop_size, $mut_rate);
#$p->dump_indivs();

do {
# breed the next generation using crossover and mutation
	$gen = $p->breed();
	printf STDERR "[%5d]\t", $gen;

# built in data dumper

	$p->dump_best();


# best_fit return the fitness of the best
} while ($p->best_fit() > 0 && $gen < $iter);

print "----------------------------------\n";
print "PMRF $PMRF_Model\n";
print "SEQ $target\n";
print "POPULATION $pop_size\n";
print "MutationRate $mut_rate\n";
print "SurviveRate $survive_rate\n";
print "Iterations $iter\n";
print "----------------------------------\n";

$p->dump_indivs();

exit(0);



# Now all we need is the fitness function that needs to understand
# the representation of an individual

sub fitness_function {
	
	my ($pop_ref) = @_ ;
	
	my $fn = &write_pop_fas($pop_ref);
	my $log = `pmrf infer $PMRF_Model $fn`;
	my @lines=split(/[\r\n]/,$log);
	
	#1     392.65       0.00      85.52       0.00 AM1_1557
	my @fit_scores;
	
	for(my $i=3;$i<@lines;$i++){
		$lines[$i]=~s/^\s+//;
		my @val=split(/\s+/,$lines[$i]);
		next unless(@val < 6 or $val[1]=~/^[-\w+]/); 
		push @fit_scores,$val[1];
		#print "$val[-1] $val[1]\n";
	}
	my $fit_len = @fit_scores;
	
	
	return \@fit_scores;
}

sub write_pop_fas {
	my ($pop_ref)=@_;
	my $fn="pop$$.fas";
	open(Out,">$fn") or die "can't open file:$fn\n";

	my @sums;
	my $pcnt=1;
	for my $i (@{$pop_ref}){
		my @indiv_tokens = @{$i} ;
		my $s = join "", @indiv_tokens;
		print Out ">p$pcnt\n$s\n";
		$pcnt++;
	}	
	close Out;
	
	return $fn;
}

