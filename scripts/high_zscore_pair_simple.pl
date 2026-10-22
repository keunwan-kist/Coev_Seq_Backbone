use strict;

use Storable;
use Getopt::Long;
use Data::Dumper;
use File::Basename;

my $usage = "Usage: $0 -pair <sif> [-pymol_script <coevpair_line.py>]\n";

sub print_usage {
    print <<"USAGE";
Usage: perl high_zscore_pair_simple.pl -pair <file> [-pymol_script <file>]

Reads a coevolution network in SIF format and prints the coevolving residue
pairs as a python list, ready to paste into PyMOL.

Options:
-pair          SIF format coev pair file (required)
-pymol_script  path to coevpair_line.py; if given, the matching 'run' line
               is printed as well (optional)
-h             Show help message

Example:
  perl high_zscore_pair_simple.pl -pair AM1_1557g2.1.0.sif \\
       -pymol_script /path/to/coevpair_line.py

USAGE
}

# ---- Special Arguments ----
my %args = (
	chain => "A",
);

GetOptions("pair:s"=>\$args{pair},
           "pymol_script:s"=>\$args{pymol_script},
					 'help:s'   => \$args{help},
) or die "Incorrect usage!\n";

my $pair = $args{pair} or (print_usage() and die $usage);

if(defined $args{h} or defined $args{help}) {
	print_usage();
	exit;
}
print "Arguments received:\n";
foreach my $key (keys %args) {
	print "$key => $args{$key}\n";
}


my (%Res,%PDB,%HET,%RI,%Pocket);

## Reading coev pair info 
my %CoevPair;
my %CoevPos;

# F5 aa T9
my @coevpair_ary;
my @coevpair_ary_long;
open(IN,$pair) or die "can't open file:$pair\n";
while(<IN>){
	chomp;
	my $line=$_;
	my ($aa1,$type,$aa2)=split;

	my $pos1 = $aa1;
	my $pos2 = $aa2;
	$pos1=~s/^[A-Z]//;
	$pos2=~s/^[A-Z]//;
	my $diff = abs($pos1-$pos2);
	push @coevpair_ary,"($pos1,$pos2)";
	push @coevpair_ary_long,"($pos1,$pos2)" if($diff > 5);

}
close IN;


my $coevpair_ary_str = join(",",@coevpair_ary);
print "\ncoevolving_pairs = [$coevpair_ary_str]\n";
my $coevpair_ary_str_long = join(",",@coevpair_ary_long);
print "\ncoevolving_pairs = [$coevpair_ary_str_long]\n";
print "run $args{pymol_script}\n" if(defined $args{pymol_script});
print "[example in pymol]: draw_coevo_lines_cgo( obj_name=\"crtw_wt\", chain_id=\"A\", pairs=coevolving_pairs, line_color=cmd.get_color_tuple(\"cyan\"), line_radius=0.5, cgo_object_name=\"coevo_lines\" )\n\n";

