package EA::GA_COEV_ar;

# Original code written by Stefan Kamphausen <mail@skamphausen.de> 
# Copyright 2001 Stefan Kamphausen.
# https://www.perlmonks.org/?node_id=55453

# Modified by Keunwan Park <keunwan@kist.re.kr>, 2019.09.25 
# Ref: "Optimizing protein properties with functional constraints", Park et al.(2019)  
# redistributed under the GNU General Public License

####################################################################
##                             LICENSE
####################################################################
# This program is free software; you can redistribute it
# and/or modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
####################################################################


use strict;
use Data::Dumper;

$GA::VERSION='0.1';

sub new {
  my ($class,$fitness,$tokenref,$refseq,$allowed_pos_ref,$restrict_pos_ref) = @_;
  my $self = {};
  # An array which contains all the individuals as arrays of tokens
  $self->{pop} = ();
  # A Hash which contains all the fitness-values adressed by
  # the concatenated tokens
  $self->{fitvals} = ();
  # the 'alphabet' of allowed symbols
  $self->{tokens} = $tokenref;
  # user provided fitness function, gets array of tokens as arg
  $self->{fitness} = $fitness;
  # default value for mutation probability; may be overridden
  $self->{mut_prob} = 0.05;
  # a counter for the generation
  $self->{generation} = 0;
  # reference query sequence to be mutated 
  $self->{refseq} = $refseq;
  # restricted amino acid sites, usually core sites for activity 
  $self->{allowed_pos} = $allowed_pos_ref;
  $self->{restrict_pos} = $restrict_pos_ref;
  # Survived population percentage 
  $self->{survive_rate} = 0.2;
  bless $self, $class;
  return $self;
}

sub init_pop {
  my ($self,$size,$opt_mutation_rate) = @_;
  my $pmut = $opt_mutation_rate || $self->{mut_prob};
  my ($i,$j);
  my $refseq = $self->{refseq};
  my $apos_ref = $self->{allowed_pos};
  my $rpos_ref = $self->{restrict_pos};
  my $seq_len = length($refseq);
  for ($i=0;$i<$size;$i++) {	# population size
	for ($j=0;$j<$seq_len;$j++) {	# sequence length
	  my $aa_site = $j + 1;
	  my $aa = substr($refseq,$j,1);
	  my $rtok;
	  if($aa !~ /[A-Z]/){
	  
		$rtok = $self->random_token();
	  
	  }elsif(defined $rpos_ref->{$aa_site}){	# restriction has priority
		$rtok = $aa;

	  }elsif(defined $apos_ref->{$aa_site}){	# allowed position			
			
		my $ran = rand();
    	if ($ran < $pmut){
			$rtok = $self->random_token();
		}else{
			$rtok = $aa;
		}
	  }else{	
		$rtok = $aa;
      }

	  #print "RTOK: $rtok\n";
	  push @{$self->{pop}[$i]}, $rtok;
	}
	#my $s = join("",@{$self->{pop}[$i]});
	#print "$s\n";
  }
  $self->calculate_fitness();
  $self->sort_pop();
}


sub sort_pop {
  my ($self) = @_;
  @{$self->{pop}} = sort {$self->{fitvals}{join "",@{$b}} <=> $self->{fitvals}{join "",@{$a}}} @{$self->{pop}};
}

sub calculate_fitness {
  my ($self) = @_;
  
  %{$self->{fitvals}} = ();
  my $fits_ref = &{$self->{fitness}}(\@{$self->{pop}});	# treat all elements, renew $self->{fitvals}{seq}
 
 my $lenp=@{$self->{pop}};
 my $lenf= @{$fits_ref};   
  
  for my $i (@{$self->{pop}}) {
	my $f=shift @{$fits_ref};   

	$self->{fitvals}{join "",@{$i}} = $f;
	#print "CALC: Indiv:",@{$i}," $f\n";
  }
}

### Return some values
sub best_fit {
  my $self = shift;
  return $self->{fitvals}{join("",@{$self->{pop}[0]})};
}
sub generation {
  my $self = shift;
  return $self->{generation};
}
sub random_token {
  my ($self) = @_;
  my $max = scalar(@{$self->{tokens}});
  
  my $ran = int rand $max;
  my $tok = @{$self->{tokens}}[$ran];
  #print "Ran: $ran, Max: $max ";
  #print "Token: $tok\n";
  return $tok;
}


### Dawn of The Next Generation 

# combines mutation and crossover
sub breed {
  # FIXME: optional mut_prob
  my $self = shift;
  my $opt_mutation_rate = shift;
  
  my ($p1,$p2,$c,$i);
  # mutation
  my $p_mut = $opt_mutation_rate || $self->{mut_prob};
  # prepare for roulette wheel
  my $a_sum = 0;
  my @fit = ();
  foreach (@{$self->{pop}}) {
	push @fit, $self->{fitvals}{join "",@{$_}};
  }

  @fit = sort {$b <=> $a} @fit;

  #print "BREED: Fitness\n",Dumper(@fit),"\n\n";
  my $length = scalar(@{$self->{pop}});
  
  my @new_pop = ();
  
  # N % survive
  my $survive_n = int( $length * $self->{survive_rate} );
  for(my $i=0;$i<$survive_n;$i++){
	push @new_pop, @{$self->{pop}}[$i];
  }
	

  @fit=@fit[1..$survive_n];

  # Choose Parents
  for ($i=$survive_n;$i<$length;$i++) {
	$p1 = rwheel(\@fit);
	$p2 = rwheel(\@fit);
	#print "P1: $p1 P2: $p2\n";
	
	my $new_gene;
	my $iter=0;
	while($iter<100){
		$new_gene = $self->crossover_mut($p1,$p2,$p_mut);
		my $s = join("",@{$new_gene});
		if(defined $self->{fitvals}{$s}){
			# generate redundant gene, try again for N-times
		}else{
			# uniq gene
			last;
		}
		$iter++;
	}
	
	push @new_pop, $new_gene;
  }
  @{$self->{pop}} = @new_pop;
  $self->calculate_fitness();
  $self->sort_pop();
  return ++$self->{generation};
}

sub mutate {
  my ($self,$rate) = @_;

  my ($ran,$i,$t);
  my $rpos_ref = $self->{restrict_pos};
  my $apos_ref = $self->{allowed_pos};

  my $the_rate = $rate || $self->{mut_prob};
  foreach $i (@{$self->{pop}}) {
	for ($t=0;$t<scalar(@{$i});$t++) {
	  my $aa_site = $t+1;	
	  next if(defined $rpos_ref->{$aa_site});	# next if restrict sites, priority.. 	
	  next if(not defined $apos_ref->{$aa_site});	# next if not allowed sites, only for allowed pos 	
	  $ran = rand();
	  if ($ran < $the_rate) {
		@{$i}[$t] = $self->random_token();
	  }
	}
  }
}

sub crossover_mut {
  my $self = shift;
  my $p1 = shift;
  my $p2 = shift;
  my $opt_mutation_rate = shift;
  my $apos_ref = $self->{allowed_pos};
  my $rpos_ref = $self->{restrict_pos};
  
  my ($ran,$t,$new_size);
  $ran = rand();
  my $pmut = $opt_mutation_rate || $self->{mut_prob};
  my $pp1 = (1.0-$pmut)/2.0;

  # 50:50 for the size of the new one 
  my @new = @{$self->{pop}[$p1]};
  my $new_size = @new;	

  for($t=0;$t<$new_size;$t++) {
	
	my $aa_site = $t+1;	
	
	$ran = rand();
	# 50:50 to take gene from p1 or p2 unless mutation
	if ($ran < $pmut) {
	  
	  if(defined $rpos_ref->{$aa_site}){
		$new[$t] = @{@{$self->{pop}}[$p1]}[$t];	# same aa
	  
	  }elsif(defined $apos_ref->{$aa_site}){
		$new[$t] = $self->random_token();
	  }else{	  
		$new[$t] = @{@{$self->{pop}}[$p1]}[$t];	# same aa
	  }	

	} elsif ($ran < $pp1) {
	  #print "1";
	  $new[$t] = @{@{$self->{pop}}[$p1]}[$t];	  
	} else {
	  #print "2";
	  $new[$t] = @{@{$self->{pop}}[$p2]}[$t];	  
	}
  }
  return \@new;
}

sub crossover {
  my $self = shift;
  my $p1 = shift;
  my $p2 = shift;
  
  my ($ran,$t,$new_size);
  $ran = rand();
  # 50:50 for the size of the new one 
  my @new = ();
  if ($ran < 0.5) {
	$new_size = scalar(@{$self->{pop}[$p1]})
  } else {
	$new_size = scalar(@{$self->{pop}[$p2]})	
  }
  for($t=0;$t<$new_size;$t++) {
	
	$ran = rand();
	# 50:50 to take gene from p1 or p2 unless mutation
	if ($ran < 0.5) {
	  $new[$t] = @{@{$self->{pop}}[$p1]}[$t];	  
	} else {
	  $new[$t] = @{@{$self->{pop}}[$p2]}[$t];	  
	}
  }
  return \@new;
}

### Print-Outs
sub dump_indivs {
  my $self = shift;
  my $i;
  my $len = scalar(@{$self->{pop}});
  for ($i=0;$i<$len;$i++) {
	my $s = join("",@{$self->{pop}[$i]});
	printf "%4d\t",$i;
	print $s;
	printf "\t%5f\n",$self->{fitvals}{$s};
  }
}
sub dump_best {
  my $self = shift;
  my $s = join("",@{$self->{pop}[0]});
  print STDERR $s;
  printf STDERR "\t%5f\n",$self->{fitvals}{$s};
}
sub set_mut_rate {
  my ($self,$mut_rate) = @_;
  $self->{mut_prob} = $mut_rate;
}
sub set_survive_rate {
  my ($self,$survive_rate) = @_;
  $self->{survive_rate} = $survive_rate;
}


### Random
sub rwheel {
  # random element of an array according to it's value
  # aka roulette wheel
  my ($a_ref) = @_;
  my @arr = @{$a_ref};
  
  my $a_sum = 0;
  foreach (@arr){
	$a_sum+= $_
  }
  
  my $sum = 0;
#  print "RWHEEL: length = ",scalar(@arr),"\n";
#  print "RWHEEL ARRAY: ",join(" ",@arr),"\n";
  my $ran = rand $a_sum;
#  print "RWHEEL: RAN $ran < $a_sum\n";
  for (my $i=0;$i<scalar(@arr);$i++) {
	$sum += $arr[$i];
#	print "\tSUM: $sum \$arr[$i] = $arr[$i]\n";
	if ($sum > $ran ) {
	  return $i;
	}
  }
  die "ARGH! I never should have reached this point!\n";
}

1;

