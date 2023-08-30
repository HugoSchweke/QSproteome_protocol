use strict;
use POSIX;
#require("../00_general/seq_general.pm");
require("/media/elusers/users/hugo/15_alphafold/37_revision_Cell/general.pm");
#require("../00_general/CPLX_STATIC.pl");

my @aas = ("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "ASP", "ASN", "LYS", "GLU", "GLN", "ARG", "HIS", "PHE", "CYS", "TRP", "TYR", "MET", "PRO");
my @aa2 = ("G"  , "A",   "V",   "L",   "I",   "S",   "T",   "D",   "N",   "K",   "E",   "Q",   "R",   "H",   "F",   "C",   "W",   "Y",   "M",   "P");
my %AAS; map { $AAS{$aas[$_]}=$aa2[$_] } (0..$#aas);

###
### Takes the alignment string from MMalign and returns chain-based match information, i.e., for each chain,
### The number of matching residues divided by the number of total residues in the chain.
###
#sub MMalign_chainInfo($){
#
#    my ($string) = @_ ;
#    my @sequences = split /\*/, $string;
#
#    my @query_sequences = ();
#    my @target_sequences = ();
#    my @ali = ();
#
#    my $query = 1; # flag to indicate if we are processing query or target sequences
#
#    for my $sequence (@sequences) {
#	if ($sequence =~ /[\.: ]+/) {
#	    $query = 0; # change the flag
#	    push @ali, $sequence;
#
#	} elsif($query==0){
#	    push @target_sequences, $sequence;
#	} else {
#	    push(@query_sequences, $sequence);
#	} 
#    }
#
#    #print STDERR "Q SEQS = ".join("\n",@query_sequences)."\n";
#    #print STDERR "T SEQS = ".join("\n",@target_sequences)."\n";
#
#    # now calculate the scores
#    my @query_scores = ();
#    my @target_scores = ();
#
#    my $Qmin = 100;
#    my $Tmin = 100;
#    my $Qmax = 0;
#    my $Tmax = 0;
#    
#    for my $Q (0..$#query_sequences) {
#	my $query   = $query_sequences[$Q];
#	my $matches = $ali[$Q];
#	$query  =~ s/-//g;
#	$matches =~ s/[^:]//g;
#	my $score;
#	if(length($query)==0){
#	    $score = 0;
#	} else {
#	    $score = sprintf("%.0f",100*length($matches)/length($query));
#	}
#	if($Qmin > $score){
#	    $Qmin=$score;
#	}
#	if($Qmax < $score){
#	    $Qmax=$score;
#	}
#	push @query_scores, $score;
#    }
#
#    for my $Q (0..$#target_sequences) {
#	my $target   = $target_sequences[$Q];
#	my $matches = $ali[$Q];
#	#print STDERR "Trimmed target = $target\n";
#	$target  =~ s/-//g;
#	$matches =~ s/[^:]//g;
#	#print STDERR "Trimmed target = $target\n";
#	#print STDERR "Trimmed matches = $matches\n";
#
#	my $score;
#	if(length($target)==0){
#	    $score=0;
#	} else {
#	    $score = sprintf("%.0f",100*length($matches)/length($target));
#	}
#	
#	if($Tmin > $score){
#	    $Tmin=$score;
#	}
#	if($Tmax < $score){
#	    $Tmax=$score;
#	}
#	push @target_scores, $score;
#    }
#    
#    #print STDERR "Q = @query_scores\n";
#    
#    return($Qmin, $Qmax, $Tmin, $Tmax, join(";",@query_scores), join(";",@target_scores));
#
#}
#
##
## TAKES A PDB FILE AS ARGUMENT AND RETURN A TABLE
##
#sub PDB2TABLE($){
#
#    my ($file) = @_;
#
#    my ($ch_name, $r_name, $a_name, $r_num, $a_num, $x, $y, $z, $ALTERN, $bf) ;
#    my ($line, $ch_size) ;
#    my $ch_name_tmp = "";
#    my $r_name_tmp  = "";
#    my $r_num_tmp   = "";
#    my $ALTERN_TMP  = "";
#    my $to_skip     = -1;
#
#    my %aa = (ALA=>"A", ARG=>"R", ASN=>"N", ASP=>"D", CYS=>"C",
#	      GLN=>"Q", GLU=>"E", GLY=>"G", HIS=>"H", ILE=>"I",
#	      LEU=>"L", LYS=>"K", MET=>"M", PHE=>"F", PRO=>"P",
#	      SER=>"S", THR=>"T", TRP=>"W", TYR=>"Y", VAL=>"V", MSE=>"M");
#    
#    open (PDB_FH, $file) or print STDERR "$file could not be opened in PDB2TABLE - $!\n";
#
#    $ch_size = 0 ;
#
#    ## This is the index used to reprint the PDB later on.
#    my $line_index = 0;
#
#    my @PDBTABLE = ();
#    # ---- goes through the PDB file and load the atoms & positions in a hash
#    while(<PDB_FH>){
#	    
#	$line = $_;
#	##            altern resid is undef OR is in A-form (we skip all B-forms)
#	if($line=~/^ATOM/ && (substr($line,16,1) eq " " ||  substr($line,16,1) eq "A")){
#
#	    $ch_name= substr($line,21,1);	
#		
#	    if($ch_name_tmp eq ""){
#		$ch_name_tmp=$ch_name ;
#	    }
#
#	    #print "The chain is equal to $ch_name, and ch_tmp = $ch_name_tmp --- size = $ch_size\n";
#
#	    $a_num  = substr($line,7,4);  # --- Changes at each line
#	    $a_name = substr($line,11,5);
#	    $x      = substr($line,30,8);
#	    $y      = substr($line,38,8);
#	    $z      = substr($line,46,8); # This has changed!!! BE CAREFULLLLLL AGAIN CHANGED FROM 7 TO 8
#		
#	    $r_name = substr($line,17,4); # --- Changes at each residue / chain
#	    #$ALTERN = substr($line,16,1); # --- If there is an ambiguity, ALTERN can take A,B,C ... values.
#	    $r_num  = substr($line,22,5);
#	    $bf     = substr($line,60,6);
#
#	    # --- Remove white spaces
#	    $a_num  =~ s/ //g ;
#	    $a_name =~ s/ //g ;
#	    $r_name =~ s/ //g ;
#	    $ch_name=~ s/ //g ;
#	    $r_num  =~ s/ //g ;
#	    $x      =~ s/ //g ;
#	    $y      =~ s/ //g ;
#	    $z      =~ s/ //g ;
#	    $bf      =~ s/ //g ;
#	    
#	    if($a_name =~ "CA"){          ### --- push only CA to the table
#
#		my @line = ($r_name, $aa{$r_name}, $ch_name, $r_num, $x, $y, $z, $bf);
#		push(@PDBTABLE, \@line );		
#	    }
#	}
#    }
#    return \@PDBTABLE ;
#}
#
#
#
## takes a hash representing a graph
## produces another way to represent it using
## an 2D array and 1D list:
## The 1D array stores the labels of the nodes
## the 2D array stores the connexions between the nodes.
## HERE the number of contacts is given for each connexion. 
#sub h_2_graph_wedge($$$){
#
#  my ($H, $name, $DOMDEF) = @_ ;
#  my ($nb_nodes, $node, $node_2, @nodes, @labels, @connex, %node_index);
#
#  #
#  # Nodes names
#  #
#  @nodes = sort {$a cmp $b} (keys %$H) ;
#  $nb_nodes = $#nodes+1 ;
#
#  #
#  # Nodes Labels
#  #
#  my $i = 0;
#  foreach $node (@nodes){
#    $node_index{$node} = $i++ ;
#    @labels = (@labels,$H->{$node}->{"prop"}->{$DOMDEF});
#  }
#
#  #
#  # Connectivity matrix
#  #
#  $i = 0 ;
#  foreach $node (@nodes){
#
#    foreach $node_2 (@nodes){
#
#      if(exists $H->{$node}->{$node_2}){
#	$connex[$node_index{$node}][$node_index{$node_2}] = $H->{$node}->{$node_2}->{"size"};
#      } else {
#	$connex[$node_index{$node}][$node_index{$node_2}] = 0 ;#$H->{$node}->{$node_2}->{"size"};
#      }
#    }
#  }
#
#  #
#  # Puts everything in a Hash
#  #
#  my %result ;
#  $result{"labels"} = \@labels ;
#  $result{"connex"} = \@connex ;
#  $result{"real_labs"} = \@nodes ;
#  $result{"name"} = $name ;
#  return \%result ;
#}
#
##
## Takes an EPPIC string and returns whether the interface is biological or not
##
#sub eppic_isbio($){
#    
#    my ($eppic_string) = @_;
#
#    my $number_of_xtal; my $number_of_bio;	
#    if(!defined($eppic_string) || length($eppic_string)==0 || $eppic_string eq ""){
##	print STDERR "COME HERE\n";
#	return(0);
#    }
#
#    if($eppic_string=~/X/){
#	$number_of_xtal = $eppic_string;
#    	$number_of_xtal =~ s/[^X]//g;
#	}
#
#    if($eppic_string=~/B/){	
#    $number_of_bio  = $eppic_string;
#    $number_of_bio  =~ s/[^B]//g;
#    }	
#
#    if(length($number_of_bio) >= length($number_of_xtal)){
#	return(1);
#    } else {
#	return(0);
#    }
#}
#
#
## takes a hash representing a graph
## produces another way to represent it using
## an 2D array and 1D list:
## The 1D array stores the labels of the nodes
## the 2D array stores the connexions between the nodes.
## HERE the number of contacts is given for each connexion. 
#sub h_2_graph_wedge_eppic($$){
#
#  my ($H, $name) = @_ ;
#  my ($nb_nodes, $node, $node_2, @nodes, @labels, @connex, %node_index);
#
#  #
#  # Nodes names
#  #
#  @nodes = sort {$a cmp $b} (keys %$H) ;
#  $nb_nodes = $#nodes+1 ;
#
#  #
#  # Nodes Labels
#  #
#  my $i = 0;
#  foreach $node (@nodes){
#    $node_index{$node} = $i ;
#    @labels = (@labels, $i);
#    $i++;
#  }
#
#  #
#  # Connectivity matrix
#  #
#  $i = 0 ;
#  foreach $node (@nodes){
#
#    foreach $node_2 (@nodes){
#
#      if(exists $H->{$node}->{$node_2}){
#
#	$connex[$node_index{$node}][$node_index{$node_2}] = 10*(eppic_isbio($H->{$node}->{$node_2}->{"eppic"}));
#      } else {
#	$connex[$node_index{$node}][$node_index{$node_2}] = 0 ;#$H->{$node}->{$node_2}->{"size"};
#      }
#    }
#  }
#
#  #
#  # Puts everything in a Hash
#  #
#  my %result ;
#  $result{"labels"} = \@labels ;
#  $result{"connex"} = \@connex ;
#  $result{"real_labs"} = \@nodes ;
#  $result{"name"} = $name ;
#  return \%result ;
#}
#
## Take - 2 sorted array refs
## Return the dist between them and the numb of gaps
## which equals the difference of length
##
## Exple : [6,6,1,1,1]
##         [6,5,2,1]
##
## Have a distance of 
## abs(6-6)+abs(6-5)+abs(1-2)+abs(1-1) = 2
## and 1 gap
#sub array_dist($$){
#
#  my ($array_1, $array_2);
#  ($array_1, $array_2) = @_ ;
#
#  my ($i, $dist, $gap);
#  ($dist, $gap)=(0,0);
#
#  if($#{$array_1} > $#{$array_2}){
#
#    for $i (0..$#{$array_1}){
#      if(defined $array_2->[$i]){
#	$dist+=abs($array_1->[$i]-$array_2->[$i]);
#      } else {
#	$gap += $array_1->[$i] ;
#      }
#    }
#
#  } else {
#
#    for $i (0..$#{$array_2}){
#      if(defined $array_1->[$i]){
#	$dist+=abs($array_2->[$i]-$array_1->[$i]);
#      } else {
#	$gap += $array_2->[$i] ;
#      }
#    }
#  }
#  return ($dist, $gap);
#}
#
#
## Some domains are split across two chains. Those means that each cannot be
## considered as containing an entire domain.
##
## This script assesses, given a domain architecture of a chain, whether one can
## reasonably assume that the chain contains FULL domains.
##
##
#sub split_domain($){
#
#  my ($dom_arch) = @_;
#
#  my @doms = split( ";", $dom_arch);
#
#  my ($dom, %dom);
#
#  foreach $dom (@doms){
#
#    if($dom =~ /(\d+)_f/){
#
#      if(exists($dom{$1})){
#	delete $dom{$1};
#      } else {
#	$dom{$1} = 1;
#      }
#    }
#  }
#
#  # -- each dom in the form "dom_f" that appears an ODD number of times will be there.
#  my $n_split = keys %dom ;
#
#  return($n_split);
#}
#
## Takes - a ref to the big hash (loaded by the function above)
##       - a pdb_code
##       - a ref to a homology hash (loaded with load_homologies)
##       - a boolean keep_sf
## Exple: 2 identical chains (A,E), 2 identical (C,D) that are homologous to A & E,
## and 1 isolated (E) will give something of this form :
##
##
##   /->{"0"}->["AB", "CD"]
##  /
## H
##  \
##   \->{"1"}->["E"]
##
## This allows to get an idea of the COMPOSITION of a complex
##
## The superfamilies numbers are renumbered from 0 UNLESS $keep_sf=1
##
##
## I have to be carefull about the fact that this function changes the CPLX HASH!!!
##
#sub get_compo_accurate($$$){
#
#  my ($hash, $pdb_code, $homologies) = @_ ;
#  my ($chain, $chain_tmp, $group_num, $i);
#  my ($sf, @not_found, @to_find, @groups_id, @groups_ho, $identical);
#
#  my %sf_groups=();
#  my %sf_groups_2=();
#
#  # ---- 0/ I need to over-ride the SF information in order to consider homologies in an accurate fashion.
#  foreach my $chain1 (keys %{$hash->{$pdb_code}}){
#
#      foreach my $chain2 (keys %{$hash->{$pdb_code}}){
#
#	  ## IF both chains meet some identity criteria, SFNUM is transferred from the first
#	  if(  ($hash->{$pdb_code}->{$chain1}->{"prop"}->{"sf_num2"} ne $hash->{$pdb_code}->{$chain2}->{"prop"}->{"sf_num2"}) &&
#	       (
#		($hash->{$pdb_code}->{$chain1}->{"prop"}->{"sf_num3"} eq $hash->{$pdb_code}->{$chain2}->{"prop"}->{"sf_num3"} && $hash->{$pdb_code}->{$chain2}->{"prop"}->{"sf_num3"} ne "")
#		||
#		($hash->{$pdb_code}->{$chain1}->{"prop"}->{"sf_num"} eq $hash->{$pdb_code}->{$chain2}->{"prop"}->{"sf_num"}   && $hash->{$pdb_code}->{$chain2}->{"prop"}->{"sf_num"} ne "")
#		||
#		(exists $homologies->{$pdb_code.$chain1.$chain2} &&                  ### THIS IS BASED ON ATOMIC SEQUENCE, NOT SEQRES
#		 $homologies->{$pdb_code.$chain1.$chain2}->{"overlap"} > 80 &&
#		 $homologies->{$pdb_code.$chain1.$chain2}->{"id"} > 80 &&
#		 exists $hash->{$pdb_code}->{$chain1}->{"prop"}->{"sf_num2"}
#		)  
#	       ) ){
#	      $hash->{$pdb_code}->{$chain2}->{"prop"}->{"sf_num2"} = $hash->{$pdb_code}->{$chain1}->{"prop"}->{"sf_num2"} ;
#	  }
#	  
#	  ## IF both chains meet some identity criteria, but overlap (alignment length) is not sufficient (likely due to truncated sequence), SFNUM is changed to avoid looking for symmetries that may not exist
#	  if( exists $hash->{$pdb_code}->{$chain1}->{"prop"}->{"sf_num2"} &&
#	      ($hash->{$pdb_code}->{$chain1}->{"prop"}->{"sf_num2"} eq $hash->{$pdb_code}->{$chain2}->{"prop"}->{"sf_num2"}) &&
#	      exists $homologies->{$pdb_code.$chain1.$chain2} && 
#	      ($homologies->{$pdb_code.$chain1.$chain2}->{"overlap"} < 50 ||  ## 2018 re-activated that condition, not sure why it was removed.
#	       $homologies->{$pdb_code.$chain1.$chain2}->{"ali_len"} < 20 )
#	      ){
#	      $hash->{$pdb_code}->{$chain2}->{"prop"}->{"sf_num2"} = $chain2 ;
#	  }
#      }
#  }
#  
#  # ---- 1/ group the chains into superfamilies
#  foreach $chain (keys %{$hash->{$pdb_code}}){
#
#    if( exists $sf_groups{$hash->{$pdb_code}->{$chain}->{"prop"}->{"sf_num2"}}){
#      push(@{$sf_groups{ $hash->{$pdb_code}->{$chain}->{"prop"}->{"sf_num2"} }}, $chain);
#    } else {
#	$sf_groups{ $hash->{$pdb_code}->{$chain}->{"prop"}->{"sf_num2"} } = [$chain] ;
#    }
#  }
#
#  # ---- 2/ In each superfamily: look for identical chains
#  $i=0;
#  foreach $group_num (keys %sf_groups){
#
#    @to_find = @{$sf_groups{$group_num}} ;
#    @groups_id=(); # -- store prots that are identical as one string
#
#    while($#to_find >= 0 ){
#
#      @not_found=();
#
#      # -- if it's the last one then it has no identical partner
#      if($#to_find == 0){
#	push(@groups_id,shift(@to_find));
#
#      } else {
#
#	# -- else look for it
#	$chain = shift(@to_find);
#	$identical="";
#
#	foreach $chain_tmp (@to_find){
##	    my $tmp = are_identical($pdb_code.$chain, $pdb_code.$chain_tmp, $identities) ;
##	    if($tmp == 1){
#
#	    if(exists $homologies->{$pdb_code.$chain.$chain_tmp} &&
#	       $homologies->{$pdb_code.$chain.$chain_tmp}->{"overlap"}> 80 && 
#	       $homologies->{$pdb_code.$chain.$chain_tmp}->{"id"} > 90 
#		){
#		$identical .= $chain_tmp ;
#	    } else {
#		push(@not_found,$chain_tmp );
#	    }
#	}
#
#	# -- if no ident chain found then put it in the "no identical chain" group
#	if($identical eq ""){
#	    push(@groups_id,$chain);
#	} else {
#	    # -- else look put it in the identical group with its partners
#	    $identical .= $chain ;
#	    push(@groups_id, $identical);
#	}
#	# -- continue on chains that have not been found previously
#	@to_find = @not_found ;
#      }
#    }
#
#    @{$sf_groups_2{$i}} = @groups_id ;
#    $i+=1;
# }
#  return \%sf_groups_2;
#}
#
#
## Takes - a ref to the big hash (loaded by the function above)
##       - a pdb_code
##       - a ref to a homology hash (loaded with load_homologies)
##       - a boolean keep_sf
## Exple: 2 identical chains (A,E), 2 identical (C,D) that are homologous to A & E,
## and 1 isolated (E) will give something of this form :
##
##
##   /->{"0"}->["AB", "CD"]
##  /
## H
##  \
##   \->{"1"}->["E"]
##
## This allows to get an idea of the COMPOSITION of a complex
##
## The superfamilies numbers are renumbered from 0 UNLESS $keep_sf=1
##
#sub get_compo($$$$){
#
#  my ($hash, $pdb_code, $homologies,$keep_sf) = @_ ;
#  my ($chain, $chain_tmp, $group_num, $i);
#  my ($sf, @not_found, @to_find, @groups_id, @groups_ho, $identical);
#
#  my %sf_groups=();
#  my %sf_groups_2=();
#
#  # ---- 1/ group the chains into superfamilies
#  foreach $chain (keys %{$hash->{$pdb_code}}){
#
#    if( exists $sf_groups{$hash->{$pdb_code}->{$chain}->{"prop"}->{"sf_num3"}}){
#      push(@{$sf_groups{ $hash->{$pdb_code}->{$chain}->{"prop"}->{"sf_num3"} }}, $chain);
#    } else {
#      $sf_groups{ $hash->{$pdb_code}->{$chain}->{"prop"}->{"sf_num3"} } = [$chain] ;
#    }
#  }
#
#  # ---- 2/ In each superfamily: look for identical chains
#  $i=0;
#  foreach $group_num (keys %sf_groups){
#
#    @to_find = @{$sf_groups{$group_num}} ;
#    @groups_id=(); # -- store prots that are identical as one string
#
#    while($#to_find >= 0 ){
#
#      @not_found=();
#
#      # -- if it's the last one then it has no identical partner
#      if($#to_find == 0){
#	push(@groups_id,shift(@to_find));
#
#      } else {
#
#	# -- else look for it
#	$chain = shift(@to_find);
#	$identical="";
#
#	foreach $chain_tmp (@to_find){
#	  my $tmp = are_identical($pdb_code.$chain, $pdb_code.$chain_tmp, $homologies) ;
#	  if($tmp == 1){
#	    $identical .= $chain_tmp ;
#	  } else {
#	    push(@not_found,$chain_tmp );
#	  }
#	}
#
#	# -- if no ident chain found then put it in the "no identical chain" group
#	if($identical eq ""){
#	  push(@groups_id,$chain);
#	} else {
#	  # -- else look put it in the identical group with its partners
#	  $identical .= $chain ;
#	  push(@groups_id, $identical);
#	}
#	# -- continue on chains that have not been found previously
#	@to_find = @not_found ;
#      }
#    }
#
#    if($keep_sf == 1){
#      @{$sf_groups_2{$group_num}} = @groups_id ;
#    } else {
#      @{$sf_groups_2{$i}} = @groups_id ;
#    }
#    $i+=1;
# }
#  return \%sf_groups_2;
#}
#
## Takes - a ref to the big hash (loaded by the function above)
##       - a pdb_code
##       - a ref to a homology hash (loaded with load_homologies)
##       - a boolean keep_sf
## Exple: 2 identical chains (A,E), 2 identical (C,D) that are homologous to A & E,
## and 1 isolated (E) will give something of this form :
##
##
##   /->{"0"}->["AB", "CD"]
##  /
## H
##  \
##   \->{"1"}->["E"]
##
## This allows to get an idea of the COMPOSITION of a complex
##
## The superfamilies numbers are renumbered from 0 UNLESS $keep_sf=1
##
#sub get_compo_pfam($$$$){
#
#  my ($hash, $pdb_code, $homologies,$keep_sf) = @_ ;
#  my ($chain, $chain_tmp, $group_num, $i);
#  my ($sf, @not_found, @to_find, @groups_id, @groups_ho, $identical);
#
#  my %sf_groups=();
#  my %sf_groups_2=();
#
#  # ---- 1/ group the chains into superfamilies
#  foreach $chain (keys %{$hash->{$pdb_code}}){
#
#    if( exists $sf_groups{$hash->{$pdb_code}->{$chain}->{"prop"}->{"sf_num2"}}){
#      push(@{$sf_groups{ $hash->{$pdb_code}->{$chain}->{"prop"}->{"sf_num2"} }}, $chain);
#    } else {
#      $sf_groups{ $hash->{$pdb_code}->{$chain}->{"prop"}->{"sf_num2"} } = [$chain] ;
#    }
#  }
#
#  # ---- 2/ In each superfamily: look for identical chains
#  $i=0;
#  foreach $group_num (keys %sf_groups){
#
#    @to_find = @{$sf_groups{$group_num}} ;
#    @groups_id=(); # -- store prots that are identical as one string
#
#    while($#to_find >= 0 ){
#
#      @not_found=();
#
#      # -- if it's the last one then it has no identical partner
#      if($#to_find == 0){
#	push(@groups_id,shift(@to_find));
#
#      } else {
#
#	# -- else look for it
#	$chain = shift(@to_find);
#	$identical="";
#
#	foreach $chain_tmp (@to_find){
#	  my $tmp = are_identical($pdb_code.$chain, $pdb_code.$chain_tmp, $homologies) ;
#	  if($tmp == 1){
#	    $identical .= $chain_tmp ;
#	  } else {
#	    push(@not_found,$chain_tmp );
#	  }
#	}
#
#	# -- if no ident chain found then put it in the "no identical chain" group
#	if($identical eq ""){
#	  push(@groups_id,$chain);
#	} else {
#	  # -- else look put it in the identical group with its partners
#	  $identical .= $chain ;
#	  push(@groups_id, $identical);
#	}
#	# -- continue on chains that have not been found previously
#	@to_find = @not_found ;
#      }
#    }
#
#    if($keep_sf == 1){
#      @{$sf_groups_2{$group_num}} = @groups_id ;
#    } else {
#      @{$sf_groups_2{$i}} = @groups_id ;
#    }
#    $i+=1;
# }
#  return \%sf_groups_2;
#}
#
#
#sub compare_compo_strict($$){
#
#  my ($compo_1, $compo_2)=@_;
#  my ($group_num, $group_num_2, $group_1, $group_2);
#  my ($dist, $gap);
#  my (@fams_1, @fams_2);
#
#  # ---- Looks wether the families are strictky identical.
#  @fams_1 = keys %{$compo_1};
#  @fams_2 = keys %{$compo_2};
#
#  my $common = inter_exclu_union(\@fams_1, \@fams_2) ;
#
#  if( $#fams_1 == $#fams_2 && $#{$common->{"intersect"}} == $#fams_1 ){
#
#    # ---- for each group of compo_1, looks if the same group exists in compo_2
#    foreach $group_num (@fams_1){
#      $group_1 = group_2_num($compo_1->{$group_num});
#      $group_2 = group_2_num($compo_2->{$group_num});
#      ($dist, $gap) = array_dist($group_1, $group_2) ;
#
#      if( $dist != 0 || $gap !=0){
#	return 0;
#      }
#    }
#  } else {
#    return 0 ;
#  }
#  return 1;
#}
#
#sub compare_compo_relaxed($$){
#
#  my ($compo_1, $compo_2)=@_;
#  my ($group_num, $group_num_2, $group_1, $group_2);
#  my ($dist, $gap);
#  my (@fams_1, @fams_2);
#
#  # ---- Looks wether the families are strictky identical.
#  @fams_1 = keys %{$compo_1};
#  @fams_2 = keys %{$compo_2};
#
#  my $common = inter_exclu_union(\@fams_1, \@fams_2) ;
#
#  if( $#fams_1 == $#fams_2 && $#{$common->{"intersect"}} == $#fams_1 ){
#    return 1;
#  } else {
#    return 0 ;
#  }
#}
#
#
## converts a group of chains into a group of digits
## corresponding to the number of chains
#sub group_2_num($){
#
#  my ($array) = @_;
#  my (@new_array, $elmt) ;
#
#  foreach $elmt (@$array){
#    push(@new_array, length($elmt));
#  }
#  @new_array = sort {$b <=> $a} @new_array ;
#  return \@new_array ;
#}
#
## Take   - 1 hash of a composition
## Return the number of proteins in it
#sub sum_compo($){
#
#  my ($hash) = @_ ;
#  my ($gr_num, @gr_size, $elmt, $size_tot) ;
#  my %new_hash ;
#  $size_tot=0;
#
#  foreach $gr_num (keys %{$hash}){
#    @gr_size = @{$hash->{$gr_num}} ;
#    foreach $elmt (@gr_size){
#      $size_tot+=length($elmt);
#    }
#  }
#  return $size_tot ;
#}
#
##  Collection of CRAP COMPLEX: 6ins, 1hlt, 1tq0, 1id5, 1xmn, 1tbr
##
##
##
## Load all the PDBs in a Hash
## Take the array REF to their name as argument
## - save interface = 0 / 1, if 1, it saves the PDB lines to
##   be able to print the interface later on to visualize it
#
## Note that I don't handle that case: 
##
##    ****   NORMALLY I NOW HANDLE IT WITH THE $ALTERN SCALAR. ****
##
##ATOM   1203  O   ALA 2  29      12.441  28.720  58.969  1.00 22.13           O
##ATOM   1204  CB  ALA 2  29       9.875  27.325  60.216  1.00 22.49           C
##ATOM   1205  N  AVAL 2  30      13.014  27.746  60.934  0.50 22.00           N
##ATOM   1206  CA AVAL 2  30      14.424  27.471  60.607  0.50 22.16           C
##ATOM   1207  C  AVAL 2  30      14.812  26.043  60.988  0.50 22.77           C
##ATOM   1208  O  AVAL 2  30      14.537  25.602  62.112  0.50 23.21           O
##ATOM   1209  CB AVAL 2  30      15.349  28.484  61.337  0.50 21.67           C
##ATOM   1210  CG1AVAL 2  30      16.826  28.132  61.131  0.50 21.59           C
##ATOM   1211  CG2AVAL 2  30      15.058  29.931  61.003  0.50 21.71           C
##ATOM   1212  N  BMET 2  30      12.846  27.477  60.708  0.50 20.65           N
##ATOM   1213  CA BMET 2  30      14.125  27.065  60.112  0.50 20.32           C
##ATOM   1214  C  BMET 2  30      14.615  25.745  60.722  0.50 20.96           C
##ATOM   1215  O  BMET 2  30      14.405  25.479  61.908  0.50 20.32           O
##ATOM   1216  CB BMET 2  30      15.157  28.159  60.219  0.50 22.63           C
##ATOM   1217  CG BMET 2  30      15.815  28.373  61.513  0.50 24.42           C
##ATOM   1218  SD BMET 2  30      16.689  30.047  61.346  0.50 25.92           S
##ATOM   1219  CE BMET 2  30      15.204  31.060  61.332  0.50 27.96           C
##ATOM   1220  N   HIS 2  31      15.500  25.353  60.112  1.00 22.87           N
##ATOM   4814  CG  GLN B 289      25.829  58.731  26.834  1.00  5.21
##ATOM      7  CA  TYR A   5      15.490  39.455  26.171  1.00 12.03
#
##
## And the atoms CG1 and CG2 are stored with a MET instead of a VAL
##
sub load_PDBs($$){

    my ($name_list, $save_interface) = @_;
    my %pdb ;
    my $counter=0 ;
    my $smallest_chain = 22;
    my ($name);

    #print STDERR "Loading the PDBs from " ;
    #print STDERR "\"".$name_list->[0]."\" to \"".$name_list->[-1]."\" ...";

    foreach $name (@{$name_list}){

	my ($ch_name, $r_name, $a_name, $r_num, $a_num, $x, $y, $z, $ALTERN) ;
	my ($line, $ch_size) ;
	my $ch_name_tmp = "";
	my $r_name_tmp  = "";
	my $r_num_tmp   = "";
	my $ALTERN_TMP  = "";
	my $to_skip     = -1;
	
	$name =~ s/ //g;

	if( $name =~ m{/}){
	    open (PDB_FH, $name) or print STDERR "$name could not be opened in load_PDBs - $!\n";
	} else {
	    #print STDERR "Be carreful, PDB is loaded from default directory $::BU_all_renum\n";
	    open (PDB_FH, $::BU_RES_RENUM.$name) or print STDERR "$::BU_RES_RENUM $name could not be opened in load_PDBs - $!\n";
	}

	# ---- print the progression as points
	if ($counter%100 == 0 && $counter !=0){
	    #print STDERR "\n";
	}
	if( $counter%10 == 0 && $counter !=0){
	    #print STDERR "."
	}
	$counter+=1 ;
	# -----------------------------------
	#    print STDERR $name."\n";  # -- debug

	$ch_size = 0 ;

	## This is the index used to reprint the PDB later on.
	my $line_index = 0;

	#print "The PDB = ".Dumper(\%pdb);

	# ---- goes through the PDB file and load the atoms & positions in a hash
	while(<PDB_FH>){
	    
	    $line = $_;
	                   ##            altern resid is undef OR is in A-form (we skip all B-forms)
	    if($line=~/^ATOM/ && (substr($line,16,1) eq " " ||  substr($line,16,1) eq "A")){


		$ch_name= substr($line,21,1);	
		
		if($ch_name_tmp eq ""){
		    $ch_name_tmp=$ch_name ;
		}

		#print "The chain is equal to $ch_name, and ch_tmp = $ch_name_tmp --- size = $ch_size\n";

		if($ch_name ne $ch_name_tmp){ # -- then a new chain just arrived

		    # ---- this is to look if the chain is a protein or something else
		    # ---- if it is something else (no Calpha), we remove the entry
		    if($ch_size < $smallest_chain){
			#print "The chain is equal to $ch_name, size = $ch_size -- We delete the chain\n";
			delete $pdb{$name}->{$ch_name_tmp};
			$ch_size = 0 ;
		    } else {
			$pdb{$name}->{$ch_name_tmp}->{"size"}=$ch_size ;
			$ch_size = 0 ;
			
			### 
			### I need to check that the order of chain in the original PDB is increasing - otherwise there will be bugs.
			###
			if ( $::chain_number{$ch_name_tmp} > $::chain_number{$ch_name} ){
			    print "There is a problem in $name, $ch_name_tmp ($::chain_number{$ch_name_tmp}) is after $ch_name ($::chain_number{$ch_name}) but chains should always appear in the following order in the PDB file: A->Z, a->z, 0->9\n";
			}
		    }
		    $ch_name_tmp=$ch_name ;
		}

		$a_num  = substr($line,7,4);  # --- Changes at each line
		$a_name = substr($line,11,5);
		$x      = substr($line,30,8);
		$y      = substr($line,38,8);
		$z      = substr($line,46,8); # This has changed!!! BE CAREFULLLLLL AGAIN CHANGED FROM 7 TO 8
		
		$r_name = substr($line,17,4); # --- Changes at each residue / chain
		#$ALTERN = substr($line,16,1); # --- If there is an ambiguity, ALTERN can take A,B,C ... values.
		$r_num  = substr($line,22,5);

		# --- Remove white spaces
		$a_num  =~ s/ //g ;
		$a_name =~ s/ //g ;
		$r_name =~ s/ //g ;
		$ch_name=~ s/ //g ;
		$r_num  =~ s/ //g ;
		$x      =~ s/ //g ;
		$y      =~ s/ //g ;
		$z      =~ s/ //g ;

		# --- if wanted, add the line to
		# --- print the interface
		if( $save_interface == 1){
		    if($a_name eq "CA"){
			$pdb{$name}->{$ch_name}->{"lines"}->{$r_num}->{"CA"} = $line ;			
		    }
		    $pdb{$name}->{$ch_name}->{"lines"}->{$r_num}->{$line_index++} = $line ;
		}


		# ATOM   3944  N   TYR H 100A      2.196  19.583 107.279  1.00 18.94           N
		# ---- Sometimes, people are really dumb ... and put the A,B,C ALTERN stuff AFTER the residue number.
		# ---- --> So In that case I put it back in the right format.
		if($r_num =~ /^(-{0,1}[0-9]+)([a-zA-Z])$/){
		    $r_num  = $1.$2;     # This is to keep all residues distinct
		    #$ALTERN = $2;
		} elsif(! ($r_num =~ /^-{0,1}[0-9]+$/) ){
		    print STDERR $r_num." has not been recognized in pdb $name\n";
		}

		#if($ALTERN eq " "){
		#    $ALTERN_TMP = "OK";
		#}

		if( ($r_num_tmp ne ""  && $r_num_tmp == $r_num && $r_name_tmp ne $r_name) ||    # if resid def and same resid number but different name
		    $to_skip == $r_num  ##||
		    ##(exists($pdb{$name}->{$ch_name}->{$r_num}) && $ALTERN ne " " && $ALTERN ne "A" && $ALTERN_TMP ne "")        
                    # This last condition is to skip residues with alternative defs but some
		    # prots like 2g0x are entirely alternatively def ... so I hate PDBs - so I add the conditions 
		    # that 
		    # -1/ a TRUE RESIDUE (with no ALTER DEF) must have been seen in order to make the condition true.	   	    
		    # -2/ ALTERN must be different than A --> so B/C/D etc are skipped, but A is kept.
		    #
		    # Actually, a simpler solution is to only consider to OK lines to start with
		    ){
		    # Do nothing, because this is a second def for the same resid
		    #print " $r_num_tmp == $r_num && $r_name_tmp ne $r_name ".$line;
		    $to_skip=$r_num;
		    #$ALTERN_TMP = $ALTERN;

		} else {
		    
		    $to_skip=-1;
		    ######################################
		    if($a_name =~ "CA"){          ### --- Special treatment for the Ca
			### --- we index them
			$ch_size+=1;                ### --- chain size
			######################################
			# --- if first time init the list
			if( ! exists($pdb{$name}->{$ch_name}->{"x"})){
			    $pdb{$name}->{$ch_name}->{"x"}=[[$x,$r_num]];
			    $pdb{$name}->{$ch_name}->{"y"}=[[$y,$r_num]];
			    $pdb{$name}->{$ch_name}->{"z"}=[[$z,$r_num]];
			} else {                    # --- else add residues to it
			    push(@{$pdb{$name}->{$ch_name}->{"x"}},[$x,$r_num]);
			    push(@{$pdb{$name}->{$ch_name}->{"y"}},[$y,$r_num]);
			    push(@{$pdb{$name}->{$ch_name}->{"z"}},[$z,$r_num]);
			}
		    }

		    # --- add all atoms
		    $pdb{$name}->{$ch_name}->{$r_num}->{$a_name}->{"x"}=$x;
		    $pdb{$name}->{$ch_name}->{$r_num}->{$a_name}->{"y"}=$y;
		    $pdb{$name}->{$ch_name}->{$r_num}->{$a_name}->{"z"}=$z;

		    # --- and the name of the corresponding residue
		    $pdb{$name}->{$ch_name}->{$r_num}->{"Rname"}=$r_name ;   ### XXX must be carrefull when looking at KEYS of r_num!!!

		}

		$r_num_tmp = $r_num;
		$r_name_tmp = $r_name;
	    }
	}

	# ---- this is to look if the chain is a protein or something else
	# ---- if it is something else (no Calpha), we remove the entry
	if($ch_size < $smallest_chain){
	    delete $pdb{$name}->{$ch_name};                # XXX I think it should be ch_name not tmp
	    $ch_size = 0 ;
	} else {
	    $pdb{$name}->{$ch_name}->{"size"}=$ch_size ;   # XXX I think it should be ch_name not tmp
	    $ch_size = 0 ;
	}
    }

    # --- Sort the x, y, z coordinates of each residue
    my ($pdb1, $chain, @nbChains);

    #print "The PDBs are == ".Dumper(\%pdb)."\n";

    foreach $pdb1 (keys %pdb){

	@nbChains = keys %{$pdb{$pdb1}} ; # If it is not at leat a dimer

	# if($#nbChains >=1){

	foreach $chain (@nbChains){
	    @{$pdb{$pdb1}->{$chain}->{"x"}} = sort {$a->[0]<=>$b->[0]} @{$pdb{$pdb1}->{$chain}->{"x"}} ;
	    @{$pdb{$pdb1}->{$chain}->{"y"}} = sort {$a->[0]<=>$b->[0]} @{$pdb{$pdb1}->{$chain}->{"y"}} ;
	    @{$pdb{$pdb1}->{$chain}->{"z"}} = sort {$a->[0]<=>$b->[0]} @{$pdb{$pdb1}->{$chain}->{"z"}} ;
	}
	#    } else {
	#      delete $pdb{$pdb1};
	#    }
    }
    return \%pdb ;
}

####
#### Takes a path to a PDB file and returns a HASH giving correspondances between chain names and chain order.
####
#sub map_chain_name_num($){
#
#    my ($name) = @_;
# 
#    open (PDB_FH, $name) or die($name);
#
#    my ($ch_name, $line);
#    my $ch_name_tmp = "";
#    my $ch_number = 0;
#    my %H_ch_num = ();
#
#    while(<PDB_FH>){
#	$line = $_;
#	if ($line=~/^ATOM/){
#
#	    $ch_name= substr($line,21,1);	
#	    
#	    if($ch_name_tmp eq ""){
#		$ch_name_tmp=$ch_name ;
#
#		if(! exists $H_ch_num{$ch_name}){
#		    $H_ch_num{$ch_name}=$ch_number++;
#		}
#	    }
#	    if($ch_name ne $ch_name_tmp){ # -- then a new chain just arrived
#		if(! exists $H_ch_num{$ch_name}){
#		    $H_ch_num{$ch_name}=$ch_number++;
#		}
#	    }
#	}
#    }
#    return(\%H_ch_num);
#}
#
####
#### Takes:
#### 1. a hash from Load_PDBs (must contain the full path of the loaded PDBs);
#### 2. a hash giving the new chain correspondances under each full path as a key.
####
####
####   H BDBH:$VAR1 = {
####          '/tmp/1ohz/1qm8_1_1ohz.pdb' => {
####                                           '1' => 'A',
####                                           '3' => 'C',
####                                           '2' => 'B'
####                                         }
#sub write_PDB_kpax($$){
#
#    my ($PDB_H, $chain_def) = @_;
#
#     foreach my $PDB (keys %$PDB_H){
#
#	 #print STDERR "$PDB\n";
#
#	if(! ($PDB =~ /query/) ){
#
#	    my @chain_num = sort {$a <=> $b} keys %{$chain_def->{$PDB}} ;
#
#	    my ($REF, $code_tmp, $pdb_tmp);
#	    ##             REF     CODE
#	    if($PDB =~ /\/([^\/]+)\/([^\/]+.pdb)$/){
#		$REF      = $1;
#		$pdb_tmp  = $2;
#		$code_tmp = $2;
#		
#		## to get the name of the temporary code, I delete the constant part corresponding to the ref name
#		$code_tmp =~ s/_$REF.pdb//;
#	    } else {
#		print STDERR "CODE format not recognized!!\n";
#		exit(0);
#	    }
#	    
#	    #print STDERR "REF = $REF, pdb_tmp = $pdb_tmp , code_tmp = $code_tmp \n";
#
#	    if ( ! -e $ENV{KPAX_RESULTS}."$REF/pdbs"){
#		mkdir($ENV{KPAX_RESULTS}."$REF/pdbs")
#	    }
#
#	    my $PDB2 = $ENV{KPAX_RESULTS}."$REF/pdbs/$code_tmp.pdb";
#	    
#	    open(PDBOUT, ">$PDB2") or die ("cannot open $PDB2 : $!\n");
#	    my $all_pdb = "";
#	    foreach my $num (@chain_num){
#		
#		## Writes the file.		
#		# $pdb{$name}->{$ch_name}->{"lines"}->{$r_num}->{$line_index} = $line ;
#		
#		my $current_chain = $chain_def->{$PDB}->{$num};
#
#		my @residues_index = sort {$a <=> $b} keys %{$PDB_H->{$PDB}->{$current_chain}->{"lines"}};
#
#		foreach my $resid (@residues_index){
#
#		    my @line_index = sort {$a <=> $b} keys %{$PDB_H->{$PDB}->{$current_chain}->{"lines"}->{$resid}};
#
#		    foreach my $line_i (@line_index){
#
#			if($line_i ne "CA"){
#			    $all_pdb .= $PDB_H->{$PDB}->{$current_chain}->{"lines"}->{$resid}->{$line_i};
#			}
#			delete $PDB_H->{$PDB}->{$current_chain}->{"lines"}->{$resid}->{$line_i};
#		    }
#		}
#	    }
#	    
#	    print PDBOUT $all_pdb;
#	    close(PDBOUT);
#	}
#    }
#}
#
#
#
####
#### Return a type of molecule according to residues names.
####
#sub typeof(){
#
#  my @aas   = ("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "ASP", "ASN", "LYS", "GLU", "GLN", "ARG", "HIS", "PHE", "CYS", "TRP", "TYR", "MET", "PRO", "UNK");
#  my @dna   = ("DA", "DG", "DT", "DC");
#  my @rna   = ("A","G","C","U");
#
#  my %typeof=();
#
#  map { $typeof{$_}="aa" } (@aas);
#  map { $typeof{$_}="dna" } (@dna);
#  map { $typeof{$_}="rna" } (@rna);
#  return (\%typeof);
#}
#
#
#
#
####
#### SUB FUNCTION TO LOAD ONLY A SINGLE PDB FILE INTO A HASH
####
#sub load_PDB_ASA($$){
#
#  my ($file, $mytypeof) = @_;
#  my %pdb ;
#  my ($ch_name, $r_name, $a_name, $r_num, $a_num, $type, %lookup) ;
#  my ($line, $ch_size) ;
#  my $smallest_chain = 22;
#
#  open (PDB_FH, "<$file") or die("couldn't open $file");
#
#  $ch_size = 0 ;
#
#  # --- NAME: LIGAND (HETERO) - RNA (A/C/G/U/(t)) - DNA (DA/DG/DT/DC) - PROTEIN (Amino Acid+ Chain)
#
#  # ---- goes through the PDB file and load the atoms & positions in a hash
#  while(<PDB_FH>){
#    $line = $_;
#
#    if ($line=~/^ATOM/){
#
#      $ch_name= substr($line,21,1);	
#      $r_name = substr($line,17,4); # --- Changes at each residue / chain
#      $a_name = substr($line,11,5);
#      $r_num  = substr($line,22,5);
#
#      $a_name =~ s/ //g ;
#      $r_name =~ s/ //g ;
#      $r_num =~ s/ //g ;
#
#      if(exists( $mytypeof->{$r_name} )){
#	$type = $mytypeof->{$r_name};
#      } else {
#	$type = "other";
#      }
#
#      if( !exists($lookup{$ch_name}->{$r_num}) ){
#
#	if(exists $pdb{$ch_name}->{ $type } ){
#	  $pdb{$ch_name}->{ $type }++ ;
#	} else {
#	  $pdb{$ch_name}->{ $type } =1;
#	}
#	$lookup{$ch_name}->{$r_num}=1;
#      }
#    }
#  }
#
#  my %PDBfinal;
#  my $mytype;
#
#  foreach $ch_name (keys %pdb){
#
#    my @type = keys %{$pdb{$ch_name}};
#
#    my $besttype=0;
#
#    foreach $type (@type){
#
#      if( $pdb{$ch_name}->{$type} > $besttype){
#	$mytype = $type;
#	$besttype = $pdb{$ch_name}->{$type};
#      }
#    }
#
#    if( $mytype eq "aa" && $pdb{$ch_name}->{$mytype} <= $smallest_chain){
#
#      $PDBfinal{"pep"}->{$ch_name}=1;
#
#    } else {
#
#      $PDBfinal{$mytype}->{$ch_name}=1;
#    }
#  }
#
#  return \%PDBfinal ;
#}
#
#
#
## -- Extract a sub-complex from a complex
##
#sub sub_cplx($$$){
#
#  my ($complex, $chains, $out) = @_;
#
#  #print STDERR "$complex, CH = ".Dumper($chains)."\n";
#
#  my @tmp = ($complex);
#  my $cplx = load_PDBs(\@tmp, 1);
#
#  open(OUT_SUB, ">$out") or die("cannot open $out\n");
#
#  foreach my $ch (@$chains){
#
#    #print STDERR "CH = ".$ch."\n";
#
#    foreach my $rnum (keys %{$cplx->{$complex}->{$ch}}){
#
#      # print "Rnum = $rnum\n";
#
#      foreach my $aname (keys %{$cplx->{$complex}->{$ch}->{"lines"}->{$rnum}}){
#
#	print OUT_SUB $cplx->{$complex}->{$ch}->{"lines"}->{$rnum}->{$aname};
#
#      }
#    }
#  }
#}
#
#
## There was a bug in load_PDB so here we get the name of the structures that have actually been
## affected by the bug.
#sub pdb_bad($$$){
#
#  my ($file, $path, $N) = @_ ;
#  my ($pdb, $size, $sym, $groups, @GR, @tmp, @groups);
#
#  my $counter=0;
#  my  $OK2=0;
#  open(IN_SYM, "<$file") or die "$file doesn't exists\n";   # For each structure
#  my %prob=();
#  while(<IN_SYM>){
#
#    # ---- print the prgression as points
#    if ($counter%500 == 0){
#      print STDERR "\n";
#    }
#    if( $counter%10 == 0){
#      print STDERR "."
#    }
#    $counter+=1 ;
#
#    my $line = $_ ;		#<IN_SYM>;
#    chomp($line);
#    ($pdb, $size, $sym, $groups) = split(/\t/, $line) ;
#
#    if($sym eq "S1"){                                       # IF OK
#
#      @GR = split(/ /, $groups);                            # Get the chains
#      @groups=();
#      map { @{$groups[$_]} = split(//,$GR[$_]) } (0..$#GR);
#
#      my ($each_group, $each_chain, $resid_i,$N_res, $pdbs, $chain_ref, @resids, $OK, @selected);
#      # Write into a file the C alpha coordinates of the selected residues
#
#      # Load the corresponding PDB
#      $pdbs = load_PDBs(["$pdb"],1);
#
#      # For each group of identical chains
#      # Look for N corresponding resudies (same number)
#      # and write their corresponding number.
#      foreach $each_group (@groups) {
#
#	$N_res = $N;
#	@selected=();
#	$chain_ref = $each_group->[0];
#	@resids = @{$pdbs->{$pdb}->{$chain_ref}->{"x"}} ;
#	@resids = sort {$a->[1] <=> $b->[1]} @resids ;
#	$resid_i=$OK=1;
#
#	while ($N_res > 0 && $resid_i < @resids ) {
#	  foreach $each_chain (@$each_group) {
#	    if ( !exists($pdbs->{$pdb}->{$each_chain}->{$resids[$resid_i][1]}->{"CA"})) {
#	      $OK=0 ;
#	    }
#	  }
#
#	  if ($OK==1) {
#	    push(@selected,$resids[$resid_i][1]);
#	    $N_res--;
#	  }
#	  $resid_i+=1;
#	}
#
#	foreach $each_chain (@$each_group) {
#
#	  foreach $resid_i (@selected) {
#	    if( ($pdbs->{$pdb}->{$each_chain}->{$resid_i}->{"CA"}->{"x"}<=-100) ||
#		($pdbs->{$pdb}->{$each_chain}->{$resid_i}->{"CA"}->{"y"}<=-100) ||
#		($pdbs->{$pdb}->{$each_chain}->{$resid_i}->{"CA"}->{"z"}<=-100) ){
#	      $prob{$pdb}=1;
#	    }
#	  }
#	}
#      }
#    }
#  }
#  print join("\n",(keys %prob));
#}
#
### Load the result of the fasta alignment
### i.e. the fasta.all file
###
#sub load_percent_identities($){
#
#  my ($file) = @_;
#  open(FILE_IN, "<$file") or die "$file does not exists!\n";
#
#  my ($pdb1, $pdb2, $ident, %IDs, @line);
#
#  while(<FILE_IN>){
#
#    ## 1b3oA   1jr1A   514     514     98.833  514     1-514:1-514     3285    3839.6  720.0   2.3e-207        NA
#    @line = split(/\t/, $_);
#    $pdb1 = $line[0];
#    $pdb2 = $line[1];
#    $ident= $line[4];
#    #print STDERR "$pdb1 $pdb2 $ident \n";
#    $IDs{$pdb1}->{$pdb2}=$ident;
#  }
#  return \%IDs;
#}
#
#
## Takes - a cplx hash
##       - a homologies hash
## Returns the class of the cplx according to its composition:
##    1   2   4 
## 1  I          = 1
## 2      P      = 2
## 3  I + P      = 3
## 4          D  = 4
## 5  I   +   D  = 5
## 6      P + D  = 6
## 7  I + P + D  = 7
##
## Where I, P & D mean that identical/paralogous/unrelated chains 
## exist in the complex
##
#sub get_cplx_type3($$$){
#
#  my ($hash, $pdb_code, $homologies) = (@_) ;
#  my ($ident_index, $homol_index, $diff_index)= (0,0,0);
#  my $compo = get_compo_pfam($hash, $pdb_code, $homologies,1);
#
#  # ------------- To see if identical / homologous / unrelated chains exists in the complex ------ #
#
#  # --- 1/ Different chains ?
#  my @compo = keys %$compo ;
#  if($#compo > 0){
#    $diff_index = 1 ;
#  }
#
#  foreach my $each_sf (@compo){
#  # --- 2/ Identical chains ?
#    map { if(length($_) > 1){ $ident_index = 1 } } @{$compo->{$each_sf}} ;
#
#  # --- 3/ Paralogous chains ?
#    if( $#{$compo->{$each_sf}} > 0){
#      $homol_index = 1 ;
#    }
#  }
#  return ($ident_index + 2*$homol_index + 4*$diff_index) ;
#}
#
#sub load_cluster($){
#
#  my ($file) = @_;
#
#  open(CLUST_IN, "<$file") or die ("$file does not exists!\n");
#
#  my (@line, @codes, $ref, %clust);
#
#  while(<CLUST_IN>){
#
#    chomp($_);
#
#    @line  = split(/\t/,$_);
#    $ref = $line[0];
#
#    @codes = split(/\s/,$line[1]);
#
#    map {$clust{$codes[0]}->{$_}=1  } (@codes);
#  }
#  return(\%clust);
#}
#
#
#sub load_cluster_by_num($){
#
#  my ($file) = @_;
#
#  open(CLUST_IN, "<$file") or die ("$file does not exists!\n");
#
#  my (@line, @codes, $ref, %clust);
#
#  while(<CLUST_IN>){
#
#    chomp($_);
#
#    @line  = split(/\t/,$_);
#    $ref = $line[0];
#
#    @codes = split(/\s/,$line[1]);
#
#    map {$clust{$_}=$ref  } (@codes);
#  }
#  return(\%clust);
#}
#
#
##sub SF_2_CPLX($){
#
##  my ($hash) = @_;
##  my $homologies = load_identities2();
#
##  my ($cplx, $sf, %SFS, $compo);
#
##  foreach $cplx (keys %$hash){
#
##    $compo = get_compo($hash, $cplx, $homologies, 1);
#
##    foreach $sf (keys %$compo){
##      $SFS{$sf}->{$cplx}=1;
##    }
##  }
##  return(\%SFS);
##}
#
#sub get_uniq_chains($){
#
#  my ($grps) = @_;
#
#  my (@uniq_chains, $sf, $g);
#
#  foreach $sf (keys %{$grps}){
#
#    foreach $g (@{$grps->{$sf}}){
#      push(@uniq_chains, substr($g,0,1) );
#    }
#  }
#  return(\@uniq_chains);
#}
#
## FORMAT OF A SCOP ENTRY:
##
## d1jqga	1jqg	A:4P-100P	d.58.3.1	71792	cl=53931,cf=54861,sf=54897,fa=54898,dm=54899,sp=75429,px=71792
##
## <code>     <code>  <chain>:<start>-<end>,<chain>:<start>-<end> ... <rest>
##
## Now <start> and <end> can display unexpected features:
##
## 1/ they can end with a letter --> if in the PDB file they actually end with a letter ...
##
## 2/ they can start with a dash --> if the structure comes from NMR, or for other reasons beyond my understanding.
##
## The apprach is thus to:
##
## ---- 1. split according to the comma
## ---- 2. for each domain fragment (if any) --> extract the chain, the start, the end.
## ---- 3. fill a hash of the form:
##
##
##      {H1} --> {pdb} --> {chain} --> {position} --> {id}  = domain_id ; an "_f" is appended at the end of domains_id if it is a fragment.
##                                               \-> {num} = internal number
## ---- 4. to have a complete definition of fragments, a second hash is filled
##
##      {H2} --> {num} = "full" or "num1;num2" etc ...
##
#
## ---- Takes a scop file
## ---- and load a hash of the following structure :
#sub load_scop_file($){
#
#  my ($scop_file) = @_;
#
#  open(LOG,">../log/load_scop_file.log") or die ("couldn't open ../log/load_scop_file.log\n") ;
#
#  my ($line, $pdb_code, $sf_num, @line_sc, @dom_defs, $dom, $chain, $start, $end, $fragment, @fragments);
#  my %H1=();
#  my %H2=();
#
#  my $counter=0;
#  my $dom_count = 1;
#
#  print STDERR "\nLoading the SCOP file \"$scop_file\" \nStarted ".(scalar localtime)."\n";
#  open(IN_SCOP, "<".$scop_file);
#
#  while(<IN_SCOP>){
#
#    # ---- print the progression as points
#    if ($counter%20000 == 0 && $counter !=0){
#      print STDERR " $counter domains parsed\n";
#    }
#    if( $counter%1000 == 0){
#      print STDERR "."
#    }
#    $counter+=1 ;
#
#    $line=$_;
#
#    if (! ($line =~ /^\#/)){
#
#      # ---- splits the line on spaces first -----------------------------------------------------
#      @line_sc = split(/\s+/,$line);
#
#      # ---- get the pdb code --------------------------------------------------------------------
#      $pdb_code = $line_sc[1] ;
#
#      # ---- gets the superfamily number ---------------------------------------------------------
#      if( $line_sc[5] =~ /.*,sf=(\d+),/){
#	$sf_num = $1;
#      } else {
#	print STDERR "PROBLEM : BAD PARSING OF THE SCOP FILE, IT IS NOT THE SUPERFAMILY THAT IS LOADED\n";
#      }   
#
#      print LOG "$pdb_code\t$line_sc[2]\t$sf_num\tstatus=";
#
#      if($sf_num ne "310607"){
#
#      if($line_sc[2] eq "-"){
#
#	$chain = "A";
#	$start = 0;
#	$end = 0;
#	$dom_count++;
#	$H1{$pdb_code}->{$chain}->{"1"}->{"id"} = $sf_num;
#      } else {
#
#	@dom_defs = split(/,/, $line_sc[2]);
#
#	$fragment = "";
#	@fragments = ();
#	if($#dom_defs >= 1){
#	  $fragment .= "_f";
#	}
#
#	my $frag_num = 0;
#
#	# print STDERR "\n\n ----------- DOM: $line_sc[2] \n\n";
#	foreach $dom (@dom_defs){
#
#	  if($dom =~ /(\w{1}):-{0,1}(\d{0,4})\w{0,1}-{0,1}(\d{0,4})\w{0,1}/){
#
#	    $frag_num++;
#
#	    $chain = $1;
#
#	    if(defined $2 && $2 ne ""){
#	      $start = $2;
#	    } else {
#	      $start = 0;
#	    }
#
#	    if(defined $3 && $3 ne ""){
#	      $end = $3;
#	    } else {
#	      $end = 0;
#	    }
#
#	    # print STDERR "$pdb_code - $chain - $start - $end - SF NUM = ".$sf_num.$fragment."\n";
#
#	    if( exists $H1{$pdb_code}->{$chain}->{$start}){
#	      print STDERR "PDB $pdb_code - CHAIN $chain - START $start - ALREADY DEFINED!\n";
#	    }
#
#	    if($fragment ne ""){
#	      $H1{$pdb_code}->{$chain}->{$start}->{"id"}  = $sf_num.$fragment.$frag_num ;
#	    } else {
#	      $H1{$pdb_code}->{$chain}->{$start}->{"id"}  = $sf_num.$fragment ;
#	    }
#	    $H1{$pdb_code}->{$chain}->{$start}->{"num"} = $dom_count ;
#	    push(@fragments, $dom_count);
#
#	    # print STDERR "PDB $line_sc[1] - DOM $dom_count - CHAIN $chain - START $start - END = $end \n";
#	  } else {
#	    print STDERR "NOT RECOGNIZED: $pdb_code $dom\n";
#	  }
#	  $dom_count++;
#	}
#      
#	foreach $dom (@fragments){
#	  $H2{$dom} = join(";", @fragments);
#	}
#      }
#      }
#    }
#  }
#  print STDERR "\nEnded ".(scalar localtime)."\n\n";
#  return(\%H1);
#}
#
#
#
## FORMAT OF A SCOP ENTRY:
##
## uid	ecod_domain_id	manual_rep	f_id	pdb	chain	pdb_range	seqid_range	arch_name	x_name	h_name	t_name	f_name	asm_status	ligand
## 000000011	e1htr.1	MANUAL_REP	1.1.1.1	1htr	.	P:1-43,B:1-329	P:1-43,B:1-329	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	EF00082,EF00710	NOT_DOMAIN_ASSEMBLY	NO_LIGANDS_4A
## 000020647	e1avf.1	AUTO_NONREP	1.1.1.1	1avf	.	Q:1-22,J:2-329	Q:1-22,J:2-329	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	EF00082,EF00710	NOT_DOMAIN_ASSEMBLY	NO_LIGANDS_4A
## 000020649	e1avf.2	AUTO_NONREP	1.1.1.1	1avf	.	P:1-21,A:2-329	P:1-21,A:2-329	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	EF00082,EF00710	NOT_DOMAIN_ASSEMBLY	NA
## 000020650	e3psgA1	AUTO_NONREP	1.1.1.1	3psg	A	A:1P-326	A:1-370	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	EF00082,EF00710	NOT_DOMAIN_ASSEMBLY	NO_LIGANDS_4A
##
## This compares to SCOP:
## d1jqga	1jqg	A:4P-100P	d.58.3.1	71792	cl=53931,cf=54861,sf=54897,fa=54898,dm=54899,sp=75429,px=71792
##
## <code>     <code>  <chain>:<start>-<end>,<chain>:<start>-<end> ... <rest>
##
## Now <start> and <end> can display unexpected features:
##
## 1/ they can end with a letter --> if in the PDB file they actually end with a letter ...
##
## 2/ they can start with a dash --> if the structure comes from NMR, or for other reasons beyond my understanding.
##
## The apprach is thus to:
##
## ---- 1. split according to the comma
## ---- 2. for each domain fragment (if any) --> extract the chain, the start, the end.
## ---- 3. fill a hash of the form:
##
##
##      {H1} --> {pdb} --> {chain} --> {position} --> {id}  = domain_id ; an "_f" is appended at the end of domains_id if it is a fragment.
##                                               \-> {num} = internal number
## ---- 4. to have a complete definition of fragments, a second hash is filled
##
##      {H2} --> {num} = "full" or "num1;num2" etc ...
##
#
## ---- Takes a scop file
## ---- and load a hash of the following structure :
#sub load_ecod_file($){
#
#    my ($ecod_file) = @_;
#
#    open(LOG,">../log/load_ecod_file_".$::TIME_OF_CURRENT_UPDATE.".log") or die ("couldn't open ../log/load_ecod_file_".$::TIME_OF_CURRENT_UPDATE.".log\n") ;
#    open(LOGSMALL,">../log/load_ecod_file_SMALL_".$::TIME_OF_CURRENT_UPDATE.".log") or die ("couldn't open ../log/load_ecod_file_SMALL_".$::TIME_OF_CURRENT_UPDATE.".log\n") ;
#
#    my ($line, $pdb_code, $sf_num, @line_sc, @dom_defs, $dom, $chain, $start, $end, $fragment, @fragments);
#    my %H1=();
#    my %H2=();
#
#    my $counter=0;
#    my $dom_count = 1;
#
#    print STDERR "\nLoading the ECOD file \"$ecod_file\" \nStarted ".(scalar localtime)."\n";
#    open(IN_SCOP, "<".$ecod_file);
#
#    my $NB_SMALL_DOMAINS = 0;
#
#    while(<IN_SCOP>){
#
#	# ---- print the progression as points
#	if ($counter%20000 == 0 && $counter !=0){
#	    print STDERR " $counter domains parsed\n";
#	}
#	if( $counter%1000 == 0){
#	    print STDERR "."
#	}
#	$counter+=1 ;
#
#	$line=$_;
#
#	if (! ($line =~ /^\#/)){
#
#	    ## 000020649	e1avf.2	AUTO_NONREP	1.1.1.1	1avf	.	P:1-21,A:2-329	P:1-21,A:2-329	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	EF00082,EF00710	NOT_DOMAIN_ASSEMBLY	NA
#	    ## 000020650	e3psgA1	AUTO_NONREP	1.1.1.1	3psg	A	A:1P-326	A:1-370	beta barrels	"cradle loop barrel"	"RIFT-related"	"acid protease"	EF00082,EF00710	NOT_DOMAIN_ASSEMBLY	NO_LIGANDS_4A
#
#	    ##001609779       e4v8sAR7        AUTO_NONREP     1.1.2.27        4v8s    AR      AR:716-752,AR:876-998   AR:716-752,AR:876-998   beta barrels    "cradle loop barrel"    "RIFT-related"  "double psi"    EF19303 NOT_DOMAIN_ASSEMBLY     NO_L
#            ## d1jqga	1jqg	A:4P-100P	d.58.3.1	71792	cl=53931,cf=54861,sf=54897,fa=54898,dm=54899,sp=75429,px=71792
#
#	    # ---- splits the line on spaces first -----------------------------------------------------
#	    @line_sc = split(/\t/,$line);
#
#	    # ---- get the pdb code --------------------------------------------------------------------
#	    $pdb_code = $line_sc[4] ;
#
#	    # ---- gets the superfamily number ---------------------------------------------------------
#	    if( $line_sc[3] =~ /^([^\.]+\.[^\.]+)\.[^\.]+/){
#		$sf_num = $1;
#	    } else {
#		print STDERR "PROBLEM : BAD PARSING OF $line_sc[3] IN THE ECOD FILE IT IS NOT THE SUPERFAMILY THAT IS LOADED\n";
#	    }
#
#	    # print LOG "$pdb_code\t$line_sc[2]\t$sf_num\tstatus=";
#
#	    if(! ($line_sc[5] =~ /[A-Za-z0-9\.]+/) ){
#		
#		print STDERR "Problem with lines $line -- chain not recognized\n";
#
#	    } else {
#		
#		@dom_defs = split(/,/, $line_sc[7]);
#
#		my $frag_num = 0;
#		
#		## First I need to check if there are fragments or not, i.e., two pieces over 19 AA each
#		foreach $dom (@dom_defs){
#		    ##            CH $1          START $2         STOP $3
#		    if($dom =~ /(\w{1,2}):-{0,1}(\d{1,4})\w{0,1}--{0,1}(\d{1,4})\w{0,1}/){
#			$chain = $1;
#
#			if(defined $2 && $2 ne ""){
#			    $start = $2;
#			} else {
#			    $start = 0;
#			}
#
#			if(defined $3 && $3 ne ""){
#			    $end = $3;
#			} else {
#			    $end = 0;
#			}
#			if($end-$start > 10){
#			    $frag_num++;
#			}
#		    }
#		}
#
#		@fragments = ();
#		if($frag_num > 1){
#		    $fragment = "_f";
#		} else {
#		    $fragment = "";
#		}
#
#		$frag_num = 0;		
#		
#		# print STDERR "\n\n ----------- DOM: $line_sc[2] \n\n";
#		foreach $dom (@dom_defs){
#
#		    ##            CH $1          START $2         STOP $3
#		    if($dom =~ /(\w{1,2}):-{0,1}(\d{1,4})\w{0,1}--{0,1}(\d{1,4})\w{0,1}/){
#
#			$chain = $1;
#
#			if(defined $2 && $2 ne ""){
#			    $start = $2;
#			} else {
#			    $start = 0;
#			}
#
#			if(defined $3 && $3 ne ""){
#			    $end = $3;
#			} else {
#			    $end = 0;
#			}
#
#			### We consider the domain / or fragment only if it's 20 or more amino acids long
#			if($end-$start > 10){
#
#			    $frag_num++;
#
#			    # print STDERR "$pdb_code - $chain - $start - $end - SF NUM = ".$sf_num.$fragment."\n";
#			    
#			    if( exists $H1{$pdb_code}->{$chain}->{$start}){
#				print LOG "PDB $pdb_code - CHAIN $chain - START $start - ALREADY DEFINED!\n";
#			    }
#			    
#			    if($fragment ne ""){
#				$H1{$pdb_code}->{$chain}->{$start}->{"id"}  = $sf_num.$fragment.$frag_num ;
#			    } else {
#				$H1{$pdb_code}->{$chain}->{$start}->{"id"}  = $sf_num.$fragment ;
#			    }
#			    $H1{$pdb_code}->{$chain}->{$start}->{"num"} = $dom_count ;
#			    push(@fragments, $dom_count);
#
#			} else {
#			    print LOG "Domain $dom was ignored because too short $line\n";
#			    print LOGSMALL $line;
#			    $NB_SMALL_DOMAINS++;
#			}
#
#			# print STDERR "PDB $line_sc[1] - DOM $dom_count - CHAIN $chain - START $start - END = $end \n";
#		    } else {
#			print STDERR "NOT RECOGNIZED: $pdb_code $dom\n";
#		    }
#		    $dom_count++;
#		}
#
#		foreach $dom (@fragments){
#		    $H2{$dom} = join(";", @fragments);
#		}
#	    }
#	}
#    }
#    print STDERR "There are $NB_SMALL_DOMAINS small domains/segments that have been ignored\n";
#    print STDERR "\nEnded ".(scalar localtime)."\n\n";
#    return(\%H1);
#}
#
#
## - Takes a H from load_scop_file
##
## - Returns a H that contains the dom_arch:
##
## - H->{code.chain}->dom_arch
##
#sub get_dom_archS_ECOD($){
#
#    open(LOG_DOMARCH,">../log/load_ecod_file_domarch_$::TIME_OF_CURRENT_UPDATE.log") or die ("couldn't open ../log/load_ecod_file_domarch_$::TIME_OF_CURRENT_UPDATE.log\n") ;
#
#    my ($H) = @_;
#
#    my ($code, $chain, $start, %DOMS, $dom_arch, @doms);
#
#    foreach $code (keys %$H){
#
#	foreach $chain (keys %{$H->{$code}}){
#
#	    @doms = ();
#
#	    foreach $start (keys %{$H->{$code}->{$chain}}){
#		push(@doms, [ $H->{$code}->{$chain}->{$start}->{"id"}, $start]);
#	    }
#
#	    $dom_arch = get_dom_arch(\@doms);
#	    
#	    my @dom_simplified = @{$dom_arch};
#
#	    my $dom_simplified = simplify_doms(\@dom_simplified);
#	    
#	    $dom_arch = join(";", @$dom_arch);
#	    $dom_simplified = join(";", @$dom_simplified);
#	    
#	    if($dom_arch ne $dom_simplified){
#		print LOG_DOMARCH "DIFF - $code - $chain - $dom_arch ==> $dom_simplified\n";
#	    } else {
#		print LOG_DOMARCH "SAME - $code - $chain - $dom_arch ==> $dom_simplified\n";
#	    }
#
#	    $DOMS{$code.$chain} = $dom_simplified;
#	}
#    }    
#    return(\%DOMS);
#}
#
### Takes a string of domain architure and simplifies it
### E.g. if two consecutive domains appear as fragments, it merges them into one.
#sub simplify_doms($){
#
#    my ($doms) = @_;
#    my @doms_dup = @$doms;
#    #my @new_doms = ();
#    
#    my ($dom_current, $dom_next, $d1, $d2, $f1, $f2, $index, $d2_tmp, $f2_tmp, $K);
#    $d2 = "";
#    $f2 = "";
#
#    for( $index = 0; $index <= $#{$doms} ; $index++){
#
#	#print STDERR "Current dom = ".join(",",@$doms)."\n";
#	$dom_current = $doms->[$index];
#
#	if($dom_current =~ /^([^_]+)(_{0,1}f{0,1}[0-9]{0,2})$/){
#		
#	    $d1 = $1;
#	    $f1 = $2;
#
#	    my $new_fragment_num = $f1;
#	    $new_fragment_num =~ s/_f//;
#
#	    my $BROKEN=0;
#	    my $BROKEN_SAME=0;
#	    my $NEW_INSTANCE=0;
#	    my $remove_f1=0;
#	    ## only if there is a "_1" do we scan the rest of the domains to either remove /  renumber them.
#	    for( $K=$index+1 ; $K <= $#{$doms} ; $K++){
#
#		$dom_next = $doms->[$K];
#
#		if($dom_next =~ /^([^_]+)(_{0,1}f{0,1}[0-9]{0,2})$/){
#		
#		    $d2 = $1;
#		    $f2 = $2;
#
#		    if($d1 != $d2){
#			$BROKEN=1;
#		    }
#		    
#		    ## If we have a  full domain, even if it's the same (e.g., ccc in 'ccc_f1;ccc;ccc_f2' --> then we got a break.
#		    if($f2 eq ""){
#			$BROKEN=1;
#			$BROKEN_SAME=1;
#		    }
#
#		    ### This takes care of a special case like this:
#		    ### 4ud2 - S - 244.4_f1;244.4_f1;4025.1_f1;4046.1;4025.1_f2;244.4_f2;244.4_f2 ==> 244.4;244.4_f1;4025.1_f1;4046.1;4025.1_f2;244.4_f2;NA
#		    ###
#		    ### In that case I don't touch it.
#		    if($d1 eq $d2 && $f1 eq "_f1" && $f2 eq "_f1"){
#			return \@doms_dup;
#		    }
#
#		    #print STDERR "$BROKEN == $BROKEN_SAME == $NEW_INSTANCE == $d1 == $d2 \n";
#
#		    ## Another instance of the same domain before a break with another domain
#		    if($d1 eq $d2 && !$BROKEN && $f2 ne "" && $f1 ne ""){
#
#			my $new_fragment_num2 = $f2;
#			$new_fragment_num2 =~ s/_f//;
#
#			if($new_fragment_num2 == 1){
#			    $NEW_INSTANCE=1;
#			}
#
#			if(!$NEW_INSTANCE){
#			    $doms->[$K]="NA";
#			    $remove_f1=1;
#			}
#
#			##  Another instance of the same domain after a break with another domain
#		    } elsif($d1 eq $d2 && $BROKEN && !$NEW_INSTANCE && $f2 ne "" && $f1 ne "") {
#			
#			my $new_fragment_num2 = $f2;
#			$new_fragment_num2 =~ s/_f//;
#
#			if($new_fragment_num2 == 1){
#			    $NEW_INSTANCE=1;
#			} else {
#			    $new_fragment_num++;
#			    $doms->[$K]=$d1."_f".$new_fragment_num;			
#			    $remove_f1 = 0;
#			}
#		    }		    
#		} else {
#		    print STDERR "Domain dom_next $dom_current not recognized\n";
#		}		
#	    }
#
#	    ## Here I need to remove the _f1 only if there is not other matching domain left
#	    ## i.e., if 
#	    if($remove_f1){
#		## If we're on another f we dont remove it
#		if($doms->[$index] =~ /f1/){
#		    $doms->[$index] =~ s/_f[0-9]+//;
#		}
#	    }
#	
#	} elsif($dom_current eq "NA"){
#	    # DO nothing	    
#	}else {
#	    print STDERR "Domain dom_current $dom_current not recognized\n";
#	}
#    }
#
#    my @new_doms = ();
#    foreach my $each_dom (@$doms){
#	if($each_dom ne "NA"){
#	    push @new_doms, $each_dom;
#	}
#    }
#
#    return(\@new_doms);
#		  }
#    	      
#
## - Takes a H from load_scop_file
##
## - Returns a H that contains the dom_arch:
##
## - H->{code.chain}->dom_arch
##
#sub get_dom_archS($){
#
#  my ($H) = @_;
#
#  my ($code, $chain, $start, %DOMS, $dom_arch, @doms);
#
#  foreach $code (keys %$H){
#
#    foreach $chain (keys %{$H->{$code}}){
#
#      @doms = ();
#
#      foreach $start (keys %{$H->{$code}->{$chain}}){
#	push(@doms, [ $H->{$code}->{$chain}->{$start}->{"id"}, $start]);
#      }
#
#      $dom_arch = get_dom_arch(\@doms);
#
#      $dom_arch = join(";", @$dom_arch);
#
#      # print STDERR "$code - $chain - dom_arch = $dom_arch\n";
#
#      $DOMS{$code.$chain} = $dom_arch;
#
#    }
#  }
#  return(\%DOMS);
#}
#
##
## Takes a table of the form [ [id1, start1, end1], [id2, start2, end2] ]
##
## Returns a domain architecture.
#
#sub get_pfam_dom_arch($){
#
#  my ($array) = @_;
#
#  my (@pos, %Hstart, %Hend, $elmt, @dom_arch);
#
#  foreach $elmt (@$array){
#
#    push(@pos, $elmt->[1]);
#    push(@pos, $elmt->[2]);
#
#    $Hstart{$elmt->[1]}->{"end"} = $elmt->[2];
#    $Hstart{$elmt->[1]}->{"dom"} = $elmt->[0];
#
#    $Hend{$elmt->[2]} -> {"end"} = $elmt->[1];
#    $Hend{$elmt->[2]} -> {"dom"} = $elmt->[0];
#  }
#
#  @pos = sort {$a <=> $b } @pos ;
#
#  my $i=0;
#
#  foreach $i (0..$#pos-1){
#
#    # Then it is a complete domain
#    if(exists $Hstart{$pos[$i]} && exists $Hend{$pos[$i+1]}){
#
#      push @dom_arch, $Hstart{$pos[$i]}->{"dom"};
#
#      # Then there is only the begining of the domain
#    } elsif(exists $Hstart{$pos[$i]}){
#
#      push @dom_arch, $Hstart{$pos[$i]}->{"dom"}."_f1";
#
#      # Then there is only the end of the domain
#    } elsif(exists $Hend{$pos[$i+1]}){
#
#      push @dom_arch, $Hend{$pos[$i+1]}->{"dom"}."_f2";
#
#      # If it's not the begining or the end, it should be the end of a first fragment, or the begining of a last fragment
#    } else {
#      # print STDERR "Problem $pos[$i]: neither full, nor begining, nor end!\n"; THAT'S ACTUALLY OK
#    }
#  }
#  return(\@dom_arch);
#}
#
#sub get_pfam_dom_archS($){
#
#  my ($file) = @_;
#
#  open(IN, "<$file") or die(" Could not open $file \n");
#
#  my (%H, %H2, @line, $pdb, $chain, $dom, $start, $end, $dom_arch);
#
#
#  while(<IN>){
#
#    chomp($_);
#    @line = split(/\s+/, $_);
#
#    if($#line != 14) {# && !($#line == 16 && $line[16] eq "(nested)")){
#
#	#print STDERR "Line: --> @line <-- is strange!! - length = $#line\n";
#
#    } else {
#
#      if(length($line[0]) <= 3){
#	print STDERR "There is a problem with $line[0] ! The name is too short\n";
#      }
#
#      ($pdb, $chain, $dom, $start, $end) = (substr($line[0], 0, length($line[0])-1), substr($line[0], -1, 1), $line[5], $line[1], $line[2] );
#
#      #print STDERR "@line\n$pdb - $chain - $dom - $start - $end\n";
#      if(exists $H{$pdb}->{$chain}){
#
#	push @{$H{$pdb}->{$chain}}, [$dom, $start, $end] ;
#
#      } else {
#
#	$H{$pdb}->{$chain} = [] ;
#	push @{$H{$pdb}->{$chain}}, [$dom, $start, $end] ;
#      }
#    }
#  }
#
#  #print Dumper($H{"1b01"});
#
#  foreach $pdb (keys %H){
#
#    foreach $chain (keys %{$H{$pdb}}){
#
#      #print STDERR "$pdb $chain\n";
#      $dom_arch = get_pfam_dom_arch($H{$pdb}->{$chain});
#      $H2{$pdb.$chain} = join(";",@$dom_arch);
#    }
#  }
#  return(\%H2);
#}
#
#
#sub get_pfam_dom_archS_UNIPROT($){
#
#    my ($file) = @_;
#
#    open(IN, "<$file") or die(" Could not open $file \n");
#
#    my (%H, %H2, @line, $pdb, $dom, $start, $end, $dom_arch);
#
#
#    while(<IN>){
#
#	chomp($_);
#	@line = split(/\s+/, $_);
#
#	if($#line != 14) {# && !($#line == 16 && $line[16] eq "(nested)")){
#
#	    print STDERR "Line: --> @line <-- is strange!! - length = $#line\n";
#
#	} else {
#
#	    if(length($line[0]) <= 3){
#		print STDERR "There is a problem with $line[0] ! The name is too short\n";
#	    }
#
#	    ($pdb, $dom, $start, $end) = ($line[0], $line[5], $line[1], $line[2] );
#
#	    ## PROBLEM: FORMAT NOT RECOGNIZED FOR PROTEIN tr|Q5SMH7|Q5SMH7_THET8
#	    
#	    if($pdb =~ /^[^\|]+\|([^\|]+)\|/){
#		$pdb = $1;
#	    } else {
#		print STDERR "PROBLEM: FORMAT NOT RECOGNIZED FOR PROTEIN $pdb\n";
#	    }
#	    
#	    if(exists $H{$pdb}){
#
#		push @{$H{$pdb}}, [$dom, $start, $end] ;
#
#	    } else {
#
#		$H{$pdb} = [] ;
#		push @{$H{$pdb}}, [$dom, $start, $end] ;
#	    }
#	}
#    }
#    #print Dumper($H{"1b01"});
#
#    foreach $pdb (keys %H){
#	$dom_arch = get_pfam_dom_arch($H{$pdb});
#	$H2{$pdb} = join(";",@$dom_arch);
#    }
#    return(\%H2);
#}
#
#
#sub get_pfam_dom_archS_ATOM($){
#
#  my ($file) = @_;
#
#  open(IN, "<$file") or die(" Could not open $file \n");
#
#  my (%H, %H2, @line, $pdb, $chain, $dom, $start, $end, $dom_arch);
#
#  while(<IN>){
#
#    chomp($_);
#    @line = split(/\s+/, $_);
#
#    if($#line != 14) {# && !($#line == 16 && $line[16] eq "(nested)")){
#
#	#print STDERR "Line: --> @line <-- is strange!! - length = $#line\n";
#
#    } else {
#
#      if(length($line[0]) <= 3){
#	  print STDERR "There is a problem with $line[0] ! The name is too short\n";
#      }
#
#      if($line[0] =~ /^([a-zA-Z0-9_]+)-([a-zA-Z0-9]+)$/){
#	  $pdb = $1;
#	  $chain = $2;
#      } else {
#	  print STDERR "Problem with $line[0], not regonized in REGEX\n";
#      }
#      ($dom, $start, $end) = ($line[5], $line[1], $line[2] );
#
#      #print STDERR "@line\n$pdb - $chain - $dom - $start - $end\n";
#      if(exists $H{$pdb}->{$chain}){
#
#	push @{$H{$pdb}->{$chain}}, [$dom, $start, $end] ;
#
#      } else {
#
#	$H{$pdb}->{$chain} = [] ;
#	push @{$H{$pdb}->{$chain}}, [$dom, $start, $end] ;
#      }
#    }
#  }
#
#  #print Dumper($H{"1b01"});
#
#  foreach $pdb (keys %H){
#
#    foreach $chain (keys %{$H{$pdb}}){
#
#      #print STDERR "$pdb $chain\n";
#      $dom_arch = get_pfam_dom_arch($H{$pdb}->{$chain});
#      $H2{$pdb.$chain} = join(";",@$dom_arch);
#    }
#  }
#  return(\%H2);
#}
#
#
#
## Take an array of the following form
## [ [SF1, p1], [SF2, p2], [SF3, p3] ....]
## where SF1/2/3 correspond to 3 superfam. found in the corresponding chain
## and p1/2/3 correspond to the first position where the domain appears in the chain.
## Return the SFs in the right order.
#sub get_dom_arch($){
#  my ($array) = @_;
#  my (@sorted, @dom_arch);
#  @sorted = sort {$a->[1] <=> $b->[1]} @$array ;
#  map {push(@dom_arch, $_->[0])} @sorted ;
#  return \@dom_arch;
#}
#
####
#### Functions used for the NEW ASA calculation.
####
#sub aa_areas(){
#  my %VALs = ("ALA" => 108.0,
#	      "ARG" => 238.7,
#	      "ASN" => 143.9,
#	      "ASP" => 140.4,
#	      "CYS" => 134.0,
#	      "GLN" => 178.4,
#	      "GLU" => 172.2,
#	      "GLY" => 80.1,
#	      "HIS" => 183.0,
#	      "ILE" => 175.1,
#	      "LEU" => 178.0,
#	      "LYS" => 200.9,
#	      "MET" => 194.0,
#	      "PHE" => 199.5,
#	      "PRO" => 136.2,
#	      "SER" => 116.5,
#	      "THR" => 139.2,
#	      "TRP" => 249.0,
#	      "TYR" => 212.7,
#	      "VAL" => 151.3
#	     );
#  return(\%VALs);
#}
#
#sub aa_cor(){
#  my %VALs = ("ALA" => "A",
#	      "ARG" => "R",
#	      "ASN" => "N",
#	      "ASP" => "D",
#	      "CYS" => "C",
#	      "GLN" => "Q",
#	      "GLU" => "E",
#	      "GLY" => "G",
#	      "HIS" => "H",
#	      "ILE" => "I",
#	      "LEU" => "L",
#	      "LYS" => "K",
#	      "MET" => "M",
#	      "PHE" => "F",
#	      "PRO" => "P",
#	      "SER" => "S",
#	      "THR" => "T",
#	      "TRP" => "W",
#	      "TYR" => "Y",
#	      "VAL" => "V"
#	     );
#  return(\%VALs);
#}
#
#sub load_HETATM($){
#
#    my ($file) = @_;
#  
#    open(ASA_IN, "<$file") or die ("File  $file does not exists\n");
#
#    my ($line, $ch_name, $r_num, $r_name, $i);
#    my $H = {};
#
#    while(<ASA_IN>){
#
#	$line = $_;
#
#	if ($line=~/^HETATM/){
#
#	    #ATOM    176  CB AVAL A  22      44.028  16.418   8.487  0.77 15.30           C  
#	    #ATOM   1211  CB  UNK D   4       9.456   5.533  55.562  1.00 28.20
#	    $ch_name= substr($line,21,1);	
#	    $r_name = substr($line,17,3); # --- Changes at each residue / chain
#	    $r_num  = substr($line,22,4);
#	    
#	    $r_name =~ s/ //g ;
#	    $r_num  =~ s/ //g ;
#	    
#	    $H->{$ch_name.$r_num.$r_name}=1;
#	}
#    }
#    close(ASA_IN);
#    return($H);
#}
#
#
#sub load_ASA($$$$){
#
#  my ($HASH, $file, $type, $aa_names, $HETATOMS) = @_;
#
#  if( -e $file){
#
#      open(ASA_IN, "<$file");
#  } else {
#      #or  ("File  $file does not exists\n");
#      return -1;
#  } 
#
#  my ($line, $ch_name, $r_num, $ASA, $r_name, $x, $y, $z, $a_num, $i, $a_name);
#
#  while(<ASA_IN>){
#
#    $line = $_;
#
#    if ($line=~/^ATOM/){
#
#      #ATOM    176  CB AVAL A  22      44.028  16.418   8.487  0.77 15.30           C    PDB FORMAT
#      #ATOM      2   CA VAL A   1       7.854  18.800   3.718  1.00   7.3           C   AREAIMOL FORMAT
#
#      $ch_name= substr($line,21,1);	
#      $r_name = substr($line,17,3); # --- Changes at each residue / chain
#      $r_num  = substr($line,22,4);
#      $ASA    = substr($line,60,6);
#      $a_num  = substr($line,6,5);
#      $a_name  = substr($line,14,2);
#
#      $x      = substr($line,30,8);
#      $y      = substr($line,38,8);
#      $z      = substr($line,46,8);
#
#      $x      =~ s/ //g ;
#      $y      =~ s/ //g ;
#      $z      =~ s/ //g ;
#      $r_name =~ s/ //g ;
#      $r_num  =~ s/ //g ;
#      $a_num  =~ s/ //g ;
#      $ASA    =~ s/ //g ;
#
#      $ASA =  sprintf "%.2f", $ASA;
#
#      ### We do not consider the changes for molecules other than Amino Acids.
#      ### Here I need the HET_indexes --> I need to get this from before.
#      if(! exists $HETATOMS->{$ch_name.$r_num.$r_name} ){ 
#
#	if ( $type eq "all") {
#
#	  if (exists $HASH->{$ch_name}->{$r_num}->{$type}) {
#
#	      ## ASA is added only if the name corresponds (because otherwise the ASA of HETATOM can
#	      ## be added to that of the residue on the same chain/resid number
#	      if($HASH->{$ch_name}->{$r_num}->{"name"} eq $r_name){
#		  $HASH->{$ch_name}->{$r_num}->{$type} += $ASA;
#	      }
#
#	  } else {
#
#	      $HASH->{$ch_name}->{$r_num}->{$type} = $ASA;
#	  }
#
#	  $HASH->{$ch_name}->{$r_num}->{"burried"}->{$type}->{$a_name} = $ASA;
#
#	  ## Usually HETATOMS ARE AFTER THE ATOM so here it gets defined only if it does not overwite an ATOM
#	  if( ! exists $HASH->{$ch_name}->{$r_num}->{"name"} ){
#	      $HASH->{$ch_name}->{$r_num}->{"name"} = $r_name;
#	  }
#
#	  #print STDERR "THE ATOM is $a_name /  in $line\n";
#
#	  ### OK NOW TAKES CA WHEN IT CAN
#	  if($a_name eq "CA"){
#	    $HASH->{$ch_name}->{$r_num}->{"posx"} = $x;
#	    $HASH->{$ch_name}->{$r_num}->{"posy"} = $y;
#	    $HASH->{$ch_name}->{$r_num}->{"posz"} = $z;
#	    $HASH->{$ch_name}->{$r_num}->{"CA"} = 1;
#	    #print STDERR "LoadASA enters in 'CA' \n";
#	  } elsif(! exists $HASH->{$ch_name}->{$r_num}->{"CA"}){
#	    $HASH->{$ch_name}->{$r_num}->{"posx"} = $x;
#	    $HASH->{$ch_name}->{$r_num}->{"posy"} = $y;
#	    $HASH->{$ch_name}->{$r_num}->{"posz"} = $z;
#	  }
#
#	} elsif ( ! exists $HASH->{$ch_name}->{$r_num} ) {
#	
#	  #die("Problem with ${file}, $r_num does not exist in the 'all' but does in the '${type}'\n");
#	
#	} elsif ($HASH->{$ch_name}->{$r_num}->{"name"} ne $r_name) {
#
#	    ## A problem here is that HETATOMS can be amino acids!! In that case it will get confused because the same residue number of supposedly the same chain 
#	    ## corresponds to two different r_name
#	    #die("Problem with ${file}, resnum $r_num (ch=$ch_name) ".($HASH->{$ch_name}->{$r_num}->{"name"})."' in the 'all' is not the same as in no'${type}' where it is '${r_name}'\n");
#
#	} else {
#
#	  if ($ASA > 2) {
#	    print STDERR "($file) Weird: ASA > 2 for type $type at resid $r_num \n";
#	  }
#	
#	  $ASA = abs($ASA);
#
#	  if (exists $HASH->{$ch_name}->{$r_num}->{$type}) {
#	    $HASH->{$ch_name}->{$r_num}->{$type} += $ASA;
#	  } else {
#	    $HASH->{$ch_name}->{$r_num}->{$type} = $ASA;
#	  }
#
#	  $HASH->{$ch_name}->{$r_num}->{"burried"}->{$type}->{$a_name} = $ASA;
#	}
#	
#      } else {
#
#	## IF the AA name does not exists and it is the HET stuff, I'll try to keep the largest HET EFFECT
#	##
#	## Now the problem is that some HETATM, e.g. 1eb3, are badly given!!! Only the coordinates change and all the rest is the same.
#	## ----> so in order to make the difference between different ligand, the only way to go is to check that their co-factors are close enough.
#	## ----> I CANNOT EVEN USE THE ATOM INFO BECAUSE AREAIMOL RENUMBERS THE HETATOMSSS!!!!!!!!!
#	##
#	## HETATM 2660  C1  DSB A1341      21.829  26.895  13.592  1.00 15.86           C
#
#	if( $type eq "oHET"){
#
#	  $ASA = abs($ASA);
#
#	  if(exists ( $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name} ) ){
#
#	    $a_num = 100;
#	    foreach $i (0 .. $#{$HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"x"}}){
#	      if( $a_num > distance($x, $y, $z,
#				    $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"x"}->[$i],
#				    $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"y"}->[$i],
#				    $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"z"}->[$i])){
#
#		$a_num = distance($x, $y, $z,
#				  $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"x"}->[$i],
#				  $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"y"}->[$i],
#				  $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"z"}->[$i]);
#	      }
#	    }
#
#	    if( $a_num > 7 ){
#	      ### DO NOTHING
#	    } else {
#
#	      $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"ASA"} += $ASA;
#	      push(@{$HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"x"}}, $x);
#	      push(@{$HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"y"}}, $y);
#	      push(@{$HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"z"}}, $z);
#	    }
#	  } else {
#	    $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"ASA"} = $ASA;
#	    $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"x"} = [$x];
#	    $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"y"} = [$y];
#	    $HASH->{"HET"}->{$r_name."-".$r_num."-".$ch_name}->{"z"} = [$z];
#	  }
#	}
#      }
#    }
#  }
##  return($HASH);
#  return 0;     ## Now the HASH is passed and modified directly into the function, so no need to return
#}

# Process contacts between pdbs chains
sub process_contacts($$$){

  ### my $CODE = ["file_path.pdb"];
  my ($to_proceed_file, $out, $Num) = @_ ;
  print "x: $to_proceed_file, y: $out, z: $Num\n";

  my $VDW = vdw();

  my $DO_GNM = 0;          # Calculate the contacts using the GNM view
  my $DO_INTER_ONLY = 0;   # Do not calculate the contacts WITHIN each chain
                           # ---- Some constants
  my $cut_long  = 20;      # ---- cut-off to look for a contact between 2 residues
  my $cut_GNM  = 8;        # ---- cut-off to count a contact between 2 residues (GNM view)
  my $smallest_chain=22;   # ---- size of the smallest chain to consider (smaller chains are discarded)
  my $save_interface=0;    # ---- If you want to have a look at the interface, setting this option to 1
                           # ---- will save it in a file.
  
  print("Beginning\n"); ## Hugo
  # ----- General

  # Output will be placed here
  # my $out = "../../results/03_contacts/chainsContacts_$num"; # ---- OUPUT FILE
  #my $out = "chainsContacts_$num"."_TMP"; # ---- OUTPUT FILE

  # --- Clears the files.
  if($DO_GNM){
    open(OUT_GNM, ">$out"."_res_GNM.txt");
  }
  open(OUT_FULL, ">$out"."_FULL.txt");
  #open(OUT_ATOM, ">$out"."_ATOM.txt");

  my ($pdb_name, $pdb_short, $res_n1, $res_n2, @chains, $ch_i, $ch_j, $x, $y, $z, $resid_num, $contact_list_big_rough, $contact_list_big, $contact_list_small, $contact, @pdb_array, $pdbs) ;
  my ($x_1, $y_1, $z_1, $atom_1, $atom_2, $x_2, $y_2, $z_2, $VDW1, $VDW2, $res_name, $res_name2, $cut_short, %pretty_print, %seen, @Ds, $d_min, $d_max, $d_avg);
  my $counter=0 ;
  my (@resid_I_inter_m1, @resid_J_inter_m1, @resid_I_inter_m2, @resid_J_inter_m2, %tmp);

  my ($B4, $AFTER) = format_stderr($Num);
  


  print("Before loop\n"); ## Hugo
  print("array pdb files: $to_proceed_file\n");

  my @to_proceed = ($to_proceed_file);
  print("to_proceed: @to_proceed\n");

  foreach $pdb_name (@to_proceed){
      print("Inside loop\n"); ## Hugo

      print "pdb name: $pdb_name\n";

      print STDERR $B4."Processing contacts for ".(sprintf "%02d", $Num)." .. ".(sprintf "%06d", $counter)."/$#{to_proceed}\r".$AFTER;

      #print STDERR "Processing: PDB = $pdb_name\n" ;

    %seen=();
    %pretty_print=();
    if($pdb_name =~ /([^\/]+)\.pdb$/){
      $pdb_short = $1;
    } else {
      $pdb_short = $pdb_name;
      $pdb_short =~ s/\.pdb$// ;
    }
    print("pdb_short: $pdb_short\n");
    #print STDERR "PDB SHORT = $pdb_short\n" ;

    ## SHOULD LOAD HERE !!!
    @pdb_array = ($pdb_name);
    $pdbs = load_PDBs(\@pdb_array, $save_interface) ;

      my $PDB_CONTACTS = "";

    # -----------------------------------
    #print STDERR $pdb_name."\n";  # -- debug

    # ----------------------------------------------------- INTERF BETWEEN 2 CHAINS ----------------------------------------------------------
    # -
    @chains = keys %{$pdbs -> {$pdb_name}} ;
    for $ch_i (0...$#{chains}){

      for $ch_j (($ch_i+$DO_INTER_ONLY)...$#{chains}){

	  #print STDERR $pdb_name." $ch_i $ch_j\n";  # -- debug

	my @table = @{$pdbs->{$pdb_name}->{$chains[$ch_i]}->{"x"}};
	@table = sort {$a->[1] <=> $b->[1]} @table ;
	
	foreach $resid_num (@table){

	  $x = $pdbs->{$pdb_name}->{$chains[$ch_i]}->{$resid_num->[1]}->{"CA"}->{"x"} ;
	  $y = $pdbs->{$pdb_name}->{$chains[$ch_i]}->{$resid_num->[1]}->{"CA"}->{"y"} ;
	  $z = $pdbs->{$pdb_name}->{$chains[$ch_i]}->{$resid_num->[1]}->{"CA"}->{"z"} ;

	  # ---------------------- FIRST ROUGHLY ----------------------
	  if($DO_GNM==1)
	    {
	      #
	      # ---- Here we look at residues which Calphas are closer than cut_GNM angstrom
	      $contact_list_big_rough = get_close($pdbs, $pdb_name, $x, $y, $z, $cut_GNM, $chains[$ch_j]);
	      # print STDERR "CONTACTS = ".join(" ",@$contact_list_big_rough)."\n";
	      foreach $contact (@{$contact_list_big_rough}){

		$res_n1 = $pdbs->{$pdb_name}->{$chains[$ch_i]}->{$resid_num->[1]}->{"Rname"};
		$res_n2 = $pdbs->{$pdb_name}->{$chains[$ch_i]}->{$contact}->{"Rname"};
		if(exists( $AAS{$res_n1} )) {$res_n1 = $AAS{$res_n1}; } else { $res_n1 = "X" ;}
		if(exists( $AAS{$res_n2} )) {$res_n2 = $AAS{$res_n2}; } else { $res_n2 = "X" ;}

		if( ($ch_i != $ch_j) || (abs($resid_num->[1]-$contact)>=3)){
		  # ---- code ch1 ch2 resnum1 resnum2 resname1 resname2 N dmin dmax davg
		  print OUT_GNM "pdb_short $chains[$ch_i] $chains[$ch_j] ".($resid_num->[1])." $contact $res_n1 $res_n2 8 8 8\n";
		}
	      }
	    }

	  #
	  # --------------------- THEN ACCURATELY ---------------------
	  #
	  # ---- Here we look at residues which Calphas are closer than cut_long angstrom
	  # ---- this is the major calculation                    ---##
	  $contact_list_big  = get_close($pdbs, $pdb_name, $x, $y, $z, $cut_long, $chains[$ch_j]);

	  foreach $contact (@{$contact_list_big}){

	    @Ds = ();

	    if( ($ch_i != $ch_j) || (abs($resid_num->[1]-$contact)>=3) ){

	      # Tries to calculate a contact only if the atomic gorup is defined --> get the residue name and check
	      $res_name  = $pdbs->{$pdb_name}->{$chains[$ch_i]}->{$resid_num->[1]}->{"Rname"} ;
	      $res_name2 = $pdbs->{$pdb_name}->{$chains[$ch_j]}->{$contact}->{"Rname"} ;
	      if(exists( $AAS{$res_name} )) {$res_n1 = $AAS{$res_name}; } else { $res_n1 = "X" ;}
	      if(exists( $AAS{$res_name2} )) {$res_n2 = $AAS{$res_name2}; } else { $res_n2 = "X" ;}

	      $VDW1 = 1.75; # Defaut value		
	      foreach $atom_1 (keys %{$pdbs->{$pdb_name}->{$chains[$ch_i]}->{$resid_num->[1]}}){

		if($atom_1 ne "Rname"){

		  $x_1 = $pdbs->{$pdb_name}->{$chains[$ch_i]}->{$resid_num->[1]}->{$atom_1}->{"x"} ;
		  $y_1 = $pdbs->{$pdb_name}->{$chains[$ch_i]}->{$resid_num->[1]}->{$atom_1}->{"y"} ;
		  $z_1 = $pdbs->{$pdb_name}->{$chains[$ch_i]}->{$resid_num->[1]}->{$atom_1}->{"z"} ;

		  if(exists $VDW->{$res_name}->{$atom_1}){
		    $VDW1 = $VDW->{$res_name}->{$atom_1} ;
		  }

		  $VDW2 = 1.75;
		  foreach $atom_2 (keys %{$pdbs->{$pdb_name}->{$chains[$ch_j]}->{$contact}}){

		    if($atom_2 ne "Rname"){

		      $x_2 = $pdbs->{$pdb_name}->{$chains[$ch_j]}->{$contact}->{$atom_2}->{"x"} ;
		      $y_2 = $pdbs->{$pdb_name}->{$chains[$ch_j]}->{$contact}->{$atom_2}->{"y"} ;
		      $z_2 = $pdbs->{$pdb_name}->{$chains[$ch_j]}->{$contact}->{$atom_2}->{"z"} ;

		      if(exists $VDW->{$res_name2}->{$atom_2}){
			$VDW2 = $VDW->{$res_name2}->{$atom_2} ;
		      }

		      # CUT OFF TO COUNT A CONTACT BETWEEN 2 ATOMS
		      $cut_short = $VDW1 + $VDW2 + 0.5 ;

		      # ---- if the distance between both atoms is smaller than cut_short, then print a contact line
		      my $d = distance($x_1, $y_1, $z_1, $x_2, $y_2, $z_2) ;

		      if ( $d < $cut_short){
			#print OUT_ATOM "$pdb_short $chains[$ch_i] $chains[$ch_j] ".($resid_num->[1])." $contact $res_n1 $res_n2 $atom_1 $atom_2 $d\n";
			push(@Ds, $d);
		      }
		    }
		  }
		}
	      }

	      if( $#Ds >= 0){
		$d_min = sprintf "%.3f", MIN(\@Ds);
		$d_max = sprintf "%.3f", MAX(\@Ds);
		$d_avg = sprintf "%.3f", MEAN(\@Ds);
		$PDB_CONTACTS .= "$pdb_short $chains[$ch_i] $chains[$ch_j] ".($resid_num->[1])." $contact $res_n1 $res_n2 ".($#Ds+1)." $d_min $d_max $d_avg\n";
	      }

	    } # if($ch_i =! $ch_j ...
	  }
	}
      }
    }

      print OUT_FULL $PDB_CONTACTS ;

      if($counter == $#{to_proceed}){
	  print STDERR $B4."Processing contacts for ".(sprintf "%02d", $Num)." .. ".(sprintf "%06d", $counter)."/$#{to_proceed} --- Finnished at ".localtime()."\b".$AFTER;
	  #print STDERR "Processing contacts .. ".(sprintf "%06d", $counter)."/$#{$to_proceed}\n";
      }
      $counter++;
  }
  
  #print STDERR "\nFinnished processing PDBs for interChain contacts...".."\n" ;
}

# Takes - a ref to the PDB Hash structure,
#       - a pdb name
#       - a point (x
#                  y
#                  z)
#       - a radius Rc
#       - a 2d pdb chain name to compare to CH2
# RETURNS : the residues numbers which are comprised in the bowl of
#           center (x,y,z) and radius rc
sub get_close($$$$$$$){

  my ($pdbs, $name, $x_1, $y_1, $z_1, $rc, $CH2) = @_ ;
  #print STDERR "GET close : $pdbs, $name, $x_1, $y_1, $z_1, $rc, $CH2 \n";
  my ($d, $ar, $elmt, $inf, $sup, $pos_inf, $pos_sup) ;
  my %x_ok=() ;
  my %y_ok=() ;
  my %z_ok=() ;

  # --- First look wether ...
  if ( $pdbs->{$name}->{$CH2}->{"x"}->[0][0] - $rc > $x_1 ||
       $pdbs->{$name}->{$CH2}->{"x"}->[-1][0] + $rc < $x_1 ||
       $pdbs->{$name}->{$CH2}->{"y"}->[0][0] - $rc > $y_1 ||
       $pdbs->{$name}->{$CH2}->{"y"}->[-1][0] + $rc < $y_1 ||
       $pdbs->{$name}->{$CH2}->{"z"}->[0][0] - $rc > $z_1 ||
       $pdbs->{$name}->{$CH2}->{"z"}->[-1][0] + $rc < $z_1
     ) {
    # --- if one of these conditions is true then $R1 can't be in contact with 
    # --- any residue of the chain 2
    my @empty ;
    return \@empty ;
  } else {
    # --- else, we look for the potential residues in contact (included in the bowl)

    ## --- look for wich residues have their X coord inside the bowl
    $inf = $x_1 - $rc ;
    $sup = $x_1 + $rc ;

    $pos_inf = dichoSearch($inf, $pdbs->{$name}->{$CH2}->{"x"});
    $pos_sup = dichoSearch($sup, $pdbs->{$name}->{$CH2}->{"x"});

    for $elmt ($pos_inf..$pos_sup){
      $ar = $pdbs->{$name}->{$CH2}->{"x"}->[$elmt] ; # Those "CA" that are OK on X
      $x_ok{$ar->[1]}= ($x_1 - $ar->[0])**(2);
    }

    ## --- look for wich residue having their x coord that is inside the bowl,
    ## --- have also their Y coord inside the bowl
    $inf = $y_1 - $rc ;
    $sup = $y_1 + $rc ;

    $pos_inf = dichoSearch($inf, $pdbs->{$name}->{$CH2}->{"y"});
    $pos_sup = dichoSearch($sup, $pdbs->{$name}->{$CH2}->{"y"});

    for $elmt ($pos_inf..$pos_sup){
      $ar = $pdbs->{$name}->{$CH2}->{"y"}->[$elmt] ;
      if( exists $x_ok{$ar->[1]} ){
	$y_ok{$ar->[1]}= $x_ok{$ar->[1]} + ($y_1-$ar->[0])**(2);
      }
    }

    ## -- same for Z
    $inf = $z_1 - $rc ;
    $sup = $z_1 + $rc ;

    $pos_inf = dichoSearch($inf, $pdbs->{$name}->{$CH2}->{"z"});
    $pos_sup = dichoSearch($sup, $pdbs->{$name}->{$CH2}->{"z"});

    for $elmt ($pos_inf..$pos_sup){
      $ar = $pdbs->{$name}->{$CH2}->{"z"}->[$elmt] ;
      if( exists $y_ok{$ar->[1]}){
	$d = sqrt( $y_ok{$ar->[1]} + ($z_1 - $ar->[0])**(2) ) ;
	if ( $d < $rc){
	  $z_ok{$ar->[1]}= $d;
	  # print STDERR "PDB = $name, ORIG = $x_1, $y_1, $z_1, Found=Resid number ".$ar->[1].", cut = $d\n";
	}
      }
    }

    my @good_res  = keys(%z_ok) ;
    return \@good_res ;
  }
}

# The diff with get_close is that this one looks at all chains together.
# Takes - a ref to the PDB Hash structure,
#       - a pdb name
#       - a point (x
#                  y
#                  z)
#       - a radius Rc
# RETURNS : the residues numbers which are comprised in the bowl of
#           center (x,y,z) and radius rc
sub get_close_all($$$$$){

  my ($pdbs, $x_1, $y_1, $z_1, $rc) = @_ ;
#  print STDERR "GET close : $pdb, $x_1, $y_1, $z_1, $rc, $CH2 \n";
  my ($d, $ar, $elmt, $inf, $sup, $pos_inf, $pos_sup) ;
  my %x_ok=() ;
  my %y_ok=() ;
  my %z_ok=() ;

  ## --- X -------------------------------------------------------------------
  ## --- look for wich residues have their X coord inside the bowl
  $inf = $x_1 - $rc ;
  $sup = $x_1 + $rc ;

  $pos_inf = dichoSearch($inf, $pdbs->{"all"}->{"x"});
  $pos_sup = dichoSearch($sup, $pdbs->{"all"}->{"x"});

  for $elmt ($pos_inf..$pos_sup){
    $ar = $pdbs->{"all"}->{"x"}->[$elmt] ; # Those "CA" that are OK on X
    $x_ok{$ar->[1]} = ($x_1 - $ar->[0])**(2);
  }

  ## --- Y -------------------------------------------------------------------
  ## --- look for wich residue having their x coord that is inside the bowl,
  ## --- have also their Y coord inside the bowl
  $inf = $y_1 - $rc ;
  $sup = $y_1 + $rc ;

  $pos_inf = dichoSearch($inf, $pdbs->{"all"}->{"y"});
  $pos_sup = dichoSearch($sup, $pdbs->{"all"}->{"y"});

  for $elmt ($pos_inf..$pos_sup){

    $ar = $pdbs->{"all"}->{"y"}->[$elmt] ;
    if( exists $x_ok{$ar->[1]} ){
      #print STDERR "YOYOMA\n";
      $y_ok{$ar->[1]} = $x_ok{$ar->[1]} + ($y_1-$ar->[0])**(2);
    }
  }

  ## --- Z -------------------------------------------------------------------
  ## -- same for Z
  $inf = $z_1 - $rc ;
  $sup = $z_1 + $rc ;

  $pos_inf = dichoSearch($inf, $pdbs->{"all"}->{"z"});
  $pos_sup = dichoSearch($sup, $pdbs->{"all"}->{"z"});

  for $elmt ($pos_inf..$pos_sup){
    $ar = $pdbs->{"all"}->{"z"}->[$elmt] ;
    if( exists $y_ok{$ar->[1]}){
      $d = sqrt( $y_ok{$ar->[1]} + ($z_1 - $ar->[0])**(2) ) ;
      # print STDERR "ORIG = $x_1, $y_1, $z_1, Found=Resid number ".$ar->[1].", cut = $d\n";
      if ( $d < $rc){
	$z_ok{$ar->[1]}= $d;
	# print STDERR "ORIG = $x_1, $y_1, $z_1, Found=Resid number ".$ar->[1].", cut = $d\n";
      }
    }
  }

  my @good_res  = keys(%z_ok) ;
  return \@good_res ;
}


# euclidian distance beetween 2 points
# x1, y1, z1, x2, y2, z2
sub distance($$$$$$){
  my ($x_1, $y_1, $z_1, $x_2, $y_2, $z_2) = @_ ;
  my $dist=sqrt(($x_1 - $x_2)**(2) + ($y_1-$y_2)**(2) + ($z_1-$z_2)**(2));
  return $dist;
}

#
# Takes a SORTED Array and an Integer as argument,
# Return the index of the smaller value in the array
# which is bigger or equal to the Integer.
#
sub dichoSearch($$){
  my ($toFind, $array) = @_ ;
  my ($oldIndex, $newIndex) =($#{$array},0);
  my ($temp,$direction)=(0,1) ;
  my $i=0;

  while( $i< $#{$array}+1){
    $i++;
    #print STDERR "INDEX : $newIndex, array : ".join(", ",@{$array})."-> array[0] = $array->[0]\n";
    #print STDERR "toFind ; $toFind\n";
    if($array->[$newIndex][0] == $toFind){
      return $newIndex ;
    }

    $temp = $newIndex ;
    $newIndex = int( ($newIndex+$oldIndex)/2 ) ;

    if($newIndex == $temp){

      if($#{$array} > 0){
	
	if( $array->[$newIndex+1][0] > $toFind){
	  return $newIndex;
	} else {
	  return $newIndex+1;
	}

      } else {
	
	if($toFind > $array->[$newIndex][0]){
	  return 1 ;
	} else {
	  return 0;
	}
      }
    }

    if($direction == 1){ # Then we come from the left

      if( $array->[$newIndex][0] > $toFind ){ # Then have to go to the left
	# So the direction changes to opposite
	$oldIndex = $temp ;
	$direction = -1 ;
      }
    }
    elsif($direction == -1){ # Then we come from the right

      if( $array->[$newIndex][0] < $toFind ){ # Then have to go to the right
	# So the direction changes
	$oldIndex = $temp ;
	$direction = 1 ;
      }
    }
  }
  print STDERR "ERROR IN DICHO SEARCH\n";
}

sub vdw(){

  my $radii ;
  my $aa ;

  # ALL AMINO ACIDS
  my @aas = ("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "ASP", "ASN", "LYS", "GLU", "GLN", "ARG", "HIS", "PHE", "CYS", "TRP", "TYR", "MET", "PRO");

  # THESE ARE THE VALUES USED FOR THE DIFFERENT AMINO ACIDS
  foreach $aa (@aas){
    $radii->{$aa}->{"N"}=1.64 ;
    $radii->{$aa}->{"CA"}=1.88 ;
    $radii->{$aa}->{"C"}=1.64 ;
    $radii->{$aa}->{"O"}=1.42 ;
    $radii->{$aa}->{"CB"}=1.88;
  }

  $aa= "ALA";

  $aa= "VAL";
  $radii->{$aa}->{"CG1"}=1.88;
  $radii->{$aa}->{"CG2"}=1.88;

  $aa= "LEU";
  $radii->{$aa}->{"CG"}=1.88;
  $radii->{$aa}->{"CD1"}=1.88;
  $radii->{$aa}->{"CD2"}=1.88;

  $aa= "ILE";
  $radii->{$aa}->{"CG1"}=1.88;
  $radii->{$aa}->{"CG2"}=1.88;
  $radii->{$aa}->{"CD1"}=1.88;

  $aa = "SER";
  $radii->{$aa}->{"OG"}=1.46;

  $aa = "THR";
  $radii->{$aa}->{"CG2"}=1.88;
  $radii->{$aa}->{"OG1"}=1.46;

  $aa = "ASP";
  $radii->{$aa}->{"CG"}=1.61;
  $radii->{$aa}->{"OD1"}=1.42;
  $radii->{$aa}->{"OD2"}=1.42;

  $aa = "ASN";
  $radii->{$aa}->{"CG"}=1.61;
  $radii->{$aa}->{"OD1"}=1.42;
  $radii->{$aa}->{"ND2"}=1.64;

  $aa = "LYS";
  $radii->{$aa}->{"CG"}=1.88;
  $radii->{$aa}->{"CD"}=1.88;
  $radii->{$aa}->{"CE"}=1.88;
  $radii->{$aa}->{"NZ"}=1.64;

  $aa = "GLU";
  $radii->{$aa}->{"CG"}=1.88;
  $radii->{$aa}->{"CD"}=1.61;
  $radii->{$aa}->{"OE1"}=1.42;
  $radii->{$aa}->{"OE2"}=1.42;

  $aa = "GLN";
  $radii->{$aa}->{"CG"}=1.88;
  $radii->{$aa}->{"CD"}=1.61;
  $radii->{$aa}->{"OE1"}=1.42;
  $radii->{$aa}->{"NE2"}=1.64;

  $aa = "ARG";
  $radii->{$aa}->{"CG"}=1.88;
  $radii->{$aa}->{"CD"}=1.88;
  $radii->{$aa}->{"NE"}=1.64;
  $radii->{$aa}->{"CZ"}=1.61;
  $radii->{$aa}->{"NH1"}=1.64;
  $radii->{$aa}->{"NH2"}=1.64;

  $aa = "HIS";
  $radii->{$aa}->{"CG"}=1.61;
  $radii->{$aa}->{"CD2"}=1.76;
  $radii->{$aa}->{"ND1"}=1.64;
  $radii->{$aa}->{"NE2"}=1.64;
  $radii->{$aa}->{"CE1"}=1.76;

  $aa = "PHE";
  $radii->{$aa}->{"CG"}=1.61;
  $radii->{$aa}->{"CD1"}=1.76;
  $radii->{$aa}->{"CD2"}=1.76;
  $radii->{$aa}->{"CE1"}=1.76;
  $radii->{$aa}->{"CE2"}=1.76;
  $radii->{$aa}->{"CZ"}=1.76;

  $aa = "CYS";
  $radii->{$aa}->{"SG"}=1.77;
  $radii->{$aa}->{"SH"}=1.77; #?

  $aa = "TRP";
  $radii->{$aa}->{"CG"}=1.61;
  $radii->{$aa}->{"CD1"}=1.76;
  $radii->{$aa}->{"CD2"}=1.61;
  $radii->{$aa}->{"NE1"}=1.64;
  $radii->{$aa}->{"CE2"}=1.88;
  $radii->{$aa}->{"CE3"}=1.76;
  $radii->{$aa}->{"CZ2"}=1.76;
  $radii->{$aa}->{"CZ3"}=1.76;
  $radii->{$aa}->{"CH2"}=1.76;

  $aa = "TYR";
  $radii->{$aa}->{"CG"}=1.61;
  $radii->{$aa}->{"CD1"}=1.76;
  $radii->{$aa}->{"CD2"}=1.76;
  $radii->{$aa}->{"CE1"}=1.76;
  $radii->{$aa}->{"CE2"}=1.76;
  $radii->{$aa}->{"CZ"}=1.61;
  $radii->{$aa}->{"OH"}=1.46;

  $aa = "MET";
  $radii->{$aa}->{"CG"}=1.88;
  $radii->{$aa}->{"SD"}=1.77;
  $radii->{$aa}->{"CE"}=1.88;

  $aa = "PRO";
  $radii->{$aa}->{"CG"}=1.88;
  $radii->{$aa}->{"CD"}=1.88;

  return $radii;
}



# euclidian distance beetween 2 points
# x1, y1, z1, x2, y2, z2
#sub distance($$$$$$){
#  my ($x_1, $y_1, $z_1, $x_2, $y_2, $z_2) = @_ ;
#  my $dist=sqrt(($x_1 - $x_2)**(2) + ($y_1-$y_2)**(2) + ($z_1-$z_2)**(2));
#  return $dist;
#}

######
######  NACCESS - ASA related function
######
sub load_all_feed_DB($$$$){

  my ($list, $dbh, $path_all, $path_chain) = @_;
  my ($pdb, $h, $h2, $h3, $h4, $chain, $resnum, $resname, $asa_1, $asa_2, $cat_1, $cat_2, @data, $rel_1, $rel_2, $ratio, $chains, $count, $x, $y, $z);

  my ($sth_one, %present);
  my $present = oneColQuery("SELECT distinct code FROM residue", $dbh);
  map { $present{$_}=1 } @{$present} ; # 2013
  ###while (my $row = $sth_one->fetchrow_arrayref) { $present{substr($row->[0],0,4)}=1;} $sth_one->finish; OLD

  # -- Prepare query
  my $sql = "INSERT INTO residue (code, chain, resnum, resname, pos_x, pos_y, pos_z, letter, asa_in_BU, asa_alone, delta_asa, asa_rel_in_BU, asa_rel_alone, ratio_rel, category_complex, category_alone) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
  my $sth_two = $dbh->prepare($sql);

  if($sth_two){
    # Ok
  } else {
    die("INSERT QUERY WAS BAD  $!");
  }

  $count=1;

  foreach $pdb (@$list){

    if(!exists($present{$pdb})){

	print LOG $pdb."\n";
	
	if($count % 10 == 0){
	    print STDERR "Filling the residue table ... ".(sprintf "%06d", $count)."/".(sprintf "%06d", $#{$list}+1)."\r";	
	}
	$count++;
	
	print LOG "$pdb \n";

	# print STDERR "Come here $pdb \n";
	($h, $chains) = load_one_rsa_BU($path_all.$pdb);
	$h2 = load_one_rsa_split($path_chain.$pdb, $chains);
	$h3 = load_x_y_z($path_all.$pdb);
	$h4 = merge($h, $h2, $h3);

      #print "H=".Data::Dumper::Dumper($h);
      #print "\n----------------------------------------------------------\nH2=".Data::Dumper::Dumper($h2);
      #print "\n----------------------------------------------------------\nH4=".Data::Dumper::Dumper($h4);

      # H  -> {chain} --> {resnum} -> { resname / asa / category }
      #  \           \--> {total}  -> float
      #   \-> {tot}   \-> float

      foreach $chain (keys %{$h4->{"all"}}){

	if($chain ne "tot" && $chain ne " "){

	  foreach $resnum (keys %{$h4->{"all"}->{$chain}}){

	    if($resnum ne "tot"){

	      $resname = $h4->{"all"}->{$chain}->{$resnum}->{"name"};

	      $asa_1   = $h4-> {"all"} ->{$chain}->{$resnum}->{"asa"};
	      $asa_2   = $h4->{"split"}->{$chain}->{$resnum}->{"asa"};

	      $cat_1   = $h4-> {"all"} ->{$chain}->{$resnum}->{"category"};
	      $cat_2   = $h4->{"split"}->{$chain}->{$resnum}->{"category"};

	      $rel_1   = $h4-> {"all"} ->{$chain}->{$resnum}->{"asa_rel"};
	      $rel_2   = $h4->{"split"}->{$chain}->{$resnum}->{"asa_rel"};

	      if($rel_2 == 0){

		if($rel_1 !=0){
		  print LOG "Problem with $pdb, $chain, $resnum, ASA_REL in BU is $rel_1 and it is $rel_2 in the split\n";
		}
		$ratio=0;
	      } else {
		$ratio = $rel_1 / $rel_2 ;
	      }

	      if(exists($h4->{"xyz"}->{$chain}->{$resnum}->{"x"})){
		$x = $h4->{"xyz"}->{$chain}->{$resnum}->{"x"};
	      } else {
		$x=0;
	      }
	      if(exists($h4->{"xyz"}->{$chain}->{$resnum}->{"y"})){
		$y = $h4->{"xyz"}->{$chain}->{$resnum}->{"y"};
	      } else {
		$y=0;
	      }
	      if(exists($h4->{"xyz"}->{$chain}->{$resnum}->{"z"})){
		$z = $h4->{"xyz"}->{$chain}->{$resnum}->{"z"};
	      } else {
		$z=0;
	      }

	      @data = ($pdb, $chain, $resnum, $resname, $x, $y, $z, $AAS{$resname}, $asa_1, $asa_2, ($asa_2-$asa_1), $rel_1, $rel_2, $ratio, $cat_1, $cat_2);
	      $sth_two->execute(@data);
	      print $!;
	    }
	  }
	}
      }
    }
  }
  print STDERR "Filling the residue table ... ".(sprintf "%06d", $count)."/".(sprintf "%06d", $#{$list}+1)." -- finished ".localtime()."\n";
}



######
######  NACCESS - ASA related function
######
sub load_all_feed_DB_freeASA($$$$){

  my ($list, $dbh, $path_all, $path_chain) = @_;
  my ($pdb, $h, $h2, $h3, $h4, $chain, $resnum, $resname, $asa_1, $asa_2, $cat_1, $cat_2, @data, $rel_1, $rel_2, $ratio, $chains, $count, $x, $y, $z);

  my ($sth_one, %present);
  my $present = oneColQuery("SELECT distinct code FROM residue", $dbh);
  map { $present{$_}=1 } @{$present} ; # 2013
  ###while (my $row = $sth_one->fetchrow_arrayref) { $present{substr($row->[0],0,4)}=1;} $sth_one->finish; OLD

  # -- Prepare query
  my $sql = "UPDATE residue SET asa_in_BU='?', asa_alone='?', asa_rel_in_BU='?', asa_rel_alone='?' WHERE code='?' AND chain='?' AND resnum='?'";
  my $sth_two = $dbh->prepare($sql);

  if($sth_two){
    # Ok
  } else {
    die("INSERT QUERY WAS BAD  $!");
  }

  $count=1;

  foreach $pdb (@$list){

    if(!exists($present{$pdb})){

	print LOG $pdb."\n";
	
	if($count % 10 == 0){
	    print STDERR "Filling the residue table ... ".(sprintf "%06d", $count)."/".(sprintf "%06d", $#{$list}+1)."\r";	
	}
	$count++;
	
	print LOG "$pdb \n";

	# print STDERR "Come here $pdb \n";
	($h, $chains) = load_one_rsa_BU($path_all.$pdb);
	$h2 = load_one_rsa_split($path_chain.$pdb, $chains);
	$h3 = load_x_y_z($path_all.$pdb);
	$h4 = merge($h, $h2, $h3);

      #print "H=".Data::Dumper::Dumper($h);
      #print "\n----------------------------------------------------------\nH2=".Data::Dumper::Dumper($h2);
      #print "\n----------------------------------------------------------\nH4=".Data::Dumper::Dumper($h4);

      # H  -> {chain} --> {resnum} -> { resname / asa / category }
      #  \           \--> {total}  -> float
      #   \-> {tot}   \-> float

      foreach $chain (keys %{$h4->{"all"}}){

	if($chain ne "tot" && $chain ne " "){

	  foreach $resnum (keys %{$h4->{"all"}->{$chain}}){

	    if($resnum ne "tot"){

	      $resname = $h4->{"all"}->{$chain}->{$resnum}->{"name"};

	      $asa_1   = $h4-> {"all"} ->{$chain}->{$resnum}->{"asa"};
	      $asa_2   = $h4->{"split"}->{$chain}->{$resnum}->{"asa"};

	      $cat_1   = $h4-> {"all"} ->{$chain}->{$resnum}->{"category"};
	      $cat_2   = $h4->{"split"}->{$chain}->{$resnum}->{"category"};

	      $rel_1   = $h4-> {"all"} ->{$chain}->{$resnum}->{"asa_rel"};
	      $rel_2   = $h4->{"split"}->{$chain}->{$resnum}->{"asa_rel"};

	      if($rel_2 == 0){

		if($rel_1 !=0){
		  print LOG "Problem with $pdb, $chain, $resnum, ASA_REL in BU is $rel_1 and it is $rel_2 in the split\n";
		}
		$ratio=0;
	      } else {
		$ratio = $rel_1 / $rel_2 ;
	      }

	      if(exists($h4->{"xyz"}->{$chain}->{$resnum}->{"x"})){
		$x = $h4->{"xyz"}->{$chain}->{$resnum}->{"x"};
	      } else {
		$x=0;
	      }
	      if(exists($h4->{"xyz"}->{$chain}->{$resnum}->{"y"})){
		$y = $h4->{"xyz"}->{$chain}->{$resnum}->{"y"};
	      } else {
		$y=0;
	      }
	      if(exists($h4->{"xyz"}->{$chain}->{$resnum}->{"z"})){
		$z = $h4->{"xyz"}->{$chain}->{$resnum}->{"z"};
	      } else {
		$z=0;
	      }
	      ### my $sql = "UPDATE residue SET asa_in_BU='?', asa_alone='?', asa_rel_in_BU='?', asa_rel_alone='?' WHERE code='?' AND chain='?' AND resnum='?'";
	      @data = ($asa_1, $asa_2, $rel_1, $rel_2,$pdb, $chain, $resnum);
	      $sth_two->execute(@data);
	      print $!;
	    }
	  }
	}
      }
	$dbh->commit();
    }    
  }
  print STDERR "Filling the residue table ... ".(sprintf "%06d", $count)."/".(sprintf "%06d", $#{$list}+1)." -- finished ".localtime()."\n";
}



#
# Creates a H storing the ASA info for a ---** BU **---
#
# H  -> {chain} --> {resnum} -> { resname / asa / category }
#  \           \--> {total}  -> float
#   \-> {tot}   \-> float
#
#
sub load_one_rsa_BU($){

  my ($file) = @_;
  #print STDERR "$file loaded\n";

  open(RSA_IN, "<$file".".rsa") or die("Couldn't open $file".".rsa");

  my (%h, %chains, $line, @line, $chain, $resnum, $resname, $abs_acc, $rel_acc);

  while(<RSA_IN>){

    $line = $_;
    # print STDERR "$_";

    if($line =~ /^RES/){

      # --------- RES GLU C 315   190.68 110.7 147.74 109.6  42.94 114.5  73.31 121.6 117.37 104.8
      @line = split(/\s+/,$line);
      if(@line != 14){
	#print LOG "RES field didn't have 14 fields in file $file, line:\n$line";
      }
      #RES HIS A1000    17.29   9.5  17.29  11.8   0.00   0.0  13.85  14.3   3.44   4.0
      ($chain, $resnum, $resname, $abs_acc, $rel_acc) = (substr($line,8,1),substr($line,9,4),substr($line,4,3),substr($line,15,7),substr($line,22,6));
      $chain =~ s/ //g;
      $resnum =~ s/ //g;
      $resname =~ s/ //g;
      $abs_acc =~ s/ //g;
      $rel_acc =~ s/ //g;

      $h{$chain}->{$resnum}->{"name"}=$resname;
      $h{$chain}->{$resnum}->{"asa"}=$abs_acc;
      $h{$chain}->{$resnum}->{"asa_rel"}=$rel_acc;
      $chains{$chain}=1;
      if($abs_acc == 0){
	$h{$chain}->{$resnum}->{"category"}=1;
      } elsif($abs_acc > 0 && $abs_acc <= 20){
	$h{$chain}->{$resnum}->{"category"}=2;
      } elsif($abs_acc > 20 && $abs_acc <= 60){
	$h{$chain}->{$resnum}->{"category"}=3;
      } elsif($abs_acc > 60 && $abs_acc <= 100){
	$h{$chain}->{$resnum}->{"category"}=4;
      } elsif($abs_acc > 100){
	$h{$chain}->{$resnum}->{"category"}=5;
      } else {
	print LOG "Error, ASA not recognized: $abs_acc\n";
      }

    } elsif ($line =~ /^CHAIN/){
      # --------- CHAIN  1 C    13179.2      11216.0       1963.2       7955.1       5224.1
      @line = split(/\s+/,$line);
      if(@line != 8){
	print LOG "CHAIN field didn't have 8 fields in line\n$line";
      }
      ($chain, $abs_acc) = (substr($line,9,1),substr($line,10,11));
      $chain =~ s/ //g;
      $abs_acc =~ s/ //g;
      $h{$chain}->{"tot"}=$abs_acc;

    } elsif ($line =~ /^TOTAL/){
      # --------- TOTAL         13179.2      11218.0       1965.2       7957.2       5222.1
      @line = split(/\s+/,$line);
      if(@line != 6){
	print LOG "TOTAL field didn't have 6 fields in file $file, line:\n$line";
      }
      $abs_acc = substr($line,10,11);
      $abs_acc =~ s/ //g;
      $h{"tot"}=$abs_acc;

    } elsif ($line =~ /^REM/){
      #print $line;
      # --------- skip
      # --------- REM RES _ NUM      All-atoms   Total-Side   Main-Chain    Non-polar    All polar
      # --------- REM                ABS   REL    ABS   REL    ABS   REL    ABS   REL    ABS   REL
    } elsif ($line =~ /^END/){
      #print $line;
      # --------- skip
      # --------- END  Absolute sums over all chains
    } else {
      print LOG "This line wasn't expected in file $file, line:\n$line ";
    }
  }
  return(\%h, \%chains);
}

#
# Creates a H storing the ASA info for a ---** SPLIT CHAINS **---
#
# H  -> {chain} --> {resnum} -> { resname / asa / category }
#  \           \--> {total}  -> float
#   \-> {tot}   \-> float
#
#
sub load_one_rsa_split($$){

  my ($file, $chains) = @_;
  #print STDERR "$file loaded\n";

  my (%h, $line, @line, $chain, $resnum, $resname, $abs_acc, $each_chain, $rel_acc);

  foreach $each_chain (keys %$chains){

    if( -e $file."_$each_chain".".rsa"){

      open(RSA_IN, "<".$file."_$each_chain".".rsa");

      while(<RSA_IN>){

	$line = $_;
	# print STDERR "$_";

	if($line =~ /^RES/){

	  # --------- RES GLU C 315   190.68 110.7 147.74 109.6  42.94 114.5  73.31 121.6 117.37 104.8
	  @line = split(/\s+/,$line);
	  if(@line != 14){
	    #print LOG "RES field didn't have 14 fields in file $file, line:\n$line";
	  }

	  ($chain, $resnum, $resname, $abs_acc, $rel_acc) = (substr($line,8,1),substr($line,9,4),substr($line,4,3),substr($line,15,7),substr($line,22,6));
	  $chain =~ s/ //g;
	  $resnum =~ s/ //g;
	  $resname =~ s/ //g;
	  $abs_acc =~ s/ //g;
	  $rel_acc =~ s/ //g;

	  if($each_chain ne $chain){
	    print LOG "There is a problem in file $file!! chain $each_chain is different from $chain\n";
	  }

	  $h{$chain}->{$resnum}->{"name"}=$resname;
	  $h{$chain}->{$resnum}->{"asa"}=$abs_acc;
	  $h{$chain}->{$resnum}->{"asa_rel"}=$rel_acc;

	  if($abs_acc == 0){
	    $h{$chain}->{$resnum}->{"category"}=1;
	  } elsif($abs_acc > 0 && $abs_acc <= 20){
	    $h{$chain}->{$resnum}->{"category"}=2;
	  } elsif($abs_acc > 20 && $abs_acc <= 60){
	    $h{$chain}->{$resnum}->{"category"}=3;
	  } elsif($abs_acc > 60 && $abs_acc <= 100){
	    $h{$chain}->{$resnum}->{"category"}=4;
	  } elsif($abs_acc > 100){
	    $h{$chain}->{$resnum}->{"category"}=5;
	  } else {
	    print LOG "Error, ASA not recognized in file $file: $abs_acc\n";
	  }

	} elsif ($line =~ /^CHAIN/){
	  # --------- CHAIN  1 C    13179.2      11216.0       1963.2       7955.1       5224.1
	  @line = split(/\s+/,$line);
	  if(@line != 8){
	    print LOG "CHAIN field didn't have 8 fields in file $file, line:\n$line";
	  }
	  ($chain, $abs_acc) = (substr($line,9,1),substr($line,10,11));
	  $chain =~ s/ //g;
	  $abs_acc =~ s/ //g;
	  $h{$chain}->{"tot"}=$abs_acc;

	} elsif ($line =~ /^TOTAL/){
	  # --------- skip
	} elsif ($line =~ /^REM/){
	  #print $line;
	  # --------- skip
	} elsif ($line =~ /^END/){
	  #print $line;
	  # --------- skip
	} else {
	  print LOG "This line wasn't expected in file $file, line:\n$line ";
	}
      }
    } else {
      print LOG "Beware, one chain hasn't been found in: ".$file."_$chain".".rsa\n";
    }
  }
  return(\%h);
}

#
# Creates a H storing the coordinates info
#
# H  -> {chain} --> {resnum} -> {resname} -> { x / y / z }
#
sub load_x_y_z($){

  my ($file) = @_;
  #print STDERR "$file loaded\n";

  open(RSA_IN, "<$file".".asa") or die("Cannot open $file".".asa");

  my (%h, $line, @line, $chain, $resnum, $resname, $atom, $x, $y, $z);

  while(<RSA_IN>){
    #print STDERR $line ;
    $line = $_;
    # print STDERR "$_";

    if($line =~ /^ATOM/){

      #print STDERR $line ;

      # ATOM   4378  CA  HIS D 146       4.753 -12.229 -14.962   3.816  1.87
      #@line = split(/\s+/,$line);
      #if(@line != 14){
      #	print LOG "RES field didn't have 14 fields in file $file, line:\n$line";
      #}

      ($chain, $resnum, $resname, $atom, $x, $y, $z) = (substr($line,21,1),
							substr($line,22,4),
							substr($line,17,3),
							substr($line,11,4),
							substr($line,30,8),
							substr($line,38,8),
							substr($line,46,8));
      $chain =~ s/ //g;
      $resnum =~ s/ //g;
      $resname =~ s/ //g;
      $atom =~ s/ //g;
      $x =~ s/ //g;
      $y =~ s/ //g;
      $z =~ s/ //g;

      if($atom eq "CA"){

	$h{$chain}->{$resnum}->{"x"}=$x;
	$h{$chain}->{$resnum}->{"y"}=$y;
	$h{$chain}->{$resnum}->{"z"}=$z;
      }

    } else {
      print LOG "This line wasn't expected in file $file".".asa, line:\n$line ";
    }
  }
  return(\%h);
}


# Takes the Hash containing infos about the ASA of a complex
# in two forms:
# - complete
# - split by chain
#
# Return a HASH of the form:
#
#   |--> {all}    -> h
# H |--> {split}  -> h2
#   |--> {merged} -> {chain} -> {resnum} -> { name / category_1 / category_2 / asa_1 / asa_2 } = value
#   |           \           \-> {tot} = delta
#   |            \-> {tot} = delta
#   |
#   |  NOT INTERFACE
#   |--> {1}      ->    # Those that do not change category and for which the ASA does NOT vary more than 10 A
#   |--> {2}      ->    # OR those that change category BUT for which the ASA does NOT vary more than 10 A
#   |--> {3}      ->    # in the latter case, keep the original category.
#   |--> {4}      ->
#   |--> {5}      ->
#   |  INTERFACE
#   |--> {55}     ->    # Those that do not change category BUT for which the ASA varies more than 10 A
#   |--> {44}     ->
#   |--> {33}     ->
#   |--> {22}     ->
#   |--> {11}     ->
#   |--> {54}     ->    # Those that do change category BUT AND for which the ASA varies more than 10 A
#   |--> {53}     ->
#   |--> {52}     ->
#   |--> {51}     ->
#   |--> {43}     ->
#   |--> {42}     -> {chain} -> {resnum}  = delta
#   |--> {41}     ->
#   |--> {32}     ->
#   |--> {31}     ->
#   |--> {21}     ->
#
#
#
sub merge($$$){

  my($h, $h2, $h3) = @_;

  my (%r, $chain, $resnum, $tot, $asa_1, $asa_2, $cat_1, $cat_2);
  my @c2 = keys %$h2;

  foreach $chain (keys %$h){

    if($chain ne "tot" &&  !exists($h2->{$chain})){
      print LOG "Couldn't compare those two structures, they don't have the same chains ($chain)\n";
      print LOG "Chains in h2 = ".join("\t", @c2)."\n";
      return \%r;
    }
  }

  $tot=0;

  foreach $chain (keys %$h){

    if($chain ne "tot"){

      foreach $resnum (keys %{$h->{$chain}}){

	if($resnum ne "tot"){

	  $asa_1 = $h->{$chain}->{$resnum}->{"asa"};
	  $asa_2 = $h2->{$chain}->{$resnum}->{"asa"};

	  if($asa_2 - $asa_1 < -0.1){
	    print LOG "Error, Delta ASA < 0, ".($asa_2 - $asa_1)." at chain $chain, res $resnum!!!\n";
	  } elsif($asa_2 - $asa_1 > 10) {

	    $cat_1 = $h->{$chain}->{$resnum}->{"category"};
	    $cat_2 = $h2->{$chain}->{$resnum}->{"category"};

	    $r{$cat_2.$cat_1}->{$chain}->{$resnum} = $asa_2 - $asa_1 ;

	  } else {

	    $cat_2 = $h2->{$chain}->{$resnum}->{"category"};
	    $r{$cat_2}->{$chain}->{$resnum} = $asa_2 - $asa_1 ;

	  }

	} else {
	  $tot += $h2->{$chain}->{"tot"};
	}
      }
    }
  }

  $r{"merged"}->{"tot"}= $tot - $h->{"tot"} ;
  $r{"split"} = $h2 ;
  $r{"all"} = $h ;
  $r{"xyz"} = $h3 ;
  $r{"split"}->{"tot"} = $tot ;

  return(\%r);
}

######
######  END OF NACCESS ASA RELATED FUNCTIONS
######

######
###### GRAPHVIZ FUNCTIONS
######
#
# takes a hash and a pdb code
# prints the code to generate a graph of a complex using graphviz
sub perl2graphviz($$$$$){

  my ($hash, $pdb, $homologues, $file_out, $homologies) = @_ ;
  my $compo = get_compo_accurate($hash,$pdb, $homologies);

  #print STDERR Dumper($compo);

  my @colors = ("chocolate ", "brown ", "coral3", "darkorange","gold2 ","gold4", "aquamarine4","azure3","blue3","blue4","blueviolet","chartreuse4","darkolivegreen3","darkolivegreen4", "aliceblue", "antiquewhite","antiquewhite3","aquamarine3");
  my @shapes = ("circle","pentagon","hexagon", "septagon","octagon","rect", "diamond","trapezium","parallelogram");
  my ($current_color, $current_shape) = (0,0);
  my ($chain, $chain_2, @sfs, $each_sf, $each_grp, @chains);

  my $size = int( ((keys %{$hash->{$pdb}})+1)/4  );
  if( $size > 3){
    $size=3;
  }
  $size=3;
  open(OUT, ">".$file_out."/".$pdb.".viz");
  print OUT "graph G\n";
  print OUT "{\n";
  print OUT "\tsize=\"$size,$size\"\n";
  print OUT "\tmargin=0\n";
  ## print OUT "\tmode=hier\n";
  ## print OUT "\tdim=6\n";

  ## Determines the colors and shapes in the graph
  @sfs = keys %{$compo} ;
  for $each_sf (0..$#sfs){

    if($current_shape > $#shapes){
      print STDERR "Be carefull: not enought shapes for $pdb\n";
      $current_shape=0;
    }

    for $each_grp (0..$#{$compo->{$sfs[$each_sf]}}){

      if($current_color > $#colors){
	print STDERR "Be carefull: not enought colors for $pdb\n";
	$current_color = 0;
      }

      @chains = split(//,$compo->{$sfs[$each_sf]}->[$each_grp]);
      # print OUT "\t{node [width=.01, height=.01, style=filled, shape=$shapes[$current_shape], color=$colors[$current_color], fontsize=20] ".join(" ",@chains)."};\n";

      print OUT "\t{node [fixedsize=true, width=.3, height=.3, style=filled, shape=$shapes[$current_shape], color=$colors[$current_color], fontname=Helvetica, fontsize=15, URL=\"$pdb$chains[0]\" ] ".join(" ",@chains)." }\n";

      #{node [width=.01, height=.01,style=filled, shape=circle, color=red, fontsize=16] V1 V2};
      $current_color++ ;
    }
    $current_shape++ ;
  }

  ## gives the wiring
  my %node_seen=();
  foreach $chain (keys %{$hash->{$pdb}}){
    $node_seen{$chain}=0 ;

    foreach $chain_2 (keys %{$hash->{$pdb}->{$chain}}){
      if($chain_2 ne "prop"){

	if( ! defined $node_seen{$chain_2}){
	  print OUT "\t$chain -- $chain_2 "."[label=\"".$hash->{$pdb}->{$chain}->{$chain_2}->{"size"}."\",fontsize=10, fontname=Helvetica];\n";
	}
      }
    }
  }
  print OUT "}\n";
  close(OUT);
}


# takes a hash
#       a pdb code
#       a path
# prints the code in PATH to generate a graph of the given 
# complex using graphviz
sub write_viz_topo($$$){
  my ($hash, $pdb, $path) = @_ ;
  my ($chain, $chain_2, $monomer);

  open(OUT, ">".$path.$pdb.".viz");
  print OUT "graph G\n";
  print OUT "{\n";
  print OUT "\tsize=\"1,1\"\n";
  print OUT "\tmargin=0\n";
  print OUT "\tmode=KK\n";
  print OUT "\tdim=2\n";

  ## Determines the colors and shapes in the graph
  print OUT "\tnode [fixedsize=true, width=.2, height=.2, style=filled, shape=circle, color=black, fontsize=0.001, fillcolor=black]\n";

  ## gives the wiring
  $monomer=1;
  my %node_seen=();
  foreach $chain (keys %{$hash->{$pdb}}){

    $node_seen{$chain}=0 ;

    foreach $chain_2 (keys %{$hash->{$pdb}->{$chain}}){
      if($chain_2 ne "prop"){

	if($chain_2 eq $chain){
	  print OUT "\t$chain -- $chain_2 \n";
	  $monomer=0;
	}
	if( ! defined $node_seen{$chain_2}){
	  print OUT "\t$chain -- $chain_2 \n";
	  $monomer=0;
	}
      }
    }
    if($monomer){
      print OUT "\t$chain\n";
    }
  }
  print OUT "}\n";
  close(OUT);
}


# Returns a HASH allowing to acces the hydrophobicity of a residue
sub aa_hydrophobicity()
  {
    my %hydro=();
    my @aas =    ("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "ASP", "ASN", "LYS", "GLU", "GLN", "ARG", "HIS", "PHE", "CYS", "TRP", "TYR", "MET", "PRO");

    # Scale taken from "Pacios, J. Chem. Inf. Comput. Sci., vol41, N.5, 2001"
    my @scale1 = (0.35,0.45,0.70,0.61,1.00,0.25,0.23,0.12,0.15,0.00,0.11,0.11,0.08,0.16,0.41,0.40,0.39,0.25,0.18,0.23);

    # Kyte & Doolittle Scale
    my @scale2 = (-0.4,   1.8,  4.2,    3.8,   4.5,   -0.8,  -0.7,  -3.5,  -3.5,  -3.9,  -3.5,  -3.5,  -4.5,  -3.2,  2.8,  2.5,  -0.9,  -1.3,  1.9,  1.6);

    my $i;

    foreach $i (0..$#aas)
      {
	$hydro{$aas[$i]}->{"1"}=$scale1[$i];
	$hydro{$aas[$i]}->{"2"}=$scale2[$i];
      }
    return(\%hydro);
  }


# Returns a HASH allowing to acces the hydrophobicity of a residue
sub aa_interfacity($$)
  {

    my ($org, $cut) = @_;

    my @aas =  ("GLY",  "ALA",  "VAL",  "LEU",  "ILE",  "SER",  "THR",  "ASP",  "ASN",  "LYS",  "GLU",  "GLN",  "ARG",  "HIS",  "PHE",  "CYS",  "TRP",  "TYR",  "MET",  "PRO");
    my @scale;

    if($org eq "sc" && $cut == 25){
      @scale = (0,0.34,0.71,0.83,0.74,-0.02,0.09,-0.8,-0.51,-1.4,-0.76,-0.12,0.34,0.23,1.09,1.35,1.33,0.94,0.88,-0.2);
    }
    if($org eq "sc" && $cut == 35){
      @scale = (0.01,0.19,0.68,0.98,0.94,-0.13,0.14,-0.71,-0.36,-1.2,-0.61,-0.16,0.42,0.31,1.28,1.6,1.47,1.17,0.74,-0.18);
    }
    if($org eq "ec" && $cut == 25){
	@scale =   (-0.1771,0.0062, 0.7599, 0.9138, 1.1109, 0.1376, 0.1031,-0.7485,-0.2693,-1.1806,-0.7893,-0.4114,-0.0876, 0.1214, 1.2727, 1.0372, 0.7925, 0.8806, 1.0124,-0.1799);  ## PUBLISHED in PNAS 2012          
        # @scale = (-0.17,  -0.03 , 0.82,   0.96,   1.14,   0.14,   0.12,  -0.75,  -0.26,  -1.15,  -0.8,   -0.47,  -0.1,    0.16,   1.4,    1.09,   0.95,   0.93,   1.06,  -0.19);  ## OBTAINED AFTER REMOVAL OF MEMBRANE PROTS	
    }
    if($org eq "ec" && $cut == 35){
	@scale = (-0.12,-0.1,0.67,0.91,1.05,0.09,0.02,-0.75,-0.19,-0.7,-0.81,-0.38,0.24,0.31,1.35,1.07,1.09,1.16,1.14,-0.18);
    }
    if($org eq "hs" && $cut == 25){
	@scale = (-0.06,0.25,0.61,0.82,0.87,-0.09,0.06,-0.63,-0.36,-1.09,-0.81,-0.34,-0.02,0.13,1.14,0.88,1.11,0.96,0.92,-0.37);
    }
    if($org eq "hs" && $cut == 35){
	@scale = (-0.09,0.19,0.62,0.9,0.98,-0.13,0.02,-0.58,-0.3,-0.9,-0.73,-0.25,0.13,0.25,1.25,0.9,1.29,1.13,0.99,-0.35);
    }
    if($org eq "exp"){
	@scale = (-0.15, 0.05,0.03,0.44,0.1 ,0.16,-0.06,-0.19,-0.12,-0.57,-0.28, 0.4,  0.29,0.49,0.24,0.69,0.81,0.22,-0.18, 0.39);
    }

    ## Before AREAIMOL
    ## scale = (-0.14,-0.04,0.65,0.85,0.94,0.12,0.04,-0.78,-0.27,-0.91,-0.83,-0.46,0.02,0.25,1.22,1.01,0.8,0.94,1.06,-0.17);

    # ABSOLUTE DEFINITION YEAST (50A)
    #   A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S 
    #1.45 3.74 0.54 0.53 2.81 1.29 1.19 2.14 0.25 2.29 2.21 0.70 0.96 0.89 1.02 0.98 
    #   T    V    W    Y 
    #1.24 2.04 2.68 2.16

    # RELATIVE DEFINITION YEAST (30%)
    #    A    C    D    E    F    G    H    I    K    L    M    N    P    Q    R    S 
    # 1.32 3.19 0.46 0.54 3.36 0.97 1.25 2.30 0.31 2.45 2.44 0.58 0.82 0.93 1.36 0.86 
    #    T    V    W    Y
    # 1.08 2.04 4.66 2.96

    # RELATIVE FREQUENCY YEAST (LOG SCALE)
    #    A     C     D     E     F     G     H     I     K     L     M     N     P 
    # 0.28  1.16 -0.78 -0.61  1.21 -0.03  0.22  0.83 -1.17  0.90  0.89 -0.54 -0.20 
    #    Q     R     S     T     V     W     Y 
    #-0.07  0.31 -0.15  0.07  0.71  1.54  1.09 

    my %hydro=();

    # Scale obtained by looking at the ratio of AA presence at surfaces and interfaces.
    my $i;

    foreach $i (0..$#aas)
      {
	$hydro{$aas[$i]}->{"1"}=$scale[$i];
	$hydro{$aas[$i]}->{"2"}=$scale[$i];
      }
    return(\%hydro);
  }

sub get_pax_sizeV2($){
    
    my ($file) = (@_);

    my %DEF = ();

    open(PAX_RES, "<$file") or die("Cannot open $file\n" );

    while(<PAX_RES>){
	
	chomp($_);
	my $line = $_;

	if($line =~ /^([a-zA-Z0-9_]+) +\[0\.0\.0\.0\].*CA=(\d+),.*$/){

	    $DEF{$1}=$2;		    
	}
    }
    return(\%DEF);
}

sub get_pax_sizeV3($){

    my ($file) = (@_);

    my %DEF = ();

    if( -e $file){
	open(PAX_RES, "<$file");
    } else {
	return(\%DEF);
    }

    my $line;

    while($line = <PAX_RES>){

	chomp($line);
	if(! ($line =~ /^#/)){
	    $line =~ s/^ +//;
	    my @elemts = split(/ +/,$line);	
	    $DEF{ $elemts[12] }= $elemts[10];	
	}	    
    }
    return(\%DEF);
}

sub get_pax_sizeV5($){

    my ($file) = (@_);

#Rank  K-Score  G-Score J-Score M-Score T-Score   RMSD   N/*   P/$   I/@   I/!   Len  Seg  TP   Match [Family]
#====  =======  ======= ======= ======= =======   ====   ===   ===   ===   ===   ===  ===  ==  ==================
#    1   164.00   164.00  1.0000  1.0000  1.0000   0.00   164   164   164 100.0   164    1   1  104l_2  0.0.0.0
#    2   164.00   164.00  1.0000  1.0000  1.0000   0.00   164   164   164 100.0   164    1   1  104l_p2  0.0.0.0

    my %DEF = ();

    if( -e $file){
	open(PAX_RES, "<$file");
    } else {
	return(\%DEF);
    }

    my $line;

    while($line = <PAX_RES>){

	chomp($line);
	if(! ($line =~ /^#/)){
	    $line =~ s/^ +//;
	    my @elemts = split(/ +/,$line);	
	    $DEF{ $elemts[14] }= $elemts[11];	
	}	    
    }
    return(\%DEF);
}

sub get_pax_sizeV3_pair($){

    my ($file) = (@_);

    my %DEF = ();

    if( -e $file){
	open(PAX_RES, "<$file");
    } else {
	return(0);
    }
    # open(PAX_RES, "<$file") or die("Cannot open $file\n" );

    my $line;
    my $ind=1;

#Rank  K-Score  G-Score J-Score H-Score T-Score  RMSD   M/*   N/$   I/@   Len  TP   Match [Family]
#====  =======  ======= ======= ======= =======  ====   ===   ===   ===   ===  ==  ==================
#    1   165.00   165.00  1.0000  1.0000  1.0000   0.00   165   165   165   165   1  1nmk_1  0.0.0.0


    while($line = <PAX_RES>){

	chomp($line);
	if(! ($line =~ /^#/)){
	    $line =~ s/^ +//;
	    my @elemts = split(/ +/,$line);	
	    $DEF{ $elemts[12]."_".$ind }->{"len"} = $elemts[10];	
	    $DEF{ $elemts[12]."_".$ind }->{"K"} = $elemts[1];	
	    $DEF{ $elemts[12]."_".$ind }->{"G"} = $elemts[2];	
	    $DEF{ $elemts[12]."_".$ind }->{"T"} = $elemts[5];	
	    $ind++;
	}	    
    }
    return(\%DEF);
}

sub get_pax_chainMatch($){
		       
    my ($file) = (@_);
    my $chainMatch = "";
    my ($Q, $T, $Npairs, $Ksum, $Gsum, $Kscore, $Gscore, $TMscore);
    my $BEST = {};
    my $minTM = 200;

    if( -e $file){			    
	
	open(IN_CHAINS, "<$file");
	
	while(<IN_CHAINS>){
##	    #===============================================================================
##          #  Q  T  Npairs       K-Sum       G-Sum   K-Score   G-Score  TM-Score
##          #===============================================================================
##          A  A     250    227.0353    177.5716    0.9045    0.7075    0.9201
##          A  D       0      0.0000      0.0000    0.0000    0.0000    0.0000

	    #print STDERR $_;

	    if( $_ =~ /^[ \t]+([^ ]+)[ \t]+([^ ]+)[ \t]+([^ ]+)[ \t]+([^ ]+)[ \t]+([^ ]+)[ \t]+([^ ]+)[ \t]+([^ ]+)[ \t]+([^ ]+)$/){

		($Q, $T, $Npairs, $Ksum, $Gsum, $Kscore, $Gscore, $TMscore) = ($1,$2,$3,$4,$5,$6,$7,floor(100*$8));

		if( !exists $BEST->{$Q} ){
		    $BEST->{$Q}->{"T"}=$T;
		    $BEST->{$Q}->{"Npairs"}=$Npairs;
		    $BEST->{$Q}->{"Ksum"}=$Ksum;
		    $BEST->{$Q}->{"Gsum"}=$Gsum;
		    $BEST->{$Q}->{"Kscore"}=$Kscore;
		    $BEST->{$Q}->{"Gscore"}=$Gscore;
		    $BEST->{$Q}->{"TMscore"}=$TMscore;
		} elsif ( $BEST->{$Q}->{"Npairs"} < $Npairs  ){
		    $BEST->{$Q}->{"T"}=$T;
		    $BEST->{$Q}->{"Npairs"}=$Npairs;
		    $BEST->{$Q}->{"Ksum"}=$Ksum;
		    $BEST->{$Q}->{"Gsum"}=$Gsum;
		    $BEST->{$Q}->{"Kscore"}=$Kscore;
		    $BEST->{$Q}->{"Gscore"}=$Gscore;
		    $BEST->{$Q}->{"TMscore"}=$TMscore;
		}
	    }
	}
	
    } else {
	print STDERR "$file does not exists (in get_pax_chainMatch)\n";
    }
    
    foreach my $chain_pax (keys %$BEST){
	if( $BEST->{$chain_pax}->{"TMscore"} < $minTM){
	    $minTM = $BEST->{$chain_pax}->{"TMscore"};
	}
	$chainMatch .= $chain_pax.":".$BEST->{$chain_pax}->{"T"}.":".$BEST->{$chain_pax}->{"Npairs"}.":".$BEST->{$chain_pax}->{"TMscore"}.";"
    }
    return ($chainMatch, $minTM);
}		          

## Check the consistensy of a complex considering it's number of subunits, symmetry, and whether it's a homomer or not
##
sub check_consistency($$$){
		      
    my ($NSUB, $SYM, $HOMO) = @_;

    my ($MULTIPLY, $TYPE);

    if($SYM =~ /^C(\d+)$/){
	$MULTIPLY = 1 ;
	$TYPE     = $1;

    } elsif($SYM =~ /^NPS$/){
	$MULTIPLY = 1 ;
	$TYPE     =1  ;
	
    } elsif($SYM =~ /^D(\d+)$/){
	$MULTIPLY = 2 ;
	$TYPE     =$1 ;
	
    } elsif($SYM =~ /^Tet/){
	$MULTIPLY = 12;
	$TYPE     = 1;

    } elsif($SYM =~ /^Oct/){
	$MULTIPLY = 24;
	$TYPE     = 1;
    } elsif($SYM =~ /^Ico/){
	$MULTIPLY = 60;
	$TYPE     = 1;
    } else {

	$MULTIPLY = 1000;
	$TYPE     = 1000;	
    }

    if($HOMO == 1){

	if( $NSUB == $MULTIPLY * $TYPE ){
	    return 1;
	} else {
	    return 0;
	}

    } else {

	if( $NSUB % ($MULTIPLY*$TYPE) == 0 ){
	    return 1;
	} else {
	    return 0;
	}
    }
}


return 1;

