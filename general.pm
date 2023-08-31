use strict ;
use Cwd;
use FileHandle;
use File::stat;
use File::Copy;
use File::Basename;
use Term::Cap;
use POSIX;

sub format_stderr($){

    my ($Num) = @_;

    my $termios = POSIX::Termios->new;
    $termios->getattr;
    my $ospeed = $termios->getospeed;
    my $terminal = Term::Cap->Tgetent({ TERM => undef, OSPEED => $ospeed });
    $terminal->Trequire("up");  # move cursor up
    my $UP = $terminal->Tputs("up");
    my ($B4, $AFTER) = ("","");
    for(my $k=0; $k<$Num; $k++){
	$B4 .= "\n";
	$AFTER .= $UP;
    }

    return ($B4, $AFTER);
}

sub time_stamp($$) {
  my ($d,$t);
  my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);

        $year += 1900;
        $mon++;
        $d = sprintf("%4d-%2.2d-%2.2d",$year,$mon,$mday);
        $t = sprintf("%2.2d:%2.2d:%2.2d",$hour,$min,$sec);
        return($d,$t);
}


# Takes any hash.
# H -> {id1} -> {field1} = value
#           \-> {field2} = value
#
# Take an array of fields ["field1", "field2"]
#
# Return a list of ids so that values in field 1 are
# ordered first, and then value if field 2.
#
sub order($$$){

#  my ($hash, $field1, $field2) = @_;

#  my (@list_tmp, @list_final);

  # -- first makes a list of lists according to field 1
#  foreach my $id (keys %$hash){

 #   if( exists( $list_tmp[$hash->{$id}->{$field1}]) ){
  #    push(@{$list_tmp[$hash->{$id}->{$field1}]},$id);
   # } else {
    #  $list_tmp[$hash->{$id}->{$field1}] = [$id];
    #}

  #}

}

sub cmp_3($$){

  my ($a, $b) = @_;

  if($a->[1] != $b->[1]){
    return $a->[1] < $b->[1];
  } else{
    return $a->[2] < $b->[2];
  }
}

sub rep($$){
  my ($string, $N) = @_;
  my $res = "";
  map { $res .= $string } (1..$N);
  return $res;
}

sub rep_array($$){
  my ($val, $N) = @_;
  my @res;
  map { push(@res, $val) } (1..$N);
  return \@res;
}


# Look if value1 is within X% of value2
#
sub within($$$){

  my ($value1, $value2, $X) = @_;

  if( ($value1 >= ($value2-($X*$value2))) && ($value2 <= ($value2+($X*$value2)))){
    return 1;
  } else {
    return 0;
  }
}

# Compares 2 arrays elmt by elmt
# return the number of mismatchs
sub array_eq($$){

  my ($array1, $array2) = @_ ;
  if($#{$array1} != $#{$array2}){
    return -1 ;
  }

  my ($i);
  my $nb_mismatch = 0;

  for $i (0..$#{$array1}){
    if( $array1->[$i] != $array2->[$i]){
      $nb_mismatch++;
    }
  }
  return $nb_mismatch ;
}


# Takes an array, look if all the values are equal in it
sub is_array_uniform($){

  my $array = shift ;

  # If the array is of size < 2 there is no point
  if( $#{$array} < 0){
    die("there is no point using the function is_array_uniform because the array is of size 0\n");
  } else {
    my ($elmt);
    my ($equal, $i) = (1,1);

    $elmt=$array->[0];

    while($equal==1 && $i <= $#{$array}){

      if($elmt != $array->[$i]){
	$equal=0;
      }
      $i++;
    }
    return $equal ;
  }
}

# Takes 2 arrays as argument
# Returns the number of elements that are not common
sub get_exclusion($$){

  my ($array1, $array2) = @_ ;
  my (%hash, $elmt) ;
  my $exclu = 0 ;

  foreach $elmt (@$array1){
    if(exists $hash{$elmt}){
      $hash{$elmt}++;
    } else {
      $hash{$elmt} = 1;
    }
  }

  foreach $elmt (@$array2){
    if (exists $hash{$elmt}){
      $hash{$elmt}-- ;
    } else {
      $hash{$elmt} = -1 ;
    }
  }

  foreach $elmt (keys %hash){
    $exclu += abs($hash{$elmt}) ;
  }
  return $exclu ;
}

# Takes 2 arrays as argument
# Returns a hash of array
# H->{"union"}
# \->{"intersect"}
# \->{"exclu1"}
# \->{"exclu2"}
sub inter_exclu_union($$){

  my ($array1, $array2) = @_ ;
  my (%result, %hash, $elmt, @union) ;

  foreach $elmt (@$array1){
    $hash{$elmt} = "exclu1" ;
  }

  foreach $elmt (@$array2){
    if ($hash{$elmt} eq "exclu1"){
      $hash{$elmt} = "intersect" ;
    } else {
      $hash{$elmt} = "exclu2" ;
    }
  }

  @union = keys %hash ;
  @{$result{"union"}} = @union ;
  $result{"exclu1"}=[];
  $result{"exclu2"}=[];
  $result{"intersect"}=[];
  foreach $elmt (@union){
    push(@{$result{$hash{$elmt}}}, $elmt)
  }
  return \%result ;
}

# Takes 2 arrays references
# Return a ref to an array that is the intersection of the 2 arrays
sub intersect($$){
  my ($array1, $array2);
  ($array1, $array2) = @_ ;

  my ($elmt, %defined, @intersection);

  foreach $elmt (@{$array1}){
    $defined{$elmt}=1;
  }

  foreach $elmt (@{$array2}){
    if( exists($defined{$elmt})){
      push(@intersection, $elmt);
    }
  }
  @intersection = sort {$a <=> $b} @intersection ;
  return \@intersection ;
}


# Takes 2 arrays as argument
# Returns 1 if they have the same content
#         0 otehrwise
#
sub ar_same_content($$){

  my ($array1, $array2) = @_ ;
  my (%hash, $elmt) ;

  foreach $elmt (@$array1){
    if(exists $hash{$elmt}){
      $hash{$elmt}++;
    } else {
      $hash{$elmt}=1;
    }
  }

  foreach $elmt (@$array2){
    if (exists $hash{$elmt}){
      $hash{$elmt}--;
    } else {
      return 0;
    }
  }

  foreach $elmt (keys %hash){

    if( $hash{$elmt} != 0 ){
      return 0;
    }
  }
  return 1 ;
}

# return the standard deviation
# (Sum_i( (x_i - X)^2/n-1 ))^0.5
#
sub SD($){

  my ($array) = @_;
  my $mean = MEAN($array);
  my $sum = 0;
  map {$sum += ($_ - $mean)**2 } @$array ;

  if($sum == 0){
    return "0";
  } elsif($#{$array}+1 == 0 || $#{$array} == 0){
    return "NaN";
  } else {
    return sqrt($sum / $#{$array}) ;
  }
}

# return the Variance
# (Sum_i( (x_i - X)^2/n-1 ))
#
sub VAR($){

  my ($array) = @_;
  my $mean = MEAN($array);
  my $sum = 0;
  map {$sum += ($_ - $mean)**2 } @$array ;

  if($sum == 0){
    return "0";
  } elsif($#{$array}+1 == 0){
    return "NaN";
  } else {
    return $sum / $#{$array} ;
  }
}

# Returns the mean of an array
sub MEAN($){

  my ($array) = @_ ;
  my $sum =0;
  map { $sum += $_ } @$array ;
  if($sum == 0){
    return "0";
  } elsif($#{$array}+1 == 0){
    return "NaN";
  } else {
    return $sum/($#{$array}+1);
  }
}

# Returns the sum of an array
sub SUM($){

  my ($array) = @_ ;
  my $sum =0;
  map { $sum += $_ } @$array ;
  return $sum;
}

# Take 2 values, return the smallest one
sub min($$){
  my ($a, $b)= @_;
  if($a < $b){
    return $a;
  } else {
    return $b;
  }
}

# Take 2 values, return the biggest one
sub max($$){
  my ($a, $b) = @_ ;
  if($a > $b){
    return $a ;
  } else {
    return $b ;
  }
}

# Returns the sum of an array
sub MAX($){

  my ($array) = @_ ;
  my $max = $array->[0];
  map { $max = max($max, $_) } @$array ;
  return $max;
}

# Returns the sum of an array
sub MIN($){

  my ($array) = @_ ;
  my $min = $array->[0];
  map { $min = min($min, $_) } @$array ;
  return $min;
}


# Load a table of the form:
# X     key1   key2     .... keyN
# ref1  val1.1  val1.2  ....
# ref2  val2.1  val2.2  ....
# ref3  val3.1  val3.2  ....
# .     .       .
# .     .       .
# refn  valn.1  valn.2  .... valn.N
#
# returns a hash of the form:
#
# H -->  {"ref1"} -> {"key1"} = val1.1
#  \             \-> {"key2"} = val1.2
#   \
#    \-> {"ref2"} -> {"key1"} = val2.1
#                \-> {"key2"} = val2.2
#  etc..
#
# Obviously there MUST not be 2 keys having the same value!!
sub T2H($){

  my ($path) = @_;

  open(TABLE, "<$path") or die("The file $path does not exists!");

  my ($line, %hash, @keys, @vals, $i);

  $line = "##";
  while($line =~ /^#/){
      $line = <TABLE>;
  }
  #print STDERR $line."YOYOYO\n";
  chomp($line);
  @keys = split(/\t/,$line);

  while(<TABLE>){
    $line = $_;
    chomp($line);
    @vals = split(/\t/,$line);

    for $i (1..$#vals){
      $hash{$vals[0]}->{$keys[$i]}=$vals[$i] ;
    }
  }

  close(TABLE);
  return \%hash ;
}

### Same as T2H, but the keys are given as argument
sub T2H_nokey($$){

  my ($path, $keys) = @_;

  open(TABLE, "<$path") or die("The file $path does not exists!");

  my ($line, %hash, @keys, @vals, $i);

  $line = "##";
  while($line =~ /^#/){
      $line = <TABLE>;
  }
  chomp($line);
  @keys = @$keys;

  while(<TABLE>){
    $line = $_;
    chomp($line);
    @vals = split(/\t/,$line);

    for $i (1..$#vals){
      $hash{$vals[0]}->{$keys[$i]}=$vals[$i] ;
    }
  }

  close(TABLE);
  return \%hash ;
}




##
## Same as T2H except that the resulting H uses the first two coloumns as KEYS
##
sub T2H_2keys($$){

  my ($path, $reverse) = @_;

  open(TABLE, "<$path") or die("The file $path does not exists!");

  my ($line, %hash, @keys, @vals, $i);

  $line = "##";
  while($line =~ /^#/){
      $line = <TABLE>;
  }
  #print STDERR $line."YOYOYO\n";
  chomp($line);
  @keys = split(/\t/,$line);

  print STDERR "KEYS: @keys\n";

  while(<TABLE>){
    $line = $_;
    chomp($line);
    @vals = split(/\t/,$line);

    if(exists $hash{$vals[1]}->{$vals[0]} && $hash{$vals[1]}->{$vals[0]}->{$keys[$#vals]} < $vals[$#vals]){
	print STDERR "Be carreful in T2H_2keys: a value is being replaced and should not ($line)\n";	
    } else {

	for $i (2..$#vals){
	    if($reverse == 0){
		$hash{$vals[0]}->{$vals[1]}->{$keys[$i]}=$vals[$i] ;
	    } elsif($reverse == 1) {
		$hash{$vals[1]}->{$vals[0]}->{$keys[$i]}=$vals[$i] ;	    
	    }
	}
    }
  }

  close(TABLE);
  return \%hash ;
}



# SAME EXCEPT THAT THE MATRIX CONTAINS SPACES NOT TABS
sub T2H_s($){

  my ($path) = @_;

  open(TABLE, "<$path") or die("The file $path does not exists!");

  my ($line, %hash, @keys, @vals, $i);

  $line = <TABLE>;
  chomp($line);
  @keys = split(/\s+/,$line);

  while(<TABLE>){
    $line = $_;
    chomp($line);
    @vals = split(/\s+/,$line);

    for $i (1..$#vals){
      $hash{$vals[0]}->{$keys[$i]}=$vals[$i] ;
    }
  }

  close(TABLE);
  return \%hash ;
}


# Loads a CSV file with a comma as delimiter,
# and a given field is the key.
sub T2H_ensembl($$){

  my ($path, $key) = @_;

  open(TABLE, "<$path") or die("The file $path does not exists!");

  my ($line, %hash, @keys, @vals, $i, $key_index, $passed);

  $line = <TABLE>;
  chomp($line);
  @keys = split(/,/,$line);
  $passed = 0;

  foreach $i (0..$#keys){

    if($keys[$i] == $key){
      $key_index = $i ;
      $passed = 1;
    }
  }

  if($passed == 0){
    die("The key in T2H_ensembl has not been found\n");
  }

  while(<TABLE>){
    $line = $_;
    chomp($line);
    @vals = split(/,/,$line);

    for $i (0..$#vals){
      if($i != $key_index){
	$hash{$vals[$key_index]}->{$keys[$i]}=$vals[$i] ;
      }
    }
  }

  close(TABLE);
  return \%hash ;
}


# Open a file, reads it, and return the text.
sub include_text($){
  my ($file) = @_;
  my $text = "";
  open(IN_FILE, "<$file") or die ("Cannot open $file");

  while(<IN_FILE>){
    $text .= $_;
  }
  return $text;
}

# To use with T2H:
# Takes a HASH from T2H, a list of KEYS and a property
# --> returns the list of properties corresponding to eack KEY
sub T2H_retrieve($$$)
  {
    my ($T2H, $KEYS, $property) = @_;
    my @result;
    map { push(@result, $T2H->{$_}->{$property}) } @$KEYS;
    return \@result ;
  }


# euclidian distance beetween 2 points
# x1, y1, z1, x2, y2, z2
sub distance($$$$$$){
  my ($x_1, $y_1, $z_1, $x_2, $y_2, $z_2) = @_ ;
  my $dist=sqrt(($x_1 - $x_2)**(2) + ($y_1-$y_2)**(2) + ($z_1-$z_2)**(2));
  return $dist;
}

sub fact($){
  my $x = shift;
  my $tot=1;
  while($x>1){
    $tot *= $x;
    $x--;
  }
  return $tot;
}

# Takes a matrix of connectivity
# Returns 2 references:
#  1 - table of numbers correponding to the size of the clusters found
#  2 - numbers corresponding to these clusters
#
sub singleLinkage($$){

  my ($array, $cutOff) = @_ ;
  my ($j, $seed);
  my (%passed, %passed_loc) ;
  my (@list,@groups, $tmp);

  # ----------------------------------------------------------
  # -----------  While we don't have all nodes, --------------
  # ----------------------------------------------------------
  while( (keys %passed) != ($#{$array}+1) ){

    %passed_loc=();

    # --------- Look for a node we haven't yet ---
    $j=0;
    while(exists $passed{$j}){
      $j++;

      if($j > $#{$array}){
	die("Problem");
      }
    }
    @list = ($j);

    # --------- Look for the cluster ------------
    while( $#list > -1){
      $seed = shift(@list);
      $passed{$seed}=1;
      $passed_loc{$seed}=1;

      foreach $j (0..$#{$array}){
	if($array->[$seed][$j] >= $cutOff && !exists $passed{$j} ){
	  push(@list,$j) ;
	}
      }
    }
    $tmp = keys %passed_loc;
    push(@groups, $tmp);
    # -------------------------------------------
  }
  return \@groups;
}

# Same except that returns a list of labels insteat of the size of the groups
# e.g. ['AB','CD'];
sub singleLinkage_LAB($$$){

  my ($array, $labels, $cutOff) = @_ ;
  my ($j, $seed);
  my (%passed, %passed_loc) ;
  my (@list,@groups, $tmp, %list_lab, @groups_lab, @tmp);

  # ----------------------------------------------------------
  # -----------  While we don't have all nodes, --------------
  # ----------------------------------------------------------
  while( (keys %passed) != ($#{$array}+1) ){

    %passed_loc=();
    %list_lab = ();

    # --------- Look for a node we haven't yet ---
    $j=0;
    while(exists $passed{$j}){
      $j++;

      if($j > $#{$array}){
	die("Problem");
      }
    }
    @list = ($j);

    # --------- Look for the cluster ------------
    while( $#list > -1){
      $seed = shift(@list);
      $passed{$seed}=1;
      $passed_loc{$seed}=1;
      $list_lab{$labels->[$seed]}=1;
      foreach $j (0..$#{$array}){
	if($array->[$seed][$j] >= $cutOff && !exists $passed{$j} ){
	  push(@list,$j) ;
	}
      }
    }
    $tmp = keys %passed_loc;
    push(@groups, $tmp);
    @tmp = keys %list_lab;
    push(@groups_lab, join("",@tmp));
    # -------------------------------------------
  }
  return \@groups_lab;
}

# Take a hash as argument instead of a matrix:
#
# --> H -> {A} -> {B} = X   MEANS A is linked to B with link of weight X
#
#
# Returns a list of labels
# e.g. ['AB','CD'];
sub singleLinkage_LAB_SPARSE($$$){

  my ($sparse, $labels, $cutOff) = @_ ;
  my ($j, $seed);
  my (%passed) ;
  my (@list,@groups, $tmp, %list_lab, @groups_lab, @tmp);

  # ----------------------------------------------------------
  # -----------  While we don't have all nodes, --------------
  # ----------------------------------------------------------
  while( (keys %passed) != ($#{$labels}+1) ){

    %list_lab = ();

    # --------- Look for a node we haven't yet ---
    $j=0;
    while(exists $passed{$labels->[$j]}){
      $j++;

      if($j > $#{$labels}){
	die("Problem");
      }
    }
    @list = ($labels->[$j]);

    # --------- Look for the cluster ------------
    while( $#list > -1){
      $seed = shift(@list);
      $passed{$seed}=1;
      $list_lab{$seed}=1;

      foreach $j (0..$#{$labels}){

	if(!exists($passed{$labels->[$j]}) && exists($sparse->{$seed}) && exists($sparse->{$seed}->{$labels->[$j]}) &&  ($sparse->{$seed}->{$labels->[$j]} >= $cutOff) ){
	  push(@list,$labels->[$j]) ;
	}
      }
    }
    @tmp = keys %list_lab;
    push(@groups_lab, join(",",@tmp));
    # -------------------------------------------
  }
  return \@groups_lab;
}


##
## Takes an array, return an array of couples,
## The 1st elemts of the couples contain all the different elements
## The 2d contain the number of occurences
sub abundances($)
  {
    my ($array) = @_;

    my (%elemts, @elemts);

    map
      {
	if( exists($elemts{$_}) )
	  {
	    $elemts{$_}++
	  }
	else
	  {
	    $elemts{$_}=1
	  }
      } @$array ;

    map { push(@elemts, [$_, $elemts{$_}] )  } (keys %elemts) ;

    @elemts = sort {$b->[1] <=> $a->[1] } @elemts ;

    return \@elemts ;
  }

##
## Takes an array of scalars,
## Returns it but without the redundancy.
##
sub uniq($)
  {
    my ($array) = @_;

    my (%elemts, @elemts);

    map
      {
	  if( exists($elemts{$_}) )
	  {
	    $elemts{$_}++
	  }
	else
	  {
	    $elemts{$_}=1
	  }
      } @$array ;

    @elemts = keys %elemts;

    @elemts = sort {$a <=> $b } @elemts;

    return \@elemts ;
  }

##
## Load list
##
sub load_list($){

  my ($file) = @_;
  open(IN, "<$file");
  my @list;
  while(<IN>){
    chomp($_);
    push(@list, $_);
  }
  return(\@list);
}

##
## Takes an array, a value, if it finds it then it return the place, otherwise returns -1.
##
sub ar_find($$){

  my ($ar, $to_find) = @_;
  my ($elmt,$i);
  $i=0;
  foreach $elmt (@$ar){
    if($elmt eq $to_find){
      return $i;
    }
    $i++;
  }
  return -1;
}

##
## Takes an array, a bunch number and a bunch size, and returns a sub-array with the bunch
##
sub get_bunch($$$){

    my ($array, $num, $N_BUNCH) = @_;

    my @to_proceed;

    my $BUNCH= int( ($#{$array}+1)/$N_BUNCH + 0.99 );

    if($num*$BUNCH < $#{$array}){
	@to_proceed = @$array[(($num-1)*$BUNCH)..(($num*$BUNCH)-1)];

    } elsif(($num-1)*$BUNCH > $#{$array}){
	die("The number provided ($num) is too big: not enough in the array !!");

    } else {

	@to_proceed = @$array[(($num-1)*$BUNCH)..$#{$array}];
    }

    return(\@to_proceed);
}

sub fisher_yates_shuffle {
    my $deck = shift; # $deck is a reference to an array
    my $i = @$deck;
    while (--$i) {
	my $j = int rand ($i+1);
	@$deck[$i,$j] = @$deck[$j,$i];
    }
    return($deck);
}

###
### Takes 2 hashes defining a sparse matrix and returns a HASH defining BI-directional best hits.
### as well as non defined values.
###
### HASH1 = ROWS
### HASH2 = COLS
###
sub BDBH($$){

    my ($H1, $H2) = @_ ;
    
    #print "H1=".Dumper($H1);
    #print "H2=".Dumper($H2);

    my $H1_rank = H_2_rank($H1);
    my $H2_rank = H_2_rank($H2);

    #print "H1 rank=".Dumper($H1_rank);
    #print "H2 rank=".Dumper($H2_rank);

    my $H_bdbh;
    my $H_bdbh_rev;

    foreach my $ch1 (keys %$H1_rank) {

	foreach my $ch2 (keys %{$H1_rank->{$ch1}}) {
	    
	    if ($H1_rank->{$ch1}->{$ch2} == 1 && $H2_rank->{$ch2}->{$ch1} == 1 ){
		
		$H_bdbh->{$ch1}=$ch2;
		$H_bdbh_rev->{$ch2}=$ch1;
	    
	    } else {

		#$H_bdbh->{$ch1}->{$ch2}=0;
	
	    }
	}
    }

    return ($H_bdbh,$H_bdbh_rev);
}

###
### Takes a hash of the form H ->{A}->{v1}->score1
###                                \->{v2}->score2
###                           \->{B}->{v1}->score1 etc ...
###
### Transforms the scores into ranks.
###
sub H_2_rank($){

    my ($Hs) = @_;

    my $H_res;

    foreach my $k1 (keys %$Hs){

	my @keys2 = keys %{$Hs->{$k1}};
	my @pairs = ();

	foreach my $k2 (@keys2){

	    push(@pairs, [$k2, $Hs->{$k1}->{$k2}]);
	}

	my @sorted_pairs = sort { $b->[1] <=> $a->[1] } @pairs;

	my $rank=1;

	foreach my $pair (@sorted_pairs){

	    $H_res->{$k1}->{$pair->[0]} = $rank;
	    $rank++;
	}
    }
    return $H_res;
}


1;



#foreach $fam_num (@NR2){

#  if($i++ % 4 == 0){

 #   print WEB "</tr><tr>\n";
 #   print WEB $image_string ;
 #   $image_string="";
 #   print WEB "</tr><tr>\n";
 # }

  # print PDB CODE + description
#  print WEB "<th><font size=3>".$fam_num->[0].", ".$fam_num->[1]." same topol </th>\n";
#  $image_string .= "<th> <img src=\"/gn23/elevy/pqs/results/gviz/TYPES/".$fam_num->[0].".jpg\"></img> </th> \n";
#}
#print WEB "</tr><tr>\n";
#print WEB $image_string ;
#print WEB "</table>\n";
#print WEB HtmlBot();
#close(WEB);

#print join("\n",@NR) ;
