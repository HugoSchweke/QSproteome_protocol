use strict ;
use Cwd;
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
