use strict;
use POSIX;
require("/media/elusers/users/hugo/15_alphafold/37_revision_Cell/general.pm");


my @aas = ("GLY", "ALA", "VAL", "LEU", "ILE", "SER", "THR", "ASP", "ASN", "LYS", "GLU", "GLN", "ARG", "HIS", "PHE", "CYS", "TRP", "TYR", "MET", "PRO");
my @aa2 = ("G"  , "A",   "V",   "L",   "I",   "S",   "T",   "D",   "N",   "K",   "E",   "Q",   "R",   "H",   "F",   "C",   "W",   "Y",   "M",   "P");
my %AAS; map { $AAS{$aas[$_]}=$aa2[$_] } (0..$#aas);

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
  open(OUT_FULL, ">$out"."_FULL_CONTACT.txt");
  #open(OUT_ATOM, ">$out"."_ATOM.txt");

  my ($pdb_name, $pdb_short, $res_n1, $res_n2, @chains, $ch_i, $ch_j, $x, $y, $z, $resid_num, $contact_list_big_rough, $contact_list_big, $contact_list_small, $contact, @pdb_array, $pdbs) ;
  my ($x_1, $y_1, $z_1, $atom_1, $atom_2, $x_2, $y_2, $z_2, $VDW1, $VDW2, $res_name, $res_name2, $cut_short, %pretty_print, %seen, @Ds, $d_min, $d_max, $d_avg);
  my $counter=0 ;
  my (@resid_I_inter_m1, @resid_J_inter_m1, @resid_I_inter_m2, @resid_J_inter_m2, %tmp);

  my ($B4, $AFTER) = format_stderr($Num);
  


  #print("array pdb files: $to_proceed_file\n");

  my @to_proceed = ($to_proceed_file);
  #print("to_proceed: @to_proceed\n");

  foreach $pdb_name (@to_proceed){
      #print "pdb name: $pdb_name\n";

      #print STDERR $B4."Processing contacts for ".(sprintf "%02d", $Num)." .. ".(sprintf "%06d", $counter)."/$#{to_proceed}\r".$AFTER;

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
