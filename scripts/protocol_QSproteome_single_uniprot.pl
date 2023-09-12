use strict;
use warnings;
use File::Basename;
use Cwd;
use Getopt::Long;
use File::Path;

# Get the name of the script
my $script_name = $0;

# Get the full path to the script
my $script_path = Cwd::abs_path($script_name);
my $script_dir = dirname($script_path);

require("$script_dir/functions_get_contacts.pm");

my $PDBFILE;
my $JSON;
my $OUTPATH;
my $reconstruct;

# Define the command-line options
GetOptions(
    'pdbfile=s'    => \$PDBFILE,    # Input file (string)
    'json=s'   => \$JSON,   # Output file (string)
    'outpath=s'    => \$OUTPATH,       # Verbose mode (flag)
    'reconstruct!' => \$reconstruct,
);

# Check for required options
unless (defined $PDBFILE) {
    die "Error: Input AlphaFold model pdb file is required. Use --pdb <file>\n";
}

unless (defined $JSON) {
    die "Error: Input Alphafold model json file is required. Use --json <file>\n";
}

unless (defined $OUTPATH) {
    die "Error: Path where output files will be written is required. Use --outpath <file>\n";
}


# Process the options
print "pdb file: $PDBFILE\n" if $PDBFILE;
print "json file: $JSON\n" if $JSON;
print "out path: $OUTPATH\n" if $OUTPATH;

# Extract pdb id
my $CODE = (fileparse($PDBFILE, qr/\.[^.]*/))[0];
print "code: $CODE\n";

# Check if the output directory exists
if (!-d $OUTPATH) {
    # Create the directory
    eval {
        mkpath($OUTPATH);
    };

    if ($@) {
        die "Failed to create directory: $@";
    } else {
        print "Directory created successfully.\n";
    }
} else {
    print "Directory already exists.\n";
}

### STEP 1
## Calculate the contacts in the model
## Calculate clashes
my @PDBFILES = ($PDBFILE);             
process_contacts(@PDBFILES, "$OUTPATH/$CODE", 1);

#### STEP 2
### Calculate nodiso1, 2, 3 
my $CONTACTFILE = "$OUTPATH/$CODE"."_FULL_CONTACT.txt";
#print "PDBFILE: $PDBFILE\n";
#print "JSON: $JSON\n";
#print "CONTACTFILE: $CONTACTFILE\n";
#print "OUTPATH: $OUTPATH\n";
system("Rscript $script_dir/process_and_analyze_AF_model.R $PDBFILE $JSON $CONTACTFILE $OUTPATH");

# Test if the option is specified
if ($reconstruct) {
    print "The user specified the --reconstruct option.\n";
    system("bash $script_dir/launch_reconstruct_ananas_cN.sh $OUTPATH/${CODE}_nodiso3.pdb $OUTPATH");
    #print("$OUTPATH/${CODE}_nodiso3_all_csym.dat\n");
    system("Rscript $script_dir/select_best_rmsd_clashes_byfile.R $OUTPATH/${CODE}_nodiso3_all_csym.dat");
} else {
    print "The user did not specify the --reconstruct option => full size complex not reconstructed\n";
}
