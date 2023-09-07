use strict;
use warnings;
use File::Basename;
use Cwd;
use Getopt::Long;


# Get the name of the script
my $script_name = $0;

# Get the full path to the script
my $script_path = Cwd::abs_path($script_name);
my $script_dir = dirname($script_path);

require("$script_dir/functions_get_contacts.pm");

my $PDBFILE;
my $JSON;
my $OUTPATH;

# Define the command-line options
GetOptions(
    'pdbfile=s'    => \$PDBFILE,    # Input file (string)
    'json=s'   => \$JSON,   # Output file (string)
    'outpath=s'    => \$OUTPATH,       # Verbose mode (flag)
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

#my $PDBFILE = $ARGV[0];
#my $JSON = $ARGV[1];
#my $OUTPATH = $ARGV[2];

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
my $CONTACTFILE = "$OUTPATH/$CODE"."_FULL.txt";
print "PDBFILE: $PDBFILE\n";
print "JSON: $JSON\n";
print "CONTACTFILE: $CONTACTFILE\n";
print "OUTPATH: $OUTPATH\n";
system("Rscript process_and_analyze_AF_model.R $CODE $JSON $CONTACTFILE $OUTPATH")

