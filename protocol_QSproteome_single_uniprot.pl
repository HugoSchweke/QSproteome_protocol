use File::Basename;
use strict;
use warnings;

my $PDBFILE = $ARGV[0];
my $JSON = $ARGV[1];
my $OUTPATH = $ARGV[2];

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

require("/media/elusers/users/hugo/15_alphafold/37_revision_Cell/QSproteome_protocol/functions_get_contacts.pm");

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

