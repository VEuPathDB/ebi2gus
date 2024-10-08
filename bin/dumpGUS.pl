#!/usr/bin/env perl

use strict;
use warnings;

use Getopt::Std;

use Bio::EnsEMBL::Registry;

use GUSTableDefinitionParser;
use Organism;

use EBIParser;

use Data::Dumper;

my $REGISTRY_CONF_FILE = "/usr/local/etc/ensembl_registry.conf";
my $TABLE_DEFINITIONS_XML_FILE = "/usr/local/etc/gusSchemaDefinitions.xml";


my $OUTPUT_DIRECTORY = "/tmp";

sub HELP_MESSAGE {
    print STDERR <<"EOM";

USAGE:
$0 [-w] [-h] -e FILE -t FILE	-d STRING -v STRING -n NUMBER -c FILE -p STRING -r STRING -a STRING -o DIRECTORY

OPTIONS:

-w     skip validation
-h     print this message
-e     ensemble_registry_file
-t     table_definitions_xml_file
-d     genome_database_name
-v     genome_database_version 
-n     ncbi_tax_id 
-c     chromosome_map_file
-p     project_name
-r     project_release
-a     organism_abbrev
-o     output_directory
-g     GO ExternalDatabase Specification (GO_RSRC|someversion)
-s     SO ExternalDatabase Specification (SO_RSRC|someversion)
-l     GO Evidence ExternalDatabase Specification (GOevid_RSRC|someversion)

EOM
    exit -1;
}
our ($opt_r, $opt_a, $opt_o, $opt_t, $opt_h, $opt_v, $opt_n, $opt_c, $opt_d, $opt_p, $opt_e, $opt_g, $opt_s, $opt_l, $opt_w);
getopts('hwl:g:s:p:r:o:t:v:n:c:d:e:a:') or HELP_MESSAGE();

HELP_MESSAGE() if($opt_h);
$REGISTRY_CONF_FILE = $opt_e if($opt_e);
$TABLE_DEFINITIONS_XML_FILE = $opt_t if($opt_t);
$OUTPUT_DIRECTORY = $opt_o if($opt_o);

my $skipValidation = $opt_w;

my $ncbiTaxId = $opt_n;
my $genomeDatabaseName = $opt_d;
my $genomeDatabaseVersion = $opt_v;
my $chromosomeMapFile = $opt_c;

my $organismAbbrev = $opt_a;

my $projectName = $opt_p;
my $projectRelease = $opt_r;

my $goSpec = $opt_g;
my $goEvidSpec = $opt_l;
my $soSpec = $opt_s;

if($opt_h || 
   !-e $REGISTRY_CONF_FILE || !-e $TABLE_DEFINITIONS_XML_FILE || 
   !defined($ncbiTaxId) || !$genomeDatabaseVersion || 
   !$genomeDatabaseName || !$projectName || !$projectRelease ||
   !$goSpec || !$soSpec || !$goEvidSpec || !$organismAbbrev) { 
    HELP_MESSAGE();
}

if($chromosomeMapFile && !-e $chromosomeMapFile) {
    HELP_MESSAGE();
}

my $gusTableDefinitions = GUSTableDefinitionParser->new($TABLE_DEFINITIONS_XML_FILE);

my $organism = Organism->new($ncbiTaxId, $genomeDatabaseName, $genomeDatabaseVersion, $chromosomeMapFile, $organismAbbrev);

my $registry = 'Bio::EnsEMBL::Registry';

my $count = $registry->load_all($REGISTRY_CONF_FILE, 1);

my $sliceAdaptor = $registry->get_adaptor('default', 'Core', 'Slice' );

my (@slices, %seenSeqRegions);
foreach my $topLevelSlice (@{$sliceAdaptor->fetch_all('toplevel')}) {
    my $seqRegionName = $topLevelSlice->seq_region_name();
    
    unless($seenSeqRegions{$seqRegionName}) {
	push @slices, $sliceAdaptor->fetch_by_region('toplevel', $seqRegionName);
    }
    $seenSeqRegions{$seqRegionName}++;
}


my $ebiParser = EBIParser->new($sliceAdaptor, \@slices, $gusTableDefinitions, $OUTPUT_DIRECTORY, $organism, $registry, $projectName, $projectRelease, $goSpec, $soSpec, $goEvidSpec, $skipValidation);
$ebiParser->parse();
# exit;


# my $goa = $registry->get_adaptor( 'default', 'Core', 'OntologyTerm' );
# print Dumper $goa;
# exit;


# my $analysisAdaptor = $registry->get_adaptor('default', 'Core', 'Analysis' );
# my $analysis = $analysisAdaptor->fetch_by_logic_name("interpro2go");
# print Dumper $analysis;      
# exit;


# my $analysis = $analysisAdaptor->fetch_by_logic_name('superfamily');

# print Dumper $analysis;


#TODO  sres.ontology ??
  # keep a running list of seen so_terms.  Dump them out at the end to "SRes.OntologyTerm" file




1;



