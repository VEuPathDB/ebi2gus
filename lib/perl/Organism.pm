package Organism;

use strict;

use XML::Simple;

use Data::Dumper;

sub getNcbiTaxonId { $_[0]->{_ncbi_tax_id} }
sub setNcbiTaxonId { $_[0]->{_ncbi_tax_id}  = $_[1] }

sub getAbbrev { $_[0]->{_abbrev} }
sub setAbbrev { $_[0]->{_abbrev} = $_[1] }

sub getGenomeDatabaseVersion { $_[0]->{_genome_database_version} }
sub setGenomeDatabaseVersion { $_[0]->{_genome_database_version} = $_[1] }

sub getGenomeDatabaseName { $_[0]->{_genome_database_name} }
sub setGenomeDatabaseName { $_[0]->{_genome_database_name} = $_[1] }

sub getChromosomeMap { $_[0]->{_chromosome_map} }
sub setChromosomeMap { $_[0]->{_chromosome_map} = $_[1] }


sub new {
    my ($class, $ncbiTaxId, $orgAbbrev, $genomeDatabaseName, $genomeDatabaseVersion, $chromosomeMapFile) = @_;

    unless(defined $ncbiTaxId) {
	die "Required ncbi_taxon_id missing from organism xml";
    }

    unless($orgAbbrev) {
	die "Required abbrev missing from organism xml";
    }

    unless($genomeDatabaseName && $genomeDatabaseVersion) {
	die "GenomeDatabase requires name and version. Found [$genomeDatabaseName] and [$genomeDatabaseVersion]"
    }

    open(MAP, $chromosomeMapFile) or die "Cannot open $chromosomeMapFile for reading: $!";

    my %chromosomeMap;
    while(<MAP>) {
	chomp;

	my @a = split(/\t/, $_);
	$chromosomeMap{$a[0]} = {chromosome => $a[1], chromosome_order_num => $a[2]};
    }
    close MAP;
    
    my $self = bless {}, $class;

    $self->setNcbiTaxonId($ncbiTaxId);
    $self->setAbbrev($orgAbbrev);
    $self->setGenomeDatabaseName($genomeDatabaseName);
    $self->setGenomeDatabaseVersion($genomeDatabaseVersion);
    $self->setChromosomeMap(\%chromosomeMap);

    return $self;
}

1;