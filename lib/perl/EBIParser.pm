package EBIParser;

use strict;

use GUSTableWriter;
use OutputFile;

use Data::Dumper;

use Bio::SeqIO;
use Bio::Seq;

use GUS::SRes::OntologyTerm qw(%seenOntologyTerms);
use GUS::SRes::ExternalDatabase qw(%seenExternalDatabases);
use GUS::SRes::ExternalDatabaseRelease qw(%seenExternalDatabaseReleases);
use GUS::SRes::DbRef qw(%seenDBRefs);
use GUS::SRes::EnzymeClass qw(%seenEnzymeClasses);
use GUS::Core::DatabaseInfo qw(%seenDatabases);
use GUS::Core::TableInfo qw(%seenTables);
use GUS::DoTS::GOAssociationInstanceLOE qw(%seenGOEvidences);
use GUS::DoTS::GOAssociation qw(%seenGOAssociations);

my %INTERPRO_LOGICS = ('pfam' => 1,
		       'pirsf' => 1,
		       'cdd' => 1,
		       'hamap' => 1,
		       'hmmpanther' => 1,
		       'prints' => 1,
		       'scanprosite' => 1,
		       'sfld' => 1,
		       'smart' => 1,
		       'superfamily' => 1,
		       'tigrfam' => 1,
		       'pfscan' => 1,
    );

my %SKIP_LOGICS = ('mobidblite' => 1,
		   'ncoils' => 1,
		   'blastprodom' => 1,
		   'gene3d' => 1,
    );

my $INTERPRO2GO_LOGIC = "interpro2go";

# these are required objects
# a gus table definition object will be made for each of them
sub getTables {
    return ([
'GUS::DoTS::GeneFeature',
	     'GUS::DoTS::ExternalNASequence',
	     'GUS::DoTS::NALocation',
	     'GUS::DoTS::Transcript',
	     'GUS::DoTS::RNAFeature',
	     'GUS::DoTS::SplicedNASequence',
	     'GUS::DoTS::ExonFeature',
	     'GUS::DoTS::RNAFeatureExon',
	     'GUS::DoTS::TranslatedAASequence',	     
 	     'GUS::DoTS::TranslatedAAFeature',
	     'GUS::SRes::ExternalDatabase',
	     'GUS::SRes::ExternalDatabaseRelease',
	     'GUS::DoTS::AAFeatureExon',
	     'GUS::SRes::Taxon',
	     'GUS::SRes::OntologyTerm',
	     'GUS::SRes::EnzymeClass',	
	     'GUS::DoTS::AASequenceEnzymeClass',
             'GUS::ApiDB::SeqEdit',
             'GUS::ApiDB::AaSequenceAttribute',
             'GUS::ApiDB::GeneFeatureName',
	     'GUS::DoTS::LowComplexityNAFeature',
	     'GUS::DoTS::TransposableElement',
	     'GUS::DoTS::TandemRepeatFeature',
	     'GUS::DoTS::Repeats',
	     'GUS::DoTS::SignalPeptideFeature',
 	     'GUS::DoTS::TransMembraneAAFeature',
 	     'GUS::DoTS::LowComplexityAAFeature',
	     'GUS::DoTS::AALocation',
	     'GUS::SRes::DbRef',
	     'GUS::DoTS::DbRefAAFeature',
	     'GUS::DoTS::DbRefNAFeature',
             'GUS::DoTS::DbRefNASequence',
	     'GUS::DoTS::DomainFeature',
	     'GUS::Core::ProjectInfo',
	     'GUS::Core::DatabaseInfo',
	     'GUS::Core::TableInfo',
	     'GUS::DoTS::GOAssociation',
	     'GUS::DoTS::GOAssociationInstance',
	     'GUS::DoTS::GOAssociationInstanceLOE',
             'GUS::DoTS::GOAssocInstEvidCode',
	     'GUS::DoTS::Miscellaneous',
#	     'GUS::DoTS::NAComment',
	    ]);
}

#-----------------------------------------------------------
#
sub getRepeatMaskedIO { $_[0]->{_repeat_masked_io} }
sub setRepeatMaskedIO { $_[0]->{_repeat_masked_io} = $_[1] }

sub getRepeatsIO { $_[0]->{_repeats_io} }
sub setRepeatsIO { $_[0]->{_repeats_io} = $_[1] }

sub getTrfIO { $_[0]->{_trf_io} }
sub setTrfIO { $_[0]->{_trf_io} = $_[1] }

sub getDustIO { $_[0]->{_dust_io} }
sub setDustIO { $_[0]->{_dust_io} = $_[1] }

sub getSegIO { $_[0]->{_seg_io} }
sub setSegIO { $_[0]->{_seg_io} = $_[1] }


sub getRepeatMaskedFile { $_[0]->{_repeat_masked_file} }
sub setRepeatMaskedFile { $_[0]->{_repeat_masked_file} = $_[1] }

sub getRepeatsFile { $_[0]->{_repeats_file} }
sub setRepeatsFile { $_[0]->{_repeats_file} = $_[1] }

sub getTrfFile { $_[0]->{_trf_file} }
sub setTrfFile { $_[0]->{_trf_file} = $_[1] }

sub getDustFile { $_[0]->{_dust_file} }
sub setDustFile { $_[0]->{_dust_file} = $_[1] }

sub getSegFile { $_[0]->{_seg_file} }
sub setSegFile { $_[0]->{_seg_file} = $_[1] }



#-----------------------------------------------------------

sub getSlices { $_[0]->{_slices} }
sub setSlices { $_[0]->{_slices} = $_[1] }

sub getOrganism { $_[0]->{_organism} }
sub setOrganism { $_[0]->{_organism} = $_[1] }

sub getGOSpec { $_[0]->{_go_spec} }
sub setGOSpec { $_[0]->{_go_spec} = $_[1] }

sub getGOEvidSpec { $_[0]->{_go_evid_spec} }
sub setGOEvidSpec { $_[0]->{_go_evid_spec} = $_[1] }

sub getSOSpec { $_[0]->{_so_spec} }
sub setSOSpec { $_[0]->{_so_spec} = $_[1] }

sub getGOExtDbRlsId { $_[0]->{_go_ext_db_rls_id} }
sub setGOExtDbRlsId { $_[0]->{_go_ext_db_rls_id} = $_[1] }

sub getGOEvidExtDbRlsId { $_[0]->{_go_evid_ext_db_rls_id} }
sub setGOEvidExtDbRlsId { $_[0]->{_go_evid_ext_db_rls_id} = $_[1] }


sub getSOExtDbRlsId { $_[0]->{_so_ext_db_rls_id} }
sub setSOExtDbRlsId { $_[0]->{_so_ext_db_rls_id} = $_[1] }

sub getProjectName { $_[0]->{_project_name} }
sub setProjectName { $_[0]->{_project_name} = $_[1] }

sub getProjectRelease { $_[0]->{_project_release} }
sub setProjectRelease { $_[0]->{_project_release} = $_[1] }

sub getRegistry { $_[0]->{_registry} }
sub setRegistry { $_[0]->{_registry} = $_[1] }

# sub getGlobalSeqRegionMappings { $_[0]->{_global_seq_region_mappings} }
# sub setGlobalSeqRegionMappings {
#     my ($self) = @_;

#     # THIS IS MECHANISM FOR MAPPING SEQ REGION NAMES IS TEMPORARY
#     my $genomeSummary = "Arthropod genome summary.csv";

#     my $organismAbbrev = $self->getOrganism()->getOrganismAbbrev();
#     my $seqRegionDirectory = "/usr/local/etc/seq_region_maps";

#     my $genomeSummaryFile = $seqRegionDirectory . "/" . $genomeSummary;
#     open(SUMMARY, $genomeSummaryFile) or die "Could not open file $genomeSummaryFile for reading: $!";

#     my $rv = {};
#     my $foundRow;
#     <SUMMARY>; 
#     while(<SUMMARY>) {
# 	chomp;
# 	my ($o, $hasMapping, $fn) = split(/,/, $_);

# 	if($o eq $organismAbbrev) {
# 	    $foundRow = 1;

# 	    if($hasMapping) {
# 		my $fileName = $seqRegionDirectory . "/" . $fn;
# 		$fileName =~ s/\.xlsx/.csv/;
		
# 		open(FILE, $fileName) or die "Cannot open file $fileName for reading: $!";

# 		<FILE>;
# 		while(<FILE>) {
# 		    chomp;
# 		    my ($seqRegionName, $newSeqRegionName) = split(/,/, $_);
# 		    $rv->{$seqRegionName} = $newSeqRegionName;
# 		}
# 		close FILE;
# 	    }
# 	}
# 	last if $foundRow;
#     }
#     close SUMMARY;

#     $self->{_global_seq_region_mappings} = $rv;
# }

sub addRowToGenomicBedFile {
    my ($self, $sequenceSourceId, $feature, $fh) = @_;

    print $fh ($sequenceSourceId,
              "\t", $feature->seq_region_start(),
              "\t", $feature->seq_region_end(),
              "\n");
}

sub addRowToProteinBedFile {
    my ($self, $sequenceSourceId, $feature, $fh) = @_;

    print $fh ($sequenceSourceId,
              "\t", $feature->start(),
              "\t", $feature->end(),
              "\n");
}


sub finishBedFiles {
    my ($self) = @_;

    foreach($self->getRepeatsIO(),
            $self->getTrfIO(),
            $self->getDustIO(),
            $self->getSegIO()) {
        close $_;
    }

    foreach($self->getRepeatsFile(),
            $self->getTrfFile(),
            $self->getDustFile(),
            $self->getSegFile()) {
        $self->bgzipAndTabix($_);
    }
}


sub bgzipAndTabix {
    my ($self, $file) = @_;

    system("sort -k1,1 -k2,2n $file -o $file") == 0
        or die "sort failed: $?";

    system("bgzip $file") == 0
        or die "bgzip failed: $?";

    system("tabix", "-p", "bed", "${file}.gz") == 0
        or die "tabix failed: $?";
}



sub getPreviousIdentifiersFromPatchBuild { $_[0]->{_previous_identifiers_from_patch_build} }
sub setPreviousIdentifiersFromPatchBuild {
    my ($self, $sliceAdaptor) = @_;

    my $dbc = $sliceAdaptor->dbc();
    my $sth = $dbc->prepare("select old_stable_id,new_stable_id from stable_id_event where type = 'gene' and mapping_session_id in (select max(mapping_session_id) from mapping_session)");
    $sth->execute();
    
    my $rv = {};
    while(my ($oldStableId, $stableId) = $sth->fetchrow_array()) {
	push @{$rv->{$stableId}}, $oldStableId if ($stableId && $oldStableId && $stableId ne $oldStableId);
    }
    $sth->finish();

    $self->{_previous_identifiers_from_patch_build} = $rv;
}

sub getGlobalSeqRegionMappings { $_[0]->{_global_seq_region_mappings} }
sub setGlobalSeqRegionMappings {
    my ($self) = @_;

    # THIS IS MECHANISM FOR MAPPING SEQ REGION NAMES IS TEMPORARY
    my $genomeSummary = "Arthropod genome summary.csv";

    my $organismAbbrev = $self->getOrganism()->getOrganismAbbrev();
    my $seqRegionDirectory = "/usr/local/etc/seq_region_maps";

    my $genomeSummaryFile = $seqRegionDirectory . "/" . $genomeSummary;
    open(SUMMARY, $genomeSummaryFile) or die "Could not open file $genomeSummaryFile for reading: $!";

    my $rv = {};
    my $foundRow;
    <SUMMARY>; 
    while(<SUMMARY>) {
	chomp;
	my ($o, $hasMapping, $fn) = split(/,/, $_);

	if($o eq $organismAbbrev) {
	    $foundRow = 1;

	    if($hasMapping) {
		my $fileName = $seqRegionDirectory . "/" . $fn;
		$fileName =~ s/\.xlsx/.csv/;
		
		open(FILE, $fileName) or die "Cannot open file $fileName for reading: $!";

		<FILE>;
		while(<FILE>) {
		    chomp;
		    my ($seqRegionName, $newSeqRegionName) = split(/,/, $_);
		    $rv->{$seqRegionName} = $newSeqRegionName;
		}
		close FILE;
	    }
	}
	last if $foundRow;
    }
    close SUMMARY;

    $self->{_global_seq_region_mappings} = $rv;
}

sub getGUSTableWriters { $_[0]->{_gus_table_writers} }
sub setGUSTableWriters { 
    my ($self, $gusTableDefinitionsParser, $outputDirectory, $skipValidation) = @_;

    my $tables = $self->getTables();

    my $gusTableWriters = {};

    my %outputFiles;
    
    my %allFields;

    foreach my $className (@$tables) {
	$className =~ s/^GUS:://;
	$className =~ s/::/./;

	$className = uc $className;
	
	my $gusTableDefinition = $gusTableDefinitionsParser->makeTableDefinition($className);
	my $realTableName = $gusTableDefinition->getRealTableName();

	my @impFields = $gusTableDefinition->isView() ? keys %{$gusTableDefinition->getImpToViewFieldMap()} : @{$gusTableDefinition->getFields()};

	foreach(@impFields) {
	    $allFields{$realTableName}{$_}++;
	}

	my $outputFile; #Only one fileName/FileHandle/Counter Per RealTableName
	if($outputFiles{$realTableName}) {
	    $outputFile = $outputFiles{$realTableName};
	}
	else {
#	    my $headerFields = $gusTableDefinition->getFields();

	    $outputFile = OutputFile->new($realTableName, $outputDirectory);
	    $outputFiles{$realTableName} = $outputFile;
	}
	
	$gusTableWriters->{$className} = GUSTableWriter->new($gusTableDefinition, $outputFile, $skipValidation);
    }

    foreach my $realTableName (keys %outputFiles) {
	my $outputFile = $outputFiles{$realTableName};
	my $impFields = $allFields{$realTableName};
	my @headerFields = keys %$impFields;
	$outputFile->setHeaderFields(\@headerFields);

	$outputFile->writeHeader();
    }

    $self->{_gus_table_writers} = $gusTableWriters;
}

sub new {
    my ($class, $sliceAdaptor, $slices, $gusTableDefinitions, $outputDirectory, $organism, $registry, $projectName, $projectRelease, $goSpec, $soSpec, $goEvidSpec, $skipValidation) = @_;
    
    my $self = bless {}, $class;

    $self->setSlices($slices);

    $self->setGUSTableWriters($gusTableDefinitions, $outputDirectory, $skipValidation);
    $self->setOrganism($organism);

    $self->setProjectName($projectName);
    $self->setProjectRelease($projectRelease);

    $self->setRegistry($registry);
    
    $self->importTableModules();

    $self->setGOSpec($goSpec);
    $self->setGOEvidSpec($goEvidSpec);
    $self->setSOSpec($soSpec);

    my $repeatMaskedFile = "$outputDirectory/blocked.seq";

    my $repeatMaskedIO = Bio::SeqIO->new(-file   => ">$repeatMaskedFile",
					 -format => 'fasta' );

    my ($repeatsIO, $trfIO, $dustIO, $segIO);
    my $repeatsFile = "$outputDirectory/repeatmask.bed";
    my $trfFile = "$outputDirectory/trf.bed";
    my $dustFile = "$outputDirectory/dust.bed";
    my $segFile = "$outputDirectory/seg.bed";

    open($repeatsIO, ">$repeatsFile") or die "Cannot open file $repeatsFile for writing: $!";
    open($trfIO, ">$trfFile") or die "Cannot open file $trfFile for writing: $!";
    open($dustIO, ">$dustFile") or die "Cannot open file $dustFile for writing: $!";
    open($segIO, ">$segFile") or die "Cannot open file $segFile for writing: $!";

    $self->setRepeatMaskedIO($repeatMaskedIO); # this one is fasta
    $self->setRepeatsIO($repeatsIO);
    $self->setTrfIO($trfIO);
    $self->setDustIO($dustIO);
    $self->setSegIO($segIO);


    $self->setRepeatMaskedFile($repeatMaskedFile); # this one is fasta
    $self->setRepeatsFile($repeatsFile);
    $self->setTrfFile($trfFile);
    $self->setDustFile($dustFile);
    $self->setSegFile($segFile);


    $self->setGlobalSeqRegionMappings();

    $self->setPreviousIdentifiersFromPatchBuild($sliceAdaptor);

    return $self;
}



sub importTableModules {
    my ($self) = @_;

    my $tables = $self->getTables();

    foreach(@$tables) {
	eval "require $_";
	if($@) {
	    die $@;
	}
    }
}


sub getExternalDatabaseRelease {
    my ($self, $gusTableWriters, $organism) = @_;

    my $genomeDatabaseName = $organism->getGenomeDatabaseName();
    my $genomeDatabaseVersion = $organism->getGenomeDatabaseVersion();
    
    my $gusExternalDatabase = GUS::SRes::ExternalDatabase->new($gusTableWriters, $genomeDatabaseName);
    my $gusExternalDatabaseRelease = GUS::SRes::ExternalDatabaseRelease->new($gusTableWriters, $genomeDatabaseVersion, $gusExternalDatabase->getPrimaryKey());

    return $gusExternalDatabaseRelease;
}


sub ontologyTermForSlice {
    my ($self, $slice, $gusTableWriters) = @_;

    # fallback to coord_system_name
    my $name = $slice->coord_system_name();

    # location contains an SO term if the sequence is not nuclear e.g. "apicoplast_chromosome", "mitochondrial_chromosome"
    my $location = $slice->get_all_Attributes('sequence_location');
    if ($location && scalar @$location == 1) {
      $name = $location->[0]->value(); 
    }else{
     # location would be undef for nuclear sequences. 
     my $coordSystemTags = $slice->get_all_Attributes("coord_system_tag");
     if($coordSystemTags && scalar @$coordSystemTags == 1) {
        $name = $coordSystemTags->[0]->value();
     }
    }
    return $self->ontologyTermFromName($name, $gusTableWriters);
}

sub ontologyTermFromBiotype {
    my ($self, $biotype, $gusTableWriters) = @_;

    my $name = $biotype->name();
    my $sourceId = $biotype->so_acc();
    my $soExtDbRlsId =  $self->getSOExtDbRlsId();
    
    if($seenOntologyTerms{"$sourceId|$soExtDbRlsId"}) {
	return $seenOntologyTerms{"$sourceId|$soExtDbRlsId"};
    }

    return GUS::SRes::OntologyTerm->new($gusTableWriters, $sourceId, $name, $soExtDbRlsId)->getPrimaryKey();
}

sub ontologyTermFromName {
    my ($self, $name, $gusTableWriters) = @_;

    my $soExtDbRlsId =  $self->getSOExtDbRlsId();

    if($seenOntologyTerms{"$name|$soExtDbRlsId"}) {
	return $seenOntologyTerms{"$name|$soExtDbRlsId"};
    }

    return GUS::SRes::OntologyTerm->new($gusTableWriters, undef, $name, $soExtDbRlsId)->getPrimaryKey();
}


sub ontologyTermFromEvidenceCode {
    my ($self, $evidenceCode, $gusTableWriters) = @_;

    my $goEvidExtDbRlsId =  $self->getGOEvidExtDbRlsId();
    
    if($seenOntologyTerms{"$evidenceCode|$goEvidExtDbRlsId"}) {
	return $seenOntologyTerms{"$evidenceCode|$goEvidExtDbRlsId"};
    }

    return GUS::SRes::OntologyTerm->new($gusTableWriters, $evidenceCode, $evidenceCode, $goEvidExtDbRlsId)->getPrimaryKey();
}


sub ontologyTermFromGOTerm {
    my ($self, $goTerm, $gusTableWriters) = @_;

    $goTerm =~ s/:/_/;

    my $goExtDbRlsId =  $self->getGOExtDbRlsId();
    
    if($seenOntologyTerms{"$goTerm|$goExtDbRlsId"}) {
	return $seenOntologyTerms{"$goTerm|$goExtDbRlsId"};
    }

    return GUS::SRes::OntologyTerm->new($gusTableWriters, $goTerm, $goTerm, $goExtDbRlsId)->getPrimaryKey();
}


sub dumpRepeatMaskedSeq {
    my ($self, $slice, $externalNaSeq) = @_;

    my $io = $self->getRepeatMaskedIO();

    my $sequenceSourceId = $externalNaSeq->getGUSRowAsHash()->{source_id};

    my $repeatMaskedSlice = $slice->get_repeatmasked_seq();
    $repeatMaskedSlice->soft_mask(1);

    my $repeatMaskedSeq = $repeatMaskedSlice->seq();

    my $seq = Bio::Seq->new(-display_id => $sequenceSourceId,
			    -seq => $repeatMaskedSeq);
    
    $io->write_seq($seq);
}


sub parseSlice {
    my ($self, $slice, $gusExternalDatabaseRelease, $gusTaxon) = @_;

    my $gusTableWriters = $self->getGUSTableWriters();

    my $gusSequenceOntologyId = $self->ontologyTermForSlice($slice, $gusTableWriters);
    
    my $organismAbbrev = $self->getOrganism()->getOrganismAbbrev();

    my $insdcSynonym;

    foreach my $sliceSynonym (@{$slice->get_all_synonyms()}) {
	$insdcSynonym = $sliceSynonym->name() if($sliceSynonym->dbname() eq "INSDC");
    }


    my $seqRegionMap = $self->getGlobalSeqRegionMappings();
    my $gusExternalNASequence = GUS::DoTS::ExternalNASequence->new($gusTableWriters, $slice, $gusTaxon, $gusExternalDatabaseRelease, $gusSequenceOntologyId, $insdcSynonym, $organismAbbrev, $seqRegionMap);

    # TODO: are there other dna align features we want?
    foreach my $dnaAlignFeature (@{$slice->get_all_DnaAlignFeatures("trnascan_align")}) {
	$self->parseTRNAFeature($dnaAlignFeature, $gusExternalNASequence);
    }
    
    $self->dumpRepeatMaskedSeq($slice, $gusExternalNASequence);

    foreach my $sliceSynonym (@{$slice->get_all_synonyms()}) {
	my $databaseName = $sliceSynonym->dbname();
	my $databaseVersion = 1;
	my $primaryId = $sliceSynonym->name();

	my ($dbRefId, $externalDatabaseReleaseId) = $self->getDbRefAndExternalDatabaseReleaseIds($databaseName, $databaseVersion, $primaryId, undef, undef);

	GUS::DoTS::DbRefNASequence->new($gusTableWriters, $dbRefId, $gusExternalNASequence->getPrimaryKey());
    }

    my %transcriptXrefsLogics;
    my %geneXrefsLogics;
    my %translationXrefsLogics;    

    foreach my $gene (@{$slice->get_all_Genes()}) {
	$self->parseGene($gene, $gusExternalDatabaseRelease, $gusTaxon, $gusExternalNASequence);
    }

    foreach my $repeatFeature (@{$slice->get_all_RepeatFeatures()} ) {
	$self->parseRepeatFeature($repeatFeature, $gusExternalDatabaseRelease, $gusExternalNASequence);
    }

    my $registry = $self->getRegistry();
    my $karyAdaptor = $registry->get_adaptor('default', 'Core', 'KaryotypeBand' );
    
    foreach my $band ( @{ $karyAdaptor->fetch_all_by_Slice($slice) } ) {
	if($band->stain() eq 'ACEN') {
	    $self->parseCentromere($band, $gusExternalNASequence, $gusExternalDatabaseRelease);
	}
    }
}

sub parseCentromere {
    my ($self, $band, $gusExternalNASequence, $gusExternalDatabaseRelease) = @_;

    my $gusTableWriters = $self->getGUSTableWriters();

    my $centromereSequenceOntologyId = $self->ontologyTermFromName('centromere', $gusTableWriters);
    
    my $name = $band->name();

    my $feature = GUS::DoTS::Miscellaneous->new($gusTableWriters, $name, $gusExternalNASequence, $gusExternalDatabaseRelease, $centromereSequenceOntologyId, 'centromere');
    my $gusCentromereLocation = GUS::DoTS::NALocation->new($gusTableWriters, $band, $feature);
}


sub parseTRNAFeature {
    my ($self, $tRNAFeature, $gusExternalNASequence) = @_;

    my $gusTableWriters = $self->getGUSTableWriters();
    
    my $tRNASequenceOntologyId = $self->ontologyTermFromName('tRNA', $gusTableWriters);

    my $analysis = $tRNAFeature->analysis();
    my $externalDatabaseReleaseId = $self->getExternalDatabaseReleaseFromNameVersion($analysis->program(), $analysis->program_version());    

    my $gusTRNAFeature = GUS::DoTS::RNAFeature->new($gusTableWriters, $tRNAFeature, $gusExternalNASequence, $externalDatabaseReleaseId, $tRNASequenceOntologyId);


    my $gusTRNANALocation = GUS::DoTS::NALocation->new($gusTableWriters, $tRNAFeature, $gusTRNAFeature);

    my $exonSequenceOntologyId = $self->ontologyTermFromName("exon", $gusTableWriters);
    my $eCt = 1;
    foreach my $trnaExon ( $tRNAFeature->ungapped_features() ) {
	my $trnaExonSourceId = $gusTRNAFeature->getGUSRowAsHash()->{source_id} . ".$eCt";
	my $gusTRNAExonFeature = GUS::DoTS::ExonFeature->new($gusTableWriters, $trnaExonSourceId, $gusTRNAFeature, $externalDatabaseReleaseId, $exonSequenceOntologyId);
	my $gusTRNAExonNALocation = GUS::DoTS::NALocation->new($gusTableWriters, $trnaExon, $gusTRNAExonFeature);

	$eCt++;
    }
}


sub parseRepeatFeature {
    my ($self, $repeatFeature, $gusExternalDatabaseRelease, $gusExternalNASequence) = @_;

    my $gusTableWriters = $self->getGUSTableWriters();
    my $logicName = $repeatFeature->analysis()->logic_name();

    my $gusFeature;

    my $sequenceSourceId = $gusExternalNASequence->getGUSRowAsHash()->{source_id};

    if($logicName eq "dust") {
        $self->addRowToGenomicBedFile($sequenceSourceId, $repeatFeature, $self->getDustIO());
    }
    if($logicName eq "tefam") {
        $gusFeature = GUS::DoTS::TransposableElement->new($gusTableWriters, $repeatFeature, $gusExternalNASequence, $gusExternalDatabaseRelease);
        GUS::DoTS::NALocation->new($gusTableWriters, $repeatFeature, $gusFeature);
    }
    if($logicName eq "trf") {
        $self->addRowToGenomicBedFile($sequenceSourceId, $repeatFeature, $self->getTrfIO());
    }
    if($logicName =~ /^repeatmask/) {
        $self->addRowToGenomicBedFile($sequenceSourceId, $repeatFeature, $self->getRepeatsIO());
    }
}

sub parseGene {
    my ($self, $gene, $gusExternalDatabaseRelease, $gusTaxon, $gusExternalNASequence) = @_;

    my $gusTableWriters = $self->getGUSTableWriters();

    my $isProteinCoding;
    foreach(@{ $gene->get_all_Transcripts() }) {
	if($_->translation()) {
	    $isProteinCoding = 1;
	    last;
	}
    }
    
    my $geneSequenceOntologyId = $self->ontologyTermFromBiotype($gene->get_Biotype(), $gusTableWriters);
    # TODO: product name
    my $gusGeneFeature = GUS::DoTS::GeneFeature->new($gusTableWriters, $gene, $gusExternalNASequence, $gusExternalDatabaseRelease, $geneSequenceOntologyId);
    my $gusGeneNALocation = GUS::DoTS::NALocation->new($gusTableWriters, $gene, $gusGeneFeature);

    if(my $geneName = $gene->external_name()) {
	GUS::ApiDB::GeneFeatureName->new($gusTableWriters, $geneName, $gusGeneFeature->getPrimaryKey(),0, $gusExternalDatabaseRelease->getPrimaryKey());	
    }
    
    my $taxonId = $gusTaxon->getPrimaryKey();
    
    my %exonMap;
    my $exonSequenceOntologyId = $self->ontologyTermFromName("exon", $gusTableWriters);
    foreach my $exon ( @{ $gene->get_all_Exons() } ) {
	my $gusExonFeature = GUS::DoTS::ExonFeature->new($gusTableWriters, $exon->stable_id(), $gusGeneFeature, $gusExternalDatabaseRelease->getPrimaryKey(), $exonSequenceOntologyId);
	my $gusExonNALocation = GUS::DoTS::NALocation->new($gusTableWriters, $exon, $gusExonFeature);

	my $exonKey = $exon->start()."-".$exon->end()."-".$exon->strand()."-".$exon->phase()."-".$exon->end_phase();
	$exonMap{$exonKey} = $gusExonFeature->getPrimaryKey();
    }
    
    foreach my $transcript ( @{ $gene->get_all_Transcripts() } ) {
	$self->parseTranscript($transcript, $gene, $gusGeneFeature, $taxonId, \%exonMap, $gusExternalDatabaseRelease);
    }

    foreach my $xref (@{$gene->get_all_object_xrefs()}) {
	my $databaseName = $xref->dbname();
  	
    next unless (defined $xref->analysis());
    my $databaseVersion = $xref->analysis()->logic_name();
 
	my $primaryId = $xref->primary_id();

	my @synonyms = @{$xref->get_all_synonyms()};

	push @synonyms, $primaryId;

	foreach my $synOrPrimary (@synonyms) {
	    my ($dbRefId, $externalDatabaseReleaseId) = $self->getDbRefAndExternalDatabaseReleaseIds($databaseName, $databaseVersion, $synOrPrimary, undef, undef);
	    GUS::DoTS::DbRefNAFeature->new($gusTableWriters, $dbRefId, $gusGeneFeature->getPrimaryKey());
	}
    }

    my $geneStableId = $gene->stable_id();

    my $organismAbbrev = $self->getOrganism()->getOrganismAbbrev();	
    my $databaseName = "${organismAbbrev}_PreviousGeneIDs_aliases";
    my $databaseVersion = $self->getOrganism()->getGenomeDatabaseVersion();
    foreach my $previousStableId (@{$self->getPreviousIdentifiersFromPatchBuild()->{$geneStableId}}) {
	my ($dbRefId, $externalDatabaseReleaseId) = $self->getDbRefAndExternalDatabaseReleaseIds($databaseName, $databaseVersion, $previousStableId, undef, undef, "previous id");
	GUS::DoTS::DbRefNAFeature->new($gusTableWriters, $dbRefId, $gusGeneFeature->getPrimaryKey());
    }
}



sub parseTranscript {
    my ($self, $transcript, $gene, $gusGeneFeature, $taxonId, $exonMap, $gusExternalDatabaseRelease) = @_;

    my $gusTableWriters = $self->getGUSTableWriters();

    my $splicedNASequenceOntologyId = $self->ontologyTermFromName("mature_transcript", $gusTableWriters);
    my $gusSplicedNASequence = GUS::DoTS::SplicedNASequence->new($gusTableWriters, $transcript, $taxonId, $gusExternalDatabaseRelease, $splicedNASequenceOntologyId);

    my $translation = $transcript->translation();

    my $isProteinCoding = $transcript ? 1 : 0;
    
    my $transcriptSequenceOntologyId = $self->ontologyTermFromBiotype($transcript->get_Biotype(), $gusTableWriters);
    
    # add Product name
    my $gusTranscript = GUS::DoTS::Transcript->new($gusTableWriters, $transcript, $gusGeneFeature, $gusSplicedNASequence, $gusExternalDatabaseRelease, $transcriptSequenceOntologyId);
    my $gusTranscriptNALocation = GUS::DoTS::NALocation::Transcript->new($gusTableWriters, $gusTranscript, $gusSplicedNASequence);

    my $geneType = $gene->get_Biotype()->name();

    my ($gusTranslatedAAFeature, $gusTranslatedAASequence);

    if($translation) {
	($gusTranslatedAAFeature, $gusTranslatedAASequence) = $self->parseTranslation($translation, $transcript, $gusTranscript, $taxonId, $gusExternalDatabaseRelease);
    }

    my $exonOrderNum = 1;
    foreach my $exon ( @{ $transcript->get_all_Exons() } ) {
	my $exonKey = $exon->start()."-".$exon->end()."-".$exon->strand()."-".$exon->phase()."-".$exon->end_phase();
	
	my $gusExonId = $exonMap->{$exonKey};
	    
	GUS::DoTS::RNAFeatureExon->new($gusTableWriters, $gusExonId, $gusTranscript, $exonOrderNum);

	if($translation) {
	    my $codingRegionStart;
	    my $codingRegionEnd;
	    if($exon->strand() == -1) {
		$codingRegionStart = $exon->coding_region_end($transcript);
		$codingRegionEnd = $exon->coding_region_start($transcript);
	    }
	    else {
		$codingRegionStart = $exon->coding_region_start($transcript);
		$codingRegionEnd = $exon->coding_region_end($transcript);
	    }

	    # doesn't make a lot of sense but we add a row here even if the cds locations are null
	    GUS::DoTS::AAFeatureExon->new($gusTableWriters, $gusTranslatedAAFeature, $gusExonId, $codingRegionStart, $codingRegionEnd);		    
	}
	$exonOrderNum++;
    }


    foreach my $xref(@{$transcript->get_all_object_xrefs()}) {
	my $databaseName = $xref->dbname();
	if($databaseName eq "GO") {
	    $self->parseGOAssociation($gusTranslatedAASequence, $xref)
	}
	else {
        next unless (defined $xref->analysis()); 
        my $databaseVersion = $xref->analysis()->logic_name();

	    my $primaryId = $xref->primary_id();

	    my ($dbRefId, $externalDatabaseReleaseId) = $self->getDbRefAndExternalDatabaseReleaseIds($databaseName, $databaseVersion, $primaryId, undef, undef);
	    GUS::DoTS::DbRefNAFeature->new($gusTableWriters, $dbRefId, $gusTranscript->getPrimaryKey());
	}
    }


    foreach my $seqEdit (@{$transcript->get_all_SeqEdits()}) {
	my $code = $seqEdit->code();
	my $seqEditOntologyId = $self->ontologyTermFromName($code, $gusTableWriters);	
	GUS::ApiDB::SeqEdit->new($gusTableWriters, $seqEdit, 'transcript', $seqEditOntologyId, $transcript, $transcript->stable_id());
    }
}

sub parseGOAssociation {
    my ($self, $gusTranslatedAASequence, $xref) = @_;


    my $proteinDatabaseName = "DoTS";
    my $proteinTableName = "TranslatedAASequence";

    my $loeName = ref($xref) eq 'Bio::EnsEMBL::OntologyXref' ? $xref->analysis()->logic_name() : $xref->db();

    # NOTE:  We now skip interpro go associations
    return if($loeName eq $INTERPRO2GO_LOGIC);

    my $gusTableWriters = $self->getGUSTableWriters();

    my $proteinDatabaseId = $seenDatabases{$proteinDatabaseName} ? $seenDatabases{$proteinDatabaseName} : GUS::Core::DatabaseInfo->new($gusTableWriters, $proteinDatabaseName)->getPrimaryKey();

    my $proteinTableKey = "$proteinTableName|$proteinDatabaseId";
    my $proteinTableId = $seenTables{$proteinTableKey} ? $seenTables{$proteinTableKey} : GUS::Core::TableInfo->new($gusTableWriters, $proteinTableName, $proteinDatabaseId)->getPrimaryKey();
    
    my $goId = $xref->display_id();
    my $gusGOTermId = $self->ontologyTermFromGOTerm($goId, $gusTableWriters);

    my $goAssociationKey = $proteinTableId . "|" . $gusTranslatedAASequence->getPrimaryKey() . "|" . $gusGOTermId;
    my $goAssociationId = $seenGOAssociations{$goAssociationKey};
    unless($goAssociationId) {
	$goAssociationId = GUS::DoTS::GOAssociation->new($gusTableWriters, $proteinTableId, $gusTranslatedAASequence->getPrimaryKey(), $gusGOTermId)->getPrimaryKey();
    }

    my $goEvidenceLoeId = $seenGOEvidences{$loeName};
    unless($goEvidenceLoeId) {
	$goEvidenceLoeId = GUS::DoTS::GOAssociationInstanceLOE->new($gusTableWriters, $loeName)->getPrimaryKey();
    }

    my $gusGoAssociationInstance = GUS::DoTS::GOAssociationInstance->new($gusTableWriters, $goAssociationId, $goEvidenceLoeId);

    foreach my $linkageType (@{$xref->get_all_linkage_types()}) {
	my $evidenceCodeId = $self->ontologyTermFromEvidenceCode($linkageType, $gusTableWriters);
	GUS::DoTS::GOAssocInstEvidCode->new($gusTableWriters, $evidenceCodeId, $gusGoAssociationInstance->getPrimaryKey());
    }
}


sub parseTranslation {
    my ($self, $translation, $transcript, $gusTranscript, $taxonId, $gusExternalDatabaseRelease) = @_;

    
    my $gusTableWriters = $self->getGUSTableWriters();
    
    my $translatedAASequenceOntologyId = $self->ontologyTermFromName("polypeptide", $gusTableWriters);
    my $gusTranslatedAASequence = GUS::DoTS::TranslatedAASequence->new($gusTableWriters, $transcript, $taxonId, $gusExternalDatabaseRelease, $translatedAASequenceOntologyId);

    #GUS::ApiDB::AaSequenceAttribute->new($gusTableWriters, $translation, $gusTranslatedAASequence);

    my $gusTranslatedAAFeature = GUS::DoTS::TranslatedAAFeature->new($gusTableWriters, $transcript, $gusTranslatedAASequence, $gusTranscript, $gusExternalDatabaseRelease);
    
    my %seenDomains;
    foreach my $proteinFeature (@{$translation->get_all_ProteinFeatures()}) {
	$self->parseProteinFeature($proteinFeature, $translation, $gusTranslatedAAFeature, $gusTranslatedAASequence, \%seenDomains);
    }

    foreach my $xref (@{$translation->get_all_object_xrefs()}) {
	my $databaseName = $xref->dbname();
	
    next unless (defined $xref->analysis());
    my $databaseVersion = $xref->analysis()->logic_name();
    
	my $primaryId = $xref->primary_id();
	
	my ($dbRefId, $externalDatabaseReleaseId) = $self->getDbRefAndExternalDatabaseReleaseIds($databaseName, $databaseVersion, $primaryId, undef, undef);
	GUS::DoTS::DbRefAAFeature->new($gusTableWriters, $dbRefId, $gusTranslatedAAFeature->getPrimaryKey());

	if($databaseName eq 'KEGG_Enzyme') {
	    $self->parseKeggEnzyme($primaryId, $gusTranslatedAASequence->getPrimaryKey(), $databaseName, $externalDatabaseReleaseId);
	}
	
    }

    foreach my $seqEdit (@{$translation->get_all_SeqEdits()}) {
	my $code = $seqEdit->code();
	my $seqEditOntologyId = $self->ontologyTermFromName($code, $gusTableWriters);	
	GUS::ApiDB::SeqEdit->new($gusTableWriters, $seqEdit, 'translation', $seqEditOntologyId, $transcript, $translation->stable_id());
    }
    
    return($gusTranslatedAAFeature, $gusTranslatedAASequence);
}


sub parseKeggEnzyme {
    my ($self, $keggEnzyme, $gusAASequenceId, $databaseName, $externalDatabaseReleaseId) = @_;

    my $gusTableWriters = $self->getGUSTableWriters();
    
    my @ecNumbers = split(/\+/, $keggEnzyme);
    shift @ecNumbers; # remove the first bit which is not an ec number

    foreach my $ec (@ecNumbers) {
	my $gusEnzymeClassId = $seenEnzymeClasses{$ec};
	unless($gusEnzymeClassId) {
	    $gusEnzymeClassId = GUS::SRes::EnzymeClass->new($gusTableWriters, $ec, $externalDatabaseReleaseId)->getPrimaryKey();
	}

	GUS::DoTS::AASequenceEnzymeClass->new($gusTableWriters, $gusAASequenceId, $gusEnzymeClassId, $databaseName)->getPrimaryKey();
    }
}



sub parseProteinFeature {
    my ($self, $proteinFeature, $translation, $gusTranslatedAAFeature, $gusTranslatedAASequence, $seenDomains) = @_;

    my $gusTableWriters = $self->getGUSTableWriters();
    my $logicName = $proteinFeature->analysis()->logic_name();

    if($logicName eq 'seg') {
        my $proteinSourceId = $gusTranslatedAASequence->getGUSRowAsHash()->{source_id};
        $self->addRowToProteinBedFile($proteinSourceId, $proteinFeature, $self->getSegIO());
    }

    # NOTE:  we won't take interpro,signalp or tmhmm from here

    #     my @gusFeatures;
#     if($logicName eq 'signalp') {
# #	my $f = GUS::DoTS::SignalPeptideFeature->new($gusTableWriters, $gusTranslatedAASequence);
# #	push @gusFeatures, $f;
#     }
#     elsif($logicName eq 'tmhmm') {
# #	my $f = GUS::DoTS::TransMembraneAAFeature->new($gusTableWriters, $gusTranslatedAASequence);
# #	push @gusFeatures, $f;
#     }
#     elsif($logicName =~ /^ms_/) {
# 	# TODO: mass spec peptides
#     }
#     elsif($INTERPRO_LOGICS{$logicName}) {
# 	# my $id = $proteinFeature->display_id();
# 	# if($seenDomains->{$id}) {
# 	#     return; # seen before
# 	# }
# 	# $seenDomains->{$id} = 1;
# 	# my @f = $self->parseInterpro($proteinFeature, $gusTranslatedAASequence, $gusTranslatedAAFeature);
# 	# push @gusFeatures, @f;
#     }
#     elsif($SKIP_LOGICS{$logicName}) { }
#     else {
# 	die "unrecognized logic $logicName";
#     }

#     foreach my $gusFeature (@gusFeatures) {
# 	GUS::DoTS::AALocation->new($gusTableWriters, $proteinFeature, $gusFeature);
#     }

}

# sub parseInterpro {
#     my ($self, $interproFeature, $gusTranslatedAASequence, $gusTranslatedAAFeature) = @_;

#     my $gusTableWriters = $self->getGUSTableWriters();

#     # first make the interpro rows in dbref and domainfeature
#     my $interproSecondaryId = $interproFeature->ilabel();
#     my $interproPrimaryId = $interproFeature->interpro_ac();

#     my $remark = $interproFeature->idesc(); #this is the interpro description used as the dbref remark for both

#     # next make the rows in dbref and domain feature for the domaindb
#     my $domainPrimaryId = $interproFeature->display_id();
#     my $domainSecondaryId = $interproFeature->hdescription();
#     my $evalue = $interproFeature->p_value(); # documentation says e value is gotten from p_value methohd

#     my $analysis = $interproFeature->analysis();
#     my $name = $analysis->display_label() ? $analysis->display_label() : $analysis->logic_name();

#     my $version = $analysis->db_version();

#     my $interproName = $analysis->program();
#     my $interproVersion = $analysis->program_version();

#     my ($interproDbRefId, $interproExternalDatabaseReleaseId) = $self->getDbRefAndExternalDatabaseReleaseIds($interproName, $interproVersion, $interproPrimaryId, $interproSecondaryId, $remark);
#     my ($domainDbRefId, $domainExternalDatabaseReleaseId) = $self->getDbRefAndExternalDatabaseReleaseIds($name, $version, $domainPrimaryId, $domainSecondaryId, $remark);

#     my $interproDomainFeature = GUS::DoTS::DomainFeature->new($gusTableWriters, $gusTranslatedAASequence, undef, $interproExternalDatabaseReleaseId, $interproPrimaryId, undef);
#     my $domainFeature = GUS::DoTS::DomainFeature->new($gusTableWriters, $gusTranslatedAASequence, $interproDomainFeature, $domainExternalDatabaseReleaseId, $domainPrimaryId, $evalue);

#     GUS::DoTS::DbRefAAFeature->new($gusTableWriters, $interproDbRefId, $interproDomainFeature->getPrimaryKey());
#     GUS::DoTS::DbRefAAFeature->new($gusTableWriters, $domainDbRefId, $domainFeature->getPrimaryKey());

#     return($interproDomainFeature, $domainFeature);
# }


sub getDbRefAndExternalDatabaseReleaseIds {
    my ($self, $databaseName, $databaseVersion, $primaryId, $secondaryId, $remark, $idType) = @_;

    my $gusTableWriters = $self->getGUSTableWriters();
    
    my $externalDatabaseReleaseId = $self->getExternalDatabaseReleaseFromNameVersion($databaseName, $databaseVersion);    

    my $dbRefNaturalKey = "$primaryId|$secondaryId|$externalDatabaseReleaseId";

    my $dbRefId = $seenDBRefs{$dbRefNaturalKey};

    unless($dbRefId) {
	my $gusDbRef = GUS::SRes::DbRef->new($gusTableWriters, $primaryId, $secondaryId, $remark, $externalDatabaseReleaseId);
	$dbRefId = $gusDbRef->getPrimaryKey();
    }

    return($dbRefId, $externalDatabaseReleaseId);
}


sub getExternalDatabaseReleaseFromNameVersion {
    my ($self, $name, $version, $idType) = @_;

    my $gusTableWriters = $self->getGUSTableWriters();
    
    my $externalDatabaseId = $seenExternalDatabases{$name};
    unless($externalDatabaseId) {
	my $externalDatabase = GUS::SRes::ExternalDatabase->new($gusTableWriters, $name);
	$externalDatabaseId = $externalDatabase->getPrimaryKey();
    }

    my $extDbRlsSpec = "$externalDatabaseId|" . $version;
    my $externalDatabaseReleaseId = $seenExternalDatabaseReleases{$extDbRlsSpec};

    unless($externalDatabaseReleaseId) {
	$externalDatabaseReleaseId = GUS::SRes::ExternalDatabaseRelease->new($gusTableWriters, $version, $externalDatabaseId, $idType)->getPrimaryKey();
    }

    return $externalDatabaseReleaseId;
}


sub parse {
    my ($self) = @_;

    my $gusTableWriters = $self->getGUSTableWriters();
    
    my $projectName = $self->getProjectName();
    my $projectRelease = $self->getProjectRelease();
    my $gusProjectInfo = GUS::Core::ProjectInfo->new($gusTableWriters, $projectName, $projectRelease);

    my $goSpec = $self->getGOSpec();
    my $goEvidSpec = $self->getGOEvidSpec();
    my $soSpec = $self->getSOSpec();

    my $goExtDbRlsId = $self->getExternalDatabaseReleaseFromSpec($goSpec);
    my $goEvidExtDbRlsId = $self->getExternalDatabaseReleaseFromSpec($goEvidSpec);
    my $soExtDbRlsId = $self->getExternalDatabaseReleaseFromSpec($soSpec);

    $self->setGOExtDbRlsId($goExtDbRlsId);
    $self->setGOEvidExtDbRlsId($goEvidExtDbRlsId);
    $self->setSOExtDbRlsId($soExtDbRlsId);
    
    my $topLevelSlices = $self->getSlices();
    my $organism = $self->getOrganism();
    
    my $gusTaxon = GUS::SRes::Taxon->new($gusTableWriters, $organism);

    my $gusExternalDatabaseRelease = $self->getExternalDatabaseRelease($gusTableWriters, $organism);
    
    foreach my $slice (@$topLevelSlices) {
	$self->parseSlice($slice, $gusExternalDatabaseRelease, $gusTaxon);
    }

    $self->finishBedFiles();
}

sub getExternalDatabaseReleaseFromSpec {
    my ($self, $spec) = @_;
    my ($name, $version) = split(/\|/, $spec);
    return $self->getExternalDatabaseReleaseFromNameVersion($name, $version);
}



1;
