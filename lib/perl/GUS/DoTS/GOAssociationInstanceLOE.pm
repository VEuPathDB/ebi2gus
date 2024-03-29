package GUS::DoTS::GOAssociationInstanceLOE;
use base qw(GUSRow Exporter);

use strict;

our @EXPORT = qw(%seenGOEvidences);

our %seenGOEvidences;

sub new {
    my $class = shift;

    # this bit calls init
    my $self = $class->SUPER::new(@_);

    my $naturalKey = $self->getNaturalKey();
    $seenGOEvidences{$naturalKey} = $self->getPrimaryKey();

    return $self;
}


sub init {
    my ($self, $name) = @_;

    $self->setNaturalKey($name);
    
    return {name => $name};
}

1;

