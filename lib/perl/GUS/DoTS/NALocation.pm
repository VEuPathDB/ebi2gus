package GUS::DoTS::NALocation;
use base qw(GUSRow);

use strict;



sub init {
    my ($self, $feature, $gusParentFeature) = @_;

    return {na_feature_id => $gusParentFeature->getPrimaryKey(),
	    start_min => $feature->seq_region_start(),
	    start_max => $feature->seq_region_start(),
	    end_min => $feature->seq_region_end(),
	    end_max => $feature->seq_region_end(),
	    is_reversed => $feature->seq_region_strand() == -1 ? 1 : 0,
    };
}

1;

package GUS::DoTS::NALocation::Transcript;
use base qw(GUSRow);

sub init {
    my ($self, $gusParentFeature, $gusSplicedNASequence) = @_;

    my $gusSplicedNASequenceAsHash = $gusSplicedNASequence->getGUSRowAsHash();
    
    return {na_feature_id => $gusParentFeature->getPrimaryKey(),
	    start_min => 1,
	    start_max => 1,
	    end_min => $gusSplicedNASequenceAsHash->{length},
	    end_max => $gusSplicedNASequenceAsHash->{length},
	    is_reversed => 0
    };
}


1;
