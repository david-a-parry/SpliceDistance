=head1 LICENSE

Copyright (c) 2018 David A. Parry

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

=head1 CONTACT

 David Parry <david.parry@ed.ac.uk>
   
=cut

=head1 NAME

 SpliceDistance

=head1 SYNOPSIS

 mv SpliceDistance.pm ~/.vep/Plugins
 vep -i variations.vcf --plugin SpliceDistance

=head1 DESCRIPTION

 A VEP plugin that adds an annotation field indicating distance from the
 nearest splice site. If a variant overlaps both an exon and an intron no
 annotation is added.

=cut

package SpliceDistance;

use strict;
use warnings;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub feature_types {
    return ['Transcript'];
}

sub get_header_info {
    return {SpliceDistance => "Distance from nearest splice site"};
}

sub run {
    my ($self, $tva) = @_;
    my $tr = $tva->transcript;
    if ($tr->end_Exon() == 1){ #no introns
        return {};
    }
    my $tv = $tva->transcript_variation;
    my $vf = $tva->variation_feature;
    my $intron_donor = $tr->strand > 0 ? 'start' : 'end';
    my $intron_acc = $tr->strand > 0 ? 'end' : 'start';
    my $exon_donor = $tr->strand > 0 ? 'end' : 'start';
    my $exon_acc = $tr->strand > 0 ? 'start' : 'end';
    my ($vf_start, $vf_end) = sort {$a <=> $b} ($vf->start, $vf->end);
    my @introns = @{$tv->_overlapped_introns($vf_start, $vf_end)};
    my @exons = @{$tv->_overlapped_exons($vf_start, $vf_end)};
    if (@introns > 1 or @exons > 1){
        #if overlap more than one intron or exon do nothing -
        #already alters canonical splice site
        return {};
    }
    my $start;
    my $end;
    if (@exons){
        my $rank = $exons[0]->rank($tr);
        my $last = $tr->end_Exon() == $rank;
        my $a_start = abs($vf_start - $exons[0]->{$exon_acc}) + 1;
        my $a_end = abs($vf_end - $exons[0]->{$exon_acc}) + 1;
        my $d_start = abs($exons[0]->{$exon_donor} - $vf_start) + 1;
        my $d_end = abs($exons[0]->{$exon_donor} - $vf_end) + 1;
        if (($d_start < $a_start and not $last) or $rank == 1){
            $start = "donor-$d_start";
        }else{
            $start = "acceptor+$a_start";
        }
        if ($d_end != $d_start){
            if ($d_end < $a_end){
                $end = "donor-$d_end";
            }else{
                $end = "acceptor+$a_end";
            }
        }
    }elsif (@introns){
        my $d_start = abs($introns[0]->{$intron_donor} - $vf_start ) + 1;
        my $d_end = abs($introns[0]->{$intron_donor} - $vf_end ) + 1;
        my $a_start = abs($introns[0]->{$intron_acc} - $vf_start) + 1;
        my $a_end = abs($introns[0]->{$intron_acc} - $vf_end) + 1;
        if ($d_start <= $a_start){
            $start = "donor+$d_start";
        }else{
            $start = "acceptor-$a_start";
        }
        if ($d_end != $d_start){
            if ($d_end <= $a_end){
                $end = "donor+$d_end";
            }else{
                $end = "acceptor-$a_end";
            }
        }
    }else{
        return {}; #not intronic or exonic
    }
    my @return = ($start);
    push @return, $end if defined $end;
    return {SpliceDistance => \@return};
}

    1;
