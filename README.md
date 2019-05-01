# SpliceDistance
VEP plugin to annotate distance to the nearest splice junction

## SYNOPSIS
     
     mv SpliceDistance.pm ~/.vep/Plugins
     vep -i variations.vcf --plugin SpliceDistance

## DESCRIPTION

A VEP plugin that adds an annotation field indicating distance from the
nearest splice site. For insertions the coordinate before and after the 
insertion are reported. For MNPs or deletions spanning more than one
nucleotide the first and last altered/deleted nucleotides are reported.
No annotation is added if a large variant spans multiple exons or introns
or if a transcript has no introns. 

Note that the acceptor and donor nucleotides (first and last nucleotides of an
exon) are labelled acceptor+0 and donor+0 respectively.
