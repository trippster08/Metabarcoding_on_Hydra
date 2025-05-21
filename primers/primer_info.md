### Primer File Information
These are the fasta files called in Cutadapt to search for primers to cut. For most primers, the sequence is preceded by an ^. This ^ indicates that the primer is anchored, which means it is found at the 5' end of the sequence and is required to be found in full (i.e. the read will not be kept if the primer is not found).  All these primers (except RC primers, see below) have spacers added to the 5' end to increase basepair heterogeneity. These spacers need to be included in the primer sequence to be found, since these primers are anchored. These spacers follow those used and supplied by the Travis Glenn lab.

For primer files that end in "RC.fas", these are the reverse complement of their respective primers. These primers are not anchored (therefore the read is kept regardless of whether the primer is found). Cutadapt will search for these primers on the 3' end of the complementary read (i.e. the RC Forward primer will be found on the 3' end of the R2 (reverse) read). These are used to remove primers when there is read-through (i.e. when the amplicon is shorter than sequencing length).


We currently have 5 primer-pairs available:

COImlIntF/jgCOI2198 - General COI primers, located at the 3' end of the Folmer region. 

MiFish_12SF/MiFish12SR - Fish 12S rRNA primers

18S_V4F/18S_V4R - 18S rRNA primers, amplifying the V4 region

16S_F515/16S_R806 - 16S bacterial rRNA primers, amplifying the V4 region

