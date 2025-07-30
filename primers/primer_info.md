### Primer File Information
These are the fasta files called in Cutadapt to search for primers to cut. For 5' primers, the sequence is preceded by an ^. This ^ indicates that the primer is anchored, which means it is found at the 5' end of the sequence and is required to be found in full (i.e. the read will not be kept if the primer is not found).  Most of these primers (except RC primers, see below) have spacers added to the 5' end to increase basepair heterogeneity. These spacers need to be included in the primer sequence to be found, since these primers are anchored. These spacers follow those used and supplied by the Travis Glenn lab.

For primer files that end in "RC.fas", these are the reverse complement of their respective primers. These primers are not anchored (therefore the read is kept regardless of whether the primer is found). Cutadapt will search for these primers on the 3' end of the complementary read (i.e. the RC Forward primer will be found on the 3' end of the R2 (reverse) read). These are used to remove primers when there is read-through (i.e. when the amplicon is shorter than sequencing length).


We currently have 5 primer-pairs available:

COImlIntF/jgCOI2198 - General COI primers, located at the 3' end of the Folmer region. 

MiFish_12SF/MiFish12SR - Fish 12S rRNA primers

18S_V4F/18S_V4R - 18S rRNA primers, amplifying the V4 region

16SbacF/16SbacR: 16S_F515/16S_R806 - 16S bacterial rRNA primers, amplifying the V4 region

28SAnthF/28SAnthR: 28S primers developed for certain anthozoans.

### Custom Primers
To make your own primer pairs for cutadapt to remove, follow the format of existing primers. The primer pair files should be named simply, ending in "F.fas" and "R.fas". In the file, the primers need to have the same name as the primer filename preceeding "F.fas". If you are using primers that have spacers, your file needs to include the separate primer sequences for each spacer (and include the primer without spacers). All primer sequences need to have the same name. For example, COIF.fas has five primer sequences, and all are named COI. Each sequence should be preceded by a "^", because this is a 5' sequence that is required.

If you are using a short gene region that may have read-through (i.e. the sequencer will read the primer/adapter on the 3' end), you should also include a reverse-complement primer file. This primer is not required, so no need for a "^" before the sequence, and you don't need to include the sequences of any spacers, because it does not need to be found at the absolute end of the read. Only the actual primer sequence is needed, but you need to remember to reverse complement this sequence. These files names should be the same as the 5' primer files, but end with "F_RC.fas" and "R_RC.fas". The name of the primers inside these files is not as important as in the 5' primer files.

Custom primers should be saved in the primer folder, along with the rest of the primers.

