##BRILIA path notice for version 2.0.1

A lot of changes have been made to make BRILIA more tolerable to a variety of input sequences. Below is a list of major changes, divided by categories.

###Major changes to algorithm
- N region score is modified to penalize long sequence. Equation changed from (P_TDT * L)^2 to (2*P_TDT - 1)*L^2. 
- Sequence alignment now relies partially on seed-based alignment that looks for the conserved C and W, prior to finding the V D J gene.
- Alignment score calculations have been adjusted to allow for lower alignment score for flanking non-match regions IF desired.

###General changes
- The Help comments are updated on all codes in Main_Codes.
- Error handling is improved to prevent a single entry from stopping whole annotation process. Added try and catch statements.
- BRILIA no longer asks users to choose strain when using human databases.
- Unprocessed sequences are removed from the main file and set aside in an "Unprocessed" file.

###Specific changes
BRILIA.m
- Can process multiple files, and replaces BRILIAbatch.m
- Can process variable-length sequences by padding sequences within same cluster with 'X'. Important for conforming groups and finding clusters.
- New parameter-value pair (SettingsFile,SettingsFileName) allow users to use settings stored in a text file. See SettingExamples.txt for an example of what a setting file looks like.
- New parameter-value pair (CheckSeqDir,'n') specifies whether or not to check input sequence and its rev-complement direction. 'n' is faster, but use 'y' if input sequence have complement sequences.
- parallel process enabled only for sequence files with more t han 200 sequences, since starting parallel processing takes some time.

getCurrentDatabase
- Better switching between host species databases

filterRefSeq
- Better filtering of database genes based on host strain, direction of D gene, and V gene functionality
