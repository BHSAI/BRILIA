---------------------------------------------------------------------------
##Patch notice for version 2.0.6

Generage changes
- made a simple GUI for running BRILIA and drawing lineage tree per cluster (BRILIA/GUI/BRILIAgui.m). Still in early stages of GUI, so more features will be added.
- made changes to BRILIA.m to allow for GUI operation and to allow user to set the processor numbers while computing.
- renamed many coding variables into structure, for future proofing code developments.
- removed many matlab scripts that either were temporarily made for the BRILIA paper, or will no longer work with Version 2 of BRILIA. These codes are being reworked.
- moved some folders around for organizational purposes. 

---------------------------------------------------------------------------
##Patch notice for version 2.0.5

General changes
- renamed Data_Plotting folder to Data_Plotting(Reworking) folder to let users know that I am reworking all codes here.
- new folder created called Plot_Codes that contains the reworked codes.
- fixed the tree plotting codes, which are in the BRILIA/Plot_Codes/DrawLineageTree folder. Look at the plotTreeData.m file.
- removed some unused codes in the SHMtree folder.

---------------------------------------------------------------------------
##Patch notice for version 2.0.4

General changes
- reworked the simulate VDJ codes in the BRILIA/SimulateVDJ/ folder.
- provided new examples files for BRILIA for full-length VDJ sequence processing.
- moved the older example files into separate folder in the BRILIA/Example_Files folder.

---------------------------------------------------------------------------
##Patch notice for version 2.0.1

A lot of changes have been made to make BRILIA more tolerable to a variety of input sequences. Below is a list of major changes, divided by categories.

###General changes
- The Help comments are updated on all codes in Main_Codes.
- Error handling is improved to prevent a single entry from stopping whole annotation process. Added try and catch statements.
- BRILIA no longer asks users to choose strain when using human databases.
- Unprocessed sequences are removed from the main file and set aside in an "Unprocessed" file.
- Removed dependency on using Excel file format. Now relies on semicolon-delimited file formats mostly.

###Updates to algorithm
- N region score is modified to penalize long sequence. Equation changed from (P_TDT * L)^2 to (2*P_TDT - 1)*L^2. 
- Sequence alignment now relies partially on seed-based alignment that looks for the conserved C and W, prior to finding the V D J gene.
- Alignment score calculations have been adjusted to allow for lower alignment score for flanking non-match regions IF desired.

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