--------------------------------------------------------------------------
## Patch notice for version 3.1.5
 - Fixed a bug in plotGeneUsage.m 

--------------------------------------------------------------------------
## Patch notice for version 3.1.4
 - Enabled inputs of file path or name with a space if encapsulated in quotes. EX: "C:\Temp Dir\"
 - BRILIA now removes sequences with 40% hamming distance error rate in V and J.

--------------------------------------------------------------------------
## Patch notice for version 3.1.3
 - Added Zebrafish to the database
 - Reworked codes for important database VDJ genes from IMGT

---------------------------------------------------------------------------
## Patch notice for version 3.1.2
 - Added error messages and error-skipping features around the padtrimSeqGroup.m code.

---------------------------------------------------------------------------
## Patch notice for version 3.1.1
 - Fixed a bug caused where sequences with start index going beyond the 3' end caused trimGeneEdge to error. Encountering this sequence will generate a WARNING message, and continue the annotation.
 - Enabled BRILIA to be run within it's own environment to prevent having to reinitialize multiple processors per run. BRILIA can now be started via "BRILIA.exe" (Windows) or "run_BRILIA.sh [mcr library path]" (Linux). 

---------------------------------------------------------------------------
## Patch notice for version 3.1.0
 - Version 3.1 introduces MEX routines for improved speeds.
   - alignSeqMEX replaces alignSeq.
   - alignSeqMEX does NOT compute textual alignment results unless asked for by the 4th output argument. This contrasts with alignSeq, which always generate text even when not needed.
   - calcPairDistMEX replaces calcSHMHAMdist, calcHAMdist, and calcPairDist.
   - calcAlignScoreMEX replaces calcAlignScore.
   - trimMatchResultsMEX replaces trimMatchResults.
 - MEX routines are stored under Root/Src/MEX folder
 - Improved code for accessing VDJdata via VDJheader index locations. Uses a map structure instead of parsing header every time.
 
---------------------------------------------------------------------------
## Patch notice for version 3.0.14
 - Allows for annotation-only usage of BRILIA via the 'annotonly' 'y' command. Default is 'n' which would do lineage-based annotation. If 'y', will only do VDJ annotation without correction (good for non-repertoire data).

---------------------------------------------------------------------------
## Patch notice for version 3.0.13
 - Fixed issues with nt2aa function now correctly processing ambiguous nt letters. It now uses wrapper function convNT2AA to handle ambiguous characters.

---------------------------------------------------------------------------
## Patch notice for version 3.0.12
 - Fixed problem with plotTree creating an infinite recursion if multiple copies of a germline sequence exist.
 - Improved speed of readFasta.m for reading fasta files.
 - Moved BRILIA.m into the Src folder so that all m files required for BRILIA are together. The test script m files in the Examples folder are left there as these are specific to the example file.
 - Cleaned up the status text outputs in BRILIA.
 - Parallel processing no longer automatically shuts down after being idle for 60 min.
 - If BRILIA is stopped abruptly for a non-logical error (example: job was killed), then jobs can be resumed using the 'ResumeFrom' setting that points to the TempDir, or 'Resume' 'y' command option that will attempt to find the temp dir in the default location.
 - Removed usage of evals in BRILIA source codes.

---------------------------------------------------------------------------
## Patch notice for version 3.0.8
 - Fixed tree cluster cutoff to work properly as a percentage hamming distance. Before, it was doing half of that.
 - Fixed parallel processing of pairwise distance to be faster by splitting ONLY large jobs amongst workers, whereas small jobs are done in a single core.
 - Code cleanup in the Src/LineageTools folder

---------------------------------------------------------------------------
## Patch notice for version 3.0.7
 - Fixed issues when lineage tree plots were cut off, caused when matlab opengl renderer saves high-res images. Using painters renderer.
 - BRILIA now summons plotTree, runAnalysis, etc... within itself. Ex: BRILIA('plotTree', [BRILIA_file])
 - Binary files only contains one main binary file, BRILIA(.exe) . See above note.
 - Removed dependencies of getSeqCount on fastainfo and fastqinfo.
 - Fixed lineage tree drawing to prevent dots and legend text to go off the the axes edges.

---------------------------------------------------------------------------
## Patch notice for version 3.0.5

Major changes
- Large sequence files will be processed in batches to prevent memory overload. (default is 1000 per batch, changeable using the `BatchSize N` command input).
- Linux and Windows binary files are provided for use in command line environments.
- Lineage trees now starts form the distance to a germline V(D)J to the first sequence in the data set. Before, it started from the most ancestral sequence.
- Lineage trees now use hamming distance (%) instead of BRILIA's internal shm distance metric.
- BRILIA can process paired heavy-light chain sequences. The input files must be delimited and have the "H-Seq" and "L-Seq" column headers.
- Data columns are rearranged so that numeric values are grouped at the end and annotatoins are grouped closer to the beginning.
- Data headers are labeled with either H- or L- to specify annotation for heavy or light chain.
- Reworked the way BRILIA can process fasta files downloaded from IMGT reference gene database, which includes "." gaps.
- User can add custom IG fasta files into the Databases folder. The fasta files must be the same style as IMGT's reference sequence files, including the IMGT gap notation.
- Added CDR1 and CDR2 information.
- Fixed the way sequence alignment computed the "AllowedMiss" point mutations when calculating alignment scores. It now properly uses "MissRate" percentage instead of absolute number of nts.

* Due to these changes, version 3 will not work with lower versions. Please uninstall old version and use the latest version.

---------------------------------------------------------------------------
## Patch notice for version 2.1.0

Important changes
- CDR1 and CDR2 sequences are now provided. If the sequence does not cover the CDR1 and 2 regions, BRILIA will return the germline CDR1 and 2 sequences.
- Fixed the CDR3 sequences to match with IMGT's defintion of CDR3 region, which excludes the 104C and 118W.
- Features for importing IMGT reference sequence fasta files into BRILIA have been improved.
- If the input sequences extend beyond the V gene, the extraneous sequences will be moved to the TrimmedLeft column.
- If the input sequences extend beyond the J gene, the extraneous sequences will be moved to the TrimmedRight column.

General changes
- More bug fixes in the tree plotting scripts.

---------------------------------------------------------------------------
## Patch notice for version 2.0.7

General changes
- Fixed an issue where template counts in CSV files were not being converted to double values correctly.
- Fixed logic issues in findVDJmatch.m and findGeneMatch.m in which the best gene match was incorrect if there were no match. 
- Cleaned up the example files to make it easier to see what input format BRILIA will take. 

---------------------------------------------------------------------------
## Patch notice for version 2.0.6

Generage changes
- made a simple GUI for running BRILIA and drawing lineage tree per cluster (BRILIA/GUI/BRILIAgui.m). Still in early stages of GUI, so more features will be added.
- made changes to BRILIA.m to allow for GUI operation and to allow user to set the processor numbers while computing.
- renamed many coding variables into structure, for future proofing code developments.
- removed many matlab scripts that either were temporarily made for the BRILIA paper, or will no longer work with Version 2 of BRILIA. These codes are being reworked.
- moved some folders around for organizational purposes. 

---------------------------------------------------------------------------
## Patch notice for version 2.0.5

General changes
- renamed Data_Plotting folder to Data_Plotting(Reworking) folder to let users know that I am reworking all codes here.
- new folder created called Plot_Codes that contains the reworked codes.
- fixed the tree plotting codes, which are in the BRILIA/Plot_Codes/DrawLineageTree folder. Look at the plotTreeData.m file.
- removed some unused codes in the SHMtree folder.

---------------------------------------------------------------------------
## Patch notice for version 2.0.4

General changes
- reworked the simulate VDJ codes in the BRILIA/SimulateVDJ/ folder.
- provided new examples files for BRILIA for full-length VDJ sequence processing.
- moved the older example files into separate folder in the BRILIA/Example_Files folder.

---------------------------------------------------------------------------
## Patch notice for version 2.0.1

A lot of changes have been made to make BRILIA more tolerable to a variety of input sequences. Below is a list of major changes, divided by categories.

## #General changes
- The Help comments are updated on all codes in Main_Codes.
- Error handling is improved to prevent a single entry from stopping whole annotation process. Added try and catch statements.
- BRILIA no longer asks users to choose strain when using human databases.
- Unprocessed sequences are removed from the main file and set aside in an "Unprocessed" file.
- Removed dependency on using Excel file format. Now relies on semicolon-delimited file formats mostly.

## #Updates to algorithm
- N region score is modified to penalize long sequence. Equation changed from (P_TDT * L)^2 to (2*P_TDT - 1)*L^2. 
- Sequence alignment now relies partially on seed-based alignment that looks for the conserved C and W, prior to finding the V D J gene.
- Alignment score calculations have been adjusted to allow for lower alignment score for flanking non-match regions IF desired.

## #Specific changes
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
