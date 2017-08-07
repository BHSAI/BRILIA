# BRILIA  v2.0.7
## (B-cell repertoire inductive lineage and immunosequence annotator)

## REFERENCE:
[Lee, D.W., I. Khavrutskii, A. Wallqvist, S. Bavari, C. Cooper, and S. Chaudhury. BRILIA: Integrated Tool for High-Throughput Annotation and Lineage Tree Assembly of B-Cell Repertoires. Frontiers in Immunology, 2017. 7(681).](http://journal.frontiersin.org/article/10.3389/fimmu.2016.00681/full)

## CONTACT:
  *  Donald for questions about BRILIA software, bugs, etc: dlee@bhsai.org  
  *  Sid for questions about the immunosequencing research we do: schaudhury@bhsai.org

## PURPOSE:

BRILIA is designed to annotate a repertoire of B-cell receptor, heavy-chain sequences across the VDJ junction. It returns the CDR3 regions, VDJ germline gene predictions, and also phylogeny relationships among B cells. For more information on how BRILIA works, please read the methods article cited above.
  
## INPUT FILES (see Example_Files folder): 

Takes fasta, fastq, csv, xlsx, or xlsx file containing nucleotide sequences of VDJ junctions. Users should ensure that:
  *  the 1st row of tabulated data have data labels "SeqName" and "Seq". "TemplateCount" is optional.
  *  ambiguous non-nucleotide letters are minimal, otherwise BRILIA will treat them as wildcard nucleotides.

## OUTPUT FILE (see Example_Files folder): 

Returns a semicolon-delimited CSV file listing the annotation results and phylogeny relationships among sequences. This file is saved in a new folder called BRILIA and the version number is added to the file name (Ex: FileName.BRILIAv1.2.3.csv). Sequences that could not be processed are saved in a separate file ending with Err.csv.   

See output file column definitions [here](https://github.com/BHSAI/BRILIA/blob/master/Support_Files/DataHeaderInfo.csv).

## Running BRILIA without MATLAB
1. Download the MATLAB runtime library 9.0.1 (specifically this version) from the [MathWorks website] (https://www.mathworks.com/products/compiler/mcr.html).
2. Install the MRC library on your computer.
3. Download the BRILIA exectuable file. [BRILIAv2.0.7.zip](https://github.com/BHSAI/BRILIA/files/767682/BRILIAv2.0.7.zip)
4. Unzip the exe file and run the BRILIA.exe. It may ask for permission to read/write files.
5. Follow the top-down flow of the GUI to setup the parameters, process sequence files, and draw lineage trees. The GUI itself has help popup texts if you hover your cursor over the setting name and buttons.
6. To see lineage trees, use "Tree Tool". You can also use Tree Tool directly if you have other BRILIA-processed files.

## Running BRILIA in MATLAB (requires bioinformatics toolbox):

1. Copy all codes into a folder called BRILIA (Delete older BRILIA files as this could cause code conflicts).
2. Open MATLAB and set the working directory to the BRILIA folder.
3. In the command line, invoke the addAllPaths function to add all BRILIA folders into the matlab path.  
   > addAllPaths
4. To run BRILIA GUI, run the BRILIAgui in the command line:  
   > BRILIAgui

   OR

   To run BRILIA by code, use either of the following example commands (type help BRILIA to learn more):

   EX1) Will ask user to select the input sequence file, host species, and host strain.
   > BRILIA  
   
   EX2) Will ask user to select file while using all BRILIA parameters defined in a txt file (see SettingExample.txt)
   > BRILIA( [], 'SettingFile', 'SettingExample.txt' )    

   EX3) Will process sequence file named Seqfile.fasta using the Human VDJ gene database and other settings specfied by ParamName-Value pairs.
   > BRILIA( 'Seqfile.fasta' , 'Species' , 'human' , ParamName , Value, ... )  

   HINT) Try processing the example input files in the Example_Files folder.
   > BRILIA( 'ExMouseSeq_Semicolon.csv', 'SettingFile' , 'ExMouseSeq_SettingFile.txt' );  
   > BRILIA( 'ExMouseSeq_Tabulated.csv', 'Delimiter', '\t', 'Species', 'Mouse', 'Strain', 'C57BL');

NOTE: BRILIA should create a new folder called BRILIA and save the output results in that folder.

The program is distributed under the [GNU General Public License] (http://www.gnu.org/licenses/gpl.html).  

See BRILIA patch info at [here] (https://github.com/BHSAI/BRILIA/blob/master/PatchInfo.md).  

## UPCOMING UPDATES (2017-08-07)
  The next version of BRILIA (v3) is on its way with some MAJOR changes!
  *  More error handling of non-VDJ sequences
  *  Heavy and Light chain annotations
  *  CDR1, 2, and 3 annotations according to IMGT definitions
  *  Segments annotation jobs by CDR3 lengths to reduce memory overload
  *  Binary files for using BRILIA in a command-line interface
  *  Data plotting tools for making publication-quality figures in Matlab
  *  Analysis plots to show repertoire properties such as SHMs, etc.
  *  Cleaner codes and folder names
  *  Cleaner annotation files with alignment information
  
