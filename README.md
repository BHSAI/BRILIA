# BRILIA  v2.0.5
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

## MATLAB BASIC USAGE (requires bioinformatics toolbox):

1. Copy all codes into a folder called BRILIA (Delete older BRILIA files as this could cause code conflicts).
2. Open MATLAB and set the working directory to the BRILIA folder.
3. In the command line, type "help BRILIA" to look at how it handles inputs and outputs.
4. To run BRILIA, use either of the following example commands:

   EX1) Will ask user to select the input sequence file, host species, and host strain.
   > 'BRILIA'
   
   EX2) Will ask user to select file while using all BRILIA parameters defined in a txt file (see SettingExample.txt)
   > 'BRILIA( [], 'SettingFile', 'SettingExample.txt' )'

   EX3) Will process sequence file named Seqfile.fasta using the Human VDJ gene database and other settings specfied by ParamName-Value pairs.
   > 'BRILIA( 'Seqfile.fasta' , 'Species' , 'human' , ParamName , Value, ... )'

   HINT) Try processing the example input files in the Example_Files folder.
   > 'BRILIA( 'Ex4_SimMouseBCR_FullLength.fa' , 'SettingFile' , 'Ex4_SettingFile.txt' );'
   > 'BRILIA( 'Ex5_SimHumanBCR_FullLength.fa' , 'SettingFile' , 'Ex5_SettingFile.txt' );'  

5. BRILIA should create a new folder called BRILIA and save the output results in that folder.

The program is distributed under the [GNU General Public License] (http://www.gnu.org/licenses/gpl.html).  

See BRILIA patch info at [here] (https://github.com/BHSAI/BRILIA/blob/master/PatchInfo.md).  

## UPCOMING UPDATES (Proposed on 2017-01-24)
  *  Will update the data plotting functions. The Data_Plotting folder is being all reworked, as there's unused / temporary scripts there. New codes will be placed in Plot_Codes folder.
  *  Will add CDR1 and CDR2 into the outputs
  *  Will enforce quality control to ensure nonsense CDR3 or VDJ annotations are removed and placed into [FileName]Err.csv
  *  Will replace try/catch with validation codes, since former method is generally slower.
