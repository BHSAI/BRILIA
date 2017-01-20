# BRILIA  v2.0.1
## (B-cell repertoire inductive lineage and immunosequence annotator)

## REFERENCE:

BRILIA: Integrated tool for high-throughput annotation and lineage tree assembly of B-cell repertoires

Donald W. Lee-1, Ilja Khavrutskii-1, Anders Wallqvist-1, Sina Bavari-2, Christopher L. Cooper-2, and Sidhartha Chaudhury-1*

1-Biotechnology HPC Software Applications Institute (BHSAI), Telemedicine and Advanced Technology Research Center, U.S. Army Medical  Research and Materiel Command, Fort Detrick, MD, USA

2-Molecular and Translational Sciences, U.S. Army Medical Research Institute of Infectious Diseases, Frederick, MD, USA

*Corresponding author: schaudhury@bhsai.org

## PURPOSE:

BRILIA is designed to annotate a repertoire of B-cell receptor sequences across the VDJ junction. It returns the CDR3 regions, VDJ germline gene predictions, and also phylogeny relationships among B cells. For more information on how BRILIA works, please read the methods article cited above.
  
## INPUT FILES (see Example_Files folder): 

Takes fasta, fastq, csv, xlsx, or xlsx file containing, hopefully, the full VDJ junction. Users must ensure some basic sequence rules to ensure BRILIA works properly:
  *  Sequence has ambiguous characters outside of "acgtunx" where n and x are both wildcard nucleotides
  *  Sequence contains V and J at least
  *  Tabulated data has a column header names, labeled as "SeqName", "Seq", and/or "TemplateCount". Template count is optional.
  *  NOT raw paired end sequence data. BRILIA assumes all pair end joining has been completed.
  *  PREFER +sense strand direction read to cut down alignment time. If there's a mixtures, BRILIA will perform extra steps to check directionality of sequence.

## OUTPUT FILE (see Example_Files folder): 

Returns a semicolon-delimited CSV file listing the annotation results and phylogeny relationships among sequences. Will create and save to a new folder called BRILIA, and append BRILIAvXX at the end of the file name, where "XX" is the version number (currently 13).

See output file column definitions [here](https://github.com/BHSAI/BRILIA/blob/Dev14/Support_Files/Headers_BRILIA.xlsx).

## MATLAB SOURCE CODE USAGE (requires bioinformatics toolbox):

1. Copy all codes into a folder called BRILIA.
2. Open MATLAB, find the BRILIA folder, and run in the command line "addAllPaths". This will add all folders and subfolders to the MATLAB path. 
3. Type "BRILIA" in the command line and follow instructions. You can also use "BRILIAbatch" instead to process multiple files.
4. Select the input file to process (try the example input files shown in the Example folder). 
5. BRILIA should create a new folder called BRILIA and save the output results in that folder.

The program is distributed under the [GNU General Public License] (http://www.gnu.org/licenses/gpl.html).
