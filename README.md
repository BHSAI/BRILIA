# BRILIA  v3.0.5
## (B-cell repertoire inductive lineage and immunosequence annotator)

## REFERENCE:
[Lee, D.W., I. Khavrutskii, A. Wallqvist, S. Bavari, C. Cooper, and S. Chaudhury. BRILIA: Integrated Tool for High-Throughput Annotation and Lineage Tree Assembly of B-Cell Repertoires. Frontiers in Immunology, 2017. 7(681).](http://journal.frontiersin.org/article/10.3389/fimmu.2016.00681/full)

## CONTACT:
  *  Donald for questions about BRILIA software, bugs, etc: dlee@bhsai.org  
  *  Sid for questions about the immunosequencing research: schaudhury@bhsai.org

## PURPOSE:

BRILIA annotates VDJ and/or VJ junctions from a repertoire of B-cells, returning information about CDR regions, phylogeny relationships, and somatic hypermutations. More information about BRILIA can be found [here.](http://journal.frontiersin.org/article/10.3389/fimmu.2016.00681/full)
  
## INPUT FILES (See [example input files](https://github.com/BHSAI/BRILIA/tree/master/Examples/MouseH)): 
 
  * Accepted sequence file formats: 
    * fasta (.fa or .fasta)
    * fastq (.fastq)
    * comma-delimited (.csv)
    * semicolon-delimited (.ssv)
    * tab-delimited (.tsv)
  * Does not accept paired-end reads, so please assemble them using a 3rd party software.
  * Non-nucleotide letters are treated as wildcard "X" or "N" nucleotides.
  * For delimited files, make sure that:
    * 1st row is the data header: "SeqName,H-Seq,L-Seq,TemplateCount"
    * data is stored in subsequent rows in the order of the data header. 
      * SeqName: the name of the sequence. Required info.
      * H-Seq: heavy chain sequence. Optional IF doing Light chain only.
      * L-Seq: light chain sequence. Optional IF doing Heavy chain only.
      * TemplateCount: the number of sequence copies. Optional if DNA template count or copy number is unknown.
  * Delimiter symbols cannot be used in places where they are not delimiters.
  * NOTE: BRILIA will auto-detect the delimiter and heavy/light chain if formatted properly.

## OUTPUT FILES (See [example output files](https://github.com/BHSAI/BRILIA/tree/master/Examples/MouseH/MouseH_Fasta)): 

  * Returns 3 delimited csv file:
    * [output_file_name].BRILIAvX.csv : stores final annotation and phylogeny data of productive V(D)J sequences
    * [output_file_name].BRILIAvX.Raw.csv : stores initial annotation of V(D)J sequences without lineage-base annotation correction. This is used as a save point for reclustering sequences. 
    * [output_file_name].BRILIAvX.Err.csv : stores non-productive VDJ sequences and any sequences that could not be annotated fully.
  * If the output file is not specified, results will be stored in a subfolder with the same name as the input file. 
  * See output file header details [here](https://github.com/BHSAI/BRILIA/blob/master/Tables/DataHeaderInfo.csv).

## Running BRILIA without MATLAB
### General preparations and downloads
1. Download the MATLAB Runtime Version R2017a (9.2) from the [MathWorks website](https://www.mathworks.com/products/compiler/mcr.html).
2. Install the MRC library on your computer.
3. Download the folder containing your operating system's [BRILIA exectuable files](https://github.com/BHSAI/BRILIA/tree/master/Bin).
4. Unzip the executable files into a folder called BRILIA.

### For WINDOWS
5. Open windows' command prompt.
6. Go to the BRILIA folder.
```
  \BRILIA>  cd [some_path]/BRILIA
```
7. To process data, type:
```
  \BRILIA>  BRILIA.exe [input_file] Species [species_name] Chain H
```
8. To plot lineage trees, type the following:
```
  \BRILIA>  plotTree.exe [output_file]
```

### For LINUX
5. Open up the linux terminal.
6. Go to the BRILIA folder.
```
  BRILIA]$  cd [some_path]/BRILIA
```
7. To process data, type:
```
  BRILIA]$  ./run_BRILIA.sh ../mcr/v92/ [input_file] Species [species_name] Chain H
```
8. To plot lineage tree, type:
```
  BRILIA]$  ./run_plotTree.sh ../mcr/v92/ [output_file]
```

## Running BRILIA in MATLAB 

Note: Requires bioinformatics toolbox. Use `>> ver` command in Matlab to check.
1. Copy all codes into a folder called BRILIA, deleting any older BRILIA codes to prevent conflicts.
2. Open MATLAB and set the working directory to the */BRILIA* folder.
3. In the command line, add all BRILIA folders into the Matlab path using addAllPaths.  
4. To get details of the inputs and Param-Value pairs:
```
   >> help BRILIA 
```
4. To run BRILIA, use one of the following commands:

   EX1) Will ask user to select the input file, host species, and host strain.
```
   >> BRILIA  
```
   
   EX2) Will ask user to select the input file while using all BRILIA parameters defined in seetting txt file [(see SettingExample.txt.)]()
```
   >> BRILIA( 'InputFile.fa', 'SettingFile', 'SettingExample.txt' )    
```

   EX3) Will process a specified input file using the Human VDJ gene database and other settings specfied by Param-Value pairs.
```
   >> BRILIA( 'InputFile.fa', 'Species', 'human', Param, Value, ... )  
```

   HINT) Look at the mouse heavy chain BRILIA test script found [here](https://github.com/BHSAI/BRILIA/blob/master/Examples/MouseH/testMouseH.m).
```
   >> testMouseH
```


The program is distributed under the [GNU General Public License](http://www.gnu.org/licenses/gpl.html).  

See BRILIA patch info at [here](https://github.com/BHSAI/BRILIA/blob/master/PatchInfo.md).  

## UPCOMING UPDATES (Proposed on 2017-01-24)
  *  Will update the data plotting functions. The Data_Plotting folder is being all reworked, as there's unused / temporary scripts there. New codes will be placed in Plot_Codes folder.
  *  Will add CDR1 and CDR2 into the outputs
  *  Will enforce quality control to ensure nonsense CDR3 or VDJ annotations are removed and placed into [FileName]Err.csv
  *  Will replace try/catch with validation codes, since former method is generally slower.
