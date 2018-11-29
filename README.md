# BRILIA  v3.5.1
## (B-cell repertoire inductive lineage and immunosequence annotator)

## REFERENCE:
[Lee, D.W., I. Khavrutskii, A. Wallqvist, S. Bavari, C. Cooper, and S. Chaudhury. BRILIA: Integrated Tool for High-Throughput Annotation and Lineage Tree Assembly of B-Cell Repertoires. Frontiers in Immunology, 2017. 7(681).](http://journal.frontiersin.org/article/10.3389/fimmu.2016.00681/full)

## CONTACT:
  *  brilia@bhsai.org (M-F, 9am-5pm EST, off on federal holidays)

## PURPOSE:

BRILIA is an all-in-one software to process sequencing data of B-cell receptors. It annotates V(D)J junctions, assigns B-cell lineages/clonotypes, and characerizes B-cell repertoires. [reference article](http://journal.frontiersin.org/article/10.3389/fimmu.2016.00681/full).
  
## INPUT FILES (See [example input files](https://github.com/BHSAI/BRILIA/tree/master/Examples/)): 
 
  * Accepted sequence file formats: 
    * .fasta or .fa
    * .fastq (See MinQuality input parameter to treat low-quality base reads as "N" bases)
    * .csv (comma-delimited)
    * .ssv (semicolon-delimited)
    * .tsv (tab-delimited)
  * Does not accept paired-end reads, so please assemble them using a 3rd party software.
  * Non-nucleotide letters are treated as wildcard "N" bases. Special letters are not used yet (Ex: R = A/G is unused)
  * For delimited files, make sure that:
    * delimiter symbols are not used in places where they are not delimiters.
    * first row is the data header. Ex: "SeqName, H-Seq, L-Seq, TemplateCount"
      * SeqName: the name of the sequence
      * H-Seq: heavy chain sequence
      * L-Seq: light chain sequence
      * TemplateCount: the number of sequence copies (Optional)

## OUTPUT FILES (See [example output files](https://github.com/BHSAI/BRILIA/tree/master/Examples/MouseH/MouseH_Fasta)): 

  * Returns 3 delimited csv file:
    * [output_file_name].BRILIAv3.csv : stores final annotation and phylogeny data of productive V(D)J sequences
    * [output_file_name].BRILIAv3.Raw.csv : stores initial annotation of V(D)J sequences without lineage-base annotation correction.
    * [output_file_name].BRILIAv3.Err.csv : stores non-productive VDJ sequences and any sequences that could not be annotated fully.
  * If the output file is not specified, results will be stored in a subfolder with the same name as the input file. 
  * See [output file header info](https://github.com/BHSAI/BRILIA/blob/master/Tables/DataHeaderInfo.csv).

## Running BRILIA without MATLAB license
### General preparations and downloads
1. Download the MATLAB Runtime Version R2017a (9.2) from the [MathWorks website](https://www.mathworks.com/products/compiler/mcr.html).
2. Install the MRC library on your computer.
3. Download the appropriate [BRILIA binary files](https://github.com/BHSAI/BRILIA/releases/).
4. Unzip the executable files into a folder called BRILIA.


## Running BRILIA 
### Brilia can be run either through the OS terminal or via BRILIA's own terminal. Whatever command that can be run via the BRILIA's terminal can also be run in the OS terminal by stating the BRILIA executable first before the inputs. All instructions are for within the BRILIA terminal.
1. Open the BRILIA terminal 
Windows
```BRILIA.exe```
Linux
```run_BRILIA.sh [path_to_matlabMCR]```
2. Run the commands within the BRILIA terminal 
```BRILIA> [commands]
3. To run WITHOUT opening the BRILIA terimnal, simply retype the commands above WITH the inputs you want BRILIA to process right away.
```BRILIA.exe [commands]...```
```run_BRILIA.sh [path_to_matlabMCR] [commands]...```

## Standard Usage Example
### Assuming fastq files are in the folder "/home/user/My Files", perform the following:
1. Annotate the VDJ junctions:
```BRILIA> "/home/user/My Files/*.fa*" species mouse chain h cutoff 1```
2. Do basic data analysis
### Process files as such:
```BRILIA> runAnalysis 





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
8. To plot lineage trees and other plots, type:
```
  \BRILIA>  BRILIA.exe runAnalysis [BRILIA_output_file]
```
  In the same folder where the output file is, There will be a "Trees" and "Analysis" folders.

### For LINUX
5. Open up the linux terminal.
6. Go to the BRILIA folder.
```
  BRILIA]$  cd [some_path]/BRILIA
```

   Make sure BRILIA and run_BRILIA.sh have read (r) + execute (x) permissions. Add permissions if needed.
   ```
      BRILIA]$ chmod +rx BRILIA
      BRILIA]$ chmod +rx run_BRILIA.sh
      BRILIA]$ ls -l
          total 28773
          -r-xr-xr-x. 1 username domain users 29641403 Aug 23 13:10 BRILIA
          -r-xr-xr-x. 1 username domain users      874 Aug 23 13:10 run_BRILIA.sh
   ```
7. To process data, type:
```
  BRILIA]$  ./run_BRILIA.sh [mcr_library_path] [input_file] Species [species_name] Chain H
```
8. To plot lineage trees and other plots, type:
```
  BRILIA]$  ./run_BRILIA.sh [mcr_library_path] runAnalysis [BRILIA_output_file]
```
  In the same folder where the output file is, There will be a "Trees" and "Analysis" folders.

## Running BRILIA in MATLAB 

Note: Requires bioinformatics toolbox. Use `>> ver` command in MATLAB to check.
1. Copy all codes into a folder called BRILIA, deleting any older BRILIA codes to prevent conflicts.
2. Open MATLAB and set the working directory to the */BRILIA* folder.
3. In the command window, add all BRILIA sub folders into the MATLAB path using `>> addAllPaths`.  
4. Get details about BRILIA inputs and Param-Value pairs using `>> help BRILIA`.
5. To run BRILIA, use one of the following commands:

   EX1) Will ask user to select the input file, host species, and host strain.
   ```
      >> BRILIA  
   ```
   EX2) Will process "InputFile.fa" while using parameters stored in a setting txt file (see [example setting file](https://github.com/BHSAI/BRILIA/blob/master/SettingFile.txt)).
   ```
      >> BRILIA( 'InputFile.fa', 'SettingFile', 'SettingFile.txt' )    
   ```
   EX3) Will process "InputFile.fa" using the Human gene database and settings specfied by Param-Value pairs.
   ```
      >> BRILIA( 'InputFile.fa', 'Species', 'human', Param, Value, ... )  
   ```
   HINT) Look at the mouse heavy chain [BRILIA test script](https://github.com/BHSAI/BRILIA/blob/master/Examples/MouseH/testMouseH.m).
   ```
      >> testMouseH
   ```

## Skipping the lineage-based correction

If you only want to do only the annotation without lineage determination (eg, when annotating independent sequences), add this setting:
  ``` 
      > BRILIA InputFile.fa ... AnnotOnly y
  ```


## Resuming BRILIA after a crash

If your job stopped before completion (ex: due to a server crash), then you could recover some of the work. BRILIA saves most of the time-consuming steps in a Temp directory located within the output folder name, which by default is the same as the sequence file name. To resume, follow these steps:

0. (Optional) Make a copy of the Temp dir, which should be where the sequence file is with the same name of that file. EX: if the file is /user/data/MouseH.fa, then the temp dir is /user/data/MouseH/Temp.
1. Retype the ENTIRE command you used to call this job as shown above, BUT ADD **`Resume y`** to the end. If you used a non-default output folder, then use **`ResumeFrom  [output_dir_path]`** instead.
   
   EX: when using the default output folder that stores the Temp dir
   ```
      > BRILIA InputFile.fa ... Resume y 
   ```

   EX: when using a non-default output folder that stores the Temp dir elsewhere
   ```
      > BRILIA InputFile.fa ... ResumeFrom /output_folder/Temp 
   ```
   

The program is distributed under the [GNU General Public License](http://www.gnu.org/licenses/gpl.html).  

See [BRILIA patch info](https://github.com/BHSAI/BRILIA/blob/master/PatchInfo.md). 

See a list of ongoing [developments and open issues](https://github.com/BHSAI/BRILIA/blob/master/OpenIssues.md).
