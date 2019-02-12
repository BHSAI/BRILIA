# BRILIA v3.5.9

**B-cell Repertoire Inductive Lineage and Immunosequence Annotator**

BRILIA is designed to be a single-platform software capable of processing B-cell receptor sequences, annotating VDJ junctions, assigning B-cell lineages/clonotypes, and characerizing B-cell repertoires. 

**REFERENCE**  
[Lee, D.W., I. Khavrutskii, A. Wallqvist, S. Bavari, C. Cooper, and S. Chaudhury. BRILIA: Integrated Tool for High-Throughput Annotation and Lineage Tree Assembly of B-Cell Repertoires. Frontiers in Immunology, 2017. 7(681).](http://journal.frontiersin.org/article/10.3389/fimmu.2016.00681/full)

**CONTACT**  
  E-mail: brilia@bhsai.org (address to Donald Lee for BRILIA help, or Sid Chuadhury for all other questions)  
  Hours: M-F, 9am-5pm EST, off on federal holidays  
  
## Inputs
 
  * Accepts DNA/RNA sequence files with the following extensions:
    * .fasta or .fa
    * .fastq (see MinQuality input parameter to treat low-quality base reads as "N" bases)
    * .csv (comma-delimited)
    * .ssv (semicolon-delimited)
    * .tsv (tab-delimited)
  * Does not accept unassembled, paired-end reads. These should be assembled using another software.
  * Non-ACGTU letters are treated as wildcard "N" bases. Special letters are not used yet (Ex: R for A/G will be treated as N)
  * For delimited files, make sure that the first row is the data header. Ex: "SeqName, H-Seq, L-Seq, TemplateCount"
      * SeqName: the name of the sequence
      * H-Seq: heavy chain sequence
      * L-Seq: light chain sequence
      * TemplateCount: the number of sequence copies (Optional)
  * See [example input files](https://github.com/BHSAI/BRILIA/tree/master/GitDocs/ExampleIO)

## Outputs

  * Returns 3 delimited csv file:
    * [output_file_name].BRILIAv3.csv : final annotation and phylogeny data of productive V(D)J sequences
    * [output_file_name].BRILIAv3.Raw.csv : initial annotation of V(D)J sequences without lineage-base annotation correction.
    * [output_file_name].BRILIAv3.Err.csv : non-productive VDJ sequences and any sequences that could not be annotated fully.
  * If the output file is not specified, files will be stored in a subfolder with the same name as the input file.
  * See [output file header info](https://github.com/BHSAI/BRILIA/blob/master/Tables/DataHeaderInfo.csv)
  * See [example output files](https://github.com/BHSAI/BRILIA/tree/master/GitDocs/ExampleIO/MouseH)  

## Installing BRILIA 

  1. Dowload & install the [MATLAB Runtime R2017a(9.2) from MathWorks](https://www.mathworks.com/products/compiler/matlab-runtime.html) into a folder, referred to as `RUNTIME_FOLDER`.
  2. Download & unzip the [BRILIA exe files](https://github.com/BHSAI/BRILIA/releases/) into a folder, referred to as `BRILIA_FOLDER`.
  3. Open the OS terminal and navigate to `BRILIA_FOLDER`.
  4. Run BRILIA. If the OS asks for permission to modify files, "cancel" should work. See [here for more details](https://github.com/BHSAI/BRILIA/blob/master/GitDocs/FilePermission.md).  
     ``` Windows> BRILIA.exe ```  
     ``` Linux$ run_BRILIA.sh RUNTIME_FOLDER  ```  
     You should end up with a prompt like this:  
     ```BRILIA>  ```  

## Using BRILIA 

  ```
  Assuming the files are structured as such:  
  /home/user/my files/file1.fastq
                      file2.fastq
                      fileN.fastq
                      
  Output files will be placed as such:
  /home/user/my files/file1/file1.BRILIAvN.csv
                      file2/file2.BRILIAvN.csv
                      fileN/fileN.BRILIAvN.csv
  

  BRILIA> ?    (see help)
  BRILIA> cd "/home/user/my files"    (go to sequece file parent folder. Use quotes "" if there is a space in the path)
  BRILIA> "./file*.fastq" -chain h -species mouse -strain balb    (process heavy chain for mouse balb/c)

   Developed at BHSAI (bhsai.org)
   Website: https://github.com/bhsai/brilia
   Status: ...
   Status: Finished in 3.3 sec.

  BRILIA> mergeSimilarSeq "./*" 0.01  (optional: parent-child sequences with < 1% Hamming distance will be merged under the parent)
  
  BRILIA> runAnalysis "./*"   (run analysis on all output files in subfolder of this current folder)

    Status: Drawing 46 trees.
    Status: ...
  
  BRILIA> exit 
  ```

## Tips

  1. Consider adding `BRILIA_FOLDER` ( and `RUNTIME_FOLDER` for Linux ) to the OS `PATH` to call BRILIA easily.  
  2. You could skip entering the `BRILIA>` prompt by adding inputs after the command in step 4. Example:  
     ``` Windows> BRILIA.exe species mouse chain h ... ```  
     ``` Linux$ run_BRILIA.sh RUNTIME_FOLDER species mouse chain h ... ```  
  3. To test if BRILIA is working, try the following:  
     ``` BRILIA> test all ```  
     A folder called `Examples` will be generated in the `BRILIA_FOLDER` and BRILIA will process these folders.  
  4. For Linux, make sure executable files have read (r) + execute (x) permissions. Add permissions if needed.
      ```
      BRILIA_FOLDER$ chmod +rx BRILIA
      BRILIA_FOLDER$ chmod +rx run_BRILIA.sh
      ```  

## License, Patches, and Ongoing Developements

  * The program is distributed under the [GNU General Public License](http://www.gnu.org/licenses/gpl.html).  
  * [BRILIA patch info](https://github.com/BHSAI/BRILIA/blob/master/PatchInfo.md). 
  * [Feature requests and open issues](https://github.com/BHSAI/BRILIA/issues).
