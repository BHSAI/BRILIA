##BRILIA path notice for version 2.0.1

A lot of changes have been made to make BRILIA more tolerable to a variety of input sequences. Below is a list of major changes.

- Source codes in Main_Codes have been cleaned and new description have been added.
- No longer asks users to select a "strain" when using the human VDJ germline database.
- Will pad sequences with "x" to ensure same length CDR3 lengths have same length sequences, in case sequences were not trimmed.
- Will flip sequences to the reverse complement if it find a better V gene match. In case inputs have a mix of forward/reverse sequences.
- Sequence alignment now relies partially on seed-based alignment that looks for the conserved C and W, prior to finding the V D J gene.
- N region score is modified to penalize long sequence. Equation changed from () to (). 
- Alignment score calculations have been adjusted to lower alignment score for flanking non-match regions IF desired.
- Try and Catch statments inserted to prevent a seingle entry from stopping the entire alignment process.
- Unprocessed sequences are removed from the main file and set aside in an "Unprocessed" file.
