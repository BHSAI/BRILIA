function testMouseH
MFileName = mfilename('fullpath');
SlashLoc = regexp(MFileName, '\\|\/');
SlashType = MFileName(SlashLoc(end));
FilePath = MFileName(1:SlashLoc(end));
disp(FilePath);
BRILIA([FilePath 'MouseH_Fasta.fa'], 'Species', 'Mouse', 'Strain', 'All', 'Chain', 'H', 'DevPerc', 5);
BRILIA([FilePath 'MouseH_Comma.csv'], 'Species', 'Mouse', 'Strain', 'All', 'Chain', 'H', 'DevPerc', 5);
runAnalysis([FilePath 'MouseH_Comma' SlashType 'MouseH_Comma.BRILIAv3.csv']);
