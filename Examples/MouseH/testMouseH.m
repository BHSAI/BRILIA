function testMouseH
MFileName = mfilename('fullpath');
SlashLoc = regexp(MFileName, filesep);
FilePath = MFileName(1:SlashLoc(end));
disp(FilePath);
BRILIA([FilePath 'MouseH_Fasta.fa'], 'Species', 'Mouse', 'Strain', 'All', 'Chain', 'H', 'DevPerc', 5);
runAnalysis([FilePath 'MouseH_Fasta' filesep 'MouseH_Fasta.BRILIAv3.csv']);
BRILIA([FilePath 'MouseH_Comma.csv'], 'Species', 'Mouse', 'Strain', 'All', 'Chain', 'H', 'DevPerc', 5);
runAnalysis([FilePath 'MouseH_Comma' filesep 'MouseH_Comma.BRILIAv3.csv']);


