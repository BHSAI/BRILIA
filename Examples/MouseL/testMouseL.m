function testMouseL
MFileName = mfilename('fullpath');
SlashLoc = regexp(MFileName, filesep);
FilePath = MFileName(1:SlashLoc(end));
disp(FilePath);
BRILIA([FilePath 'MouseL_Comma.csv'], 'Species', 'Mouse', 'Strain', 'All', 'Chain', 'L', 'DevPerc', 5);
runAnalysis([FilePath 'MouseL_Comma' filesep 'MouseL_Comma.BRILIAv3.csv']);