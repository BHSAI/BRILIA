function testMouseHL
MFileName = mfilename('fullpath');
SlashLoc = regexp(MFileName, '\\|\/');
SlashType = MFileName(SlashLoc(end));
FilePath = MFileName(1:SlashLoc(end));
disp(FilePath);
BRILIA([FilePath 'MouseHL_Comma.csv'], 'Species', 'Mouse', 'Strain', 'All', 'Chain', 'HL', 'DevPerc', 5);
runAnalysis([FilePath 'MouseHL_Comma' SlashType 'MouseHL_Comma.BRILIAv3.csv']);