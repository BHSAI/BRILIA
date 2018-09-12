function testMouseHL
MFileName = mfilename('fullpath');
SlashLoc = regexp(MFileName, filesep);
FilePath = MFileName(1:SlashLoc(end));
disp(FilePath);
BRILIA([FilePath 'MouseHL_Comma.csv'], 'Species', 'Mouse', 'Strain', 'All', 'Chain', 'HL');
runAnalysis([FilePath 'MouseHL_Comma' filesep 'MouseHL_Comma.BRILIAv4.csv']);