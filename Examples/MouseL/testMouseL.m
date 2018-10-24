function testMouseL
MFileName = mfilename('fullpath');
SlashLoc = regexp(MFileName, filesep);
FilePath = MFileName(1:SlashLoc(end));
disp(FilePath);
BRILIA([FilePath 'MouseL_Comma.csv'], 'Species', 'Mouse', 'Strain', 'All', 'Chain', 'L', 'resume', 'n');
runAnalysis([FilePath 'MouseL_Comma' filesep 'MouseL_Comma.BRILIAv4.csv']);