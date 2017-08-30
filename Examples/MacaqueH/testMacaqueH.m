function testMacaqueH
MFileName = mfilename('fullpath');
SlashLoc = regexp(MFileName, '\\|\/');
SlashType = MFileName(SlashLoc(end));
FilePath = MFileName(1:SlashLoc(end));
disp(FilePath);
BRILIA([FilePath 'MacaqueH_Comma.csv'], 'Species', 'Macaque', 'Strain', 'All', 'Chain', 'H', 'DevPerc', 5);
runAnalysis([FilePath 'MacaqueH_Comma' SlashType 'MacaqueH_Comma.BRILIAv3.csv']);
