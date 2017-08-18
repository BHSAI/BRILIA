function testHumanH
MFileName = mfilename('fullpath');
SlashLoc = regexp(MFileName, '\\|\/');
SlashType = MFileName(SlashLoc(end));
FilePath = MFileName(1:SlashLoc(end));
disp(FilePath);
BRILIA([FilePath 'HumanH_Comma.csv'], 'Species', 'Human', 'Strain', 'All', 'Chain', 'H');
BRILIA([FilePath 'HumanH_Fasta.fa'], 'Species', 'Human', 'Strain', 'All', 'Chain', 'H');
runAnalysis([FilePath 'HumanH_Comma' SlashType 'HumanH_Comma.BRILIAv3.csv']);