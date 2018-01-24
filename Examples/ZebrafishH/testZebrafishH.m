function testZebrafishH
MFileName = mfilename('fullpath');
SlashLoc = regexp(MFileName, filesep);
FilePath = MFileName(1:SlashLoc(end));
disp(FilePath);
BRILIA([FilePath 'ZebrafishH.csv'], 'Species', 'Zebrafish', 'Chain', 'H', 'DevPerc', 5);
runAnalysis(fullfile(FilePath, 'ZebrafishH', 'ZebrafishH.BRILIAv3.csv'));
