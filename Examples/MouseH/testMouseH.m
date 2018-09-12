function testMouseH(Option)
if nargin == 0
    Option = 0;
end
MFileName = mfilename('fullpath');
SlashLoc = regexp(MFileName, filesep);
FilePath = MFileName(1:SlashLoc(end));
disp(FilePath);
if Option == 0 || Option == 1
    BRILIA([FilePath 'MouseH_Fasta.fa'], 'Species', 'Mouse', 'Strain', 'C57BL', 'Chain', 'H');
    runAnalysis([FilePath 'MouseH_Fasta' filesep 'MouseH_Fasta.BRILIAv4.csv']);
end
if Option == 0 || Option == 2
    BRILIA([FilePath 'MouseH_Comma.csv'], 'Species', 'Mouse', 'Strain', 'C57BL', 'Chain', 'H');
    runAnalysis([FilePath 'MouseH_Comma' filesep 'MouseH_Comma.BRILIAv4.csv']);
end