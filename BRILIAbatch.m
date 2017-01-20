%BRILIAbatch will perform BRILIA on multiple files

%Setup the options here
Species = 'mouse';  %'human' 'mouse'.    If set to empty, '', will ask user.
Strain = 'C57BL';   %'all' 'C57BL' etc   If set to empty, '', will ask user.
Ddirection = 'all'; %'all' 'fwd' 'inv'   Set direction of D gene search
Vfunction = 'all';  %'all' 'f' 'p' 'orf' Set V gene functionality search
DevPerc = 3;        %0 to 100.           Used for cluster by X % similarity.
FileType = '';      %'fasta', 'fastaq', 'excel', 'delimited'  File type.
Delimiter = ';';    %';' ',' '\t'        Set delimiter for open/save files.
CheckSeqDir = 'n';  %'y' 'n'             If there're complement seq, 'y'.

%Select the files containing the sequences
[FileNames, FilePath] = uigetfile('*.xlsx','Selection files','multiselect','on');
if ischar(FileNames)
    FileNames = {FileNames};
end

%Begin processing files with BRILIA
TimeElasped = zeros(length(FileNames),1);
for j = 1:length(FileNames)
    FullFileName = [FilePath FileNames{j}];
    [~,~,RunTime]= BRILIA(FullFileName,'Species',Species,'Strain',Strain,'Ddirection',Ddirection,'Vfunction',Vfunction,'DevPerc',DevPerc,'FileType',FileType,'Delimiter',Delimiter,'CheckSeqDir',CheckSeqDir);
    TimeElasped(j) = RunTime;
end
