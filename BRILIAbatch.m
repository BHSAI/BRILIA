%BRILIAbatch will perform BRILIA on multiple files

%Setup the input parameters and database files for BRILIA
[Vmap, Dmap, Jmap] = getCurrentDatabase('change'); %Select the database
[Vmap, Dmap, Jmap] = filterRefGene(Vmap,Dmap,Jmap); %Select the host strain
DevPerc = input('Set the fraction of seq deviation for tree clustering. Default: 0.03');
if isempty(DevPerc)
    DevPerc = 0.03;
end

%Select the files containing the sequences
[FileNames, FilePath] = uigetfile('*.xlsx','Selection files','multiselect','on');
if ischar(FileNames)
    FileNames = {FileNames};
end

%Begin processing files with BRILIA
TimeElasped = zeros(length(FileNames),1);
for j = 1:length(FileNames)
    tic
    FullName = [FilePath FileNames{j}];
    BRILIA(FullName,Vmap,Dmap,Jmap,'DevPerc',DevPerc);    
    TimeElasped(j) = toc;
end