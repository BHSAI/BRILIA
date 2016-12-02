function [H1, Ax] = plotGeneFamilyUsage(FamData,varargin)
%Open the VDJ Gene database in .mat file, if it exist. 
DefaultFile = 'IMGT_Mouse_VDJ_Database.mat';
if exist(DefaultFile,'file') > 0
    load(DefaultFile);
else
    [DefaultFile, DefaultPath] = uigetfile('*.mat','Select VDJ gene database file');
    load([DefaultPath DefaultFile]);
end

%Establish the unique families and bins/
UnqVfam = cat(1,{'Unresolved'},unique(Vmap(:,4)));
UnqJfam = cat(1,{'Unresolved'},unique(Jmap(:,4)));
UnqDfam = cat(1,{'Unresolved'},unique(Dmap(:,4)));
UnqDfam(:) = UnqDfam([end:-2:3 1 2:2:end]); %resort.

%Establish bins and start counting
Vbin = zeros(length(UnqVfam),length(FamData));
Dbin = zeros(length(UnqDfam),length(FamData));
Jbin = zeros(length(UnqJfam),length(FamData));
Vedges = 0:length(UnqVfam);
Dedges = -(length(UnqDfam)-1)/2:1:(length(UnqDfam)-1)/2+1;
Jedges = 0:length(UnqJfam);
for j = 1:length(FamData)
    Vbin(:,j) = histcounts(FamData{j}(:,1),Vedges);
    Dbin(:,j) = histcounts(FamData{j}(:,2),Dedges); 
    Jbin(:,j) = histcounts(FamData{j}(:,3),Jedges);
end

%Determine the legend names
LegendName = [];
if length(varargin) == 1
    LegendName = varargin{1};
end
if isempty(LegendName)
    LegendName = cell(1,length(FamData));
    for k = 1:length(FamData)
        LegendName{k} = sprintf('Data %d',k);
    end
end

%Determine if normalization is needed
Normalize = 0;
if length(varargin) == 2
    if strcmpi(varargin{2},'norm');
        Normalize = 1;
    end
end

%Make the plots
figure
[Hv, AxV] = plotBarGraph(Vedges,Vbin,Normalize,UnqVfam,LegendName,'V Family Name');
set(AxV,'OuterPosition',[0 0 1 1]);

figure
[Hd, AxD] = plotBarGraph(Dedges,Dbin,Normalize,UnqDfam,LegendName,'D Family Name');
set(AxD,'OuterPosition',[0 0 1 1]);

figure
[Hj, AxJ] = plotBarGraph(Jedges,Jbin,Normalize,UnqJfam,LegendName,'J Family Name');
set(AxJ,'OuterPosition',[0 0 1 1]);

Ax = {AxV AxD AxJ};
H1 = {Hv Hd Hj};