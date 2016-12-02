[FileName1 FilePath1] = uigetfile('*.mat','Select 1st, reference mat file of mutations');
[FileName2 FilePath2] = uigetfile('*.mat','Select 2nd, sample mat file of mutations');

DotLoc1 = regexp(FileName1,'Mutations');
FileName1pre = FileName1(1:DotLoc1(1)-2);

DotLoc2 = regexp(FileName2,'Mutations');
FileName2pre = FileName2(1:DotLoc2(1)-2);

OutputFilePre = ['cmpr_' FileName1pre '_vs_' FileName2pre];

plotMutation([FilePath1 FileName1],[FilePath2 FileName2],OutputFilePre);
close all
clear