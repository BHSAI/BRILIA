%showGeneDatabase will list all the species folder in the Databases folder.
%This is used mainly for deployed application, when user wants to see what
%options there are for the species list.
%
%  showGeneDatabase
%
function showGeneDatabase
Species = getGeneDatabase('list');
fprintf('The valid Species options are:\n');
fprintf(' - %s\n', Species{:});
fprintf('To edit the gene database, add/edit the species folder at:\n"%s"\n', fullfile(findExeDir, 'Databases'));