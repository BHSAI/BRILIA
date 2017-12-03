function showCredits(Options)

CellOptions = lower(regexp(Options, '\s*,\s*', 'split'));
fprintf('\n');
if any(ismember(CellOptions, 'bhsai'))
    fprintf('  Running BRILIA v%s\n', BRILIA('version'));
    fprintf('  Developed at BHSAI\n');
    fprintf('  Website: bhsai.org\n');
    fprintf('  Written by: Donald Lee dlee@bhsai.org\n');
    fprintf('\n');
end

if any(ismember(CellOptions, 'imgt'))
    fprintf('  Germline gene databases were downloaded from http://www.imgt.org.\n');
    fprintf('  IMGT founder and director: Marie-Paule Lefranc, Montpellier, France\n');
    fprintf('\n');
end
