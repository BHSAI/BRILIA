function showCredits(varargin)
fprintf('\n');
if any(startsWith(varargin, 'bhsai', 'IgnoreCase', true))
    %fprintf('  Running BRILIA v%s\n', BRILIA('version'));
    fprintf('  Developed at BHSAI (bhsai.org)\n');
    fprintf('  Website: https://github.com/bhsai/brilia\n');
    fprintf('  Online Ref: https://www.frontiersin.org/articles/10.3389/fimu.2016.00681/full\n');
    fprintf('  Written by: Donald Lee dlee@bhsai.org\n');
    fprintf('\n');
end

if any(startsWith(varargin, 'imgt', 'IgnoreCase', true))
    fprintf('  VDJ gene database from: http://www.imgt.org\n');
    fprintf('  IMGT founder and director: Marie-Paule Lefranc, Montpellier, France\n');
    fprintf('\n');
end
