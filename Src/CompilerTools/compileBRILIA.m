%compileBRILIA will compile BRILIA and generate binary files appropriate
%for the OS (Windows or Linux) being used.

function compileBRILIA()
%Make sure this is running from root, and change current directy if needed
RootDir = findRoot();
cd(RootDir);
if ~strcmpi(RootDir(1:end-1), cd)
    error('%s: Run this from the BRILIA root dir [ %s]', mfilename, RootDir);
end

%Determine where to save the Bin files
if ispc
    TargetDir = [RootDir 'Bin' filesep 'Win' filesep];
elseif isunix
    TargetDir = [RootDir 'Bin' filesep 'Linux' filesep];
end
[Success, Msg] = mkdir(TargetDir);
if Success == 0
    fprintf('%s: Could not move exe files to [ %s ].\n%s\n', mfilename, TargetDir, Msg);
end

%Compile the help text for command-line retrieval of help info
fprintf('%s: %s\n', mfilename, 'Compiling help texts.');
compileHelpText('Dir', RootDir, 'CheckSub', 'y', 'SaveTo', [RootDir 'HelpText' filesep], 'Overwrite', 'y');

fprintf('%s: %s\n', mfilename, 'Compiling BRILIA.m.');
mcc -m ./BRILIA.m ...
    -a ./Databases/ ...
    -a ./Tables/ ...
    -a ./Src/ ...
    -a ./HelpText/
if ispc
    movefile('BRILIA.exe', TargetDir);
elseif isunix
    movefile('BRILIA', TargetDir);
    movefile('run_BRILIA.sh', TargetDir);
end

%Delete the other files
MatlabFiles = {'requiredMCRProducts.txt', 'mccExcludedFiles.log', 'readme.txt'};
for j = 1:length(MatlabFiles)
    if exist(MatlabFiles{j}, 'file')
        try
            delete(MatlabFiles{j});
        catch
        end
    end
end
