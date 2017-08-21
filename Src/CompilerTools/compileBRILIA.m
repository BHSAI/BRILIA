%compileBRILIA will compile BRILIA and generate binary files appropriate
%for the OS (Windows or Linux) being used.

function compileBRILIA()
%Make sure this is running from root, and change current directy if needed
RootDir = findRoot();
SlashType = RootDir(end);
cd(RootDir);
if ~strcmpi(RootDir(1:end-1), cd)
    error('%s: Run this from the BRILIA root dir [ %s]', mfilename, RootDir);
end

%Determine where to save the Bin files
if ispc
    TargetDir = [RootDir 'Bin' SlashType 'Win' SlashType];
else isunix
    TargetDir = [RootDir 'Bin' SlashType 'Linux' SlashType];
end
[Success, Msg] = mkdir(TargetDir);
if Success == 0
    fprintf('%s: Could not move exe files to [ %s ].\n%s\n', mfilename, TargetDir, Msg);
end

%Compile the help text for command-line retrieval of help info
fprintf('%s: %s.\n', mfilename, 'Compiling help texts.');
compileHelpText('Dir', RootDir, 'CheckSub', 'y', 'SaveTo', [RootDir 'HelpText' SlashType], 'Overwrite', 'y');
mcc -m ./Src/CompilerTools/showHelpText.m ...
    -a ./HelpText/ ...
    -a ./Src/CompilerTools/
if ispc
    movefile('showHelpText.exe', TargetDir);
elseif isunix
    movefile('showHelpText', TargetDir);
    movefile('run_showHelpText.sh', TargetDir);
end

fprintf('%s: %s.\n', mfilename, 'Compiling BRILIA.m.');
mcc -m ./BRILIA.m ...
    -a ./Databases/ ...
    -a ./Tables/ ...
    -a ./Src/AnnotationTools/ ...
    -a ./Src/CommonTools/ ...
    -a ./Src/DatabaseTools/ ...
    -a ./Src/FileIOTools/ ...
    -a ./Src/LineageTools/ ...
    -a ./Src/CompilerTools/
if ispc
    movefile('BRILIA.exe', TargetDir);
elseif isunix
    movefile('BRILIA', TargetDir);
    movefile('run_BRILIA.sh', TargetDir);
end

fprintf('%s: %s.\n', mfilename, 'Compiling plotTree.');
mcc -m ./Src/AnalysisTools/Lineage/Tree/plotTree.m ...
    -a ./Src/CommonTools/ ...
    -a ./Src/FileIOTools ...
    -a ./Src/PlotTools/ ...
    -a ./Src/AnalysisTools/ ...
    -a ./Src/LineageTools/ ...
    -a ./Src/CompilerTools/ 
if ispc
    movefile('plotTree.exe', TargetDir);
elseif isunix
    movefile('plotTree', TargetDir);
    movefile('run_plotTree.sh', TargetDir);
end

fprintf('%s: %s.\n', mfilename, 'Compiling runAnalysis.m.');
mcc -m ./Src/AnalysisTools/runAnalysis.m ...
    -a ./Src/CommonTools/ ...
    -a ./Src/FileIOTools ...
    -a ./Src/PlotTools/ ...
    -a ./Src/AnalysisTools/ ...
    -a ./Src/LineageTools/ ...
    -a ./Src/CompilerTools/ 
if ispc
    movefile('runAnalysis.exe', TargetDir);
elseif isunix
    movefile('runAnalysis', TargetDir);
    movefile('run_runAnalysis.sh', TargetDir);
end

%Move the the matlab files
MatlabFiles = {'requiredMCRProducts.txt', 'mccExcludedFiles.log', 'readme.txt'};
for j = 1:length(MatlabFiles)
    if exist(MatlabFiles{j}, 'file')
        try
            delete(MatlabFiles{j});
        catch
        end
    end
end
