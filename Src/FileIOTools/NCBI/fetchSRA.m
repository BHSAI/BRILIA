%fetchSRA will fetch the data from NCBI website using the SRA Tool Kit
%https://www.ncbi.nlm.nih.gov/sra/docs/toolkitsoft/
%
%  fetchSRA OutDir SRA_NUM1 SRA_NUM2 ...
%
%  fetchSRA(OutDir, SRA_NUM1, SRA_NUM2, ...)
%
%  INPUT
%    OutDir: directory to save the download SRA sequence files
%    SRA_NUM: BioProject Accession Number or SRA sample number
%
%  OUTPUT
%    Saves the SRA files into the specified output directory.
%    
function Success = fetchSRA(OutDir, varargin)
if ~isdir(OutDir)
    [Success, Msg] = mkdir(OutDir);
    assert(Success, '%s: Could not make output directory "%s".\n  %s', mfilename, OutDir, Msg);
end

FetchCmd = sprintf('prefetch%s', sprintf(' %s', varargin{:}));
system(FetchCmd)
for j = 1:length(varargin)
    fprintf('fastq conversion for %s\n', varargin{j});
    DumpCmd = sprintf('fastq-dump -O "%s" %s', OutDir, varargin{j});
    system(DumpCmd)
end

