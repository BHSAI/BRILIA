function IsSame = checkVersion
IsSame = false;
Url = 'https://github.com/BHSAI/BRILIA/blob/master/README.md';
try
    UrlTxt =  webread(Url);
catch ErrMsg
    disp(ErrMsg);
    return
end
VerLoc = regexpi(UrlTxt, 'BRILIA\s*v(?<Version>\d\.\d\.\d)', 'tokens');
if ~isempty(VerLoc)
    Version = VerLoc{1}{1};
    CurVersion = BRILIA('version');
    if strcmpi(Version, CurVersion)
        IsSame = true;
    end
else
    Version = 'Unknown';
end

if ~IsSame
    fprintf('Latest version is %s, but you are currently running %s.\n', Version, CurVersion);
    fprintf('Download latest version at https://github.com/BHSAI/BRILIA\n');
end