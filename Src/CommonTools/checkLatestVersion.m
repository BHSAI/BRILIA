%checkLatestVersion will check GitHub for the latest version, and compare
%it with the currently used version. If the GitHub version is newer, this
%will return 0 and print a message to download the latest BRILIA version.
%
%  IsLatest = checkLatestVersion
%
%  OUTPUT
%    IsLatest: 0 if there is a newer version in GitHub. 
%              1 if you have the latest version. 
%    Msg: Msg to display about the version mismatched. This is placed as an 
%         output to allow control over when to display the message.
%    
function [IsLatest, Msg] = checkLatestVersion
IsLatest = true;
Msg = {};
Url = 'https://github.com/BHSAI/BRILIA/blob/master/README.md';
try
    UrlTxt = webread(Url);
catch ErrMsg
    disp(ErrMsg);
    return
end

VerLoc = regexpi(UrlTxt, 'BRILIA\s*v(?<Version>\d\.\d\.\d)', 'tokens', 'once');
if ~isempty(VerLoc)
    WebVersion = VerLoc{1};
    WebVersion = cell2mat(convStr2NumMEX(strsplit(WebVersion, '.')));
else
    WebVersion = [];
end

MyVersion = BRILIA('version');
MyVersion = cell2mat(convStr2NumMEX(strsplit(MyVersion, '.')));

if isequal(MyVersion, WebVersion)
    IsLatest = true;
    Msg{1} = sprintf('You have the latest version (%d.%d.%d).', MyVersion);
elseif all(MyVersion >= WebVersion)
    IsLatest = true;
    Msg{1} = sprintf('Your version (%d.%d.%d) is newer than the latest web version (%d.%d.%d).', MyVersion, WebVersion);
else
    IsLatest = false;
    Msg{1} = sprintf('Your version (%d.%d.%d) is older than the latest web version (%d.%d.%d).', MyVersion, WebVersion);
    Msg{2} = sprintf('Download the newest version at https://github.com/BHSAI/BRILIA.');
end

fprintf('%s\n', Msg{:});
