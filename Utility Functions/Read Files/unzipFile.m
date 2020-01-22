function [status,result] = unzipFile(filename,outLocation)

% No output location specified- unzip to same directory
if nargin == 1 
   outLocation = fileparts(filename); 
end

% [status,result] = system(['"C:\Program Files\7-Zip\7z.exe" -y x ' '"' filename '"' ' -o' '"' outLocation '"']);


% if contains(filename,'.exe')
%     % if this is a .exe (as is the case for old NGA APC files), just run it
%     system(filename)
% else
    % Use 7zip to open!
    [status,result] = system(['"' cd '\Utility Functions\Third Party\7zip\7za.exe" -y x ' '"' filename '"' ' -o' '"' outLocation '"']);
% end

end