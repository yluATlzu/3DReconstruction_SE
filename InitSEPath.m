% Adding the path to the code in this directory
%
% Usage
%    initpath();
%

[pathstr, ~, ~] = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(pathstr,'Code')));
