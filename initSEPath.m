% Adding the path to the code in this directory
%
% Usage
%    initSEpath();
%

[pathstr, ~, ~] = fileparts(mfilename('fullpath'));

addpath(genpath(fullfile(pathstr,'Code')));
