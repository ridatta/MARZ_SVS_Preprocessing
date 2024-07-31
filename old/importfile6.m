function Ne1 = importfile6(filename, startRow, endRow)
%IMPORTFILE6 Import numeric data from a text file as a matrix.
%   NE1 = IMPORTFILE6(FILENAME)
%   Reads data from text file FILENAME for the default selection.
%
%   NE1 = IMPORTFILE6(FILENAME, STARTROW, ENDROW)
%   Reads data from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   Ne1 = importfile6('Ne.txt', 7, 431);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2022/07/18 14:38:16

%% Initialize variables.
if nargin<=2
    startRow = 7;
    endRow = 431;
end

%% Format for each line of text:
%   column4: text (%s)
%	column6: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%*11s%18s%*3s%11s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this code. If an error occurs for a different file, try regenerating the code from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Remove white space around all cell columns.
dataArray{1} = strtrim(dataArray{1});
dataArray{2} = strtrim(dataArray{2});

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post processing code is included. To generate code which works for unimportable data, select unimportable cells in a file and regenerate the script.

%% Create output variable
Ne1 = [dataArray{1:end-1}];

