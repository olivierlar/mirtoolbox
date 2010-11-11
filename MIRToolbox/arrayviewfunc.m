function varargout = arrayviewfunc(whichcall, varargin)
%ARRAYVIEWFUNC  Support function for the Variable Editor component

%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.30 $  $Date: 2009/03/09 19:16:59 $

varargout = cell(1, 1);
switch whichcall
    case 'getdata',
        varargout{1} = getData(varargin{1}, varargin{2});
    case 'setdata',
        varargout{1} = setData(varargin{1}, varargin{2}, varargin{3});
    case 'setdatalogic',
        [varargout{1}, varargout{2}] = setDataLogic(varargin{1}, varargin{2}, varargin{3});
    case 'setvarwidth',
        varargout{1} = setVarWidth(varargin{1}, varargin{2});
    case 'setvarwidthlogic',
        [varargout{1}, varargout{2}] = setVarWidthLogic(varargin{1}, varargin{2});
    case 'setvarheight',
        varargout{1} = setVarHeight(varargin{1}, varargin{2});
    case 'setvarheightlogic',
        [varargout{1}, varargout{2}] = setVarHeightLogic(varargin{1}, varargin{2});
    case 'removerowsorcolumns',
        varargout{1} = removeRowsOrColumns(varargin{:});
    case 'insertrowsorcolumns',
        varargout{1} = insertRowsOrColumns(varargin{:});
    case 'renamefield',
        varargout{1} = renameField(varargin{1}, varargin{2}, varargin{3});
    case 'valueHasAppropriateIndexing',
        varargout{1} = valueHasAppropriateIndexing(varargin{:});
    case 'isPossibleIndexedEntityName',
        varargout{1} = isPossibleIndexedEntityName(varargin{:});
    case 'getBaseVariableName',
        varargout{1} = getBaseVariableName(varargin{1});
    case 'assignmentPassthrough',
        varargout{1} = assignmentPassthrough(varargin{1});
    case 'createSpreadsheetValues',
        varargout{1} = createSpreadsheetValues(varargin{1});
    case 'reportValues',
        [varargout{1}, varargout{2}] = reportValues(varargin{:});
    case 'reportValuesCallback',
        reportValuesCallback(varargin{:});
    case 'reportValuesLogic',
        [varargout{1}, varargout{2}, varargout{3}, varargout{4}, ...
            varargout{5}, varargout{6}, varargout{7}, varargout{8}, ...
            varargout{9}, varargout{10}] = reportValuesLogic(varargin{:});
    case 'reportValueMetaInfo',
        varargout{1} = reportValueMetaInfo(varargin{:});
    case 'reportNonexistentValueMetaInfo',
        varargout{1} = reportNonexistentValueMetaInfo;
    case 'doHashedAssignment',
        varargout{1} = doHashedAssignment(varargin{:});
    case 'undoHashedAssignment',
        varargout{1} = undoHashedAssignment(varargin{:});
    case 'doVDSAssignment',
        varargout{1} = doVDSAssignment(varargin{:});
    case 'undoVDSAssignment',
        varargout{1} = undoVDSAssignment(varargin{:});
    case 'doMultiFieldAssignment',
        varargout{1} = doMultiFieldAssignment(varargin{:});
    case 'doMultiFieldAutoCopy',
        varargout{1} = doMultiFieldAutoCopy(varargin{:});
    case 'storeValue',
        storeValue(varargin{:});
        varargout = [];
    case 'retrieveAndClearValue',
        varargout{1} = retrieveAndClearValue(varargin{1});
    case 'getAllStorage',
        varargout{1} = getAllStorage;
    case 'clearValueSafely'
        clearValueSafely(varargin{:});
        varargout = [];
    case 'whosInformationForProperties',
        [varargout{1}, varargout{2}] = whosInformationForProperties(...
            varargin{1}, varargin{2}, varargin{3}, varargin{4});
    case 'getCurrentContextPropsAndPerms',
        [varargout{1}, varargout{2}] = getCurrentContextPropsAndPerms(varargin{1}, varargin{2});
    case 'structToWhosInformation',
        varargout{1} = structToWhosInformation(varargin{1});
    case 'createCellArrayValue',
        varargout{1} = createCellArrayValue(varargin{1});
    case 'getUnsupportedString',
        
        %%% WARNING: Due to incompatibilities with MIRtoolbox, some options
        %%% of Matlab Array Variable Editor are turned off.
        %%% To resume the standard behavior of the Array Variable Editor,
        %%% remove this copy of arrayviewfunc.m in MIRToolbox folder
        %%% (or simply uncomment the command below.)
        %%% In such case, please do not edit MIRtoolbox variable in the
        %%% Array Variable Editor, else you will obtain an infinite cascade
        %%% of figure windows that cannot be stopped unless the Array
        %%% Variable Editor window is closed.
        
        %varargout{1} = getUnsupportedString(varargin{:});
        
        %%% This above is the command you can uncomment if you prefer, at
        %%% your own risk..
        
    otherwise
        error('MATLAB:arrayviewfunc:UnknownOption', ...
            'Unknown command option %s', upper(whichcall));
end

%********************************************************************
function result = getData(x, format)
if (ischar(x))
    % First, make sure that there aren't a newline in the string.
    % If there are, we can't display it properly.
    if ~isempty(findstr(x, char(10))) || ~isempty(findstr(x, char(13)))
        result = [];
    else
        result = sprintf('%s', x);
    end
elseif iscellstr(x)
    % First, make sure that there are no newlines in the strings.
    % If there are, we can't display them properly.
    found = false;
    for i = 1:length(x)
        if ~isempty(findstr(x{i}, char(10))) || ~isempty(findstr(x{i}, char(13)))
            result = [];
            found = true;
        end
    end
    % have to pad with a space so java tokenizer will function
    % properly when a cell contains an empty string.
    if ~found
        result = sprintf('%s \n', x{:});
    end
elseif iscell(x)
    result = [];
else
    oldFormat = get(0, 'format');
    oldSpacing = get(0, 'formatspacing');

    set(0, 'format', format);
    set(0, 'formatspacing', 'compact');

    if length(x) > 1
        result = evalc('disp(x(:))');
    else
        result = evalc('disp(x)');
    end

    set(0, 'format', oldFormat);
    set(0, 'formatspacing', oldSpacing);
end;

%********************************************************************
function newRef = setData(var, coord, expr)
[var, passed] = setDataLogic(var, coord, expr);
if passed
    newRef = system_dependent(45, var);
else
    newRef = var;
end

%********************************************************************
function [var, passed] = setDataLogic(var, coord, expr)
try
    if ischar(var)
        var = expr;
    elseif iscellstr(var),
        var{coord{:}} = expr;
    else
        var(coord{:}) = eval(expr);
    end
    passed = true;
catch anError
    var = anError.message;
    passed = false;
end

%********************************************************************
function [var, pass] = setVarWidthLogic(var, width)
try
    sz = size(var);
    oldWidth = sz(2);

    if iscellstr(var),
        repl = {''};
    else
        repl = 0;
    end

    if width > oldWidth,
        var(:,end+1:width) = repl;
    elseif width < oldWidth,
        var(:,width+1:end) = [];
    end
    pass = true;
catch anError
    var = anError.message;
    pass = false;
end

%********************************************************************
function newRef = setVarWidth(var, width)
[var, pass] = setVarWidthLogic(var, width);
if pass
    newRef = system_dependent(45, var);
else
    newRef = var;
end;

%********************************************************************
function [var, pass] = setVarHeightLogic(var, height)
try
    sz = size(var);
    oldHeight = sz(1);

    if iscellstr(var),
        repl = {''};
    else
        repl = 0;
    end;

    if height > oldHeight,
        var(end+1:height,:) = repl;
    elseif height < oldHeight,
        var(height+1:end,:) = [];
    end;
    pass = true;
catch anError
    var = anError.message;
    pass = false;
end

%********************************************************************
function newRef = setVarHeight(var, height)
[var, pass] = setVarHeightLogic(var, height);
if pass
    newRef = system_dependent(45, var);
else
    newRef = var;
end;

%********************************************************************
function out = removeRowsOrColumns(orig, rowindices, colindices, direction, key)

toStore = [];
% Take care of the easy cases first
if isa(orig, 'char')
    % A char array (guaranteed to be 1xN)
    orig = '';
elseif strcmp(rowindices, ':')
    % Entire columns.  Rip them out.
    toStore = orig(:, colindices);
    orig(:, colindices) = [];
elseif strcmp(colindices, ':')
    % Entire rows.  Rip them out.
    toStore = orig(rowindices, :);
    orig(rowindices, :) = [];
else
    % User specified only CERTAIN cells.  More complicated.
    % We'll be removing the selected cells, and moving the
    % "neighbors" up or left, depending on the user's choice.
    empty = 0;
    if isa(orig, 'cell')
        empty = {[]};
    end
    [lastRow lastCol] = size(orig);
    numberOfRows = length(rowindices);
    numberOfCols = length(colindices);
    if strcmp(direction, 'up/down')
        for destRow = rowindices(1):lastRow
            sourceRow = destRow + numberOfRows;
            for colCounter = 1:numberOfCols
                destCol = colindices(colCounter);
                newValue = empty;
                if (sourceRow <= lastRow)
                    newValue = orig(sourceRow, destCol);
                end
                orig(destRow, destCol) = newValue;
            end
        end
    elseif strcmp(direction, 'left/right')
        for destCol = colindices(1):lastCol
            sourceCol = destCol + numberOfCols;
            for rowCounter = 1:numberOfRows
                destRow = rowindices(rowCounter);
                newValue = empty;
                if (sourceCol <= lastCol)
                    newValue = orig(destRow, sourceCol);
                end
                orig(destRow, destCol) = newValue;
            end
        end
    end
end
if numel(orig) == 0 && isnumeric(orig) && sum(size(orig)) > 0
    % Reduces the array to a 0x0 without changing its class.
    orig = repmat(orig, 0, 0);
end
if (~isempty(toStore))
    storeValue(key, toStore);
end
out = orig;

%********************************************************************
function out = insertRowsOrColumns(orig, rowindices, colindices, direction, key)

empty = 0;
if isa(orig, 'cell')
    empty = {[]};
else if isa(orig, 'struct')
        empty = createEmptyStruct(orig);
    end
end

[height width] = size(orig);

% Take care of the easy cases first
if isa(orig, 'char')
    % A char array (guaranteed to be 1xN)
    orig = '';
elseif strcmp(rowindices, ':')
    % Entire columns.  Shift all higher columns down and fill the selection
    % with 'empty' (whatever's appropriate for the data type).
    if strcmp(colindices, ':')
        colindices = 1:width;
    end
    numToShift = length(colindices);
    for i = (width + numToShift):-1:max([colindices numToShift+1])
        orig(:, i) = orig(:, i-numToShift);
    end
    if ~isempty(key) && keyExists(key)
        toReplace = retrieveAndClearValue(key);
        orig(:, colindices) = toReplace;
    else
        for i = colindices
            orig(:, i) = empty;
        end
    end
elseif strcmp(colindices, ':')
    % Entire rows.  Shift all higher rows to the left and fill the selection
    % with 'empty' (whatever's appropriate for the data type).
    % Dealt with the row is :, col is : case above.  Don't do it again.
    numToShift = length(rowindices);
    for i = (height + numToShift):-1:max([rowindices numToShift+1])
        orig(i, :) = orig(i-numToShift, :);
    end
    if ~isempty(key) && keyExists(key)
        toReplace = retrieveAndClearValue(key);
        orig(rowindices, :) = toReplace;
    else
        for i = rowindices
            orig(i, :) = empty;
        end
    end
else
    % User specified only CERTAIN cells.  More complicated.
    % We'll be moving the selected cells and their "neighbors"
    % down or to the right, depending on the user's choice.
    % Fill in the selected cells with 'empty' (whatever's
    % appropriate for the data type).

    % Move things around
    [lastRow lastCol] = size(orig);
    lastRowOfSelection = rowindices(end);
    lastColOfSelection = colindices(end);
    numberOfRows = length(rowindices);
    numberOfCols = length(colindices);
    if strcmp(direction, 'up/down')
        for sourceRow = lastRow:-1:rowindices(1)
            destRow = sourceRow + numberOfRows;
            for colCounter = 1:numberOfCols
                destCol = colindices(colCounter);
                newValue = empty;
                if (destRow > lastRowOfSelection)
                    newValue = orig(sourceRow, destCol);
                end
                orig(destRow, destCol) = newValue;
            end
        end
    elseif strcmp(direction, 'left/right')
        for sourceCol = lastCol:-1:colindices(1)
            destCol = sourceCol + numberOfCols;
            for rowCounter = 1:numberOfRows
                destRow = rowindices(rowCounter);
                newValue = empty;
                if (destCol > lastColOfSelection)
                    newValue = orig(destRow, sourceCol);
                end
                orig(destRow, destCol) = newValue;
            end
        end
    end

    % Zero out the selection...
    orig(rowindices, colindices) = empty;

    % "Patch up" cell arrays to preserve the char array nature of empties.
    if isa(orig, 'cell')
        [l, w] = size(orig);
        for i = 1:l
            for j = 1:w
                if isempty(orig{i, j})
                    orig(i, j) = {[]};
                end
            end
        end
    end
end
out = orig;

%********************************************************************
function in = renameField(in, oldFieldName, newFieldName)
if ~strcmp(oldFieldName, newFieldName)
    allNames = fieldnames(in);
    % Is the user renaming one field to be the name of another field?
    % Remember this.
    isOverwriting = ~isempty(find(strcmp(allNames, newFieldName), 1));
    matchingIndex = find(strcmp(allNames, oldFieldName));
    if ~isempty(matchingIndex)
        allNames{matchingIndex(1)} = newFieldName;
        in.(newFieldName) = in.(oldFieldName);
        in = rmfield(in, oldFieldName);
        if (~isOverwriting)
            % Do not attempt to reorder if we've reduced the number
            % of fields.  Bad things will result.  Let it go.
            in = orderfields(in, allNames);
        end
    end
end

%********************************************************************
function out = createEmptyStruct(in)
fields = fieldnames(in);
args = cell(1, 2*length(fields));
for inc = 1:length(fields)
    args{2*inc-1} = fields{inc};
    args{2*inc} = [];
end
out = struct(args{:});

%********************************************************************
function out = valueHasAppropriateIndexing(name, value) %#ok<INUSD>
out = false;
if length(name) >= 3
    special = getIndicesOfIndexingChars(name);
    if ~isempty(special) && special(1) ~= 1
        try
            eval(['value' name(special(1):end) ';']);
            out = true;
        catch anError %#ok<NASGU> ignore exceptions
        end
    end
end

%********************************************************************
function out = getBaseVariableName(name)
out = '';
done = false;
if isempty(name)
    done = true;
end
if isvarname(name)
    out = name;
    done = true;
end

if ~done
    special = getIndicesOfIndexingChars(name);
    if ~isempty(special)
        out = name(1:special(1)-1);
    end
end

%********************************************************************
function out = isPossibleIndexedEntityName(in)
out = false;
if length(in) >= 3
    special = getIndicesOfIndexingChars(in);
    if ~isempty(special)
        out = special(1) ~= 1 && special(end) ~= length(in);
    end
end

%********************************************************************
function out = getIndicesOfIndexingChars(in)
out = [];
if length(in) >= 2
    dots = strfind(in, '.');
    parens = strfind(in, '(');
    curleys = strfind(in, '{');
    out = sort([dots parens curleys]);
end

%********************************************************************
function result = assignmentPassthrough(in)
result = in;

%********************************************************************
function out = createSpreadsheetValues(in)
s = size(in);
out = cell(s);
for i = 1:s(1)
    for j = 1:s(2)
        try
            eval(['out{i, j} = ' in{i, j} ';']);
        catch err %#ok<NASGU>
            out{i, j} = in{i, j};
        end
    end
end

%********************************************************************
function [nameForAssign, toReport, startRow, startColumn, ...
    fullWidth, fullHeight, ...
    classID, customClassDescription, dims, recurse] = reportValuesLogic(...
    nameForAssign, value, startRow, startColumn, endRow, endColumn, recurse)
import com.mathworks.mlwidgets.array.data.FunctionHandleValue;
if nargin < 3
    startRow = 1;
    if nargin < 4
        startColumn = 1;
        if nargin < 5
            endRow = -1; %Patch it up below.
            if nargin < 6
                endColumn = -1; %Patch it up below.
                if nargin < 7
                    recurse = true;
                end
            end
        end
    end
end

s = size(value); % Treat this as read-only throughout.
fullHeight = s(1);
fullWidth  = s(2);
classID = getMLArrayRefClass(value);
customClassDescription = class(value);
dims = length(size(value));

% Calculate truncation limits
ischars = ischar(value);
trueEndRow = s(1);
if (endRow == -1) || ischars
    endRow = trueEndRow;
else
    endRow = min(endRow, trueEndRow);
end
trueEndColumn = s(2);
if (endColumn == -1) || ischars
    endColumn = trueEndColumn;
else
    endColumn = min(endColumn, trueEndColumn);
end
if startRow > endRow
    startRow = endRow;
end
if startColumn > endColumn
    startColumn = endColumn;
end

if endRow == 0 || endColumn == 0
    toReport = [];
else
    if isa(value, 'function_handle')
        toReport = func2str(value);
        if toReport(1) ~= '@'
            toReport = ['@' toReport];
        end
        toReport = FunctionHandleValue.valueOf(toReport);
    else
        toReport = value(startRow:endRow, startColumn:endColumn);
    end
end
if (ischars)
    fullWidth = 1;
    fullHeight = 1;
end

%********************************************************************
function reportValuesCallback(varargin)
% name, value, startRow, startColumn, endRow, endColumn, hashCode
import com.mathworks.mlwidgets.array.*;
import com.mathworks.widgets.spreadsheet.data.*;
import com.mathworks.mlwidgets.array.data.*;
isNewStyleObject = false;
uitableUsage = false;
if (ishandle(varargin{1}))
    uitableUsage = true;
    % Special case for uitable reuse
    % varargin{1} = uitable handle
    % varargin{2} = the property to use as the workspace value
    value = get(varargin{1}, varargin{2});
    % to match the contract, we are setting varargin{1} to "Data"
    % and varargin{2} to get(table, 'Data')
    varargin{1} = varargin{2};
    varargin{2} = value;
end
[nameForAssign, toReport, startRow, startColumn, ...
    fullWidth, fullHeight, ...
    classID, customClassDescription, dims] = ...
    reportValuesLogic(varargin{1:end-1});
signed = ~dataviewerhelper('isUnsignedIntegralType', toReport);
% Calls to ValueDataSection should take Java-isms (i.e. zero-based row and
% column specifiers)
if isempty(toReport)
    ca = ComplexArrayFactory.getEmptyInstance(...
        size(toReport, 1), size(toReport, 2), ...
        cast(0, class(toReport)));
else
    if (isnumeric(toReport))
        if isreal(toReport)
            if isscalar(toReport)
                if isfloat(toReport)
                    ca = ComplexScalarFactory.valueOf(toReport);
                else
                    ca = ComplexScalarFactory.valueOf(dataviewerhelper('upconvertIntegralType', toReport), signed);
                end
            else
                if isfloat(toReport)
                    ca = ComplexArrayFactory.valueOf(toReport);
                else
                    ca = ComplexArrayFactory.valueOf(dataviewerhelper('upconvertIntegralType', toReport), signed);
                end
            end
        else
            if isscalar(toReport)
                if isfloat(toReport)
                    ca = ComplexScalarFactory.valueOf(toReport, imag(toReport));
                else
                    converted = dataviewerhelper('upconvertIntegralType', toReport);
                    ca = ComplexScalarFactory.valueOf(real(converted), imag(converted), signed);
                end
            else
                if isfloat(toReport)
                    ca = ComplexArrayFactory.valueOf(toReport, imag(toReport));
                else
                    converted = dataviewerhelper('upconvertIntegralType', toReport);
                    ca = ComplexArrayFactory.valueOf(real(converted), imag(converted), signed);
                end
            end
        end
    elseif islogical(toReport)
        ca = toReport;
    elseif iscellstr(toReport)
        ca = toReport;
    elseif ischar(toReport)
        if size(toReport, 1) == 1
            ca = toReport;
        else
            ca = createReadOnlyValue(toReport);
        end
    elseif iscell(toReport)
        if uitableUsage
            % uitable doesn't want quoted or truncated strings
            ca = createCellArrayValue(toReport, false, Inf);
        else
            % Editor for cell arrays wants truncated, unquoted strings.
            % The caller will handle the quoting for us.
            ca = createCellArrayValue(toReport, false);
        end
    elseif isstruct(toReport) && numel(toReport) > 1
        ca = createStructArrayValue(toReport);
    else
        try
            isNewStyleObject = ~isempty(meta.class.fromName(customClassDescription));
        catch err %#ok<NASGU>
        end
        if isNewStyleObject
            ca = createNewObjectArrayValue(toReport, false);
        else
            ca = createReadOnlyValue(toReport);
        end
        
    end
end
vds = ValueDataSection(nameForAssign, startRow-1, startColumn-1, ...
    size(toReport, 1) + startRow - 1, size(toReport, 2) + startColumn - 1, ca);
vds.setUITableUsage(uitableUsage);
vmi = ValueMetaInfo(classID, getAttributes(varargin{2}), customClassDescription, isNewStyleObject, dims, ...
    fullWidth, fullHeight);
com.mathworks.mlwidgets.array.ValueTableModel.valueRequestCompleted(...
    varargin{end}, vds, vmi);

%********************************************************************
function [vds, vmi] = reportValues(varargin)
% name, value, startRow, startColumn, endRow, endColumn
import com.mathworks.mlwidgets.array.*;
import com.mathworks.widgets.spreadsheet.data.*;
import com.mathworks.mlwidgets.array.data.*;

isBracketedSummary = false;
if isa(varargin{2}, 'com.mathworks.widgets.spreadsheet.data.ValueSummary')
    obj = char(varargin{2}.toString);
    varargin{2} = obj;
    varargin{5} = size(obj, 1);
    varargin{6} = size(obj, 2);
    isBracketedSummary = true;
end    
    
[nameForAssign, toReport, startRow, startColumn, ...
    fullWidth, fullHeight, ...
    classID, customClassDescription, dims, recurse] = ...
    reportValuesLogic(varargin{:});
signed = ~dataviewerhelper('isUnsignedIntegralType', toReport);
% Calls to ValueDataSection should take Java-isms (i.e. zero-based row and
% column specifiers)
if isempty(toReport)
    ca = ComplexArrayFactory.getEmptyInstance(...
        size(toReport, 1), size(toReport, 2), ...
        cast(0, class(toReport)));
else
    if (isnumeric(toReport))
        if isreal(toReport)
            if isscalar(toReport)
                if isfloat(toReport)
                    ca = ComplexScalarFactory.valueOf(toReport);
                else
                    ca = ComplexScalarFactory.valueOf(dataviewerhelper('upconvertIntegralType', toReport), signed);
                end
            else
                if isfloat(toReport)
                    ca = ComplexArrayFactory.valueOf(toReport);
                else
                    ca = ComplexArrayFactory.valueOf(dataviewerhelper('upconvertIntegralType', toReport), signed);
                end
            end
        else
            if isscalar(toReport)
                if isfloat(toReport)
                    ca = ComplexScalarFactory.valueOf(toReport, imag(toReport));
                else
                    converted = dataviewerhelper('upconvertIntegralType', toReport);
                    ca = ComplexScalarFactory.valueOf(real(converted), imag(converted), signed);
                end
            else
                if isfloat(toReport)
                    ca = ComplexArrayFactory.valueOf(toReport, imag(toReport));
                else
                    converted = dataviewerhelper('upconvertIntegralType', toReport);
                    ca = ComplexArrayFactory.valueOf(real(converted), imag(converted), signed);
                end
            end
        end
    elseif islogical(toReport)
        ca = toReport;
    elseif iscellstr(toReport)
        ca = toReport;
    elseif ischar(toReport)
        if size(toReport, 1) == 1
            if isBracketedSummary
                ca = createReadOnlyValue(toReport);
            else
                ca = toReport;
            end
        else
            ca = createReadOnlyValue(toReport);
        end
    elseif isa(toReport, 'com.mathworks.mlwidgets.array.data.FunctionHandleValue')
        ca = toReport;
    elseif ~recurse
        ca = ReadOnlyValue.valueOf(workspacefunc('getshortvalue', toReport));
    elseif iscell(toReport)
        ca = createCellArrayValue(toReport);
    else
        ca = createReadOnlyValue(toReport);
    end
end
vds = ValueDataSection(nameForAssign, startRow-1, startColumn-1, ...
    size(toReport, 1) + startRow - 1, size(toReport, 2) + startColumn - 1, ca);
if nargout > 1
    vmi = ValueMetaInfo(classID, getAttributes(varargin{2}), customClassDescription, dims, ...
        fullWidth, fullHeight);
end

%********************************************************************
function out = reportValueMetaInfo(varargin)
% Arguments are:
% name, value, startRow, startColumn, endRow, endColumn
import com.mathworks.mlwidgets.array.*;

in = varargin{1};
clazz = class(in);
type = getMLArrayRefClass(in);
[height, width] = size(in);
attributes = getAttributes(in);
numOfDims = length(size(in));
isNewStyleObject = false;
if (type == com.mathworks.jmi.types.MLArrayRef.mxUNKNOWN_CLASS)
    try
        if ~isempty(meta.class.fromName(clazz))
            isNewStyleObject = true;
        end
    catch err %#ok<NASGU>
    end
end
out = ValueMetaInfo(type, attributes, clazz, isNewStyleObject, numOfDims, width, height);

%********************************************************************
function out = createReadOnlyValue(toReport) %#ok<INUSD>
out = com.mathworks.mlwidgets.array.data.ReadOnlyValue.valueOf(evalc('disp(toReport)'));

%********************************************************************
function type = getMLArrayRefClass(in)
import com.mathworks.jmi.types.MLArrayRef;
clazz = class(in);
switch clazz
    case 'double',
        type = MLArrayRef.mxDOUBLE_CLASS;
    case 'single',
        type = MLArrayRef.mxSINGLE_CLASS;
    case 'uint8',
        type = MLArrayRef.mxUINT8_CLASS;
    case 'int8',
        type = MLArrayRef.mxINT8_CLASS;
    case 'uint16',
        type = MLArrayRef.mxUINT16_CLASS;
    case 'int16',
        type = MLArrayRef.mxINT16_CLASS;
    case 'uint32',
        type = MLArrayRef.mxUINT32_CLASS;
    case 'int32',
        type = MLArrayRef.mxINT32_CLASS;
    case 'uint64',
        type = MLArrayRef.mxUINT64_CLASS;
    case 'int64',
        type = MLArrayRef.mxINT64_CLASS;
    case 'logical'
        type = MLArrayRef.mxLOGICAL_CLASS;
    case 'cell',
        type = MLArrayRef.mxCELL_CLASS;
    case 'struct',
        type = MLArrayRef.mxSTRUCT_CLASS;
    case 'char',
        type = MLArrayRef.mxCHAR_CLASS;
    otherwise
        type = MLArrayRef.mxUNKNOWN_CLASS;
end

%********************************************************************
function out = reportNonexistentValueMetaInfo
out = com.mathworks.mlwidgets.array.ValueMetaInfo.getNonExistentInstance;

%********************************************************************
function attributes = getAttributes(in)
% This is the code we should run.  Instead, we run a hand-optimized
% version for performance.
%import com.mathworks.jmi.types.*;
%attributes = MLArrayRef.COMPLEX * (~isreal(in)) + ...
%    MLArrayRef.SPARSE * issparse(in);
if ~isreal(in)
    if issparse(in)
        attributes = 3;
    else
        attributes = 1;
    end
else
    if issparse(in)
        attributes = 2;
    else
        attributes = 0;
    end
end

function out = doHashedAssignment(var, rhs, key)
oldValue = var;
out = rhs;
storeValue(key, oldValue);

function out = undoHashedAssignment(key)
out = retrieveAndClearValue(key);

function out = doVDSAssignment(var, rhs, key)
oldValue = var;
out = rhs;
storeValue(key, oldValue);

function out = undoVDSAssignment(key)
out = retrieveAndClearValue(key);

function var = doMultiFieldAssignment(var, fieldNames, values)
for i=1:length(fieldNames)
    var.(fieldNames{i}) = values{i};
end

function var = doMultiFieldAutoCopy(var, fieldNames)
for i = 1:length(fieldNames)
    var.(workspacefunc('getcopyname', fieldNames{i}, fields(var))) = ...
        var.(fieldNames{i});
end


%********************************************************************
function out = createCellArrayValue(in, varargin)
import com.mathworks.mlwidgets.array.data.CellArrayValue;
vds = javaArray('com.mathworks.mlwidgets.array.ValueDataSection', size(in));
optargin = size(varargin,2);
for i=1:size(in, 1)
    for j = 1:size(in, 2)
        thisIn = in{i, j};
        if (numel(thisIn) < 11)
            value = reportValues('', thisIn, 1, 1, size(thisIn, 1), size(thisIn, 2), false);
            if (optargin > 1) 
                value.setUITableUsage(true);
            end
            vds(i, j) = value;
        else
            obj = workspacefunc('getshortvalueobjectj', thisIn, varargin{:});
            if ~isa(obj, 'com.mathworks.widgets.spreadsheet.data.ValueSummary')
                obj = char(obj.toString);
            end
            value = reportValues('', obj, 1, 1, size(thisIn, 1), size(thisIn, 2), false);
            if (optargin > 1) 
                value.setUITableUsage(true);
            end
            vds(i, j) = value;
        end
    end

end
out = CellArrayValue.valueOf(vds);

%********************************************************************
function out = createNewObjectArrayValue(in, varargin)
import com.mathworks.mlwidgets.array.data.CellArrayValue;
vds = javaArray('com.mathworks.mlwidgets.array.ValueDataSection', size(in));
optargin = size(varargin,2);
for i=1:size(in, 1)
    for j = 1:size(in, 2)
        thisIn = in(i, j);
        if (numel(thisIn) < 11)
            value = reportValues('', thisIn, 1, 1, size(thisIn, 1), size(thisIn, 2), false);
            if (optargin > 1) 
                value.setUITableUsage(true);
            end
            vds(i, j) = value;
        else
            obj = workspacefunc('getshortvalueobjectj', thisIn, varargin{:});
            if ~isa(obj, 'com.mathworks.widgets.spreadsheet.data.ValueSummary')
                obj = char(obj.toString);
            end
            value = reportValues('', obj, 1, 1, size(thisIn, 1), size(thisIn, 2), false);
            if (optargin > 1) 
                value.setUITableUsage(true);
            end
            vds(i, j) = value;
        end
    end

end
out = CellArrayValue.valueOf(vds);

%********************************************************************
function out = createStructArrayValue(in)
import com.mathworks.mlwidgets.array.data.CellArrayValue;
vds = javaArray('com.mathworks.mlwidgets.array.ValueDataSection', size(in));
for i=1:size(in, 1)
    for j = 1:size(in, 2)
        thisIn = in(i, j);
        vds(i, j) = reportValues('', in(i, j), 1, 1, size(thisIn, 1), size(thisIn, 2), false);
    end
end
out = CellArrayValue.valueOf(vds);


%********************************************************************
function out = valueStorage(behavior, key, value)
% behavior 1 = retrieveAndClearValue
% behavior 2 = store
% behavior 0 = getAllStorage
% behavior 3 = clearValueSafely
persistent storage
switch (behavior)
    case 2
        if isempty (storage)
            storage.empty = [];
        end
        storage.(key) = value;
    case 1
        out = storage.(key);
        storage = rmfield(storage, key);
    case 0
        out = storage;
    case 3
        if isfield(storage, key)
            storage = rmfield(storage, key);
        end
end

%********************************************************************
function ret = storeValue(key, value)
if ~isempty(key)
    valueStorage(2, key, value);
end
ret = [];

%********************************************************************
function out = retrieveAndClearValue(key)
out = valueStorage(1, key);

%********************************************************************
function clearValueSafely(key)
valueStorage(3, key);

%********************************************************************
function out = keyExists(key)
out = isfield(getAllStorage, key);

%********************************************************************
function out = getAllStorage
out = valueStorage(0);

%********************************************************************
function [whosInformation, propertyAttributes] = ...
    whosInformationForProperties(object, values, names, writablesBasedOnContext)
import com.mathworks.mlwidgets.workspace.*;
import com.mathworks.mlwidgets.array.mcos.*;

mc = metaclass(object);
metapropsArray = mc.Properties;
if isempty(metapropsArray) || isempty(names)
    whosInformation = WhosInformation.getInstance;
    propertyAttributes = PropertyAttributes.getInstance;
    return
end

filteredMPA = cell(size(names));
for i = 1:length(names)
    thisName = names(i);
    for j = 1:length(metapropsArray)
        if strcmp(thisName, metapropsArray{j}.Name)
            filteredMPA{i} = metapropsArray{j};
            break;
        end 
    end
end

sizes = cell(size(names));
bytes = zeros(size(names));
classes = cell(size(names));
isSparse = zeros(size(names));
isComplex = zeros(size(names));
nestingFunctions = cell(size(names));

getAccess = javaArray('com.mathworks.mlwidgets.array.mcos.Access', length(names));
setAccess = javaArray('com.mathworks.mlwidgets.array.mcos.Access', length(names));
constant = false(size(names));
transient = false(size(names));
description = cell(size(names));
ddescription = cell(size(names));
for i = 1:length(names)
    thisVal = values{i};
    thisMetaProp = filteredMPA{i};
    sizes{i} = int64(size(thisVal));
    classes{i} = class(thisVal);
    isSparse(i) = issparse(thisVal);
    isComplex(i) = isnumeric(thisVal) && ~isreal(thisVal);
    nestingFunctions{i} = 'foo';
    
    getAccess(i) = getGetAccessForProperty(thisMetaProp);
    setAccess(i) = getSetAccessForProperty(thisMetaProp);
    constant(i) = thisMetaProp.Constant;
    transient(i) = thisMetaProp.Transient;
    description{i} = thisMetaProp.Description;
    ddescription{i} = thisMetaProp.DetailedDescription;
end

nestingLevels = ones(size(names));
nesting = NestingInformation.getInstances(nestingLevels, nestingFunctions);
isPersistent = zeros(size(names));
isGlobal = zeros(size(names));

whosInformation = WhosInformation(names, sizes, bytes, classes, ...
    isPersistent, isGlobal, isSparse, isComplex, nesting);
propertyAttributes = PropertyAttributes(names, getAccess, setAccess, ...
    constant, writablesBasedOnContext, transient, description, ddescription);

%********************************************************************
function [names, writables] = getCurrentContextPropsAndPerms(object, ...
    classnameBeingDebugged)
% classnameBeingDebugged should come from running
% mfilename('class')
% in the context of the debugged function.

mc = metaclass(object);
if isempty(classnameBeingDebugged)
    currentContextHasPrivateAccess = false;
    currentContextHasProtectedAccess = false;
else
    currentContextHasPrivateAccess = strcmp(classnameBeingDebugged, mc.Name);
    currentContextHasProtectedAccess = currentContextHasPrivateAccess || ...
        meta.class.fromName(classnameBeingDebugged) < mc;
end

metapropsArray = mc.Properties;
names = {};
writables = false(0);
for i = 1:length(metapropsArray)
    thisProp = metapropsArray{i};
    append = false;
    writable = false;
    if ~thisProp.Hidden;
        if currentContextHasPrivateAccess
            append = true;
            writable = true;
        else
            if currentContextHasProtectedAccess && getAccessAvailableToProtected(thisProp)
                append = true;
                writable = setAccessAvailableToProtected(thisProp);
            else
                % Only has public assess
                if getAccessAvailableToPublic(thisProp)
                    append = true;
                    writable = setAccessAvailableToPublic(thisProp);
                end
            end
        end
    end
    if append
        names{end+1} = thisProp.Name; %#ok<AGROW>
        writables(end+1) = writable; %#ok<AGROW>
    end
    names = names';
    writables = writables';
end

%********************************************************************
function available = getAccessAvailableToProtected(metaProperty)
available = ~strcmp(metaProperty.GetAccess, 'private');

%********************************************************************
function available = getAccessAvailableToPublic(metaProperty)
available = strcmp(metaProperty.GetAccess, 'public');

%********************************************************************
function available = setAccessAvailableToProtected(metaProperty)
available = ~strcmp(metaProperty.SetAccess, 'private');

%********************************************************************
function available = setAccessAvailableToPublic(metaProperty)
available = strcmp(metaProperty.SetAccess, 'public');

%********************************************************************
function access = getGetAccessForProperty(propertyMetaInfo)
import com.mathworks.mlwidgets.array.mcos.Access;
switch propertyMetaInfo.GetAccess
    case 'public'
        access = Access.PUBLIC;
    case 'protected'
        access = Access.PROTECTED;
    otherwise
        access = Access.PRIVATE;
end

%********************************************************************
function access = getSetAccessForProperty(propertyMetaInfo)
import com.mathworks.mlwidgets.array.mcos.Access;
switch propertyMetaInfo.SetAccess
    case 'public'
        access = Access.PUBLIC;
    case 'protected'
        access = Access.PROTECTED;
    otherwise
        access = Access.PRIVATE;
end

%********************************************************************
function out = structToWhosInformation(object)
import com.mathworks.mlwidgets.workspace.*;

names = fields(object);
if isempty(names)
    out = com.mathworks.mlwidgets.workspace.WhosInformation.getInstance;
    return
end
sizes = cell(size(names));
bytes = zeros(size(names));
classes = cell(size(names));
isSparse = zeros(size(names));
isComplex = zeros(size(names));
nestingFunctions = cell(size(names));
for i = 1:length(names)
    thisVal = object.(names{i});
    sizes{i} = int64(size(thisVal));
    classes{i} = class(thisVal);
    isSparse(i) = issparse(thisVal);
    isComplex(i) = isnumeric(thisVal) && ~isreal(thisVal);
    nestingFunctions{i} = 'foo';
end
nestingLevels = ones(size(names));
nesting = NestingInformation.getInstances(nestingLevels, nestingFunctions);
isPersistent = zeros(size(names));
isGlobal = zeros(size(names));

out = WhosInformation(names, sizes, bytes, classes, isPersistent, isGlobal, isSparse, isComplex, nesting);

%********************************************************************
function out = getUnsupportedString(val, maxElements, tooLargeMessage) %#ok<INUSD>

if numel(val) > maxElements
    out = evalc('disp(char(tooLargeMessage))');
else
    out = evalc('display(val)');
end