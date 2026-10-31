function Migration_DP_Ver3(dataFolder, outputFolder)
%MIGRATION_DP_VER3  Cell-migration trajectory analysis (ImageJ export).
%
%   MIGRATION_DP_VER3()                          % analyse the folder set below
%   MIGRATION_DP_VER3(dataFolder)                % analyse a given folder
%   MIGRATION_DP_VER3(dataFolder, outputFolder)  % choose where to write
%
%   INPUT
%   One folder containing one Excel file per cell (ImageJ "Manual Tracking"
%   style results export).  Four columns are used:
%       B : time-point index
%       C : X position (px)
%       D : Y position (px)
%       E : displacement between this and the previous time point
%           (ImageJ's default factor; the first time point carries the
%            -1 sentinel and is therefore never used as a weight)
%
%   OUTPUT (written into <outputFolder>)
%   List1.xlsx : file name | mean migration angle | mean migration distance |
%                distance-weighted migration angle
%   List2.xlsx : displacement of every time point relative to the FIRST time
%                point, one column per file.  Row 1 holds the header, column 1
%                holds the time-point numbers 2..N, column j+1 holds the
%                values of file j.
%   List3.xlsx : origin-corrected coordinates (um).  The corrected X values
%                of every file are stacked in column 1 (rows 2..21 for file 1,
%                22..41 for file 2, ...), the corrected Y values of file n go
%                to column n+1, the file names go to row 1 of column n+1.
%   List4.xlsx : per-time-point diagnostic table (only if writeDiagnostic),
%                so that every single angle can be checked by hand.
%
%   WHY THE "MEAN MIGRATION ANGLE" USUALLY GOES WRONG
%   A time point at which the cell did not move (dX = 0 AND dY = 0) has no
%   migration angle.  The classic mistake is to store such a time point as a
%   placeholder (0) inside a fixed-size angle array, or to keep the angle
%   array of the previous file / of the previous run, and then to take
%   mean() over the WHOLE array.  The untouched time points then enter the
%   average as 0 (or as a stale angle belonging to another cell) and drag
%   the result down.  A file in which the cell never trembles is unaffected,
%   which is why one file can come out right while all the others are wrong.
%
%   Note that column 4 (the distance-weighted angle) is NOT affected by the
%   same mistake, because a time point where the cell did not move carries a
%   weight of exactly 0 in column E and therefore cancels out of
%   sum(w.*a)/sum(w).  That asymmetry - only column 2 wrong - is the
%   fingerprint of this bug.
%
%   This version therefore
%     * analyses every file completely independently: no accumulator, no
%       array and no variable survives from one file - or one run - to the
%       next (the whole program is a function, so the base workspace can
%       never leak in),
%     * keeps an undefined angle as NaN and removes it from EVERY statistic
%       (mean angle and weighted angle) instead of counting it as 0,
%     * never uses "angle ~= 0" as a test for "this angle exists", because a
%       perfectly valid angle may legitimately be exactly 0 deg (purely
%       horizontal movement).
%

%% ======================= 1. SETTINGS ====================================
if nargin < 1 || isempty(dataFolder)
    dataFolder = 'D:\???';% !!!!!Please put the excel files of cell coordinates from ImageJ in a foler and input the path of the folder here
end
if nargin < 2 || isempty(outputFolder)
    outputFolder = fullfile(dataFolder, 'Analysis_Output');
end

colTP = 2;                          % column B : time-point index
colX  = 3;                          % column C : X position (px)
colY  = 4;                          % column D : Y position (px)
colE  = 5;                          % column E : displacement to previous point
scalePxToUm = 1.61;                 % px -> um, used for List2 and List3
useSquaredDisplacement = false;     % false : sqrt(dx^2+dy^2)  (as in the spec)
                                     % true  : (dx^2+dy^2)     (true square)
moveTol = 0;                        % |dx| <= tol AND |dy| <= tol  =>  no move
writeDiagnostic = true;             % write List4 (per time point angles)

if ~exist(dataFolder, 'dir')
    error('Migration_DP_Ver3:folder', 'Folder does not exist: %s', dataFolder);
end
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% ======================= 2. FILE LIST ==================================
listing = dir(fullfile(dataFolder, '*.csv*'));
listing = listing(~[listing.isdir]);
if isempty(listing)
    error('Migration_DP_Ver3:noFiles', 'No Excel files found in %s', dataFolder);
end
[~, sortIdx] = sort({listing.name});
listing = listing(sortIdx);
nFiles = numel(listing);

fprintf('Analysing %d file(s) in %s\n', nFiles, dataFolder);

%% ======================= 3. PER FILE ANALYSIS ==========================
% Everything is recomputed for every file, nothing is carried over.
res = cell(nFiles, 1);
for k = 1:nFiles
    res{k} = analyseTrack(fullfile(dataFolder, listing(k).name), ...
                          colTP, colX, colY, colE, scalePxToUm, ...
                          useSquaredDisplacement, moveTol);
    res{k}.name = listing(k).name;
end

fileNames = cell(nFiles, 1);
for k = 1:nFiles
    fileNames{k} = res{k}.name;
end

%% ======================= 4. LIST 1 =====================================
%  column 1 : file name
%  column 2 : mean migration angle (deg) over the time points that have one
%  column 3 : mean migration distance (mean of E over time points 2..N)
%  column 4 : distance weighted migration angle (deg)
list1 = zeros(nFiles, 3);
for k = 1:nFiles
    list1(k,1) = res{k}.meanAngle;
    list1(k,2) = res{k}.meanDistance;
    list1(k,3) = res{k}.weightedAngle;
end

%% ======================= 5. LIST 2 =====================================
%  row 1             : 'Time point' | file name | file name | ...
%  column 1, row 2+ : time point numbers 2..N
%  column j+1, row 2+: displacement of file j at that time point, measured
%                      from the FIRST time point of the same file
maxTP = 0;
for k = 1:nFiles
    maxTP = max(maxTP, res{k}.nTP);
end
nValues = maxTP - 1;                       % time points 2..N
list2 = nan(1 + nValues, nFiles + 1);
list2(2:end, 1) = (2:maxTP)';              % time point numbers
for k = 1:nFiles
    r = res{k};
    nVal = min(nValues, r.nTP - 1);
    if nVal > 0
        list2(2:1+nVal, k+1) = r.relDisp(2:1+nVal);
    end
end

%% ======================= 6. LIST 3 =====================================
%  X values of every file are stacked in column 1: block k occupies rows
%  2+(k-1)*maxTP ... 1+k*maxTP (= rows 2..21 for file 1, 22..41 for file 2).
%  Y values of file k go to column k+1, the same rows.
list3 = nan(1 + nFiles*maxTP, nFiles + 1);
for k = 1:nFiles
    r = res{k};
    rowStart = 2 + (k-1)*maxTP;
    nVal = min(maxTP, r.nTP);
    list3(rowStart:rowStart+nVal-1, 1)   = r.xCorr(1:nVal);
    list3(rowStart:rowStart+nVal-1, k+1) = r.yCorr(1:nVal);
end

%% ======================= 7. WRITE ======================================
writeCellFile(fullfile(outputFolder, 'List1.xlsx'), ...
              [fileNames, num2cell(list1)], 'Results');

writeNumericFile(fullfile(outputFolder, 'List2.xlsx'), list2, 'Results');
writeHeaderRow(fullfile(outputFolder, 'List2.xlsx'), ...
               [{'Time point'}, fileNames(:)'], 'Results');

writeNumericFile(fullfile(outputFolder, 'List3.xlsx'), list3, 'Results');
writeHeaderRow(fullfile(outputFolder, 'List3.xlsx'), ...
               [{'X (um)'}, fileNames(:)'], 'Results');

if writeDiagnostic
    writeDiagnosticFile(fullfile(outputFolder, 'List4_AngleCheck.xlsx'), res);
end

%% ======================= 8. REPORT =====================================
fprintf('Done.  Results written to %s\n\n', outputFolder);
fprintf('%-26s %12s %14s %14s %10s\n', ...
        'file', 'mean angle', 'mean dist', 'weighted', 'points used');
for k = 1:nFiles
    fprintf('%-26s %12.3f %14.3f %14.3f %10s\n', ...
            res{k}.name, res{k}.meanAngle, res{k}.meanDistance, ...
            res{k}.weightedAngle, sprintf('%d/%d', res{k}.nValid, res{k}.nTP-1));
end

end % =====================================================================


%% ======================= LOCAL FUNCTIONS ===============================


function R = analyseTrack(fileName, colTP, colX, colY, colE, scale, useSquared, tol)
%ANALYSETRACK  Analyse one trajectory file.  Everything lives inside this
%               function, so no state can survive from one file to the next.

    % ---------------- read the file ----------------
    try
        raw = readmatrix(fileName);            % MATLAB R2019a and later
    catch
        raw = xlsread(fileName);               % legacy fallback
    end
    if isempty(raw)
        error('Migration_DP_Ver3:empty', 'Could not read any numbers from %s', fileName);
    end
    raw = double(raw);
    nCol = size(raw, 2);

    tp = raw(:, min(colTP, nCol));
    x  = raw(:, min(colX,  nCol));
    y  = raw(:, min(colY,  nCol));
    e  = raw(:, min(colE,  nCol));

    % Keep the real data rows only.  Row 1 of an ImageJ export is the
    % header (readmatrix returns NaN there) and the last row of some
    % exports is empty.
    good = isfinite(tp) & tp >= 1 & isfinite(x) & isfinite(y);
    firstRow = find(good, 1, 'first');
    if isempty(firstRow)
        error('Migration_DP_Ver3:noData', 'No time points found in %s', fileName);
    end
    lastRow = find(good, 1, 'last');
    tp = tp(firstRow:lastRow);
    x  = x(firstRow:lastRow);
    y  = y(firstRow:lastRow);
    e  = e(firstRow:lastRow);
    nTP = numel(tp);

    % ---------------- step by step displacement ----------------
    dx = zeros(nTP, 1);
    dy = zeros(nTP, 1);
    dx(2:end) = diff(x);
    dy(2:end) = diff(y);
    stepDist = hypot(dx, dy);

    % ---------------- migration angle of every time point ----------------
    % angle = acute angle between the step and the horizontal axis:
    %         0 deg = parallel to X, 90 deg = perpendicular to X.
    % The first time point has no predecessor, and a time point where the
    % cell did not move has no direction at all; both stay NaN and are
    % excluded from every statistic.
    moved = (abs(dx) > tol | abs(dy) > tol) & isfinite(dx) & isfinite(dy);
    angle = nan(nTP, 1);
    idx = find(moved);
    angle(idx) = atan2(abs(dy(idx)), abs(dx(idx))) * 180/pi;

    % ---------------- statistics ----------------
    ang = angle(2:end);                    % time points 2..N
    valid = ~isnan(ang);
    nValid = sum(valid);

    if nValid > 0
        meanAngle = mean(ang(valid));
    else
        meanAngle = NaN;
    end

    % Weights: the ImageJ displacement of the SAME step (column E).
    % Negative values (the -1 sentinel of the first time point) are dropped.
    w = e(2:end);
    w(~isfinite(w)) = NaN;
    w(w < 0) = NaN;

    % Mean migration distance = mean of E over the time points 2..N.
    if any(isfinite(w))
        meanDistance = mean(w(isfinite(w)));
    else
        meanDistance = NaN;
    end

    % Distance weighted migration angle.  Time points without an angle are
    % left out of numerator AND denominator; time points where the cell
    % did not move carry a weight of 0 and drop out on their own.
    if nValid > 0
        angv = ang(valid);
        wv = w(valid);
        wv(~isfinite(wv)) = 0;
        if sum(wv) > 0
            % (:) on both sides: keeps the product a vector even if the
            % orientation of the two columns ever differs
            weightedAngle = sum(wv(:) .* angv(:)) / sum(wv);
        else
            weightedAngle = NaN;
        end
    else
        weightedAngle = NaN;
    end

    % ---------------- origin corrected coordinates ----------------
    xCorr = (x - x(1)) * scale;
    yCorr = (y - y(1)) * scale;

    % ---------------- displacement from the first time point ----------------
    % same unit as List3 (um), so that both lists can be compared directly
    relDisp = hypot(x - x(1), y - y(1)) * scale;
    if useSquared
        relDisp = relDisp .^ 2;
    end

    R = struct('name', '', 'nTP', nTP, 'nValid', nValid, ...
               'meanAngle', meanAngle, 'meanDistance', meanDistance, ...
               'weightedAngle', weightedAngle, 'tp', tp, 'x', x, 'y', y, ...
               'e', e, 'dx', dx, 'dy', dy, 'stepDist', stepDist, ...
               'angle', angle, 'xCorr', xCorr, 'yCorr', yCorr, ...
               'relDisp', relDisp);
end


function writeCellFile(fileName, C, sheetName)
%WRITECELLFILE  Write a cell array (numbers and/or text) to an xlsx file.
    if exist('writecell') ~= 0
        writecell(C, fileName, 'Sheet', sheetName);
    else
        xlswrite(fileName, C, sheetName);
    end
end


function writeNumericFile(fileName, M, sheetName)
%WRITENUMERICFILE  Write a numeric matrix (NaN becomes an empty cell).
    if exist('writematrix') ~= 0
        writematrix(M, fileName, 'Sheet', sheetName);
    else
        xlswrite(fileName, M, sheetName);
    end
end


function writeHeaderRow(fileName, headerCell, sheetName)
%WRITEHEADERROW  Overwrite row 1 of an existing sheet with text.  The
%                 numeric matrix written before has an all-NaN first row, so
%                 nothing is lost.
    if exist('writecell') ~= 0
        writecell(headerCell, fileName, 'Sheet', sheetName);
    else
        xlswrite(fileName, headerCell, sheetName);
    end
end


function writeDiagnosticFile(fileName, res)
%WRITEDIAGNOSTICFILE  One block per input file with every single time
%                      point, so the angle calculation can be verified by
%                      hand against the spreadsheet.
    blocks = {};
    for k = 1:numel(res)
        r = res{k};
        head = {'Time point', 'X (px)', 'Y (px)', 'dX', 'dY', ...
                'Step distance', 'Angle (deg)', 'Used in mean', 'Weight (E)'};
        title = repmat({' '}, 1, 9);
        title{1} = ['### ' r.name];
        body = cell(r.nTP, 9);
        for i = 1:r.nTP
            if i == 1
                used = 'no (first point)';
            elseif isnan(r.angle(i))
                used = 'no (cell did not move)';
            else
                used = 'yes';
            end
            body(i,:) = {r.tp(i), r.x(i), r.y(i), r.dx(i), r.dy(i), ...
                         r.stepDist(i), r.angle(i), used, r.e(i)};
        end
        blocks{end+1} = title;
        blocks{end+1} = head;
        blocks{end+1} = body;
        blocks{end+1} = repmat({' '}, 1, 9);
    end
    writeCellFile(fileName, vertcat(blocks{:}), 'Results');
end