%% ASEN 3112 Structures Lab 3 - computeFrequency.m
% Determines the frequency of a given data set. Also saves the indicies of
% the data set where frequency was computed, and the time values at those
% indices.
%
%   Author: Andrew Thompson
%   Created: 11/14/20 Edited: 11/14/20
%
%   Parameters:     time [N x 1] <double>
%                   data [N x 1] <double>
%   Returns:        frequency    <double>
%                   tOut         <double>

function [frequency, tOut] = computeFrequency(time, data)
    %%
    % Initialize full cycle counter
    wCount = 1;
    % Initialize half cycle counter
    crossCount = 0;
    % Initialize reference data point
    refPt = 0;
    % Initialize reference time
    time1 = time(1);
    % For all data points
    for index = 2:length(data)
        % Get data point 1
        pt1 = data(index-1);
        % Get data point 2
        pt2 = data(index);
        % If the points are on seperate sides of the reference
        if pt1 < refPt && pt2 > refPt || pt1 > refPt && pt2 < refPt
            % Increment the reference cross counter
            crossCount = crossCount + 1;
            % If it is the second time crossed
            if mod(crossCount, 2) == 0
                % Calculate the frequency
                frequency(wCount) = 1/(time(index) - time1);  %#ok<*AGROW>
                % Document current time
                tOut(wCount) = time1; 
                % Increment the frequency count
                wCount = wCount + 1;
                % Set a new time
                time1 = time(index);
            end     
        end
    end
    % Mean time of period
    tOut = mean(tOut);
    % Mean frequency of period
    frequency = mean(frequency);
end