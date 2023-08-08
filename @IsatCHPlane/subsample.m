function [Isat_diff_temp,t_temp] = subsample(obj,varargin)

Inpt = inputParser;

addRequired(Inpt, 'Time_Series', @(obj) isobject(obj))
addParameter(Inpt, 'Subsample_Dimension',8, @(s) isnumeric(s))
addParameter(Inpt, 'Start_Time',5, @(tstart) isnumeric(tstart))
addParameter(Inpt, 'Stop_Time',15, @(tstop) isnumeric(tstop))
addParameter(Inpt, 'End_Radius',0.9, @(rend) isnumeric(rend) && rend <= 2.7)

parse(Inpt, obj, varargin{:})

s = Inpt.Results.Subsample_Dimension; %Default value of 8 seems to work best as subsample dimension

t1 = find(obj.data.t >= Inpt.Results.Start_Time,1);
t2 = find(obj.data.t >= Inpt.Results.Stop_Time,1);
% t3 = length(obj.data.t);

rend = find(obj.data.r >= Inpt.Results.End_Radius,1);
%Subsample the time segment by skipping s points in the original time
%series
for p = 1:1:rend
    for n = 1:1:10
        k = 0;
        for m = t1:s:t2
            k = k+1; %Create a new index for the subsampled signal
            Isat_diff_temp(n,k,p) = obj.data.Isat_diff(n,m,p);
            t_temp(k) = obj.data.t(m); %Create a subsampled time series
        end
    end
end
end

