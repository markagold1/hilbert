function y = sosfilter(sos,x)
% Usage: y = sosfilter(sos,x)
%
% Filter through cascade of second order sections (SOS)
%
%   sos.............N-by-6 array, 1 row per SOS stage
%   x...............1-d array containing input to filter
%   y...............1-d array containing filtered output
%

    for kk = 1:size(sos,1)
        b = sos(kk,1:3);
        a = sos(kk,4:6);
        y = filter(b,a,x);
        x = y;
    end

end % function
