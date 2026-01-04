function uw = pm_unwrap(w)
% Usage: uw = pm_unwrap(w)
%
% Phase unwrapper.
%
%  w.........1d array of phase values to unwrap (radians)
%  uw........1d array of unwrapped phase values (radians)
%
    if numel(w) == 1
        uw = w;
        return
    end
    if ~isvector(w)
        error('Input must be a vector.');
    end
    [uw,ns] = shiftdim(w);
    ix = isfinite(uw);
    uw(ix) = do_unwrap(uw(ix));
    uw = shiftdim(uw,-ns);

end % function

function uw = do_unwrap(w)
    uw = w;
    dw = diff(uw);
    comp = dw ./ (2*pi);
    comp(abs(rem(comp,1)) <= 0.5) = fix(comp(abs(rem(comp,1)) <= 0.5));
    comp = round(comp);
    comp(abs(dw) < pi) = 0;
    uw(2:end) = uw(2:end) - (2*pi) * cumsum(comp(:));
end % function
