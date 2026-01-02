function H = sosfreqz(sos,f,fs)
% Usage: H = sosfreqz(sos,f,fs)

    if exist('freqz') == 2
        rspz = @(b,a,f,fs) freqz(b,a,f,fs);
    else
        rspz = @(b,a,f,fs) pm_freqz(b,a,f,fs);
    end

    B = sos(:,1:3);
    A = sos(:,4:6);
    H = ones(numel(f),1);
    for kk = 1:size(B,1)
        b = B(kk,:);
        a = A(kk,:);
        H = H .* shiftdim(rspz(b,a,f,fs));
    end

end % function
