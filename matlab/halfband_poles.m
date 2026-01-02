function [P,P0,P1] = halfband_poles(wp,As)
%
% Calculate the poles of a parallel allpass Half-Band IIR filter using
% the method presented in the Appendix of [1].
%
% wp............Passband edge frequency radians, normalized to fs=2pi
% As............Stopband attenuation dB
% P.............Array of magnitude-squared poles
%               P is length (N-1)/2 where N is the filter order, and
%               is always odd. The ith element in P is |p_i|, the
%               magnitude of the ith pole, which represents a purely 
%               imaginary pole in the case of a halfband filter and a
%               purely real pole for a halfband-formulated hilbert 
%               transform filter. P is sorted from largest to smallest.
% P0............Array of magnitude-squared poles
%               P0 is the partition of P corresponding to the upper
%               branch allpass section (real path for Hilbert filter)
% P1............Array of magnitude-squared poles
%               P1 is the partition of P corresponding to the lower
%               branch allpass section (imag path for Hilbert filter)
%
% [1] D. Harris Et al., An Infinite Impulse Response (IIR) Hilbert Transformer
%     Filter Design Technique for Audio, AES Convention 2010.
%     https://www.academia.edu/73278010/An_Infinite_Impulse_Response_IIR_Hilbert_Transformer_Filter_Design_Technique_for_Audio
%
% [2] Robby Tong source code
% https://github.com/robwasab/HalfBand
%

    % (1) stopband edge frequency
    ws = pi - wp; % eq (18)

    % (2) stopband ripple
    ds = 10^(-As/20); % eq (19)

    % (3) intermediate calculations
    r = tan(wp/2) / tan(ws/2); % eq (20)
    rp = sqrt(1 - r^2); % eq (21)
    q0 = 0.5 * (1 - sqrt(rp)) / (1 + sqrt(rp)); % eq (22)
    q = q0 + 2*q0^5 + 15*q0^9 + 150*q0^13; % eq (23)
    D = ((1 - ds^2) / ds^2)^2; % eq (24)

    % (4) required filter order (must be odd)
    N = log10(16*D) / log10(1/q); % eq (25)
    N = ceil(N);
    if rem(N,2) == 0
        N = N + 1;
    end

    % (5) poles akp for k = 1..(N-1)/2
    P = nan((N-1)/2,1);
    for kk = 1:(N-1)/2
        numsum = 0;
        densum = 0;
        % The loop in [1] is formulated as two loops as per [2]
        for ii = 0:10-1
            num = (-1)^ii * q^(ii*(ii+1)) * sin(((2*ii+1)*kk*pi)/N);
            numsum = numsum + num;
        end
        for ii = 1:10-1
            den = (-1)^ii * q^(ii*ii) * cos(2*ii*kk*pi/N);
            densum = densum + den;
        end
        lambdak = (2*q^0.25) * numsum / (1 + 2*densum); % eq (26)
        %disp(lambdak);
        bk = sqrt((1 - r*lambdak^2) * (1 - lambdak^2/r)); % eq (27)
        ck = 2 * bk / (1 + lambdak^2); % eq (28)
        P(kk) = (2 - ck) / (2 + ck); % eq (29) pole locations
    end
    P = fliplr(sort(P(:).'));
    if rem((N-1)/2, 2) == 0
        P0 = P(2:2:end);
        P1 = P(1:2:end); % imag path
    else
        P0 = P(1:2:end);
        P1 = P(2:2:end); % imag path
    end

end % function
