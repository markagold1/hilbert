function info = hilbert_iir(varargin)
% Usage: info = hilbert_iir(varargin)
%
% Hilbert IIR design using method of [1].
%
% Command options:
%   info = hilbert_iir()  (demo mode)
%   info = hilbert_iir( ___ ,'As',As)
%   info = hilbert_iir( ___ ,'Fs',fs)
%   info = hilbert_iir( ___ ,'Ftype',ftype)
%   info = hilbert_iir( ___ ,'TBWN',tbwn)
%   info = hilbert_iir( ___ ,'Verbose',verbose)
%
% Inputs
%  As..........Stopband attenuation dB (default 60 dB)
%  fs..........Sampling frequeny in Hz (default 1 Hz)
%  ftype.......Filter type: 'lpf' or 'hilbert' (default 'hilbert')
%  verbose.....Verbosity: 0,1,or 2 (default 2)
%                0 - run quietly
%                1 - display configuration, plots, and filter analysis
%  tbwn........Normalized transition bandwidth, 0<tbwn<1 (default 0.1)
%              tbw (Hz) = tbwn * fnyq where fnyq = fs/2
%              Controls the passband edge frequency of the halfband
%              prototype filer. Smaller values sharpen the transition
%              at the cost of higher filter order and increased group
%              delay.
%
% Outputs
%  info........Struct containing IIR filter info
%
% Examples
%
%  1. Run the demo
%     y = hilbert_iir;
%     Filter Analysis:
%                                 Type: Hilbert
%                                Order: 11
%              Sampling frequency (Hz): 2
%       Stopband (SB) attenuation (dB): 67
%            Transition bandwidth (Hz): 0.1000
%        Lower passband (PB) edge (Hz): 0.0500
%                   Upper PB edge (Hz): 0.9500
%                       PB ripple (dB): +/-4.487e-07
%             PB group delay (samples): Peak 10.9, Mean 3.9
%              PB diff phase (degrees): Mean 90.000
%        PB diff phase error (degrees): +/-5.209e-02
%
%  2. Design a hilbert analytic transformer with 80 dB stopband attenuation,
%     transition bandwidth 360 Hz (0.015*fnyq), and sampling frequency 48 kHz.
%     y = hilbert_iir('As',80,'Fs',48e3,'TBWN',0.015,'Verbose',1);
%     Filter Analysis:
%                                 Type: Hilbert
%                                Order: 21
%              Sampling frequency (Hz): 48000
%       Stopband (SB) attenuation (dB): 82
%            Transition bandwidth (Hz): 360
%        Lower passband (PB) edge (Hz): 180
%                   Upper PB edge (Hz): 23820
%                       PB ripple (dB): +/-1.490e-08
%             PB group delay (samples): Peak 87.8, Mean 8.3
%              PB diff phase (degrees): Mean 90.000
%        PB diff phase error (degrees): +/-9.492e-03
%
%  3. Design a halfband filter with 50 dB stopband attenuation, transition
%     bandwidth 1102.5 (0.05*fnyq), and sampling frequency 44.1 kHz.
%     y = hilbert_iir('As',50,'Fs',44.1e3,'TBWN',0.05,'Ftype','lpf','Verbose',1);
%     Filter Analysis:
%                                 Type: Halfband
%                                Order: 11
%              Sampling frequency (Hz): 44100
%       Stopband (SB) attenuation (dB): 54
%            Transition bandwidth (Hz): 1103
%        Lower passband (PB) edge (Hz): -10474
%                   Upper PB edge (Hz): 10474
%                       PB ripple (dB): +/-8.710e-06
%             PB group delay (samples): Peak 18.0, Mean 4.1
%              PB diff phase (degrees): Mean -0.000
%        PB diff phase error (degrees): +/-2.295e-01
%
% References
%
% [1] D. Harris Et al., An Infinite Impulse Response (IIR) Hilbert Transformer
%     Filter Design Technique for Audio, AES Convention 2010.
%
% [2] Signal Processing Stack Exchange Answer
% https://dsp.stackexchange.com/questions/37411/iir-hilbert-transformer
%
% [3] Robby Tong source code
% https://github.com/robwasab/HalfBand
%

    global VERBOSE

    % Defaults
    info = struct();
    As = 60;
    fs = 2;
    ftype = 'hilbert';
    verbose = 2;
    tbwn_in = 0.1;
    if nargin == 0
        verbose = 1;
    end
    
    % Inputs
    kk = 1;
    while (kk < nargin)
        switch varargin{kk}
            case 'As'
                As = varargin{kk+1};
            case 'Fs'
                fs = varargin{kk+1};
            case 'Ftype'
                ftype = varargin{kk+1};
            case 'Verbose'
                verbose = varargin{kk+1};
            case 'TBWN'
                tbwn_in = varargin{kk+1};
        end
        kk = kk + 2;
    end

    % Parameters
    VERBOSE = verbose;
    fnyq = fs/2;
    tbw_norm = min(max(tbwn_in, 1e-6),1-eps); 
    tbw = tbw_norm*fnyq;
    fpass = (fnyq - tbw) / 2;
    wpass = (2*pi) * (fpass/fs);
    wstop = pi - wpass;

    % Display configuration
    if VERBOSE
        fprintf(1,'Configuration:\n');
        fprintf(1,'              Filter type: %s\n', ftype);
        fprintf(1,'           Stopband atten (dB): %d\n', As);
        fprintf(1,'       Sampling frequency (Hz): %d\n', fs);
        fprintf(1,'     Transition bandwidth (Hz): %.2f\n', tbw);
        fprintf(1,'\n');
        if VERBOSE > 1
          fprintf(1,'wpass=%.4f (%.3e Hz) wstop=%.4f (%.3e Hz)\n', ...
             wpass, wpass/2/pi*fs, wstop, wstop/2/pi*fs);
        end
    end

    % Halfband IIR pole locations
    [P,P0,P1] = halfband_poles(wpass,As);

    % Parallel allpass network using second order sections (SOS)
    [sos0,sos1] = form_allpass_network(P0,P1,ftype);

    % Filter response
    [H,H0,H1,gain,phase,gd,dph,F] = filter_response(sos0,sos1,fs,tbw,ftype);

    % Collect outputs
    info.ftype = ftype;
    info.stopband_atten_db = As;
    info.sampling_frequency_hz = fs;
    info.transition_bandwidth_hz = tbw;
    info.fpass_hz = fpass;
    info.filter_order = 2 * numel(P) + 1;
    info.filter_coefs = P;
    info.sos0 = sos0; % upper branch (real part of analytic signal)
    info.sos1 = sos1; % lower branch (imag part of analytic signal)
    info.H = H;
    info.H0 = H0;
    info.H1 = H1;
    info.gain = gain;
    info.phase = phase;
    info.Gd = gd;
    info.diff_phase_rad = dph;
    info.F = F(:);

    % Report measurements
    if VERBOSE
        report_measurements(info);
    end

end % function

function [H,H0,H1,gain,phase,gd,dph,F] = filter_response(sos0,sos1,fs,tbw,ftype);
    global VERBOSE

    F = fs * (-1000:1000)/2000;
    H0 = sosfreqz(sos0,F,fs);
    H1 = sosfreqz(sos1,F,fs);
    if strcmpi(ftype,'hilbert')
        H = (H0 + j*H1) / 2;
    else
        H = (H0 + H1) / 2;
    end
    gain = 20*log10(abs(H));
    phase = pm_unwrap(angle(H));
    dw = 2*pi*diff(F(:));
    gd = -diff(phase) ./ dw * fs;
    gd(end+1) = gd(end); % preserve length
    gd = max(0,gd);
    ph0 = pm_unwrap(angle(H0));
    ph1 = pm_unwrap(angle(H1));
    dph = birem((ph0-ph1)/pi,2) * pi; % rad

    % Plot
    if VERBOSE
        if strcmpi(ftype,'hilbert')
            fc = fs/4;
        else
            fc = 0;
        end
        tbw2 = tbw / 2;
        fpb = F > fc - fs/4 + tbw2 & F < fc + fs/4 - tbw2;
        if strcmpi(ftype,'hilbert')
            fsb = F < fc - fs/4 - tbw2 & F > fc - 3*fs/4 + tbw2;
            ftype = 'Hilbert';
        else
            fsb = F < fc - fs/4 - tbw2 | F > fc + fs/4 + tbw2;
            ftype = 'Halfband';
        end

        subplot(411);
            plot(F,gain,F(fpb),gain(fpb));dz;
            axis([-0.5*fs 0.5*fs -150 10]);
            ylabel('|H| (dB)');
            title('Filter Response');
        subplot(412);
            plot(F,phase/pi,F(fpb),phase(fpb)/pi);dz;
            ax = axis; axis([-0.5*fs 0.5*fs ax(3:4)]);
            ylabel('Phase/\pi');
        subplot(413);
            plot(F,gd,F(fpb),gd(fpb));dz;
            ax = axis; axis([-0.5*fs 0.5*fs ax(3:4)]);
            ylabel('Delay (samples)');
        subplot(414);
            dph01 = birem((ph0-ph1)/pi,2);
            plot(F,dph01,F(fpb),dph01(fpb));dz;
            ax = axis; axis([-0.5*fs 0.5*fs ax(3:4)]);
            ylabel('Diff Phase/\pi');
            xlabel('Frequency (Hz)');
            %plot(F,(ph0-ph1)/pi);dz;
    end

end % function

function [sos0,sos1] = form_allpass_network(P0,P1,ftype)

    % convert P0 to SOS
    sos0 = nan(numel(P0),6);
    for kk = 1:numel(P0)
        ak0 = P0(kk);
        sos0(kk,:) = [ak0 0 1 1 0 ak0];
    end
    if strcmpi(ftype,'hilbert')
        sos0(:,3) = -sos0(:,3);
        sos0(:,6) = -sos0(:,6);
    end

    % convert P1 to SOS
    sos1 = nan(numel(P1),6);
    for kk = 1:numel(P1)
        ak1 = P1(kk);
        sos1(kk,:) = [ak1 0 1 1 0 ak1];
    end
    if strcmpi(ftype,'hilbert')
        sos1(:,3) = -sos1(:,3);
        sos1(:,6) = -sos1(:,6);
    end
    sos1 = [[0 1 0   1 0 0]; sos1];

end % function

function r = birem(a,w)
% modulo wrap input 'a' around +/- 1/2 wrap value 'w'
% range of output 'r' is [-W/2,+W/2)

    r = rem(a,w);
    r(r >= w/2) = r(r >= w/2) - w;
    r(r < -w/2) = r(r < -w/2) + w;

end % function

function dz(cols)
    if exist('dark') == 2
        dark('tr');
        grid on;
    end
    figure(gcf);
    grid on;
    zoom on;
end % function

function report_measurements(info)

    ftype = info.ftype;
    N = info.filter_order;
    fs = info.sampling_frequency_hz;
    fpass = info.fpass_hz;
    fstop = fs/2 - info.fpass_hz;
    F = info.F;
    tbw2 = info.transition_bandwidth_hz / 2;
    dph = info.diff_phase_rad*180/pi;
    gn = info.gain;
    if strcmpi(ftype,'hilbert')
        fc = fs/4;
    else
        fc = 0;
    end
    fpb = F > fc - fs/4 + tbw2 & F < fc + fs/4 - tbw2;
    if strcmpi(ftype,'hilbert')
        fsb = F < fc - fs/4 - tbw2 & F > fc - 3*fs/4 + tbw2;
        ftype = 'Hilbert';
    else
        fsb = F < fc - fs/4 - tbw2 | F > fc + fs/4 + tbw2;
        ftype = 'Halfband';
    end
    sb_atten = abs(max(gn(fsb)));
    gd_pk = max(info.Gd(fpb));
    gd_mu = mean(info.Gd(fpb));
    ph_pk2pk = max(dph(fpb)) - min(dph(fpb));
    gn_pk2pk = max(gn(fpb)) - min(gn(fpb));
    dph_mu = mean(dph(fpb));

    fprintf(1,'\nFilter Analysis:\n');
    fprintf(1,'                          Type: %s\n', ftype);
    fprintf(1,'                         Order: %d\n', N);
    fprintf(1,'       Sampling frequency (Hz): %d\n', fs);
    fprintf(1,'Stopband (SB) attenuation (dB): %d\n', round(sb_atten));
    if tbw2 > 50
        fprintf(1,'     Transition bandwidth (Hz): %d\n', round(2*tbw2));
        fprintf(1,' Lower passband (PB) edge (Hz): %d\n', round(fc - fpass));
        fprintf(1,'            Upper PB edge (Hz): %d\n', round(fc + fpass));
    else
        fprintf(1,'     Transition bandwidth (Hz): %.4f\n', 2*tbw2);
        fprintf(1,' Lower passband (PB) edge (Hz): %.4f\n', fc - fpass);
        fprintf(1,'            Upper PB edge (Hz): %.4f\n', fc + fpass);
    end
    if gn_pk2pk/2 < 1e-3
        fprintf(1,'                     PB ripple: less than +/-0.001 dB\n');
    else
        fprintf(1,'                PB ripple (dB): +/-%.3e\n', gn_pk2pk/2);
    end
    fprintf(1,'      PB group delay (samples): Peak %.1f, Mean %.1f\n', ...
              gd_pk, gd_mu);
    fprintf(1,'       PB diff phase (degrees): Mean %.3f\n', dph_mu);
    fprintf(1,' PB diff phase error (degrees): +/-%.3e\n', ph_pk2pk/2);

end % function
