function [PF] = PosnFreq_Shot(obj,r,st,varargin)
    warning off
    Inpt = inputParser;

    addRequired(Inpt, 'Time_Series', @(obj) isobject(obj))
    addParameter(Inpt, 'Subsample_Dimension',1, @(s) isnumeric(s))
    addParameter(Inpt, 'Start_Time',5, @(tstart) isnumeric(tstart))
    addParameter(Inpt, 'Stop_Time',15, @(tstop) isnumeric(tstop))

    parse(Inpt, obj, varargin{:})
    
    t1 = find(obj.data.t >= Inpt.Results.Start_Time,1);
    t2 = find(obj.data.t >= Inpt.Results.Stop_Time,1);
    
    s = Inpt.Results.Subsample_Dimension;
    
    k = 0;
    for m = t1:s:t2
        k = k+1; %Create a new index for the subsampled signal
        Isat_diff_temp(k) = obj.data.Isat_diff(st,m,r);
        t_temp(k) = obj.data.t(m); %Create a subsampled time series
    end
    
    Fs = 1/(t_temp(2)-t_temp(1));
    L = length(t_temp);
    V_freq = fft(Isat_diff_temp);
    P1 = abs(V_freq/L);
    f2 = Fs*(0:1/L:1/2);
    P2 = fftshift(P1);
    PF = P2(L/2:end);
%     periodogram(Isat_diff_temp,rectwin(length(Isat_diff_temp)),length(Isat_diff_temp),Fs);
%     pspectrum(Isat_diff_temp,Fs)
    plot(f2,10.*log10(PF));
    ax1 = gca;
    ax1.FontSize = 12;
    ax1.XLabel.String = '$f(Hz)$';
    ax1.XLabel.Interpreter = 'latex';
    ax1.YLabel.String = '$Power (dB)$';
    ax1.YLabel.Interpreter = 'latex';
    ax1.TickLabelInterpreter = 'latex';
    ax1.Layer = 'top';
    ax1.XMinorTick = 'on';
    ax1.YMinorTick = 'on';
    ax1.TickLength = [0.03 0.035];
    ax1.Title.String = "Saturation Current FFT, R = "+((r*0.1)-0.1)+" cm, Shot "+st;
    ax1.Title.FontSize = 14;
    ax1.Title.Interpreter = 'latex';
    % ax1.YLim = [0 1];
%     ax1.XLim = [0 1];
    ax1.Box = 'on';
