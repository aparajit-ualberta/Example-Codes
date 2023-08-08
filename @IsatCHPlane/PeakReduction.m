function PeakReduction(obj,rpos,varargin)
    Inpt = inputParser;
    
    addRequired(Inpt, 'Time_Series', @(obj) isobject(obj))
    addRequired(Inpt, 'Radial_Posn', @(rpos) isnumeric(rpos))
    addParameter(Inpt, 'Embedding_Dimension',5, @(d) isnumeric(d))
    addParameter(Inpt, 'Subsample_Dimension',8, @(s) isnumeric(s))
    addParameter(Inpt, 'Start_Time',5, @(tstart) isnumeric(tstart))
    addParameter(Inpt, 'Stop_Time',15, @(tstop) isnumeric(tstop))
    addParameter(Inpt, 'save', 'none', @(s) ismember(s,{'fig','jpg'}) )

    parse(Inpt, obj, rpos, varargin{:})
    
    close all
    
    %Fitting Parameters
    I = struct('Inpt0to5s', [138.7 158.3 156.7 157.9 248.3 216.3 175.7 117.1 79.9 81.1; +0.6 -0.2 -0.1 +0.1 +0.1 -0.2 -0.1 +0.0 +0.1 -0.0; 6.3997 6.1996 5.5996 5.5996 5.7996 12.3992 10.3993 5.9996 8.1995 12.7992],...
            'Inpt5to10s', [151.1 158.3 157.5 177.7 243.5 216.9 177.1 118.7 69.7 79.9; +0.2 -0.1 -0.0 -0.0 -0.6 -0.1 -0.0 +0.0 +0.0 -0.0; 3.3998 4.1997 6.1996 6.1996 5.7996 5.3997 8.3995 6.7996 6.7996 7.1995],...
            'Inpt10to15s', [99 140.8 131.1 138.4 166.5 138.5 118.5 101.5 101.3 99.7; -0.2 -0.5 -0.1 -0.4 -0.1 -0.0 -0.0 -0.0 -0.0 -0.0; 5.4 3.5998 5.9996 3.9997 5.9996 4.9997 15.599 4.9997 4.3997 4.5997],...
            'Inpt5to15s', [51.3 151.3 151.5 151.5 248.3 248.3 180.2 122 79.9 59.8; -0.2 -0.1 -0.1 +0.1 +0.1 -0.2 -0.1 -0.1 +0.5 -0.3; 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9 0.9]);
    %CH Plane
    d = Inpt.Results.Embedding_Dimension; %Subsample Dimension
    %  Calculating the min and max complexity curves to imptove memory
    [H_min, C_min] = MinComplexCurve(d); %Max-complexity curve
    [H_max, C_max] = MaxComplexCurve(d); %Min-complexity curve

    %Data
    load('fBm_CH_Curve','fBm_CH')
    tstart = find(obj.data.t >= Inpt.Results.Start_Time,1); % FFT Time Start
    tstop  = find(obj.data.t >= Inpt.Results.Stop_Time,1); % FFT Time Stop
    TimeLabel = "Inpt"+Inpt.Results.Start_Time+"to"+Inpt.Results.Stop_Time+"s";
    %% FFT

    fftlen = length(obj.data.t(tstart:tstop));
    FFTmov = zeros(10,fftlen); % Declare FFT for position
    FFT0   = zeros(10,fftlen);  % Declare FFT of background (edge)
    Fs = 1/3.2e-7;
    f = (Fs*(0:(fftlen-1))/fftlen)./1000;

    for ishot = 1:1:10
        FFTmov(ishot,:) = fft(obj.data.Isat_diff(ishot,tstart:tstop,rpos)); % Store Complex FFT for each shot  
        FFT0(ishot,:) = fft(obj.data.Isat_diff(ishot,tstart:tstop,28)); % Store Complex FFT for last shot  
    end

    F0   = log(mean(abs(FFT0/fftlen)));
    FFTabs   = log(mean(abs(FFTmov/fftlen))) - F0;
    %% Smooth FFT
    slen = 20;
    FFT_smooth = smoothdata(FFTabs,'movmean',slen);
    %% Fitting

    fexp1 = I.(TimeLabel)(3,rpos); %0.9 for 5-15s, 5.4s for rpos=1 (10-15s),
    fexp2 = I.(TimeLabel)(1,rpos);
    yint1 = I.(TimeLabel)(2,rpos); %Y-intercept

    f1=find(abs(f-fexp1) < 0.001,1);
    f2=find(f >= fexp2,1);

    ft1=f(f1);
    ft2=f(f2);

    FT1=FFTabs(f1);
    FT2=FFTabs(f2);

    FFT_Lorentz = zeros(1,fftlen);
    FFT_fit = polyfit([ft1 ft2],[FT1 FT2],1);
    FFT_Lorentz(f1:f2) = polyval(FFT_fit, f(f1:f2))+ yint1;
    FFT_Lorentz(f2+1:end-(f2+1)+2) = FFT_Lorentz(f2);
    FFT_fit2 = polyfit([f(end-f2+2) f(end-f1+2)],[FFTabs(end-f2+2) FFTabs(end-f1+2)],1);
    FFT_Lorentz(end-f2+2:end-f1+2) = polyval(FFT_fit2, f(end-f2+2:end-f1+2))+ yint1;

    %FFT plot
    fig1 = figure(1);
    fig1.WindowState = 'maximized';
    hold on
    plot(f,(FFTabs))
    plot(f,(FFT_smooth),'k','LineWidth',2)
    hold off
    ax = gca;
    ax.Position = [0.05 0.1 0.9 0.85];
    ax.FontSize = 14;
    ax.XLabel.String = '$Frequency (kHz)$';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.String = '$Power$';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    ax.Layer = 'top';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.TickLength = [0.03 0.035];
    ax.Title.String = "Peak Reduced FFT, Port 30, Subsampled "+obj.data.r(rpos)+" cm, "+obj.data.t(tstart)+"-"+obj.data.t(tstop)+"s";
    ax.Title.FontSize = 14;
    ax.Title.Interpreter = 'latex';
    % ax.YLim = [0 0.2];
    ax.XLim = [0 200];
    ax.Box = 'on';
    %     lgd = legend([a(1) a(2) a(3)], {'Edges','Gradient','Center'});
    %     lgd.Box = 'off';
    %% Peak Reduction

    sf = ones(1,fftlen);
    sf(f1:f2) = (exp((FFT_smooth(f1:f2))))./(exp(FFT_Lorentz(f1:f2)));
    sf(end-f2+2:end-f1+2) = exp((FFT_smooth(end-f2+2:end-f1+2)))./(exp(FFT_Lorentz(end-f2+2:end-f1+2)));

    % pow = 0:0.1:1;
    u = 0;
    for n = 1:1:11
        u = u+1;
        FFTPR   = (FFTmov)./(sf.^(n*0.1 - 0.1)); % Power gives the amplitude of reduced fft
        FFTnew   = log(mean(abs(FFTPR/fftlen))) - log(mean(abs(FFT0/fftlen)));
        FFTnew_smooth = smoothdata(FFTnew,'movmean',20);

    %     sig_Old = obj.P.Isat_diff(7,tstart:tstop,rpos);
    %     sig_Old_fft=abs(fft(sig_Old));
    %     tFFT    = obj.data.t(tstart:tstop);
        for p = 1:1:10
            sig_PR(p,:) = ifft((FFTPR(p,:)));
        end
        sig_unfil = mean(sig_PR);

        figure(1);
        hold on
        plot(f,FFTnew_smooth)
        hold off

        k = 0;
        s = Inpt.Results.Subsample_Dimension;
        for m = 1:s:length(sig_unfil)
            k = k+1;
            temp(k) = real(sig_unfil(m));
            t_temp(k) = obj.data.t(m);
        end
        %CH Plane plot (Markers)
        fig2 = figure(2);
    %     fig2.WindowState = 'maximized';
        hold on
        [H(u), C(u)] = EntropyComplexity(temp,d);
        plot(H(u),C(u),'o');
        hold off
        % Annotations (Textbox)
        an1 = annotation('textbox',[0.2 0.7 0.1 0.1]);
        an1.String = 'd = 5';
        an1.Interpreter = 'latex';
        an1.FontSize = 14;
        % Axis
        ax1 = gca;
    %     ax1.Position = [0.05 0.1 0.9 0.85];
        ax1.FontSize = 14;
        ax1.XLabel.String = '$H$';
        ax1.XLabel.Interpreter = 'latex';
        ax1.YLabel.String = '$C$';
        ax1.YLabel.Interpreter = 'latex';
        ax1.TickLabelInterpreter = 'latex';
        ax1.Layer = 'top';
        ax1.XMinorTick = 'on';
        ax1.YMinorTick = 'on';
        ax1.TickLength = [0.03 0.035];
        ax1.Title.String = "Peak Reduced CH Plane, Port 30, Subsampled "+obj.data.r(rpos)+" cm";
        ax1.Title.FontSize = 14;
        ax1.Title.Interpreter = 'latex';
        ax1.YLim = [0 0.5];
        ax1.XLim = [0 1];
        ax1.Box = 'on';
    %     lgd = legend([a(1) a(2) a(3)], {'Edges','Gradient','Center'});
    %     lgd.Box = 'off';

        % Time Series Plot
        fig3 = figure(3);
        fig3.WindowState = 'maximized';
        hold on
        plot(obj.data.t(tstart:tstop),real(sig_unfil));
        hold off
        ax2 = gca;
        ax2.Position = [0.05 0.1 0.9 0.85];
        ax2.FontSize = 14;
        ax2.XLabel.String = '$Time (s)$';
        ax2.XLabel.Interpreter = 'latex';
        ax2.YLabel.String = '$Isat$';
        ax2.YLabel.Interpreter = 'latex';
        ax2.TickLabelInterpreter = 'latex';
        ax2.Layer = 'top';
        ax2.XMinorTick = 'on';
        ax2.YMinorTick = 'on';
        ax2.TickLength = [0.03 0.035];
        ax2.Title.String = "Peak Reduced Time Series, Port 30, Subsampled "+obj.data.r(rpos)+" cm";
        ax2.Title.FontSize = 14;
        ax2.Title.Interpreter = 'latex';
        % ax.YLim = [0 0.2];
    %     ax1.XLim = [0 10];
        ax2.Box = 'on';
    end
    % FFT Plot fit 
    figure(1)
    hold on
    plot(f,(FFT_Lorentz),'r--','LineWidth',3)
    hold off

    % CH Plane plot (Markers) Complexity Curves
    figure(2)
    hold on
    plot(H_min, C_min,'k-');
    plot(H_max, C_max,'k-');
    hold off

    %Quiver plot set-up
    P = [H' C'];
    for q = 1:1:length(H)-1
        D(q,:) = P(q+1,:) - P(q,:);
    end

    %CH Plane plot (Arrows)
    fig4 = figure(4);
    % fig4.WindowState = 'maximized';
    for q = 1:1:length(D)
        hold on
        quiver(P(q,1),P(q,2),D(q,1),D(q,2),'r','MaxHeadSize',5)
    end
    plot(H,C,'r*')
    plot(fBm_CH.H,fBm_CH.C,'k--');
    plot(H_min, C_min,'k-','LineWidth',1.5);
    plot(H_max, C_max,'k-','LineWidth',1.5);
    hold off
    % Annotations (Textbox)
    an2 = annotation('textbox',[0.2 0.7 0.1 0.1]);
    an2.String = 'd = 5';
    an2.Interpreter = 'latex';
    an2.FontSize = 20;
    an2.LineStyle = 'none';
    % Axis
    ax3 = gca;
    % ax3.Position = [0.05 0.1 0.9 0.85];
    ax3.FontSize = 14;
    ax3.XLabel.String = '$H$';
    ax3.XLabel.Interpreter = 'latex';
    ax3.YLabel.String = '$C$';
    ax3.YLabel.Interpreter = 'latex';
    ax3.TickLabelInterpreter = 'latex';
    ax3.Layer = 'top';
    ax3.XMinorTick = 'on';
    ax3.YMinorTick = 'on';
    ax3.TickLength = [0.03 0.035];
    ax3.Title.String = "Peak Reduced CH Plane, Port 30, Subsampled "+obj.data.r(rpos)+" cm, "+obj.data.t(tstart)+"-"+obj.data.t(tstop)+"s";
    ax3.Title.FontSize = 14;
    ax3.Title.Interpreter = 'latex';
    ax3.YLim = [0 0.5];
    ax3.XLim = [0 1];
    ax3.Box = 'on';
    
if strcmp(Inpt.Results.save,'fig')
    fl1 = "Figure_PeakRed_"+obj.data.t(tstart)+"to"+obj.data.t(tstop)+"s_"+obj.data.r(rpos)+"cm_FFT."+Inpt.Results.save;
    saveas(fig1,fl1,'fig')
    fl3 = "Figure_PeakRed_"+obj.data.t(tstart)+"to"+obj.data.t(tstop)+"s_"+obj.data.r(rpos)+"cm_TimeSeries."+Inpt.Results.save;
    saveas(fig3,fl3,'fig')
    fl4 = "Figure_PeakRed_"+obj.data.t(tstart)+"to"+obj.data.t(tstop)+"s_"+obj.data.r(rpos)+"cm_CHPlane_Arrows."+Inpt.Results.save;
    saveas(fig4,fl4,'fig')
elseif strcmp(Inpt.Results.save,'jpg')
    fl1 = "Figure_PeakRed_"+obj.data.t(tstart)+"to"+obj.data.t(tstop)+"s_"+obj.data.r(rpos)+"cm_FFT."+Inpt.Results.save;
    export_fig(fig1,fl1,'-jpg')
    fl3 = "Figure_PeakRed_"+obj.data.t(tstart)+"to"+obj.data.t(tstop)+"s_"+obj.data.r(rpos)+"cm_TimeSeries."+Inpt.Results.save;
    export_fig(fig3,fl3,'-jpg')
    fl4 = "Figure_PeakRed_"+obj.data.t(tstart)+"to"+obj.data.t(tstop)+"s_"+obj.data.r(rpos)+"cm_CHPlane_Arrows."+Inpt.Results.save;
    export_fig(fig4,fl4,'-jpg')
end

    
end

