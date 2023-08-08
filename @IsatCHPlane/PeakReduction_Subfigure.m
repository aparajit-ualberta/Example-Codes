function [fig1,ax,ax1,ax2] = PeakReduction_Subfigure(obj,rpos,varargin)
%     warning off
    Inpt = inputParser;
    
    addRequired(Inpt, 'Time_Series', @(obj) isobject(obj))
    addRequired(Inpt, 'Radial_Posn', @(rpos) isnumeric(rpos))
    addParameter(Inpt, 'Embedding_Dimension',5, @(d) isnumeric(d))
    addParameter(Inpt, 'Subsample_Dimension',8, @(s) isnumeric(s))
    addParameter(Inpt, 'Start_Time',5, @(tstart) isnumeric(tstart))
    addParameter(Inpt, 'Stop_Time',15, @(tstop) isnumeric(tstop))
    addParameter(Inpt, 'Parameters',[], @(params) isvector(params) && length(params) == 3)
    addParameter(Inpt, 'Save','none', @(sv) ischar(sv))
    addParameter(Inpt, 'Format','none', @(fm) ismember(fm, {'fig','png','eps','jpg','none'}))

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
    if isempty(Inpt.Results.Parameters)
        TimeLabel = "Inpt"+Inpt.Results.Start_Time+"to"+Inpt.Results.Stop_Time+"s";
    end
    
    cf = (2*1.602e-19*0.1)/((2*1.673e-27)+(2*1.675e-27));
    %% Figures
    %FFT/CH-Plane subplot
    fig1 = figure(1);
    fig1.Units = 'centimeters';
    fig1.Position = [12.8852 0 18.5 21];
    fig1.Color = 'white';
    
    %FFT
    ax = subplot(2,1,2);
    ax.NextPlot = 'add';
    ax.Units = 'normalized';
    ax.Position = [0.1 0.08 0.85 0.3419];
    ax.FontSize = 14;
    ax.XLabel.Interpreter = 'latex';
    ax.XLabel.String = 'Frequency, $f$ (kHz)';  
    ax.YLabel.Interpreter = 'latex';
    ax.YLabel.String = '$\log(\mathcal{F}[\delta I_{\mathrm{sat}}])$';
    ax.TickLabelInterpreter = 'latex';
    ax.Layer = 'top';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.TickLength = [0.03 0.035];
    ax.XLim = [20 600];
    ax.Box = 'on';
    
%     axx = axes('XAxisLocation','top','YAxisLocation','right','Color','none');
%     axx.NextPlot = 'add';
%     axx.Units = 'normalized';
%     axx.Position = ax.Position;
%     axx.FontSize = 14;
%     axx.XLabel.Interpreter = 'latex';
%     axx.XLabel.String = '$\omega / {\Omega}_i$';
%     axx.TickLabelInterpreter = 'latex';
%     axx.Layer = 'top';
%     axx.XMinorTick = 'on';
%     axx.YMinorTick = 'on';
%     axx.TickLength = [0.03 0.035];
%     axx.YTick = [];
%     axx.XLim = [0 (2*pi*200/cf)];
%     ax.Box = 'on';
    
    %CH Plane
    ax1 = subplot(2,1,1);
    ax1.NextPlot = 'add';
    ax1.Units = 'normalized';
    ax1.Position = [0.1 0.5161 0.85 0.4516];
    ax1.FontSize = 14;
    ax1.XLabel.Interpreter = 'latex';
    ax1.XLabel.String = '$H$';
    ax1.YLabel.Interpreter = 'latex';
    ax1.YLabel.String = '$C$';
    ax1.TickLabelInterpreter = 'latex';
    ax1.Layer = 'top';
    ax1.XMinorTick = 'on';
    ax1.YMinorTick = 'on';
    ax1.TickLength = [0.03 0.035];
    ax1.Title.FontSize = 14;
    ax1.Title.Interpreter = 'latex';
    ax1.YLim = [0 0.5];
    ax1.XLim = [0 1];
    ax1.Box = 'on';
    
    %Time Series plot
    ax2 = axes;
    ax2.NextPlot = 'add';
    ax2.Units = 'normalized';
    ax2.Position = [0.5 0.29 0.4 0.1];
    ax2.FontSize = 10;
    ax2.XLabel.Interpreter = 'latex';
    ax2.XLabel.String = 'Time, $t$ (ms)';
    ax2.YLabel.Interpreter = 'latex';
    ax2.YLabel.String = '$\delta I_{\mathrm{sat}}$';
    ax2.TickLabelInterpreter = 'latex';
    ax2.Layer = 'top';
    ax2.XMinorTick = 'on';
    ax2.YMinorTick = 'on';
    ax2.TickLength = [0.03 0.035];
    ax2.XLim = [Inpt.Results.Start_Time,Inpt.Results.Stop_Time];
    ax2.Box = 'on';
    %% FFT

    fftlen = length(obj.data.t(tstart:tstop));
    FFTmov = zeros(10,fftlen); % Declare FFT for position
    FFT0   = zeros(10,fftlen);  % Declare FFT of background (edge)
    Fs = 1/3.2e-7;
    f = (Fs*(0:(fftlen-1))/fftlen)./1000;
    wcf = (2.*pi.*f)./(cf);

    for ishot = 1:1:10
        FFTmov(ishot,:) = fft(obj.data.Isat_diff(ishot,tstart:tstop,rpos)); % Store Complex FFT for each shot  
        FFT0(ishot,:) = fft(obj.data.Isat_diff(ishot,tstart:tstop,28)); % Store Complex FFT for last shot  
    end

    F0   = log(mean(abs(FFT0/fftlen)));
    FFTabs   = log(mean(abs(FFTmov/fftlen))) - F0;
    %% Smooth FFT
    slen = 1;
    FFT_smooth = smoothdata(FFTabs,'movmean',slen);
    %% Fitting

    if isempty(Inpt.Results.Parameters)
        fexp1 = I.(TimeLabel)(3,rpos); %0.9 for 5-15s, 5.4s for rpos=1 (10-15s),
        fexp2 = I.(TimeLabel)(1,rpos);
        yint1 = I.(TimeLabel)(2,rpos); %Y-intercept
    else
        fexp1 = Inpt.Results.Parameters(1);
        fexp2 = Inpt.Results.Parameters(2);
        yint1 = Inpt.Results.Parameters(3);
    end

%     f1=find(abs(f-fexp1) < 0.001,1);
    f1=find(f >= fexp1,1);
%     f1 = 3;
    f2=find(f >= fexp2,1);

    ft1=f(f1);
    ft2=f(f2);

    FT1=FFTabs(f1);
    FT2=FFTabs(f2);

    FFT_Lorentz = zeros(1,fftlen);
    FFT_fit = polyfit([ft1 ft2],[FT1 FT2],1);
    %Temp code begin
    FFT_Lorentz(1:f1) = polyval(FFT_fit,f(1:f1))+yint1;
    %Temp code end
    FFT_Lorentz(f1:f2) = polyval(FFT_fit, f(f1:f2))+ yint1;
    FFT_Lorentz(f2+1:end-(f2+1)+2) = FFT_Lorentz(f2);
    FFT_fit2 = polyfit([f(end-f2+2) f(end-f1+2)],[FFTabs(end-f2+2) FFTabs(end-f1+2)],1);
    FFT_Lorentz(end-f2+2:end-f1+2) = polyval(FFT_fit2, f(end-f2+2:end-f1+2))+ yint1;

    
    %% Peak Reduction

    sf = ones(1,fftlen);
    for q = f1:1:f2
        if FFT_smooth(q)<FFT_Lorentz(q)
            sf(q) = 1;
        else
            sf(q) = (exp((FFT_smooth(q))))./(exp(FFT_Lorentz(q)));
        end
    end

    for q = fftlen-f2+2:1:fftlen-f1+2
        if FFT_smooth(q)<FFT_Lorentz(q)
            sf(q) = 1;
        else
            sf(q) = (exp((FFT_smooth(q))))./(exp(FFT_Lorentz(q)));
        end
    end

    % pow = 0:0.1:1;
    u = 0;
    for n = 1:1:11
        u = u+1;
        FFTPR   = (FFTmov)./(sf.^(n*0.1 - 0.1)); % Power gives the amplitude of reduced fft
        FFTnew   = log(mean(abs(FFTPR/fftlen))) - log(mean(abs(FFT0/fftlen)));
        FFTnew_smooth = smoothdata(FFTnew,'movmean',1);

        for p = 1:1:10
            sig_PR(p,:) = ifft((FFTPR(p,:)));
        end
        sig_unfil = mean(sig_PR);

        figure(1)
        hold on
        plot(ax,f,FFTnew_smooth,'LineWidth',2)
        hold off

        k = 0;
        s = Inpt.Results.Subsample_Dimension;
        for m = 1:s:length(sig_unfil)
            k = k+1;
            temp(k) = real(sig_unfil(m));
            t_temp(k) = obj.data.t(m);
        end
        %CH Plane plot (Markers)
        figure(1)
        hold on
        [H(u), C(u)] = EntropyComplexity(temp,d);
        ln1(n) = plot(ax1,H(u),C(u),'o','LineWidth',2);
        hold off
        % Annotations (Textbox)
        an1 = annotation('textbox',[0.25 0.8 0.1 0.1]);
        an1.Interpreter = 'latex';
        an1.String = '$d = 5$';
        an1.FontSize = 20;
        an1.LineStyle = 'none';
        
        % Time Series Plot
        hold on
        ln2(n) = plot(ax2,obj.data.t(tstart:tstop),real(sig_unfil));
        if sum(n == [1,6,11]) == 0
            ln2(n).Visible = 'off';
        end
        hold off
    end
    
    % FFT Plot fit 
    figure(1)
    hold on
    plot(ax,f,(FFT_Lorentz),'r--','LineWidth',3)
%     plot(axx,wcf,FFTnew_smooth,'Visible','off')
    hold off
   
    % CH Plane plot (Markers) Complexity Curves
    figure(1)
    hold on
    plot(ax1,H_min, C_min,'k-','LineWidth',2);
    plot(ax1,H_max, C_max,'k-','Linewidth',2);
    plot(ax1,fBm_CH.H,fBm_CH.C,'k--','LineWidth',2);
    hold off
    
    an2 = annotation('textbox',[0.58 0.68 0.1 0.1]);
    an2.Interpreter = 'latex';
    an2.String = 'fBm';
    an2.FontSize = 16;
    an2.LineStyle = 'none';
    
    an3 = annotation('arrow',[0.64 0.75],[0.86 0.82],'LineWidth',2);
    
    
    lgd = legend(ax1,ln1, {'0\%','10\%','20\%','30\%','40\%','50\%','60\%','70\%','80\%','90\%','100\%'});
    lgd.Location = 'northwest';
    lgd.Box = 'off';
    lgd.Interpreter = 'latex';
    %% Saving
    if strcmp(Inpt.Results.Save,'none') == 0
        if strcmp(Inpt.Results.Format,'fig')
            savefig(Inpt.Results.Save)
        elseif strcmp(Inpt.Results.Format,'png') || strcmp(Inpt.Results.Format,'eps') || strcmp(Inpt.Results.Format,'jpg')
            export_fig(Inpt.Results.Save,'-'+Inpt.Results.Format,'-native')
        elseif endsWith(Inpt.Results.Save,'.png') || endsWith(Inpt.Results.Save,'.eps') || endsWith(Inpt.Results.Save,'.jpg')
            export_fig(Inpt.Results.Save,'-native')
        elseif endsWith(Inpt.Results.Save,'.fig')
            savefig(Inpt.Results.Save)
        else
            error('File Type not specified')
        end
    end
end

