function [H,C,CHFig] = CHPlane_Posn(obj,rad,varargin)
Inpt = inputParser;

addRequired(Inpt, 'Time_Series', @(obj) isobject(obj))
addRequired(Inpt, 'Radius', @(rad) isnumeric(rad) && rad <= 2.7)
addParameter(Inpt, 'Embedding_Dimension',5, @(d) isnumeric(d))
addParameter(Inpt, 'Subsample_Dimension',8, @(s) isnumeric(s))
addParameter(Inpt, 'Start_Time',5, @(tstart) isnumeric(tstart))
addParameter(Inpt, 'Stop_Time',15, @(tstop) isnumeric(tstop))
addParameter(Inpt, 'Shots',1:1:10, @(s) isvector(s) && sum(s>=1)==length(s) && sum(s<=10)==length(s) )
addParameter(Inpt, 'Additional','none', @(fm) ismember(fm, {'none','Mean'}))
addParameter(Inpt, 'Save','none', @(sv) ischar(sv))
addParameter(Inpt, 'Format','none', @(fm) ismember(fm, {'fig','png','eps','jpg','none'}))


parse(Inpt, obj, rad, varargin{:})

load('fBm_CH_Curve','fBm_CH')
t1 = find(obj.data.t >= Inpt.Results.Start_Time,1);
t2 = find(obj.data.t >= Inpt.Results.Stop_Time,1);
% t3 = length(obj.data.t);
r = find(obj.data.r >= Inpt.Results.Radius,1);
d = Inpt.Results.Embedding_Dimension;
s = Inpt.Results.Subsample_Dimension;
st = Inpt.Results.Shots;
lst = length(st);
[H_min, C_min] = MinComplexCurve(d); %Minimum Complexity Curve
[H_max, C_max] = MaxComplexCurve(d); %Maximum Complexity Curve

if lst == 1
    k = 0;
    for m = t1:s:t2
        k = k+1; %Create a new index for the subsampled signal
        Isat_diff_temp(k) = obj.data.Isat_diff(st,m,r);
        t_temp(k) = obj.data.t(m); %Create a subsampled time series
    end
    
    [H,C] = EntropyComplexity(Isat_diff_temp,d);
    
    CHFig = figure(1);
    hold on
    plot(H,C,'ro')
    plot(H_min, C_min,'k-','LineWidth',1.5); %Min complexity curve
    plot(H_max, C_max,'k-','LineWidth',1.5); %Max complexity curve
    plot(fBm_CH.H,fBm_CH.C,'k--'); %Brownian Motion Curve
    hold off
    ax = gca;
    ax.FontSize = 14;
    ax.XLabel.String = '$H$';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.String = '$C$';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    ax.Layer = 'top';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.TickLength = [0.03 0.035];
    ax.Title.String = "Port 30, Subsampled, "+obj.data.t(t1)+"-"+obj.data.t(t2)+" s, "+obj.data.r(r)+"cm";
    ax.Title.FontSize = 14;
    ax.Title.Interpreter = 'latex';
    % ax.YLim = [0 0.2];
    ax.XLim = [0 1.3];
    ax.Box = 'on';
    % Annotations (Textbox)
    an = annotation('textbox',[0.2 0.7 0.1 0.1]);
    an.String = "d = "+d;
    an.Interpreter = 'latex';
    an.FontSize = 20;
    an.LineStyle = 'none';
else
    for n = 1:1:lst
        k = 0;
        for m = t1:s:t2
            k = k+1; %Create a new index for the subsampled signal
            Isat_diff_temp(n,k) = obj.data.Isat_diff(st(n),m,r);
            t_temp(k) = obj.data.t(m); %Create a subsampled time series
        end
    end
    
    for n = 1:1:lst
        [H(n),C(n)] = EntropyComplexity(Isat_diff_temp(n,:),d);
    end
    
    CHFig = figure(1);
    hold on
    for n = 1:1:lst
        a(n) = plot(H(n),C(n),'o');
        l{1,n} = "Shot"+st(n);
        drawnow;
        pause(0.5)
    end
    if strcmp(Inpt.Results.Additional,'Mean')
        a(n+1) = plot(mean(H),mean(C),'s');
        l{1,n+1} = "Mean";
        drawnow;
        pause(0.5)
    end
    plot(H_min, C_min,'k-','LineWidth',1.5); %Min complexity curve
    plot(H_max, C_max,'k-','LineWidth',1.5); %Max complexity curve
    plot(fBm_CH.H,fBm_CH.C,'k--'); %Brownian Motion Curve
    hold off
    ax = gca;
    ax.FontSize = 14;
    ax.XLabel.String = '$H$';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.String = '$C$';
    ax.YLabel.Interpreter = 'latex';
    ax.TickLabelInterpreter = 'latex';
    ax.Layer = 'top';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.TickLength = [0.03 0.035];
    ax.Title.String = "Port 30, Subsampled, "+obj.data.t(t1)+"-"+obj.data.t(t2)+" s, "+obj.data.r(r)+"cm";
    ax.Title.FontSize = 14;
    ax.Title.Interpreter = 'latex';
    % ax.YLim = [0 0.2];
    ax.XLim = [0 1.3];
    ax.Box = 'on';
    % Annotations (Textbox)
    an = annotation('textbox',[0.2 0.7 0.1 0.1]);
    an.String = "d = "+d;
    an.Interpreter = 'latex';
    an.FontSize = 20;
    an.LineStyle = 'none';
    
    lgd = legend(a,l);
    lgd.Box = 'off';
    lgd.Interpreter = 'latex';
end

%% Saving
if strcmp(Inpt.Results.Save,'none') == 0
    if strcmp(Inpt.Results.Format,'fig')
        savefig(Inpt.Results.Save)
    elseif strcmp(Inpt.Results.Format,'png') || strcmp(Inpt.Results.Format,'eps') || strcmp(Inpt.Results.Format,'jpg')
        export_fig(Inpt.Results.Save,'-'+Inpt.Results.Format)
    elseif endsWith(Inpt.Results.Save,'.png') || endsWith(Inpt.Results.Save,'.eps') || endsWith(Inpt.Results.Save,'.jpg')
        export_fig(Inpt.Results.Save)
    elseif endsWith(Inpt.Results.Save,'.fig')
        savefig(Inpt.Results.Save)
    else
        error('File Type not specified')
    end
end

end