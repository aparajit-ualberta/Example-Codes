function [H,C,fig1] = CHPlane_PosnRange_Shot(obj,varargin)
    warning off
    Inpt = inputParser;

    addRequired(Inpt, 'Time_Series', @(obj) isobject(obj))
    addParameter(Inpt, 'Embedding_Dimension',5, @(d) isnumeric(d))
    addParameter(Inpt, 'Subsample_Dimension',8, @(s) isnumeric(s))
    addParameter(Inpt, 'Start_Time',5, @(tstart) isnumeric(tstart))
    addParameter(Inpt, 'Stop_Time',15, @(tstop) isnumeric(tstop))
    addParameter(Inpt, 'End_Radius',0.9, @(rend) isnumeric(rend) && rend <= 2.7)
    addParameter(Inpt, 'Save','none', @(sv) ischar(sv))
    addParameter(Inpt, 'Format','none', @(fm) ismember(fm, {'fig','png','eps','jpg'}))

    parse(Inpt, obj, varargin{:})

    load('fBm_CH_Curve','fBm_CH')
    t1 = find(obj.data.t >= Inpt.Results.Start_Time,1);
    t2 = find(obj.data.t >= Inpt.Results.Stop_Time,1);
%     t3 = length(obj.data.t);
    d = Inpt.Results.Embedding_Dimension;
    s = Inpt.Results.Subsample_Dimension;
    rend = find(obj.data.r >= Inpt.Results.End_Radius,1);
    [H_min, C_min] = MinComplexCurve(d); %Minimum Complexity Curve
    [H_max, C_max] = MaxComplexCurve(d); %Maximum Complexity Curve
    
    Fs = 1/(obj.data.t(2)-obj.data.t(1));
    Fl = Fs/s;
    if s>1
        for rm = 1:1:28
            for sn = 1:1:10
                Isat_Filt(sn,:,rm) = lowpass(obj.data.Isat_diff(sn,:,rm),Fl,Fs);
            end
        end
        Temp_Obj = IsatCHPlane(obj.data.t,obj.data.r,Isat_Filt,obj.data.Isat_DC);
    else
        Temp_Obj = obj;
    end

    [Isat_diff_temp, ~] = obj.subsample('Time_Series',Temp_Obj,'Start_Time',Inpt.Results.Start_Time,'Stop_Time',Inpt.Results.Stop_Time,'End_Radius',Inpt.Results.End_Radius,'Subsample_Dimension',s);

    mkr = ['o' '+' '*' 'x' 's' 'd' 'h' 'p'];
    clr = ['k' 'k' 'r' 'r' 'r' 'r' 'r' 'r'];
    for m = 1:1:rend
        for n = 1:1:10
            [H(n,m), C(n,m)] = EntropyComplexity(Isat_diff_temp(n,:,m),d);
        end
    end
    
    HR = mean(H,1);
    HS = std(H,1);
    CR = mean(C,1);
    CS = std(C,1);

    fig1 = figure(1);
    fig1.Units = 'centimeters';
    fig1.Position = [12.8852 0 20 18];
    fig1.Color = 'white';
    
    hold on
    for m = 1:1:rend
        if m <= 8
%             as(m) = plot(H(1:1:10,m),C(1:1:10,m),'Marker',mkr(m),'Color',clr(m),'LineStyle','none');
            as(m) = plot(HR(m),CR(m),'Marker',mkr(m),'Color',clr(m),'LineStyle','none');
%             as(m) = errorbar(HR(m),CR(m),CS(m),CS(m),HS(m),HS(m),'Marker',mkr(m),'Color',clr(m),'LineStyle','none');
            drawnow;
        end
    end

    if rend > 8
        as(9) = plot(mean(mean(H(:,9:rend))),mean(mean(C(:,9:rend))),'Marker','^','MarkerSize',10,'Color','b','LineStyle','none','LineWidth',2);
        drawnow;
    end

    plot(H_min, C_min,'k-','LineWidth',2); %Min complexity curve
    plot(H_max, C_max,'k-','LineWidth',2); %Max complexity curve
    plot(fBm_CH.H,fBm_CH.C,'k--','LineWidth',2); %Brownian Motion Curve
    hold off

    ax1 = gca;
    ax1.Position = [0.096171802054155,0.077753779697624,0.873015873015873,0.873650107991359];
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
    % ax.YLim = [0 0.2];
    ax1.XLim = [0 1];
    ax1.Box = 'on';
    lgd1 = legend(as, {'$0$ cm','$0.1$ cm','$0.2$ cm','$0.3$ cm','$0.4$ cm','$0.5$ cm','$0.6$ cm','$0.7$ cm','Mean of Edge'});
    lgd1.Box = 'off';
    lgd1.Position = [0.12 0.65 0.24 0.24];
    lgd1.Interpreter = 'latex';
    % Annotations (Textbox)
    an1 = annotation('textbox',[0.3 0.8 0.15 0.1]);
    an1.String = {'$t = 4 - 15 ms$','$d = 5$'};
    an1.Interpreter = 'latex';
    an1.FontSize = 20;
    an1.LineStyle = 'none';
    
    an2 = annotation('textbox',[0.59 0.49 0.15 0.15]);
    an2.String = 'fBm';
    an2.Interpreter = 'latex';
    an2.FontSize = 18;
    an2.LineStyle = 'none';

    subax1 = axes('Position',[0.35,0.18,0.36,0.2]);
    subax1.Box = 'on';
    IsatR_mean = mean(mean(obj.data.Isat_DC));
    yyaxis left
    hold on
    asub(1) = plot(obj.data.r(1:2),permute(IsatR_mean(1,1,1:2), [3 2 1]),'k-','LineWidth',1.5);
    asub(2) = plot(obj.data.r(2:8),permute(IsatR_mean(1,1,2:8), [3 2 1]),'r-','LineWidth',1.5);
    asub(3) = plot(obj.data.r(8:rend),permute(IsatR_mean(1,1,8:rend), [3 2 1]),'b-','LineWidth',1.5);
    area(obj.data.r(1:2),permute(IsatR_mean(1,1,1:2), [3 2 1]),'FaceColor','k','FaceAlpha',0.5,'EdgeColor','k','LineStyle','none');
    area(obj.data.r(2:8),permute(IsatR_mean(1,1,2:8), [3 2 1]),'FaceColor','r','FaceAlpha',0.5,'EdgeColor','r','LineStyle','none');
    area(obj.data.r(8:rend),permute(IsatR_mean(1,1,8:rend), [3 2 1]),'FaceColor','b','FaceAlpha',0.5,'EdgeColor','b','LineStyle','none');
    hold off
    subax1.YColor = 'black';
    subax1.FontSize = 14;
    subax1.XLabel.Interpreter = 'latex';
    subax1.XLabel.String = '$r$ (cm)';
    subax1.YLabel.Interpreter = 'latex';
    subax1.YLabel.String = '$I_{sat}$';
    subax1.TickLabelInterpreter = 'latex';
    subax1.Layer = 'top';
    subax1.XMinorTick = 'on';
    subax1.YMinorTick = 'on';
    subax1.TickLength = [0.03 0.035];
    subax1.Box = 'on';
%     subax1.XLim = [0 2.7];
    
    yyaxis right
    hold on
    asub(4) = plot(obj.data.r(1:rend),HR,'m--','LineWidth',1.5);
    asub(5) = plot(obj.data.r(1:rend),CR,'c--','LineWidth',1.5);
    hold off
    subax1.YColor = 'black';
    subax1.FontSize = 14;
    subax1.YLabel.Interpreter = 'latex';
    subax1.YLabel.String = '$H,C$';
    subax1.TickLabelInterpreter = 'latex';
    subax1.TickLabelInterpreter = 'latex';
    subax1.Layer = 'top';
    subax1.YMinorTick = 'on';
    subax1.TickLength = [0.03 0.035];
    subax1.Box = 'on';
    subax1.XLim = [0 2.7];
    
    sublgd1 = legend(asub, {'Center','Gradient','Edge','Entropy','Complexity'});
    sublgd1.Box = 'off';
    sublgd1.Interpreter = 'latex';
    
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

