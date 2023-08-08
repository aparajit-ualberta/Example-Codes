function [H,C,CH_Mean_Fig] = CHPlane_PosnRange_Mean(obj,varargin)
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
% t3 = length(obj.data.t);
d = Inpt.Results.Embedding_Dimension;
s = Inpt.Results.Subsample_Dimension;
rend = find(obj.data.r >= Inpt.Results.End_Radius,1);
[H_min, C_min] = MinComplexCurve(d); %Minimum Complexity Curve
[H_max, C_max] = MaxComplexCurve(d); %Maximum Complexity Curve

% [Isat_diff_temp, ~] = obj.subsample('Start_Time',Inpt.Results.Start_Time,'Stop_Time',Inpt.Results.Stop_Time,'End_Radius',Inpt.Results.End_Radius,'Subsample_Dimension',s);
Isat_diff_temp = obj.data.Isat_diff(:,t1:t2,:);

for m = 1:1:rend
    for n = 1:1:10
%         [H(n,m), C(n,m)] = EntropyComplexity(Isat_diff_temp(n,:,m),d);
        [H(n,m), C(n,m)] = EntropyComplexity(Isat_diff_temp(n,:,m),d,s);
    end
end
%% Plotting
CH_Mean_Fig = figure(1);
hold on
for m = 1:1:rend
    a(m) = plot(mean(H(:,m)),mean(C(:,m)),'Marker','o','LineStyle','none');%,'Color','r') %Plotting mean of all shots
    l{1,m} = ""+obj.data.r(m)+" cm";
    drawnow;
    pause(0.5)
end
plot(H_min, C_min,'k-','LineWidth',1.5); %Min complexity curve
plot(H_max, C_max,'k-','LineWidth',1.5); %Max complexity curve
plot(fBm_CH.H,fBm_CH.C,'k--'); %Brownian Motion Curve
hold off

ax2 = gca;
ax2.FontSize = 14;
ax2.XLabel.String = '$H$';
ax2.XLabel.Interpreter = 'latex';
ax2.YLabel.String = '$C$';
ax2.YLabel.Interpreter = 'latex';
ax2.TickLabelInterpreter = 'latex';
ax2.Layer = 'top';
ax2.XMinorTick = 'on';
ax2.YMinorTick = 'on';
ax2.TickLength = [0.03 0.035];
ax2.Title.String = "Port 30, Subsampled, "+obj.data.t(t1)+"-"+obj.data.t(t2)+" s";
ax2.Title.FontSize = 14;
ax2.Title.Interpreter = 'latex';
% ax.YLim = [0 0.2];
ax2.XLim = [0 1.3];
ax2.Box = 'on';
% lgd2 = legend(a, {'0 cm','0.1 cm','0.2 cm','0.3 cm','0.4 cm','0.5 cm','0.6 cm','0.7 cm','0.8 cm','0.9 cm','1 cm','1.1 cm','1.2 cm','1.3 cm','1.4 cm','1.5 cm','1.6 cm','1.7 cm','1.8 cm','1.9 cm','2 cm','2.1 cm','2.2 cm','2.3 cm','2.4 cm','2.5 cm','2.6 cm','2.7 cm'});
lgd2 = legend(a,l);
lgd2.Box = 'off';
lgd2.Interpreter = 'latex';
% Annotations (Textbox)
an2 = annotation('textbox',[0.2 0.7 0.1 0.1]);
an2.String = "d = "+d;
an2.Interpreter = 'latex';
an2.FontSize = 20;
an2.LineStyle = 'none';

if strcmp(Inpt.Results.Save,'none') == 0
    switch Inpt.Results.Format
        case 'fig'
            savefig(Inpt.Results.Save)
        case 'png'
            export_fig(Inpt.Results.Save,'-'+Inpt.Results.Format)
        case 'eps'
            export_fig(Inpt.Results.Save,'-'+Inpt.Results.Format)
        case 'jpg'
            export_fig(Inpt.Results.Save,'-'+Inpt.Results.Format)
        otherwise
            if endsWith(Inpt.Results.Save,'.png') || endsWith(Inpt.Results.Save,'.eps') || endsWith(Inpt.Results.Save,'.jpg')
                export_fig(Inpt.Results.Save)
            elseif endsWith(Inpt.Results.Save,'.fig')
                savefig(Inpt.Results.Save)
            else
                error('File Type not specified')
            end
    end
end            
end

