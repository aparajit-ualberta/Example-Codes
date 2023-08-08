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