function DR_exponential_curve_fit(a,c,date,channel,cellline,drug,doses,nr,dosetextexcel,imaginginterval,yNorm_all,yStd_all,ymin_all,ymax_all,additionaltext,outputfolder_plots,outputfolder_excel)

%Carolin Ector, 27.09.2023
%modified by CE 20.11.2023 to include different seeding densitites of the cisplatin experiment

%Function fits an exponential function to growth curves in dose-response experiments.

%input: partially stored in '[date]_DR_workspace.mat', remaining variables are defined in 'DR_analysis_pipeline.m' script
% a: loop a channel being analyzed from live-imaging: 1 = cell number (red fluorescent channel), 2 = confluency (brightfield channel)
% c: loop c celllines
% date: date of the experiment being analyzed
% channel: channel being analyzed from live-imaging (see above)
% cellline: names of the cell lines being analysed
% drug: names of drugs used for drug treatments of the dose-response experiments
% doses: tested doses (eg. 1x8 cell, with each cell containing the range of the resp. drug)
% nr: number of plots per row in the subplot
% dosetextexcel: doses of each drug administered to the cells for titles
% imaginginterval: intervals of image acquisition during the live-imaging experiment
% yNorm_all/yStd_all: normalized growth curves + standard devation (normalization to growth at initial time point of the live recordings) 
% ymin_all/ymax_all: minimum and maximum normalized growth (used to define lower and upper limit of y-axis --> to be equal for all doses)
% denschosenformanuscript: density chosen for mansucript

disp('DR_exponential_curve_fit.m is executed')

% Define remaining variables
experiment = str2num(date);
yaxisnames = {'Cell Number';'Confluency'};
dd = numel(drug);
ee = (numel(doses{1})+1); %loop doses, +1 = control

%% Fit in loop
% Set up fittype and options.
ft = fittype( 'a*exp(x*k)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Lower = [0.9 -0.1]; % [a k]
opts.Robust = 'on';
opts.StartPoint = [1 0.5]; % [a k]
opts.Upper = [1.1 0.1]; % [a k]

for d = 1:dd
    fig = figure('Visible','off');

    if experiment == 2021
        fig.Position = [1,65,1439,732];
    else
        fig.Position = [500,64,940,733];
        fig.InnerPosition = [500,64,940,733];
        fig.OuterPosition = [500,64,940,812];
    end

    %define lower and upper limit of y-axis (to be equal for all doses)
    yRange_lower = min(cell2mat(ymin_all(:,d)))-min(cell2mat(ymin_all(:,d)))*0.5;
    yRange_upper = max(cell2mat(ymax_all(:,d)))+max(cell2mat(ymax_all(:,d)))*0.15;
    
    %define values for x-axis and for subtitles
    Time = (0:imaginginterval:96);
    dosesincluding0 = [0;doses{d}];

    for e = 1:ee %loop e doses

        yNorm = cell2mat(yNorm_all(e,d));
        yStd = cell2mat(yStd_all(e,d));

        %define initial and final cell number or confluency
        ystart = median(yNorm(1:3,:)); %initial
        SE_ystart = median(yStd(1:3,:)); %corresponding error
        yend = median(yNorm(end-3:end,:)); %final
        SE_yend = median(yStd(end-3:end,:)); %corresponding error

        % calculate doubling times
        doublingtime{e,d} = 96*log(2)/log(yend/ystart); %doubling time
        doublingtime_std{e,d} = (96 * log(2) / (log(yend / ystart))^2) * sqrt((SE_yend / yend)^2 + (SE_ystart / ystart)^2); %corresponding error

        [xData, yData] = prepareCurveData(Time', yNorm);
        [fitresult, gof] = fit( xData, yData, ft, opts );

        % plot exponential fit
        h=subplot(2,nr,e);
        hold all

        x = Time;
        y = transpose(yNorm);
        err = transpose(yStd);

        patch([x, fliplr(x)], [y - err, fliplr(y + err)], 'b','FaceAlpha', 0.1, 'EdgeAlpha', 0, 'HandleVisibility', 'off');
        plot( fitresult, xData, yData );

        ylabel(h,[]);
        xlabel(h,[]);
        vectcoeff=coeffvalues(fitresult);
        parameters=strcat(['Curve Fit',newline,'\color{red}k = ',num2str(vectcoeff(2),'%6.3f'),newline,'\color{black}R^2 = ',num2str(gof.rsquare,'%6.3f')]);
        well = string(dosesincluding0(e));
        firstlegendentry = strcat(['\bf',well{1},'\bf µM']);

        % Create labels and title
        legend( h, firstlegendentry, parameters, 'Location', 'NorthWest', 'Interpreter', 'tex','color','w');

        % Set the remaining axes and box properties
        ax = gca;
        grid on;
        ylim([yRange_lower yRange_upper]);
        xticks(0:24:144);
        set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
            'XGrid','off','Box','on','Color','none');

        %save exponential growthrates (k) with corresponding confidence bounds and Rsq value of the fit
        k_all{e,d} = vectcoeff(2);
        confBounds = confint(fitresult);
        CIlow_all{e,d} = vectcoeff(2)-confBounds(1,2); %lower confidence bounds k
        CIup_all{e,d} = confBounds(2,2)-vectcoeff(2); %upper confidence bounds k
        Rsq_all{e,d} = gof.rsquare;

        valuestoclear1 = {'yNorm';'yStd';'fitresult';'gof';'vectcoeff';'parameters';'well';'firstlegendentry'};
        clear(valuestoclear1{:});

    end %loop e doses

    ylabeltext = append('Normalized ', yaxisnames{a});
    titletext = append(cellline{c},' ',drug{d},newline,'f(x) = a*exp(x*k)',newline,' ');

    % Give common xlabel, ylabel and title to your figure
    han=axes(fig,'visible','off');
    han.Title.Visible='on'; han.XLabel.Visible='on'; han.YLabel.Visible='on';
    xlabel(han,'Time post treatment (h)','FontWeight','bold','FontSize',11);
    ylabel(han,ylabeltext,'FontWeight','bold','FontSize',11);
    title(han,titletext,'FontSize',11);

    %save figure
    filetext = append(outputfolder_plots,'/',date,'_DR_',cellline{c},'_',drug{d},'_Fig4_ExpFit_combined',additionaltext,'.svg');
    saveas(fig, filetext);

end %loop d drugs

%create excel file where values of the exponential fit will be stored
outputfile=append(outputfolder_excel,'/',date,'_DR_Parameters_',cellline{c},additionaltext,'.xlsx');
valuestosave = {k_all,CIlow_all,CIup_all,Rsq_all,doublingtime,doublingtime_std};
outputsheets = {'exp_k';'exp_CIlow';'exp_CIup';'exp_Rsq';'doublingtime';'doublingtime_std'};

for i = 1:numel(valuestosave)
    output = cell2table(valuestosave{i},"VariableNames",drug);
    rownames = cell2table(dosetextexcel);
    outputtable = [rownames,output];
    outputsheet = append(outputsheets{i},'_',channel{a});
    writetable(outputtable,outputfile,'sheet',outputsheet);
end

end %function