function [parameters] = DR_response_curve_fit(a,c,date,cellline,drug,doses,additionaltext,outputfolder_plots,outputfolder_excel);

%Carolin Ector, 28.09.2023
%Adapted from the publication from Hafner, M. et al.  Growth rate inhibition metrics correct for confounders in measuring sensitivity to cancer drugs.
%Nat Methods 13, 521–527 (2016). https://doi.org/10.1038/nmeth.3853

%Time-of-Day-Cancer-Drugs Manuscript Fig. 3i, Fig. S3

%Function computes the growth rate inhibition (GR) per drug dose based on the growth rates obtained from an exponential curve fit
%Utilizes the GR values to fit a sigmoidal dose-response curve from which various drug sensitivity parameters are extracted (see below: output)

%input: partially stored in '[date]_DR_workspace.mat'
% a: loop a channel being analyzed from live-imaging: 1 = cell number (red fluorescent channel), 2 = confluency (brightfield channel)
% c: loop c celllines
% channel: channel being analyzed from live-imaging (see above)
% date: date of the experiment being analyzed
% cellline: names of the cell lines being analysed
% drug: names of drugs used for drug treatments of the dose-response experiments
% % doses: doses of each drug administered to the cells
% finalresponse: final nucleus counts (relative to control) for each drug dose
% sd_finalresponse: standard deviation of final nucleus counts (relative to control) for each drug dose
% parameters: array to store final drug sensitvity parameters

%output:
% parameters: drug sensitivty parameters computed from the sigmoidal curve fit:
% parameters(1,:) = GEC50, dose at half-maximal drug effect
% parameters(2,:) = GR50, dose GR=0.5
% parameters(3,:) = GRinf, maximal drug effect
% parameters(4,:) = Hill Coefficient, slope of the dose-response curve
% parameters(5,:) = AOC, area over the dose-response curve
% dose_GRinf: closest and smallest dose eliciting a GRinf response with corresponding index in Doserange array

disp('DR_response_curve_fit.m is executed')

%Define remaining variables
dd = numel(drug);
ee = numel(doses{1})+1; %+1 = control
channel = {'CellNr';'Conf'};
yaxisnames = {'Cell Number';'Confluency'};
measurement = append('(',yaxisnames{a},')');

%Load data
inputfile = append(outputfolder_excel,'/',date,'_DR_Parameters_',cellline{c},additionaltext,'.xlsx');

%exponential growth rate (k)
sheet1 = append('exp_k_',channel{a});
[t_k] = readtable(inputfile,'sheet',sheet1);

%lower confidence intervals of the exponential growth rate (k)
sheet2 = append('exp_CIlow_',channel{a});
[t_CIlow] = readtable(inputfile,'sheet',sheet2);

%upper confidence intervals of the exponential growth rate (k)
sheet3 = append('exp_CIup_',channel{a});
[t_CIup] = readtable(inputfile,'sheet',sheet3);

for d = 1:dd %loop d drug

    Doserange = doses{d};
    factor = 2;

    xaxislabels = [Doserange(1,:)/factor;Doserange;Doserange(end,:)*factor];
    xaxislabels = xaxislabels(1:2:end);

    k = table2array(t_k(:,d+1));
    CIlow = table2array(t_CIlow(:,d+1));
    CIup = table2array(t_CIup(:,d+1));

    %calculate normalized growth rate inhibition
    for e = 1:ee
        GR(e,:) = 2^(k(e,:)/k(1,:))-1;
        error1(e,:) = 2^(CIlow(e,:)/k(1,:))-1;
        error2(e,:) = 2^(CIup(e,:)/k(1,:))-1;
    end

    %exclude controls from sigmoidal response curve fit
    GR(1,:) = [];
    error1(1,:) = [];
    error2(1,:) = [];

    titletext = append(cellline{c},' ',drug{d});

    Response = GR;
    [xData, yData] = prepareCurveData(Doserange, Response);

    minDose=min(Doserange);
    maxDose=max(Doserange);
    midDose=mean([minDose maxDose]);

    % Set up fittype and options.
    ft = fittype( 'GRinf+((1-GRinf)/(1+(c/GEC50)^hGR))', 'independent', 'c', 'dependent', 'GR' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Lower = [minDose/100 -1  0.5]; % [GEC50,GRinf,hGR] %previously hGR was set to 1
    opts.Robust = 'Bisquare';
    %     opts.Weights = Rsq;
    opts.StartPoint = [midDose 0 5]; % [GEC50,GRinf,hGR]
    opts.Upper = [maxDose*100 1 10]; % [GEC50,GRinf,hGR]

    % Fit model to data
    try
        [fitresult, gof] = fit( xData, yData, ft, opts);
        vectcoeff=coeffvalues(fitresult);
        parameters(1,d) = vectcoeff(1); %GEC50
        parameters(5,d) = trapz(Doserange, 1-Response); %AOC
        parameters(3,d) = vectcoeff(2); %GRinf
        parameters(4,d) = vectcoeff(3); %Hill Coefficient

        % Define a function that calculates the squared difference between the fitted curve and 0.5
        squared_difference = @(Doserange) (feval(fitresult, Doserange) - 0.5)^2;

        % Use fminsearch to minimize the squared difference and find the x value
        initial_guess = vectcoeff(1); % Your initial guess for x
        parameters(2,d) = fminsearch(squared_difference, initial_guess); %GR50

        %plot results
        fig = figure('Visible','off');

        xpoints=logspace(log10(minDose/factor),log10(maxDose*factor),1000);
        yfit = feval(fitresult, xpoints);
        h = plot(xpoints, yfit); hold on; % Fitted line on the extended range
        set(h,'lineWidth',1,'Color','r','HandleVisibility','off');
        errorbar(Doserange,GR,error1,error2,'o','MarkerSize',9,'MarkerFaceColor',[0 0 0],'LineWidth',1.5,'Color',[0 0 0],'HandleVisibility','off');

        %create dummyline for legend entries
        line(xData, yData, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none'); %GEC50
        line(xData, yData, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none'); %GR50
        line(xData, yData, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none'); %GRinf
        line(xData, yData, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none'); %Hill
        line(xData, yData, 'Linestyle', 'none', 'Marker', 'none', 'Color', 'none'); %AOC
        %     hold off

        %create legend with computed drug sensitvity values
        lgd_GEC50=strcat(['GEC_{50} = ',num2str(parameters(1,d),'%6.3f')]);
        lgd_GR50=strcat(['GR_{50} = ',num2str(parameters(2,d),'%6.3f')]);
        lgd_GRinf=strcat(['GR_{inf} = ',num2str(parameters(3,d),'%6.3f')]);
        lgd_Hill=strcat(['Hill = ',num2str(parameters(4,d),'%6.3f')]);
        lgd_AOC=strcat(['AOC = ',num2str(parameters(5,d),'%6.3f')]);

        legend(lgd_GEC50,lgd_GR50,lgd_GRinf,lgd_Hill,lgd_AOC,'Interpreter','tex','color','none','EdgeColor','none','FontSize',14);

        %label axes
        title(titletext,'FontSize',18,'FontWeight','normal');
        xlabel('Dose (µM)','FontSize',18,'FontName','Helvetica Neue');
        ylabeltext = append('GR ',measurement);
        ylabel(ylabeltext,'FontSize',18,'FontName','Helvetica Neue');

        ax = gca;

        set(ax,'XScale','log','XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
            'XGrid','off','Box','on','Color','none','FontSize',18);

        xticks(xaxislabels);
        xticklabels(string(xaxislabels));

        %save figure
        filetext = append(outputfolder_plots,'/',date,'_DR_',cellline{c},'_',drug{d},'_',channel{a},'_Fig5_RC_combined',additionaltext,'.svg');
        saveas(fig, filetext);

    catch ME

        parameters(1,d) = NaN; %GEC50
        parameters(2,d) = NaN;
        parameters(3,d) = NaN; %GRinf
        parameters(4,d) = NaN; %Hill Coefficient
        parameters(5,d) = NaN; %AOC
 
    end

    valuestoclear1 = {'Doserange';'GR';'error1';'error2';'xData';'yData';'GRinf';'closest_indices';'indices'};
    clear(valuestoclear1{:});

end %loop d drug
end %loop function