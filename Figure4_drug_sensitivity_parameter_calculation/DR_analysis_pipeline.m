function DR_analysis_pipeline(date,cellline,pathtofolder,drug,doses,colorarray,nr,dosetextexcel)

%Carolin Ector, 11.09.2023
%modified by CE 20.11.2023 to include different seeding densitites of the cisplatin experiment

%Function processes growth data of dose-response (DR) experiments for multiple drugs and cell lines
%runs the following scripts:
% DR_growth_plots: Plot growth curves from DR experiments
% DR_exponential_fit: Fit an exponential function to growth curves from DR experiments
% DR_response_curve_fit: Compute growth rate inhibition parameters

%input: stored in '[date]_DR_workspace.mat'
% date: date of the experiment being analyzed
% cellline: names of the cell lines being analysed
% drug: names of drugs used for drug treatments of the dose-response experiments
% doses: doses of each drug administered to the cells
% v: column order in excel sheet to load data as follows:
% 1. solvent control
% 2. lowest dose (dose(1)),
% 3. intermediate doses (dose(2),...,dose(n-1))
% 4. highest dose (dose(n))
% colorarray: colors for different drug doses

% Define remaining variables
experiment = str2num(date);
channel = {'CellNr';'Conf'};
dd = numel(drug); %loop drugs
ee = (numel(doses{1})+1); %loop doses, +1 = control

%create folders to save the different figures and result files in
outputfolder_plots = 'DR_plots';
outputfolder_excel = 'DR_results';
mkdir(outputfolder_plots);
mkdir(outputfolder_excel); 

%% Load and normalize data
for a = 1:2 %loop a channel

    k_all = cell(1,dd);
    CIlow_all = cell(1,dd);
    CIup_all = cell(1,dd);
    Rsq_all = cell(1,dd);
    ystart = nan(4,15);

    cc =  numel(cellline);

    startingvalues = nan(4,cc);

    if a == 1
        cc = cc-2; % no cell number data for MDAMB436 & BT549
    end

    for c = 1:cc %loop c cell lines

        disp(cellline{c})

        % Load excel file where data is stored
        filename = append(date,'_DR_',cellline{c},'.xlsx');
        pathtofile = append(pathtofolder,filename);
        sheets = sheetnames(pathtofile);

        is_metric = contains(sheets, channel{a,:});
        filtered_sheets = sheets(is_metric,:);

        densities = numel(filtered_sheets);

        for dens = 1:densities

            inputsheet = append(filtered_sheets{dens,:});

            [num] = xlsread(pathtofile,inputsheet);

            for d = 1:dd %loop d drug

                %define order how to load columns (data organization is partially different for the different experiments)
                %goal: load columns serially by increasing doses starting with the control
                if c == 1 || c == 10
                    yvalues = [3,12,10,8,6,4,13,11,9,7,5];
                else
                    yvalues = [3,14,12,10,8,6,4,13,11,9,7];
                end

                %load x-values (elapsed recording time)
                xdata=num(:,2);

                for e = 1:ee %loop e concentrations

                    %load and smooth growth data
                    ysmooth=smooth(num(:,yvalues(e)),0.35,'rloess');
                    ydata_smooth(:,e) = ysmooth;
                    stdev_smooth(:,e) = smooth(num(:,(yvalues(e)+12)),0.35,'rloess');

                     % 1. normalize ydata to time of treatment (= row number)
                    timeoftreat = 1;
                    if c == 10 || c == 13
                        timeoftreat = 20;
                    elseif c == 12
                        timeoftreat = 19;
                    end

                    imaginginterval = xdata(5,1)-xdata(4,1);

                    % adjust lengths of datasets to range from time of treatment to 96h post treatment
                    row_end = timeoftreat + 96/imaginginterval; % find corresponding row number where time after treatment = 96h
                    ysmooth1 = ysmooth(timeoftreat:row_end,:);

                    %initial value cell number / confluency for controls
                    if e == 1
                        startingvalues(dens,c) = mean(ysmooth1(1:3,:));
                    end

                    ydata_normx0(:,e) = ysmooth1./ysmooth1(1,:);
                    stdev_normx0(:,e) = stdev_smooth(timeoftreat:row_end,e)./ysmooth1;

                    % save values in cell array for exponential fit
                    yNorm_all{e,d} = ydata_normx0(:,e);
                    yStd_all{e,d} = stdev_normx0(:,e);
                    ymin_all{e,d} = min(ydata_normx0(:,e));
                    ymax_all{e,d} = median(ydata_normx0(end-3:end,e));
                    
                    % 2. normalize ydata_normx0 to control over time
                    ydata_normctrl(:,e) = ydata_normx0(:,e)./ydata_normx0(:,1);
                    stdev_normctrl(:,e) = stdev_normx0(:,e)./ydata_normx0(:,1);

                    ySmooth_all{e,d} = ydata_smooth(timeoftreat:row_end,e);
                    ySmooth_std_all{e,d} = stdev_smooth(timeoftreat:row_end,e);

                    valuestoclear1 = {'ysmooth'};
                    clear(valuestoclear1{:});

                    % 3. extract final responses (for sigmoidal curve fit, traditional EC50 extraction)
                    finalresponse{e,d} = median(ydata_normctrl(end-3:end,e));
                    sd_finalresponse{e,d} = median(stdev_normctrl(end-3:end,e));

                    if e == 1
                        ystart(dens,c) = ysmooth1(1,:);
                    end

                end %loop concentration
    
                %save smoothed raw data  of controls
                controls(:,1) = xdata(timeoftreat:row_end,:);
                controls(:,dens+1) = ySmooth_all{1,d};
                controls(:,dens+1+densities) = ySmooth_std_all{1,d};

                additionaltext = append('_',inputsheet);

                %%plot growth curves (smoothed raw growth, normalized growth to t=0 and relative growth to control)
                DR_growth_plots(a,c,d,date,channel,cellline,drug,doses,colorarray,imaginginterval,xdata,ydata_smooth,stdev_smooth,ydata_normx0,stdev_normx0,ydata_normctrl,stdev_normctrl,additionaltext,outputfolder_plots)

            end %loop d drug

            %exponential fit of growth curves
            DR_exponential_curve_fit(a,c,date,channel,cellline,drug,doses,nr,dosetextexcel,imaginginterval,yNorm_all,yStd_all,ymin_all,ymax_all,additionaltext,outputfolder_plots,outputfolder_excel)

            %sigmoidal curve fit - Hafner et al. 2017 approach
            [parameters] = DR_response_curve_fit(a,c,date,cellline,drug,doses,additionaltext,outputfolder_plots,outputfolder_excel);

            close all

            % %create excel file where drug sensitivity parameters will be store
            outputsheet = append('parameters_',additionaltext);
            output = array2table(parameters,"VariableNames",drug);
            valuenames = {'GEC50';'GR50';'GRinf';'Hill';'GRAOC'};
            rownames = cell2table(valuenames);
            outputtable = [rownames,output];

            writetable(outputtable,'drug_sensitivity_parameters.xlsx','sheet',outputsheet);

            valuestoclear2 = {'xdata';'yvalues_1';'yvalues_2';'ydata_smooth';'stdev_smooth';'ydata_normx0';'stdev_normx0';
                'ydata_normx0';'y_normx0_2';'yval_normctrl_1';'yval_normctrl_2';'ydata_normctrl';'stdev_normctrl';
                'finalresponse';'sd_finalresponse';'parameters';'growth_GRinf';'growth_GRinf_std';'columnheader'};
            clear(valuestoclear2{:});

        end %loop dens densitities

        outputexcel_growth_controls = append('growth_controls_',channel{a},'_norm.xlsx');
        output2 = array2table(controls);
        writetable(output2,outputexcel_growth_controls,'sheet',cellline{c});

        clear controls
        clear output2

    end %loop c cell lines

    outputexcel_startingval = append('startingvalues.xlsx');
    output3 = array2table(startingvalues,"VariableNames",cellline);
    writetable(output3,outputexcel_startingval,'sheet',channel{a});
    clear output3

    valuestoclear3 = {'sheets';'filtered_sheets';'is_metric';};
    clear(valuestoclear3{:});

end %loop a channel

end %function