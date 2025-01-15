function f3_fit_exponential_decay_amplitude(detrended_all,envelope_all,celllinenames_file,reporters,reporter_colors,axOpt,recordingtime_all)

%Carolin Ector, 14.11.2023

%Function fits an exponential function to a circadian signal's envelope and calculates the exponential decay rate.

%Clock-TNBC Manuscript, part of Deep Circadian Phenotyping Appraoch 

%input: stored in "amplitudedecay_autocorrelation_MRA_pipeline.mat" & calculated in preceding function "f1_sort_lumicycle_and_mra_data.m"
% detrended_all: detrended_all: 2x19 cell array containing detrended Luc-data for different cell lines (columns). 
    %Row 1 contains Bmal1-Luc data; Row 2 contains Per2-Luc data.
    %each cell stores a double matrix where rows are time points and columns are replicates for the corresponding reporter and cell line.
% envelope_all: see "detrended_all", but compiled of signal's amplitude (envelope) data
% celllinenames_file: names of the cell lines being analysed, as written in the file names
% reporters: circadian gene names for the luciferase reporters
% reporter_colors: colors used in the graphs for the two circadian clock luciferase reporters
% axOpt: setting for appearance of axes of a plot
% recordingtime_all: total recording time as a time series (1-137.7 hours) stored for replicate1

disp('f3_fit_exponential_decay_amplitude.m is executed')

emptyarray = nan(8,numel(celllinenames_file));

for a = 1:numel(reporters) %loop a Luc-reporters

    Rsq = emptyarray;
    decay = emptyarray;

    for c = 1:length(celllinenames_file) %loop c celllines

        detrended = detrended_all{a,c};
        envelope = envelope_all{a,c};
        numberofsamples = size(detrended,2);

        fig = figure('Visible','off');
        fig.Position = [1921,1,1440,821];

        for r = 1:numberofsamples %loop r replicates

            %exponential decay
            fitresult = cell( r, c );
            gof = struct( 'sse', cell( r, c ), 'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] ); %gof is somehow not saved to the workspace

            [xData, yData] = prepareCurveData(recordingtime_all, envelope(:,r));

            startvalue = median(envelope(1:3,r));
            rangestart = startvalue*0.1;

            equation = 'a*exp(x*k)';

            % Set up fittype and options.
            ft = fittype( equation, 'independent', 'x', 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
            opts.Robust = 'on';
            opts.Lower = [(startvalue-rangestart),-0.1]; % [k] %exp decay
            opts.StartPoint = [startvalue,0.001]; % [k] %exp decay
            opts.Upper = [(startvalue+rangestart),0.1]; % [k]  %exp decay

            s=subplot(2,ceil(numberofsamples/2),r);
            hold all

            [fitresult{r,c}, gof(r,c)] = fit( xData, yData, ft, opts );
            % Determine the range for the fit
            xRange = min(recordingtime_all):0.01:max(recordingtime_all);
            % Evaluate the fit over this range
            yFit = feval(fitresult{r,c}, xRange);
            % Plot the fit result over the specified range
            plot(xRange, yFit, 'k', 'LineWidth', 2);
            % Plot [xData, yData]
            plot(xData, yData, 'Color', reporter_colors{a}, 'LineWidth', 2);
            vectcoeff=coeffvalues(fitresult{r,c});
            Rsq(r,c) = gof(r,c).rsquare;

            ylabel(s,[]);
            xlabel(s,[]);

            rsqval = round(Rsq(r,c),1);

            if rsqval < 0.6
                decay(r,c) = NaN;
                nanval=1;
            else
                decay(r,c) = vectcoeff(2);
                nanval=0;
            end

            %Set the remaining axes and box properties
            ax = gca;
            grid on;
            xticks(0:24:120);
            set(ax,axOpt{:});

            if nanval == 1
                parameters=strcat(['\color[rgb]{0.5,0.5,0.5}k = ',num2str(vectcoeff(2),'%6.3f'),newline,'\color{black}R^2 = ',num2str(rsqval)]);
            else
                parameters=strcat(['\color{red}k = ',num2str(vectcoeff(2),'%6.3f'),newline,'\color{black}R^2 = ',num2str(rsqval)]);
            end

            %Create labels and title
            legend( s, parameters,'Location', 'NorthEast', 'Interpreter', 'tex','color','w','FontSize',14);
            samplename = append(reporters{a},'-Luc ',num2str(r));
            title(samplename,'FontSize',14,'Interpreter','none');

            varstoclear1 = {'maxvalue','xData','yData','samplename','yFit'};
            clear(varstoclear1{:})

        end %loop r replicates

        ylabeltext = append(['Envelope (a.u.)',newline,' ']);
        xlabeltext = append([' ',newline,'Time (h)']);
        titletext = append(['f(x) = ',equation,newline,' ']);

        % Give common xlabel, ylabel and title to your figure
        han=axes(fig,'visible','off');
        han.Title.Visible='on'; han.XLabel.Visible='on'; han.YLabel.Visible='on';
        xlabel(han,xlabeltext,'FontWeight','bold','FontSize',22);
        ylabel(han,ylabeltext,'FontWeight','bold','FontSize',22);
        title(han,titletext,'FontSize',22);

        %save figure
        %filetext = append(pathtofigure,'Exponential_decay_envelope_',celllinenames_file{c},'_',reporters{a}); %exp decay
        %saveas(fig, [ filetext, '.svg']);

        %close all

    end % loop c celllines

    allfinalvaluesbyrep = {decay;Rsq};
    outputsheets = {'expdecay_envelope';'rsq_expdecay_envelope'};
    circadian_parameters_excel = append('extracted_circadian_parameters_by_replicate_',reporters{a},'.xlsx');
    replicate = {'1';'2';'3';'4';'5';'6';'7';'8'};
    t_row = cell2table(replicate);

    %save individual values by replicate
    for n = 1:numel(allfinalvaluesbyrep)
        valuestosave = allfinalvaluesbyrep{n};
        t_valuestosave = array2table(valuestosave,'VariableNames',celllinenames_file);
        t_final = [t_row,t_valuestosave];
        writetable(t_final,circadian_parameters_excel,'sheet',outputsheets{n});
        clear t_final
        clear t_valuestosave
    end

    varstoclear2 = {'allfinalvaluesbyrep','expdecay','Rsq','maxenvelope'};
    clear(varstoclear2{:})

end %loop a reporter

disp('f3_fit_exponential_decay_amplitude.m is completed')

end %function