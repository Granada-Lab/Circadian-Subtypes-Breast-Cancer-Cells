function f4_autocorrelation_calculation(detrended_all,reporters,celllinenames_file,reporter_colors,axOpt)

%Carolin Ector, 22.08.2023,
%adapted from https://www.mathworks.com/help/econ/autocorr.html

%%Function returns sample autocorrelation function (ACF) of the univariate time series y, including 95.4% ACF confidence bounds
%%For all time series, the lag 0 autocorrelation acf(1) = 1.

%Clock-TNBC Manuscript Fig. 2A and Box1

%input: stored in "amplitudedecay_autocorrelation_MRA_pipeline.mat"
% detrended_all: 2x19 cell array containing detrended Luc-data for different cell lines (columns).
%Row 1 contains Bmal1-Luc data; Row 2 contains Per2-Luc data.
%each cell stores a double matrix where rows are time points and columns are replicates for the corresponding reporter and cell line.
% celllinenames_file: names of the cell lines being analysed, as written in the file names
% reporters: circadian gene names for the luciferase reporters
% reporter_colors: colors used in the graphs for the two circadian clock luciferase reporters
% axOpt: setting for appearance of axes of a plot
% recordingtime_all: total recording time as a time series (1-137.7 hours)

disp('f4_autocorrelation_calculation.m is executed')

emptyarray = nan(8,numel(celllinenames_file));

for a = 1:2 %loop j reporter

    %create empty cell arrays to store data
    autocorr_peak = emptyarray;
    autocorr_lag = emptyarray;

    for c = 1:length(celllinenames_file) %loop c celllines

        detrended = detrended_all{a,c};
        numberofsamples = size(detrended,2);

        for r = 1:numberofsamples %loop r replicates

            data = rmmissing(detrended(:,r));
            total_lags = size(data,1);
            titletext = append(celllinenames_file{c},' ',num2str(r));

            fig = figure('Visible','off');
            hold on

            % calculate autocorrelation
            [acf,lags,bounds] = autocorr(data,NumLags=total_lags-1);
            h0 = stem(lags,acf);
            set(h0,'LineStyle','none','color',reporter_colors{a},'MarkerSize',4);

            %confidence intervals
            l1 = line(lags,bounds(1)*ones(length(acf),1));
            l2 = line(lags,bounds(2)*ones(length(acf),1));
            set(l1,'color','b','linewidth',1);
            set(l2,'color','b','linewidth',1,'HandleVisibility','off');

            %identify second peak with period >= 12h
            [pks,locs]=findpeaks(acf(72:end,:),lags(72:end,:),'MinPeakHeight',-0.2,'MinPeakProminence',0.01,'MinPeakWidth',12);

            TF = isempty(pks); %check if peak was detected

            if TF == 0
                peak = pks(1);
                if peak >= bounds(1,:) || peak <= bounds(2,:)
                    abcissa = locs(1)/6;
                    text(locs(1,:)+.03,pks(1,:),'2','FontSize',14) % to show peaks in plots by arrow
                else
                    text(locs(1,:)+.03,peak,'X','FontSize',14) % to show peaks in plots by arrow
                    peak = NaN;
                    abcissa = NaN;
                end
            else
                peak = NaN;
                abcissa = NaN;
            end

            hold off

            %save values
            autocorr_peak(r,c) = peak;
            autocorr_lag(r,c) = abcissa;

            %set graph appearance options
            xlabel('Lag (h)','FontSize',20,'FontName','Helvetica Neue');
            ylabel('Autocorrelation','FontSize',20,'FontName','Helvetica Neue');

            ax = gca;
            grid on;
            ylim([-1 1]);
            xlim([-36 900]);
            xticks(0:144:864);
            xticklabels({'0','24','48','72','96','120','144'});

            set(ax,axOpt{:});

            title(titletext,'FontSize',18,'FontName','Helvetica Neue','Interpreter', 'none');

            %figurename = append('autocorrelation_',celllinenames_file{c},'_',reporters{a},'_',num2str(r),'.svg');
            %saveas(fig,figurename);

            variabletoclear = {'data';'pks';'locs';'acf';'bounds';'lags';'total_lags'};
            clear(variabletoclear{:})

        end %loop r replicates

    end %loop c cellline

    allfinalvaluesbyrep = {autocorr_lag;autocorr_peak};
    outputsheets = {'autocorrelation_lag';'autocorrelation_peak'};
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

    varstoclear = {'replicates','autocorr_period','autocorr_peak','replicate','varnames'};
    clear(varstoclear{:})

end %loop a reporter

disp('f4_autocorrelation_calculation.m is completed')

end %function
