function f2_overlay_circadian_signals_mean(detrended_all,norm_all,celllinenames_file,reporters,reporter_colors,reporter_colors_dark,axOpt,recordingtime_all)

%Carolin Ector, 25.10.2023

%Function overlays signals from t he Bmal1 and Per2 luciferase reporters per cell line, averaged across all replicates
%Loops through detrended and normalized detrended signals (normalized to the signal envelope)

%Clock-TNBC Manuscript Box1

%input: stored in "amplitudedecay_autocorrelation_MRA_pipeline.mat" & calculated in preceding function "f1_sort_lumicycle_and_mra_data.m"
% detrended_all: detrended_all: 2x19 cell array containing detrended Luc-data for different cell lines (columns). 
    %Row 1 contains Bmal1-Luc data; Row 2 contains Per2-Luc data.
    %each cell stores a double matrix where rows are time points and columns are replicates for the corresponding reporter and cell line.
% norm_all: see "detrended_all", but compiled of amplitude (envelope)-normalized data
% celllinenames_file: names of the cell lines being analysed, as written in the file names
% reporters: circadian gene names for the luciferase reporters
% reporter_colors: colors used in the graphs for the two circadian clock luciferase reporters
% reporter_colors_dark: darker version of colorbp to use for text and thin lines
% axOpt: setting for appearance of axes of a plot
% recordingtime_all: total recording time as a time series (1-137.7 hours) stored for replicate1

disp('f2_overlay_circadian_signals_mean.m is executed')

datasets = {detrended_all;norm_all};
description = {'detrended','normalized'};

for i = 1:2 %loop i datasets

    dataset=datasets{i,1};
    
    for c = 1:length(celllinenames_file) %loop c celllines

        %create figure
        fig = figure;%('Visible','off');
        hold on

        % Iterate over each reporter
        for a = 1:numel(reporters) %loop a Luc-reporters

            dataSubset = dataset{a,c};

            if isempty(dataSubset) 
                % Handle the empty case
            else
                % Calculate the mean across replicates (columns)
                meanValues = mean(dataSubset, 2, 'omitnan'); % 'omitnan' ignores NaN values
                stdValues = std(dataSubset, [], 2, 'omitnan');

                y1 = transpose(meanValues);
                e1 = transpose(stdValues);
                x1 = recordingtime_all';

                if a == 1
                    yyaxis left
                else
                    yyaxis right
                end
                patch([x1 fliplr(x1)], [(y1-e1)  (fliplr(y1+e1))], reporter_colors{a}, 'FaceAlpha',0.1, 'EdgeAlpha',0, 'HandleVisibility','off');
                plot(x1,y1,'LineWidth',2.5,'LineStyle','-','Color',reporter_colors{a});
                ylabeltext1 = append(reporters{a},'-Luc signal');
                ylabel(ylabeltext1);
                leftAxisLim = max(abs([y1-e1, y1+e1]));
                ylim(leftAxisLim * [-1, 1]);
                leftAxis = gca;
                leftAxis.YColor = reporter_colors_dark{a};
            end
            clear dataSubset
        end %loop a reporter

        hold off

        xlabel('Time (days)','FontSize',20,'FontName','Helvetica Neue');

        titletext = append(celllinenames_file{c}, '-',description{i});
        title(titletext,'FontSize',20,'FontName','Helvetica Neue','Interpreter','none');

        ax = gca;
        grid on;
        xticks(0:24:144);
        xticklabels({'0','1','2','3','4','5','6'});

        set(ax,axOpt{:});

        % add legend
        lgd = legend(reporters);
        set(lgd,'FontSize',15,'Orientation','vertical','Location','northeast','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');

        %save figure
        %filetext = append(description{i},'_mean_signal_bmal1_per2_overlay_',celllinenames_file{c},'.svg');
        %saveas(fig, filetext);

        varstoclear = {'x1','x2','y1','y2','e1','e2','Bmal1','Bmal1_mean','Bmal1_std','Per2','Per2_mean','Per2_std','lgd'};
        clear(varstoclear{:});

        close all

    end %loop c celllines

end %loop i datasets

disp('f2_overlay_circadian_signals_mean.m is completed')

end %function
