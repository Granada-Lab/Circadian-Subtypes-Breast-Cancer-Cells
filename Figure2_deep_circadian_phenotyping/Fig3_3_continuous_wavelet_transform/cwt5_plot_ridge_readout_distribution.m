function cwt5_plot_ridge_readout_distribution(reporters,reporter_colors,celllinenames_file,reporter_colors_dark)

%Carolin Ector, 02.11.2023

%Function plots distribution of continuous amplitudes and periods per cellline and per subtype

%Clock-TNBC Manuscript Fig. 2C

%input: stored in "cwt_analysis_pipeline.mat" & defined in cwt_analysis_pipeline.m (run in sequence before this script)
% folderPath_global: folder name that contains the ridge readout files extracted with a global ridge detection threshold.
% folderPath_adaptive: folder name that contains the ridge readout files extracted with an adaptive ridge detection threshold.
% AllFiles_global: all file names in 'folderPath_global'
% AllFiles_adaptive: all file names in 'folderPath_adaptive'
% is_bmal / is_per: files that contain variants of either "Bmal1" or "Per2" in their name.
% celllinenames_file: names of the cell lines being analysed, as written in the file names
% reporters: Bmal1-Luc or Per2-Luc

%define remaining parameters
values = {'period';'amplitude'};
xaxisvalues = {'Period (hours)';'Amplitude (a.u.)'};
weightsheets = {'weight_ridgelength_norm';'weight_ridgelength_unnorm'};
finer_subtype = {'TNBC-BL1';'TNBC-BL2';'TNBC-M';'LumA';'Epithelial';'Sarcoma'};
finer_sub = [4,3,4,3,1,1]; %for loop finer subtypes

disp('cwt5_plot_ridge_readout_distribution.m is executed')

load('sorted_periods_amplitudes.mat')

for v = 1:numel(values) %loop v values

    fig = figure('Visible','off');
    fig.Position = [1,1,1440,821];

    %calculate number of subplots per figure (including all cell models together)

    div = 3;
    n = ceil([numel(celllinenames_file)/div]);

    %load sorted continuous periods or amplitdes 
    sorted_data = sorted_periods_amplitudes{v};

    for c = 1:numel(celllinenames_file) % loop c celllines

        if c < 17
            aa = 2; %no Per2 data for knockout cell lines
        else
            aa = 1;
        end

        %load cell line data
        celllinedata = sorted_data{:,c};

        for a = 1:aa %loop a reporters

            %load data of specific reporter
            data = cell2mat(celllinedata(a,:));

            %load weights (calculated from fitting a sigmoidal function done in "weighted_boxplot_circadian_values.m"
            inputfile_weights = append('extracted_circadian_parameters_by_replicate_',reporters{a},'.xlsx');
            [weights] = table2array(readtable(inputfile_weights, 'Sheet',weightsheets{v}));
            w = rmmissing(transpose(weights(:,c+1)));

            %% plot histogram per cell line model, overlay reporters

            for f = 1:size(data,1)
                meanval = mean(data(f,:),'omitnan');
                meantimeseries(f,c) = meanval;
            end

            subplot(div,n,c);
            x = meantimeseries(:,c);

            %plot bars
            h = histfit(x,3); hold on
            h(1).FaceColor = reporter_colors{a};
            h(1).FaceAlpha = 0.3;
            h(2).Color = reporter_colors_dark{a};

            varstoclear = {'data','w','weighted_sum','count','x','x_values','y_values'};
            clear(varstoclear{:})

            meandata{a,c} = meantimeseries(:,c);
            % resultexcel = append('mean_by_cellline_',values{v},'_',reporters{a},'.xlsx');
            % outputsheet = celllinenames_file{c};
            % writematrix(meantimeseries(:,c),resultexcel,'sheet',outputsheet);

        end %loop a reporters

        hold off

        ax = gca;  % Get the current axes
        title(celllinenames_file{c},'FontSize',10,'FontName','Arial','Interpreter','none');
        set(ax,'XLimitMethod','padded','linewidth',1.5,'FontSize',12,'FontName','Arial');

        if c == 1
            legendentries = {'Bmal1','Bmal1 fit','Per2','Per2 fit'};
            legend(legendentries,'Location','northeast','FontSize',10,'FontName','Arial');
        end

    end %loop c celllines

    han1=axes(fig,'visible','off');

    han1.Title.Visible='on'; han1.XLabel.Visible='on'; han1.YLabel.Visible='on';
    xlabel(han1,xaxisvalues{v},'FontWeight','bold','FontSize',15);
    ylabel(han1,'Count','FontWeight','bold','FontSize',15);

    hold off
    %figurename1 = append('histogram_',values{v},'_avgpercellline');
    %saveas(fig, figurename1, 'svg');

    %% overlay results per finer subtypes

    for a2 = 1:2

        fig = figure('visible','off');
        fig.Position = [1,1,1440,821];

        s2 = 0;
        aa = 0;

        for k = 1:numel(finer_sub) %loop k subtypes

            s2 = s2+finer_sub(:,k);
            s1 = s2-finer_sub(:,k)+1;

            subplot(3,2,k);

            legend_entries = {}; % Collect cell line names for this subtype

            for s = s1:s2 %loop celllines per subtype

                aa = aa+1;
                %load data for specific combination
                x1 = rmmissing(cell2mat(meandata(a2,s)));

                %plot histogtam
                h2 = histfit(x1,3); hold on
                h2(1).FaceAlpha = 0;
                h2(1).EdgeColor = 'none';
                h2(1).HandleVisibility = 'off';
                barColor = h2(1).FaceColor;
                h2(2).Color = barColor;

                % Collect cell line name for the legend
                legend_entries{end+1} = celllinenames_file{s};

                varstoclear1 = {'x1'};
                clear(varstoclear1{:});

            end %loop celllines per subtype

            meanbysubtype_mat = mean(cell2mat(meandata(a2,[s1:s2])),2,'omitnan');
            meanbysubtype{a2,k} = meanbysubtype_mat;
            % resultexcel2 = append('mean_by_finer_subtype_',values{v},'_',reporters{a2},'.xlsx');
            % outputsheet = finer_subtype{k};
            % writematrix(meanbysubtype_mat,resultexcel2,'sheet',outputsheet);
            % clear meanbysubtype_mat;

            %add legend and modify appearance of the polarhistogram
            ax = gca;  % Get the current axes
            legend(legend_entries,'Location','eastoutside','FontSize',12,'FontName','Arial');
            set(ax,'XLimitMethod','padded','YLimitMethod','padded','linewidth',2,'FontSize',15,'FontName','Arial');
            title(finer_subtype{k},'FontSize',15,'FontName','Arial');

        end %loop k subtypes

        han1=axes(fig,'visible','off');

        han1.Title.Visible='on'; han1.XLabel.Visible='on'; han1.YLabel.Visible='on';
        title(han1,reporters{a2},'FontWeight','bold','FontSize',18);
        xlabel(han1,xaxisvalues{v},'FontWeight','bold','FontSize',18);
        ylabel(han1,'Count','FontWeight','bold','FontSize',18);

        hold off
        clear legend_entries

        % figurename2 = append('histogramfits_',values{v},'_',reporters{a2},'_persubtype');
        % saveas(fig, figurename2, 'svg');

        clear inputdata

    end %loop a2 reporters

    %% overlay results per subtype in single histogram

    fig = figure;%('Visible','off');
    fig.Position = [1,888,758,474];

    for k2 =1:numel(finer_sub) %loop k2 subtypes

        legend_entries = {}; % Collect cell line names for this subtype

        x2_1 = cell2mat(meanbysubtype(1,k2));
        x2_2 = cell2mat(meanbysubtype(2,k2));
        x2 = mean([x2_1,x2_2],2,'omitnan');
        % resultexcel2 = append('mean_by_finer_subtype_',values{v},'.xlsx');
        % outputsheet = finer_subtype{k2};
        % writematrix(x2,resultexcel2,'sheet',outputsheet);

        %plot histogtam
        h2 = histfit(x2,3); 
        hold on
        h2(1).FaceAlpha = 0;
        h2(1).EdgeColor = 'none';
        h2(1).HandleVisibility = 'off';

        varstoclear1 = {'x2'};
        clear(varstoclear1{:});

    end %loop k2 subtypes

    %add legend and modify appearance of the histogram
    ax = gca;  % Get the current axes
    legend(finer_subtype,'Location','eastoutside','FontSize',15,'FontName','Arial');
    set(ax,'XLimitMethod','padded','YLimitMethod','padded','linewidth',2,'FontSize',18,'FontName','Arial');
    if v == 1
        ax.XLim = [20.5,34];
        ax.XTick = [21 24 27 30 33];
    elseif v == 2
        ax.XScale = 'log';
        ax.XLim = [6,1200];
        ax.XTick = [12.5000 50 200];
    end
    xlabel(xaxisvalues{v},'FontWeight','bold','FontSize',18);
    ylabel('Count','FontWeight','bold','FontSize',18);

    hold off
    clear legend_entries

    %figurename3 = append('histogramfits_',values{v},'_combined_allsubtypes');
    %saveas(fig, figurename3, 'svg');

end %loop v values

clear input
disp('cwt5_plot_ridge_readout_distribution.m is completed')

end %function
