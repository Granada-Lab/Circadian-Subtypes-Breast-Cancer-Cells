function cwt6_plot_weighted_phase_difference(celllinenames_file)

%Carolin Ector, 31.10.2023

%Function plots weighted Bmal1-Per2-phase differences per cell line and per subtype in polarhistograms

%Clock-TNBC Manuscript Fig. 2E, F

%input: stored in "cwt_analysis_pipeline.mat" & defined in cwt_analysis_pipeline.m (run in sequence before this script)
% phdiff_all: Bmal1-Per2 phase differences for each combination, cellline-by-cellline
% time_all: corresponding times of common Bmal1-Per2 phases for each combination, cellline-by-cellline
% celllinenames_file: names of the cell lines being analysed, as written in the file names

%load weights (calculated from common ridge in script "weighting_of_cwt_parameters.m")

load('extracted_phase_differences.mat')

disp('cwt6_plot_weighted_phase_difference.m is executed')

celllines = celllinenames_file(1:16,:); %exclude U2OS-KO cell lines (no Per2-Luc reporter)
inputfile_weights = 'extracted_circadian_parameters_by_replicate_Bmal1.xlsx';
weights = table2array(readtable(inputfile_weights, 'Sheet','weight_phdiff_commonridge'));

subtypes = {'TNBC-BL1';'TNBC-BL1';'TNBC-BL1';'TNBC-BL1';'TNBC-BL2';'TNBC-BL2';'TNBC-BL2';'TNBC-M';
    'TNBC-M';'TNBC-M';'TNBC-M';'LumA';'LumA';'LumA';'Epithelial';'Sarcoma';'Sarcoma';'Sarcoma';'Sarcoma'};
finer_subtype = {'TNBC-BL1';'TNBC-BL2';'TNBC-M';'LumA';'Epithelial';'Sarcoma'};
main_subtype = {'TNBC';'LumA';'Epithelial';'Sarcoma'};
sub = [4,3,4,3,1,1]; %for loop finer subtypes
sub2 = [11,3,1,1]; %for loop main subtypes
total_time = (0:0.16666667:137.7)'; %total recording time

%create figure for polarhistograms
fig1 = figure('Visible','off');
fig1.Position = [1,1,2560,1361];
hold all

for c = 1:numel(celllines)

    phdiff = phdiff_all{c};
    weights2 = rmmissing(transpose(weights(:,c+1)));
    time = time_all{c};

    row_indices = [];
    g = 0;

    for h = 1:numel(phdiff)
        test = isnan(phdiff{h});
        if test == 1
            g = g+1;
            row_indices(g,:) = h;
        end
        clear test
    end

    if ~isempty(row_indices)
        phdiff(row_indices,:) = [];
        time(row_indices,:) = [];
    end

    u = 0;
    for i = 1:numel(weights2)
        p = phdiff{i};
        t = time{i};
        p1 = isnan(p);
        if p1 ~= 1
            u = u+1;
            if t(1,:) ~= 0
                coltoadd = numel(0:0.1666666667:(t(1,:)-0.1));
                emptycols(1:coltoadd,:) = NaN;
                t = [emptycols;t];
                p = [emptycols;p];
                clear coltoadd
                clear emptycols
            end

            le_total = length(total_time);

            if length(t) < le_total
                t((end+1:le_total),:) = NaN;
                p((end+1:le_total),:) = NaN;
            end

            TF1 = isnan(t);

            x(:,u) = p;
            w((1:le_total),u) = weights2(:,u);
            w(TF1,u) = NaN;

            varstoclear6 = {'t','p'};
            clear(varstoclear6{:})
        else

        end
        clear p1
    end

    % resultexcel_x = append('continuous_phdiff_by_cellline_replicate.xlsx');
    % resultexcel_w = append('weights_continuous_phdiff_by_cellline_replicate.xlsx');
    % outputsheet = celllinenames_file{c};
    % writematrix(x,resultexcel_x,'sheet',outputsheet);
    % writematrix(w,resultexcel_w,'sheet',outputsheet);

    % Initialize variables to store the weighted sum and count
    weighted_sum = zeros(size(x, 1), 1);
    count = zeros(size(x, 1), 1);

    % Loop through each row in x
    for row = 1:size(x, 1)
        % Initialize variables to track the sum and count for the current row
        row_sum = 0;
        row_count = 0;

        % Loop through each time series
        for col = 1:size(x, 2)
            % Check if the value is not NaN
            if ~isnan(x(row, col))
                % Update the sum and count for the current row
                row_sum = row_sum + x(row, col) * w(row,col);
                row_count = row_count + 1;
            end
        end

        % Check if there are at least two non-NaN values in the current row
        if row_count >= 3
            % Update the weighted sum and count for the final result
            weighted_sum(row) = row_sum;
            count(row) = row_count;
        end
    end

    % Calculate the weighted mean time-series
    meanphasespercellline(:,c) = weighted_sum ./ count;

    TF = isnan(meanphasespercellline(:,c));
    x2 = total_time;
    x2(TF,:) = NaN;

    varstoclear3 = {'x','phdiff2','subtitles','times','commonridge1','commonridge2','w','weights2','row_indices'};
    clear(varstoclear3{:});

    %% phase difference over time per cell line - weighted mean of all replicates (weight = common ridge length)
    subplot(4,4,c);

    plot(x2,meanphasespercellline(:,c));

    ylim([-2*pi,2*pi]);
    yticks(-2*pi:pi:2*pi);
    yticklabels({'-2π','-π','0','π','2π'});
    xticks(0:24:144);
    xticklabels({'0','24','48','72','96','120','144'});
    title(celllines{c},'FontSize',10,'FontName','Arial','Interpreter','none');
    ax=gca;
    set(ax,'XLimitMethod','padded','linewidth',1.5,'FontSize',10,'FontName','Arial');

    clear x2
    clear y
    clear w

end %celllines

% resultexcel_x = append('continuous_phdiff_by_cellline_mean.xlsx');
% table = array2table(meanphasespercellline,'VariableNames',celllines);
% writetable(table,resultexcel_x,'Sheet','weighted_mean_phdiff');

han1=axes(fig1,'visible','off');

han1.Title.Visible='on'; han1.XLabel.Visible='on'; han1.YLabel.Visible='on';
xlabel(han1,'Time (h)','FontWeight','bold','FontSize',11);
ylabel(han1,'Phase difference (radians)','FontWeight','bold','FontSize',11);

hold off

%save figure
%figurename1 = append('phase_difference_over_time_wavgpercellline.svg');
%saveas(fig1, figurename1);

%% polarhistogram per subtype, TNBC subtypes separate -> Figure 2E

fig = figure('Visible','off');
fig.Position = [1,1,2560,1361];
s2 = 0;
aa = 0;

for k = 1:numel(sub) %loop k subtypes

    s2 = s2+sub(:,k);
    s1 = s2-sub(:,k)+1;

    subplot(3,2,k);

    legend_entries = {}; % Collect cell line names for this subtype

    for s = s1:s2 %loop celllines per subtype

        aa = aa+1;
        %load data for specific combination
        alpha2 = rmmissing(meanphasespercellline(:,s));

        %plot polarhistogtam
        polarhistogram(alpha2,6,'Normalization','probability','FaceAlpha',0.4); %v2
        hold all

        % Collect cell line name for the legend
        legend_entries{end+1} = celllines{s};

        varstoclear4 = {'alpha2'};
        clear(varstoclear4{:});

    end %loop celllines per subtype

    %add legend and modify appearance of the polarhistogram
    ax = gca;  % Get the current axes
    legend(legend_entries,'Location','eastoutside','FontSize',15,'FontName','Arial');
    set(ax,'linewidth',3,'FontSize',20,'FontName','Arial');
    ax.ThetaZeroLocation = 'top';  % Adjust as needed
    title(finer_subtype{k},'FontSize',20,'FontName','Arial');

    hold off
    clear legend_entries

end %loop k subtypes

%save figure
%figurename3 = append('polarhistogram_phase_difference_persubtype.svg');
%saveas(fig, figurename3);

%% polarhistogram per main subtype (TNBC combined)

fig = figure('Visible','off');
fig.Position = [1,1,2560,1361];
s2 = 0;
aa = 0;

for k2 =1:numel(sub2) %loop k subtypes

    s2 = s2+sub2(:,k2);
    s1 = s2-sub2(:,k2)+1;

    subplot(2,2,k2);

    legend_entries = {}; % Collect cell line names for this subtype

    for s = s1:s2 %loop celllines per subtype

        aa = aa+1;
        %load data for specific combination
        alpha2 = rmmissing(meanphasespercellline(:,s));

        %plot polarhistogtam
        %             polarhistogram(alpha2,10,'Normalization', 'pdf','FaceAlpha',0.4); %v1
        polarhistogram(alpha2,6,'Normalization', 'probability','FaceAlpha',0.4); %v2
        hold all

        % Collect cell line name for the legend
        legend_entries{end+1} = celllines{s};

        varstoclear4 = {'alpha2'};
        clear(varstoclear4{:});

    end %loop celllines per subtype

    %add legend and modify appearance of the polarhistogram
    ax = gca;  % Get the current axes
    legend(legend_entries,'Location','eastoutside','FontSize',15,'FontName','Arial');
    set(ax,'linewidth',3,'FontSize',20,'FontName','Arial');
    ax.ThetaZeroLocation = 'top';  % Adjust as needed
    title(main_subtype{k2},'FontSize',20,'FontName','Arial');

    hold off
    clear legend_entries

end %loop k subtypes

%save figure
%figurename4 = append('polarhistogram_phase_difference_permainsubtype.svg');
%saveas(fig, figurename4);

%% polarhistogram per subtype - all in one polarhistogram (merge cell lines per subtype) -> Figure 2F

fig = figure('Visible','off');
fig.Position = [1,1,2560,1361];
s2 = 0;

for k3 =1:numel(sub) %loop k subtypes

    s2 = s2+sub(:,k3);
    s1 = s2-sub(:,k3)+1;

    bb = 0;

    for s3 = s1:s2 %loop celllines per subtype

        bb = bb+1;
        %load data for specific combination
        alpha2(:,bb) = meanphasespercellline(:,s3);

    end %loop s3 celllines per subtype

    meanalpha = mean(alpha2,2,'omitnan');

    %         polarhistogram(alpha2,8,'Normalization','probability','FaceAlpha',0.4); %v1

    polarhistogram(meanalpha,6,'Normalization', 'probability','FaceAlpha',0.4); %v2
    hold all

    varstoclear7 = {'alpha2','meanalpha'};
    clear(varstoclear7{:});

end %loop k3 subtypes

%add legend and modify appearance of the polarhistogram
ax = gca;  % Get the current axes
legend(finer_subtype,'Location','eastoutside','FontSize',15,'FontName','Arial');
set(ax,'linewidth',3,'FontSize',20,'FontName','Arial');
ax.ThetaZeroLocation = 'top';  % Adjust as needed

hold off

%save figure
%figurename5 = append('polarhistogram_phase_difference_allsubtypescombined_norm_probability_6nbins.svg');
%saveas(fig, figurename5);

disp('cwt6_plot_weighted_phase_difference.m is completed')

end %function
