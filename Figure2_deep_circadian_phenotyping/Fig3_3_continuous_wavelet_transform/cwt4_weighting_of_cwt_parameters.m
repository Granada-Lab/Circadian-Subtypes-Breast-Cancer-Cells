function cwt4_weighting_of_cwt_parameters(celllinenames_file,reporter_colors)

%Carolin Ector, 05.12.2023

%Function calculates the relative weight for different circadian clock parameters and creates weighted boxplots to sort cell lines by the respective parameter
    %Continuous wavelet parameters are weighed based on their underlying ridgelength
        %The longer the ridge length from which the parameters is derived from, the higher the weight for the statistics (eg. median)
    %Relative weights are estimated by fitting a sigmoidal function to all input weights (ridglength or Rsq)

%Clock-TNBC Manuscript Fig. 2, Appendix Figure S1B

%input: stored in "cwt_analysis_pipeline.mat" & defined in cwt_analysis_pipeline.m (run in sequence before this script)
% pathtofolder: path to manuscript folder
% celllines_text: names of the cell lines being analysed, optimized for being used as x-labels
% reporter_colors: colors used in the graphs for the two circadian clock luciferase reporters, or the merged version (Bmal1+Per2)

%define remaining parameters
metric = {'period_median';'period_coeffvar';'amplitude_median';'amplitude_coeffvar';'phdiff_median';'phdiff_coeffvar'};
yaxisnames = {'Period (h)';'CV period';'Amplitude (a.u.)';'CV amplitude';'Phase difference';'CV phase difference'};
colors = hsv(numel(celllinenames_file));

disp('cwt4_weighting_of_cwt_parameters.m is executed')

for m = 1:numel(metric) %loop m metric

    %load data
    if m == 5 || m == 6
        file_bmal = 'extracted_circadian_parameters_by_replicate_Bmal1.xlsx';
        file_per = 'extracted_circadian_parameters_by_replicate_Bmal1.xlsx';
    else
        file_bmal = 'extracted_circadian_parameters_by_replicate_Bmal1.xlsx';
        file_per = 'extracted_circadian_parameters_by_replicate_Per2.xlsx';
    end

    [t_bmal] = readtable(file_bmal,'sheet',metric{m});
    [t_per] = readtable(file_per,'sheet',metric{m});

    bmal = table2array(t_bmal(:,2:end));
    per = table2array(t_per(:,2:end));

    %no Per2 data for U2OS-KO cell lines
    if size(bmal,2) == 19
        per(1:6,17:19) = NaN; %add NaN values to make per2 variable compatible with bmal1
    end

    %% Load data (calculated in previous scripts)

    mridge = (1:1:4);

    %specify which excel sheet should be loaded
    if ismember(m, mridge)
        if m <= 2
            normalization = 'norm';
        elseif m > 2
            normalization = 'unnorm';
        end
        weightsheet = append('ridgelength_',normalization);
    elseif m == 5 || m == 6
        weightsheet = append('phdiff_commonridge');
        normalization = 'phdiff_commonridge';
    end

    [t_weights_bmal] = readtable(file_bmal,'sheet',weightsheet);
    [t_weights_per] = readtable(file_per,'sheet',weightsheet);

    weights_bmal = table2array(t_weights_bmal(:,2:end));
    weights_per = table2array(t_weights_per(:,2:end));

    if size(weights_bmal,2) == 19
        weights_per(1:6,17:19) = NaN;
    end

    if m ~= 7
        for bb = 1:2
            if bb == 1
                input1 = bmal;
                input2 = weights_bmal;
            else
                input1 = per;
                input2 = weights_per;
            end
            for ii = 1:size(input1,1)
                for jj = 1:size(input1,2)
                    if isnan(input1(ii,jj))
                        input2(ii,jj) = NaN;
                    end
                end
            end
            if bb == 1
                weights_bmal = input2;
            else
                weights_per = input2;
            end
        end

        %% Weighting of parameters based on ridge lengths 

        maxweight = max([max(weights_bmal(:)),max(weights_per(:))]);
        minweight = min([min(weights_bmal(:)),min(weights_per(:))]);

        % Parameters for the sigmoid function
        x0 = mean([maxweight, minweight], 'all'); % Midpoint of the data
        hill = 0.1; % Control the steepness of the logistic curve

        % Calculate weights using the sigmoid function
        w_bmal = 1 ./ (1 + exp(-hill * (weights_bmal - x0)));
        w_per = 1 ./ (1 + exp(-hill * (weights_per - x0)));

        %save relative weights to excel
        weightstosave = {w_bmal,w_per};
        outputfile = {file_bmal,file_per};
        outputsheet = append('weight_',weightsheet);

        if m == 5 || m == 6
            nn = 1;
        else
            nn = 2;
        end

        for n = 1:nn
            valuestosave = weightstosave{n};
            varnames = celllinenames_file(1:size(valuestosave,2));
            t_row = array2table([1:1:size(valuestosave,1)]');
            t_valuestosave = array2table(valuestosave,'VariableNames',varnames);
            t_final = [t_row,t_valuestosave];
            writetable(t_final,outputfile{n},'sheet',outputsheet);

            clear t_final
            clear t_row
            clear varnames
        end

        %% create figure of sigmoidal fit (Appendix Figure S1B)

        fig1 = figure('Visible','off');

        hold all

        for h = 1:size(bmal,2)
            scatter(weights_bmal(:,h),w_bmal(:,h),60,colors(h,:),'filled','o');
            scatter(weights_per(:,h),w_per(:,h),60,colors(h,:),'filled','v');
        end

        legend('Bmal1','Per2','location','northeast');
        hold off

        xlabel('Ridge length (h)','FontSize',20,'FontName','Helvetica Neue');
        ylabel('Relative weight','FontSize',20,'FontName','Helvetica Neue');

        x0_str = sprintf('%.2f',x0);
        titletext = append('1 ./ (1 + exp(-0.1 * (ridgelength - ',x0_str,')');
        title(titletext,'FontSize',20,'FontName','Helvetica Neue','Interpreter','none');

        ax = gca;
        grid on;
        xticks(0:24:144);

        set(ax,'XLimitMethod','padded','YLimitMethod','padded','linewidth',1.5,'YGrid','on', ...
            'XGrid','off','Box','on','Color','none','FontSize',18,'FontName','Helvetica Neue');

        %save figure
        % filetext = append('sigmoidalfit_ridge_',normalization,'_weights.svg');
        % saveas(fig1, filetext);

    else 

        w_bmal = weights_bmal;
        w_per = weights_per;

    end 
    %% Plot weighted parameters in a boxplot -> Appendix Figure S1C

    %merge Bmal1 and Per2 data for all cell lines
    if m == 5 || m == 6
        datacombined = bmal;
        weightscombined = w_bmal;
    else
        datacombined = [bmal;per];
        weightscombined = [w_bmal;w_per];
    end

    fig = figure('Visible','off');
    fig.Position = [420,285,525,425];

    % Replicate each data point by its weight
    inputdata = datacombined;
    inputweight = weightscombined;
    replicatedData = cell(size(datacombined, 2), 1); % Preallocate a cell array to hold replicated data for each group
    for i = 1:size(datacombined, 2) % Loop over each category
        currentData = datacombined(:, i);
        currentWeights = inputweight(:, i);
        groupData = [];
        for j = 1:length(currentData) % Replicate each data point by its weight
            replicationFactor = max(1, round(currentWeights(j) * 100)); % Scale the weights to a reasonable replication factor
            groupData = [groupData; repmat(currentData(j), replicationFactor, 1)]; % Concatenate the replicated data
        end
        replicatedData{i} = groupData; % Store the replicated data for this group
        wMed(:,i) = median(groupData,'omitnan');
    end

    for c = 1:size(inputdata,2)
        n_samples(:,c) = numel(rmmissing(inputdata(:,c)));
        xlabeltext(c,:) = append(celllinenames_file(c,:),' (',num2str(n_samples(:,c)),')');
    end

    %sort cell lines by their median lag or peak value (ascending)
    if m == 2 || m == 4 || m == 6
        order = 'ascend';
    else
        order = 'descend';
    end

    [~, sortIdx] = sort(wMed,order);
    bmalSorted = bmal(:,sortIdx);
    w_bmalSorted = w_bmal(:,sortIdx);
    perSorted = per(:,sortIdx);
    w_perSorted = w_per(:,sortIdx);
    celllinesSorted = xlabeltext(sortIdx,:);
    if  m > 4 && m < 7
        dataSorted = bmalSorted;
        weightsSorted = w_bmalSorted;
    else
        dataSorted = [bmalSorted;perSorted];
        weightsSorted = [w_bmalSorted;w_perSorted];
    end

    % Replicate and sort data according to the new order
    replicatedDataSorted = cell(1, length(dataSorted(1, :)));
    for i = 1:length(dataSorted(1, :))
        currentData = dataSorted(:, i);
        currentWeights = weightsSorted(:, i);
        groupData = [];
        for j = 1:length(currentData)
            replicationFactor = max(1, round(currentWeights(j) * 100));
            groupData = [groupData; repmat(currentData(j), replicationFactor, 1)];
        end
        replicatedDataSorted{i} = groupData;
    end

    % Combine all sorted and replicated data into a single vector for the box plot
    allData = vertcat(replicatedDataSorted{:});
    groupVector = [];
    for i = 1:length(replicatedDataSorted)
        groupVector = [groupVector; repmat(i, length(replicatedDataSorted{i}), 1)];
    end

    % Create the box plot
    boxplot(allData, groupVector, 'Labels', celllinesSorted);
    boxes = findobj(gca, 'Tag', 'Box');
    for j=1:length(boxes)
        set(boxes(j), 'Color', [1 1 1],'LineWidth',0.8); % Set edge color to black
        patch(get(boxes(j), 'XData'), get(boxes(j), 'YData'), [0.4 0.4 0.4], 'FaceAlpha', 0.2); % Set face color to grey
    end

    hold all

    % Get the number of categories and data points per category
    numCategories = numel(celllinesSorted);

    % Define a base size for unweighted points
    %     baseSize = 50;

    % Loop through each category to overlay data points
    for i = 1:numCategories

        % Extract the data and weights for the current category
        data_bmal = rmmissing(bmalSorted(:, i));
        wbmal = rmmissing(w_bmalSorted(:, i));
        data_per = rmmissing(perSorted(:, i));
        wper = rmmissing(w_perSorted(:, i));

        for pp = 1:2
            if pp == 1
                inputweight = wbmal;
            else
                inputweight = wper;
            end
            sizefactor(:,1:length(inputweight)) = NaN;
            for ff = 1:length(inputweight)
                if inputweight(ff) <= 0.33
                    sizefactor(ff) = 33;
                elseif inputweight(ff) < 0.67 && inputweight(ff) > 0.33
                    sizefactor(ff) = 67;
                elseif inputweight(ff) >= 0.67
                    sizefactor(ff) = 100;
                end
            end
            if pp == 1
                sizebmal = sizefactor;
            elseif pp == 2
                sizeper = sizefactor;
            end
            clear inputweight
            clear sizefactor
        end

        % Calculate x-axis positions for the scatter plot
        xPositions_b = i * ones(size(data_bmal));
        xPositions_p = i * ones(size(data_per));

        % Calculate sizes for the scatter plot
        %         pointSizes_b = baseSize * wbmal;
        %         pointSizes_p = baseSize * wper;

        % Overlay the weighted data points on the box plot
        if m > 4 && m < 7
            scatter(xPositions_b, data_bmal, sizebmal,reporter_colors{3},'filled','MarkerFaceAlpha',0.8','jitter','on','jitterAmount',0.15);
        else
            scatter(xPositions_b, data_bmal, sizebmal,reporter_colors{1},'filled','MarkerFaceAlpha',0.8','jitter','on','jitterAmount',0.15);
            scatter(xPositions_p, data_per, sizeper, reporter_colors{2},'filled','MarkerFaceAlpha',0.8','jitter','on','jitterAmount',0.15);
        end
    end

    hold off

    % Aesthetics
    ax = gca;
    set(ax, 'XTick', 1:numCategories, 'XTickLabel', celllinesSorted, 'XTickLabelRotation', 45);
    ylabel(yaxisnames{m});
    xlim([0, numCategories + 1]);
    hold off;

    box on;
    grid on
    set(ax,'YLimitMethod','padded','LineWidth',0.9,'FontName','Helvetica Neue','FontSize',18,'XMinorGrid','off','YMinorGrid','off');

    % figurename = append('weighted_boxplot_ranking_',metric{m},'.svg');
    % saveas(fig,figurename);
    % clear figurename
    clear sortIdx  
    clear wMed

end %loop m metric

disp('cwt4_weighting_of_cwt_parameters.m is completed')

end %function