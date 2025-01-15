function f5_scatterplot_autocorrelation_mra(reporters,celllinenames_file,subtype_colors,markers)

%Carolin Ector, 14.11.2023


%Function creates 2D-scatterplots to visualize the distribution of circadian parameters of individual samples

%Clock-TNBC Manuscript Fig. 2A and B

%input: stored in "amplitudedecay_autocorrelation_MRA_pipeline.mat" & calculated in preceding functions
% celllinenames_file: names of the cell lines being analysed, as written in the file names
% markers:
% subtype_colors:
% reporters: circadian gene names for the luciferase reporters
% axOpt: setting for appearance of axes of a plot

disp('f5_scatterplot_autocorrelation_mra.m is executed')

%define excel sheet names of parameters you would like to plot together
inputsheet1 = {'autocorrelation_lag';'mra_circadian'};
inputsheet2 = {'autocorrelation_peak';'mra_noise'};

for s = 1:numel(inputsheet1) %loop s inputsheet-combination

    fig = figure;
    fig.Position = [827,469,790,584];
    hold on

    for a = 1:numel(reporters) %loop a Luc-reporters

        %load data
        circadian_parameters_excel = append('extracted_circadian_parameters_by_replicate_',reporters{a},'.xlsx');
        [data1] = readtable(circadian_parameters_excel,'sheet',inputsheet1{s});
        [data2] = readtable(circadian_parameters_excel,'sheet',inputsheet2{s});

        %exclude first column of the table, as it contains the replicate number
        x1 = table2array(data1(:,2:end));
        y1 = table2array(data2(:,2:end));

        %exclude negative indices in autocorrelation values -> convert to NaN value
        if s == 1
            negindices = y1<0;
            for b1 = 1:size(negindices,1)
                for b2 = 1:size(negindices,2)
                    test = negindices(b1,b2);
                    if test == 1
                        x1(b1,b2) = NaN;
                        y1(b1,b2) = NaN;
                    end
                end
            end
        end


        %find common rows between datasets (as sometimes a value might be missing, while the other sheet has a value for that sample).
        cols = size(x1,2);

        for c = 1:cols %loop c celllines

            % Removing rows with missing values in either x or y
            validRows = ~isnan(x1(:,c)) & ~isnan(y1(:,c));
            x2 = x1(validRows,c);
            y2 = y1(validRows,c);

            %plot replicate-averaged data for mra and individual samples for autocorrelation
            if s ==2
                x2 = mean(x1(:,c),'omitnan');
                y2 = mean(y1(:,c),'omitnan');
                x_err = std(x1(:,c),[],'omitnan'); %stdev
                err_y = std(y1(:,c),[],'omitnan'); %stdev

                %plot errorbars first so that they appear behind the scatterplot markers
                eb(1) = errorbar(x2,y2,x_err, 'horizontal', 'LineStyle', 'none','HandleVisibility','off');
                eb(2) = errorbar(x2,y2,err_y, 'vertical', 'LineStyle', 'none','HandleVisibility','off');
                set(eb, 'LineWidth', 0.75, 'Color','black');
            end

            if a == 1 %Bmal1
                scatter(x2, y2, 150, markers{c,:},'MarkerEdgeColor',subtype_colors(c,:),'MarkerFaceColor',subtype_colors(c,:));
            elseif a == 2 %Per2
                scatter(x2, y2, 150, markers{c,:},'LineWidth',1.5,'MarkerEdgeColor',subtype_colors(c,:),'MarkerFaceColor','w');
            end
            varstoclear1 = {'x2';'y2'};
            clear(varstoclear1{:})

        end %loop c celllines

    end %loop reporter

    varstoclear1 = {'x1';'y1'};
    clear(varstoclear1{:})

    legend(celllinenames_file,'location','northeastoutside');

    %set graph appearance options
    xlabel(inputsheet1{s},'FontSize',20,'FontName','Helvetica Neue','Interpreter','none');
    ylabel(inputsheet2{s},'FontSize',20,'FontName','Helvetica Neue','Interpreter','none');

    %add cutoff line on log-lin scale to mra scatterplot
    if s == 2
        Y=(0:0.0001:1);
        Y2 = flip(Y);
        plot(Y,Y2,'Handlevisibility','off','Color',[0.3 0.3 0.3],'LineWidth',1);
        hold off
        set(gca,'YScale','log')
        yticklabels({'0.0001','0.001','0.01','0.1','1'});
        xlim([-0.03 1.03]);
        xticks(0:0.2:1);
    end

    hold off

    ax = gca;
    grid on
    set(ax,'linewidth',2,'YGrid','on','XGrid','on','Box','on','Color','none','FontSize',22);

    %figurename = append('scatter_',inputsheet1{s},'_vs_',inputsheet2{s},'.svg');
    %saveas(fig, figurename);

    clear data

end %loop s inputsheet-combination

disp('f5_scatterplot_autocorrelation_mra.m is completed')

end %function
