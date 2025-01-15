function DR_growth_plots(a,c,d,date,channel,cellline,drug,doses,colorarray,imaginginterval,xdata,ydata_smooth,stdev_smooth,ydata_normx0,stdev_normx0,ydata_normctrl,stdev_normctrl,additionaltext,outputfolder_plots)

%Carolin Ector, 11.09.2023

%Function plots growth curves of cell lines treated with various drugs and drug concentrations in dose response (DR) experiments:
% 1. plot raw growth curves
% 2. plot normalized growth curves (normalized to growth at timepoint 0)
% 3. plot relative growth curves (normalized to growth at timepoint 0, relative growth of each treated condition to respective solvent control)

%input: partially stored in '[date]_DR_workspace.mat'
% a: loop a channel being analyzed from live-imaging: 1 = cell number (red fluorescent channel), 2 = confluency (brightfield channel)
% c: loop c celllines
% date: date of the experiment being analyzed
% channel: channel being analyzed from live-imaging (see above)
% cellline: names of the cell lines being analysed
% drug: names of drugs used for drug treatments of the dose-response experiments
% doses: doses of each drug administered to the cells
% colorarray: colors for different drug doses
% imaginginterval: intervals of image acquisition during the live-imaging experiment
% xdata: elapsed time of the live-recordings (from beginning to end)
% ydata_smooth/stdev_smooth: smoothed growth curves + standard devation (loess smoothed) 
% ydata_normx0/stdev_normx0: normalized growth curves + standard devation (normalization to growth at initial time point of the live recordings)
% ydata_normctrl/stdev_normctrl: relative growth curves + standard devation (relative to growth of the solvent control) 

disp('DR_growth_plots.m is executed')

% Define remaining variables
yaxisnames = {'Cell Number';'Confluency'};
font = {'Helvetica Neue'};
ee = (numel(doses{1})+1); %loop doses, +1 = control
black = [0,0,0]; %control color
color  = [black;colorarray];

%% Figure 1: Smoothed growth curves as a function of time
figure('Visible','off');

time1 = transpose(xdata);

for j = 1:ee

    y1 = transpose(ydata_smooth(:,j));
    error1 = transpose(stdev_smooth(:,j));

    if j == 1
        linestyle = '--'; %dashed line for control
    else
        linestyle = '-'; %solid line for drug-treated conditions
    end

    h1 = patch([time1 fliplr(time1)], [(y1-error1)  (fliplr(y1+error1))], color(j,:), 'FaceAlpha',0.1, 'EdgeAlpha',0, 'HandleVisibility','off'); hold all
    plot(time1,y1,'LineWidth',2,'LineStyle',linestyle,'Color',color(j,:));

end

hold off

% Create labels and title
if a == 2
    ylabeltext1 = append(yaxisnames{a},' %');
else
    ylabeltext1 = yaxisnames{a};
end
xlabel('Time (h)','FontSize',22,'FontName',string(font));
ylabel(ylabeltext1,'FontSize',22,'FontName',string(font));

titletext = append(cellline{c});
title(titletext,'FontSize',22,'FontName',string(font),'FontWeight','normal');

% Set the remaining axes and box properties
ax = gca;
grid on;
xticks(0:24:144);

if a==2
    yticks(0:20:100); %ylim([0 110]);
else
    ax.YAxis.Exponent = 3;
end

set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
    'XGrid','off','Box','on','Color','none','FontSize',22,'FontName',string(font));

% add legend
dosesincluding0 = [0;doses{d}];
legendentries = num2str(dosesincluding0);
lgd = legend(legendentries);
set(lgd,'FontSize',15,'Orientation','vertical','Location','northwest','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');

legendtitle=append(drug{d},' (µM)');
title(lgd,legendtitle,'FontWeight','normal','FontName',string(font),'FontSize',16);

%save figure1
filetext1 = append(outputfolder_plots,'/',date,'_DR_',cellline{c},'_',drug{d},'_',channel{a},'_Fig1_RawData',additionaltext,'.svg');
saveas(h1, filetext1);

%% Figure 2: Growth curves normalized to time = 0 (time of treatment)
figure('Visible','off');

time2 = (0:imaginginterval:96);

for k = 1:ee

    y2 = transpose(ydata_normx0(:,k));
    error2 = transpose(stdev_normx0(:,k));

    if k == 1
        linestyle = '--'; %dashed line for control
    else
        linestyle = '-'; %solid line for drug-treated conditions
    end

    h2 = patch([time2 fliplr(time2)], [(y2-error2)  (fliplr(y2+error2))], color(k,:), 'FaceAlpha',0.1, 'EdgeAlpha',0, 'HandleVisibility','off'); hold all
    plot(time2,y2,'LineWidth',2,'LineStyle',linestyle,'Color',color(k,:));

end

hold off

% Create labels and title
ylabeltext2 = append('Normalized ',yaxisnames{a});
xlabel('Time post treatment (h)','FontSize',22,'FontName',string(font));
ylabel(ylabeltext2,'FontSize',22,'FontName',string(font));

title(titletext,'FontSize',22,'FontName',string(font),'FontWeight','normal');

ax = gca;
grid on;
xticks(0:24:144);
%ylim([0 Inf]);
set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
    'XGrid','off','Box','on','Color','none','FontSize',22,'FontName',string(font));

% add legend
lgd = legend(legendentries);
set(lgd,'FontSize',15,'Orientation','vertical','Location','northwest','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');

title(lgd,legendtitle,'FontWeight','normal','FontName',string(font),'FontSize',16);

%save figure 2
filetext2 = append(outputfolder_plots,'/',date,'_DR_',cellline{c},'_',drug{d},'_',channel{a},'_Fig2_Normx0',additionaltext,'.svg');
saveas(h2, filetext2)

%% Figure 3: Relative normalized growth to the control
figure('Visible','off');

for l = 2:ee

    y3 = transpose(ydata_normctrl(:,l));
    error3 = transpose(stdev_normctrl(:,l));

    h3 = patch([time2 fliplr(time2)], [(y3-error3)  (fliplr(y3+error3))], color(l,:), 'FaceAlpha',0.1, 'EdgeAlpha',0, 'HandleVisibility','off'); hold all
    plot(time2,y3,'LineWidth',2,'LineStyle',linestyle,'Color', color(l,:));

end

% Create labels and title
ylabeltext3 = append('Relative ',yaxisnames{a});
xlabel('Time post treatment (h)','FontSize',22,'FontName',string(font));
ylabel(ylabeltext3,'FontSize',22,'FontName',string(font));

title(titletext,'FontSize',22,'FontName',string(font),'FontWeight','normal');

% insert line at y=1 (=control)
yline(1, 'color', color(1,:), 'LineWidth',2,'LineStyle','--','HandleVisibility','off');

hold off

ax = gca;
grid on;
xticks(0:24:144);
yticks(0:0.25:1); ylim([0 1.2]);
set(ax,'XLimitMethod', 'padded','YLimitMethod', 'padded','linewidth',1.5,'YGrid','on', ...
    'XGrid','off','Box','on','Color','none','FontSize',22,'FontName',string(font));

% add legend
legendentries3 = num2str(doses{d});
lgd = legend(legendentries3);
set(lgd,'FontSize',15,'Orientation','vertical','Location','southwest','FontWeight','normal','EdgeColor','none','Color','#f5f5f5');

title(lgd,legendtitle,'FontWeight','normal','FontName',string(font),'FontSize',16);

%save figure 3
filetext3 = append(outputfolder_plots,'/',date,'_DR_',cellline{c},'_',drug{d},'_',channel{a},'_Fig3_RelCTRL',additionaltext,'.svg');
saveas(h3, filetext3)

end %function