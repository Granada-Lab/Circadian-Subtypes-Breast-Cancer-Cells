function cwt1_extract_ridge(folderPath_global,folderPath_adaptive,AllFiles_global,AllFiles_adaptive,is_bmal,is_per,celllinenames_file,reporters)

%Carolin Ector, 31.10.2023

%Function ...
% 1. cleans the ridge by identifying discontinuity in periods.
% 2. saves ridge length values.

%Clock-TNBC Manuscript Fig. 2

%input: stored in "cwt_analysis_pipeline.mat" & defined in cwt_analysis_pipeline.m (run in sequence before this script)
% folderPath_global: folder name that contains the ridge readout files extracted with a global ridge detection threshold.
% folderPath_adaptive: folder name that contains the ridge readout files extracted with an adaptive ridge detection threshold.
% AllFiles_global: all file names in 'folderPath_global'
% AllFiles_adaptive: all file names in 'folderPath_adaptive'
% is_bmal / is_per: files that contain variants of either "Bmal1" or "Per2" in their name.
% celllinenames_file: names of the cell lines being analysed, as written in the file names
% reporters: Bmal1-Luc or Per2-Luc
% thresholds: ridge detection thresholds
% global = 1/4 of median half-maximal wavelet power across all samples(fixed threshold = 118.4/4) --> amplitude, ridge length
% adaptive = 1/4 of maximal wavelet power per sample (adaptive threshold)--> period, phase difference

%define period cut-off for analysis (i.e. only ridges below this period are considered
maxperiod = 36;

disp('cwt1_extract_ridge.m is executed')

for un = 1:2 %loop 'un' unnormalized & normalized data

    if un == 1 %unnormalized
        folderPath = folderPath_global;
        AllFileNames = AllFiles_global;
        normalization = 'unnorm';
        disp('loop unnormalized data & global ridge threshold')
    elseif un == 2 %normalized
        folderPath = folderPath_adaptive;
        AllFileNames = AllFiles_adaptive;
        normalization = 'norm';
        disp('loop normalized data & adaptive ridge threshold')
    end

    for a = 1:numel(reporters) %loop a Luc-reporters

        %create empty array with all nan's
        ridge_length = nan(8, numel(celllinenames_file));

        % Filter and sort data for reporter group
        if a == 1
            filtered_files = AllFileNames(:,is_bmal);
        elseif a == 2
            filtered_files = AllFileNames(:,is_per);
        end

        for c = 1:length(celllinenames_file) %loop c celllines

            % Use 'contains' to find which filenames include the cell line name
            % matches = logical array where 1 indicates a match
            matches = contains(filtered_files, celllinenames_file{c});
            % Extract only the filenames that match the cell line name
            matchingFiles = filtered_files(matches);

            disp(celllinenames_file{c})

            for r = 1:numel(matchingFiles) %loop r replicates

                pathtomatchingfile = append(folderPath,'/',matchingFiles{r});
                [pyboatdata] = readtable(pathtomatchingfile);
                pyboatdata(:,1) = [];

                %disp(matchingFiles{r})

                %% clean the ridge by identifying discontinuity.

                %exclude values above max period
                indices = find(abs(pyboatdata.periods)>maxperiod);
                pyboatdata(indices,:) =[];

                %identify sudden jumps in periods (>=10% difference between two consecutive datapoints)
                ipt(1,:) = NaN;
                for g = 1:size(pyboatdata,1)-1
                    reldiff(g,:) = abs(1-(pyboatdata.periods(g+1,:)/pyboatdata.periods(g,:)));
                    if reldiff(g,:) >= 0.1 %~10%
                        row = length(ipt);
                        ipt(row+1,:) = g+1;
                    end
                end

                ipt(1,:) = [];
                TF = isempty(ipt);

                % Exclude data after or before a sudden jump and keep the data array with the overall higher power
                % or if the array with the overall higher power is < 3 times shorter than the other one, keep the longer one.

                if TF == 0 %sudden jump identified.

                    % Calculate the lengths of segments before and after each change point
                    lengthridges(:,1) = ipt(1)-1;
                    ridgestartindices(:,1) = 1;
                    ridgestartpoints(:,1) = pyboatdata.time(1);
                    ipt(end+1) = size(pyboatdata,1);

                    for j = 2:length(ipt)
                        lengthridges(:,j) = ipt(j) - ipt(j-1);
                        ridgestartindices(:,j) = ipt(j-1) + 1;
                        ridgestartpoints(:,j) = pyboatdata.time(ridgestartindices(j),:);
                    end

                    % Find ridges starting within the first 36 hours
                    withinFirstXHours = find(ridgestartpoints <= 36);

                    % Find the longest ridge among those
                    [~, longestRidgeIdx] = max(lengthridges(withinFirstXHours));

                    % Keep only the longest ridge starting within the first 36 hours
                    longestRidgeLength = lengthridges(withinFirstXHours(longestRidgeIdx));
                    longestRidgeStartIndex = ridgestartindices(withinFirstXHours(longestRidgeIdx));

                    % Extract the longest ridge from the original time-series
                    pyboatdata_clean = pyboatdata(longestRidgeStartIndex:longestRidgeStartIndex + longestRidgeLength - 2,:);

                    ipt(end) = [];

                    %%plot overlay of time-resolved periods from originaland "cleaned" ridge.
                    % plot(pyboatdata_clean.time,pyboatdata_clean.periods,'r','LineStyle','-','LineWidth',1);
                    % plot(pyboatdata.time(ipt,:),pyboatdata.periods(ipt,:),'o','Color','g','MarkerSize',12);
                    % legend('raw','cleaned','changepoint','Location','best');
                    % title(matchingFiles{r},'Interpreter','none');
                    % hold off
                    %
                    % figurename1 = append('ridge_cleaning_',normalization,'_',matchingFiles{r},'_v2.svg');
                    % destination = append(pathtofolder,'Figures/Fig1/continuous_wavelet_transform/ridge_threshold_determination/',resultfolder,normalization,'_detrended/',th);
                    % savelocation1 = append(destination,figurename1);
                    % saveas(fig1,savelocation1);
                    w = 1;

                else %no sudden jump identified.

                    % No change points found, keep original data
                    pyboatdata_clean = pyboatdata;
                    w = 0;
                end

                %save cleaned (continuous) ridge readout data.
                outputfolder = append(folderPath,'_continuous/');
                if ~exist(outputfolder, 'dir')
                    % Folder does not exist
                    mkdir(outputfolder);
                end
                outputfile_clean = append(outputfolder,matchingFiles{r});
                writetable(pyboatdata_clean,outputfile_clean);

                %% extract ridge length.
                lengthnonans = rmmissing(pyboatdata_clean);

                ridge_length(r,c) = size(lengthnonans,1)/6;

                if ridge_length(r,c) == 0
                    ridge_length(r,c) = NaN; %not evaluable sample
                end

                vars = {'pyboatdata','pyboatdata_clean','ipt'};
                clear(vars{:})

            end %loop r replicates

        end %loop c celllines

        %% save ridge lengths per replicate.
        outputsheet = append('ridgelength_',normalization);
        circadian_parameters_excel = append('extracted_circadian_parameters_by_replicate_',reporters{a},'.xlsx');
        replicate = {'1';'2';'3';'4';'5';'6';'7';'8'};
        t_row = cell2table(replicate);
        t_valuestosave = array2table(ridge_length,'VariableNames',celllinenames_file);
        t_final = [t_row,t_valuestosave];
        writetable(t_final,circadian_parameters_excel,'sheet',outputsheet);

        clear ridge_length

    end %loop a Luc-reporters

end %loop un unnormalized & normalized data

disp('cwt1_extract_ridge.m is completed')

end %function