function [detrended_all, envelope_all, norm_all] = f1_sort_lumicycle_and_mra_data(lumicycle_datafile,celllinenames_file,reporters)

%Carolin Ector, 13.11.2023

% Function loads and sorts excel file where processed circadian time-series data is stored and multi-resolution values (output from python function)

%Clock-TNBC Manuscript Fig. 2

%input: stored in "amplitudedecay_autocorrelation_MRA_pipeline.mat"
% lumicycle_datafile: excel file where processed circadian time-series data is stored and multi-resolution values (output from python function)
% celllinenames_file: names of the cell lines being analysed, as written in the file names
% reporters: circadian gene names for the luciferase reporters

%output:
% detrended_all: 2x19 cell array containing detrended Luc-data for different cell lines (columns). 
    %Row 1 contains Bmal1-Luc data; Row 2 contains Per2-Luc data.
    %each cell stores a double matrix where rows are time points and columns are replicates for the corresponding reporter and cell line.
% envelope_all: see "detrended_all", but compiled of signal's amplitude (envelope) data
% norm_all: see "detrended_all", but compiled of amplitude (envelope)-normalized data

disp('f1_sort_lumicycle_and_mra_data.m is executed')

%define different names for Bmal1 and Per2 in the files
bmal = {'BLH','BMAL','Bmal','bmal'};
per = {'PER','Per','per'};
datasets = {'detrended';'envelope';'normalized';'mra'};
ii = numel(datasets);
cc = numel(celllinenames_file);

for i = 1:ii

    inputdata = datasets{i}

    % Read the Excel file into a table
    dataTable = readtable(lumicycle_datafile,'Sheet',inputdata);

    if i == ii %mra sample names are stored in first column
        sampleNames = dataTable{:, 1};
        dataTable(:,1) = [];
    else %all other sample names are stored in first row of each column
        sampleNames = dataTable.Properties.VariableNames;
    end

    if i == ii
        is_bmal = false(length(sampleNames),1);
        is_per = false(length(sampleNames),1);
        noise = nan(8,cc);
        ultradian = nan(8,cc);
        circ = nan(8,cc);
        infr = nan(8,cc);
    else
        is_bmal = false(1, length(sampleNames));
        is_per = false(1, length(sampleNames));
    end

    % check for reporter names
    for b = 1:numel(bmal)
        is_bmal = is_bmal | contains(sampleNames, bmal{b});
    end
    for p = 1:numel(per)
        is_per = is_per | contains(sampleNames, per{p});
    end

    % Iterate over each reporter
    for a = 1:numel(reporters)

        % Filter column names based on the current reporter
        if a == 1
            filtered_samples = sampleNames(is_bmal);
        elseif a == 2
            filtered_samples = sampleNames(is_per);
        end

        for c = 1:length(celllinenames_file)

            %disp(celllinenames_file{c})

            % Use 'contains' to find which filenames include the cell line name
            % matches = logical array where 1 indicates a match
            matches = contains(filtered_samples, celllinenames_file{c});
            % Extract only the filenames that match the cell line name
            matchingSamples = filtered_samples(matches);

            if i == ii %mra sample names are stored in first column
                indices = find(matches);
                noise(1:length(indices), c) = cellfun(@str2double, dataTable{indices, 1});
                ultradian(1:length(indices), c) = cellfun(@str2double, dataTable{indices, 2});
                circ(1:length(indices), c) = cellfun(@str2double, dataTable{indices, 3});
                infr(1:length(indices), c) = cellfun(@str2double, dataTable{indices, 4});
            else
                % Extract data for matching columns
                dataSubset{a,c} = dataTable{:, matchingSamples};
            end

            varstoclear1 = {'matchingSamples','indices','matches'};
            clear(varstoclear1{:})

        end %loop c celllines

        if i == ii

            % export sorted mra data to results sheet
            mra_all = {noise;ultradian;circ;infr};
            sheetnames = {'mra_noise';'mra_ultradian';'mra_circadian';'mra_infradian'};
            circadian_parameters_excel = append('extracted_circadian_parameters_by_replicate_',reporters{a},'.xlsx');
            replicate = {'1';'2';'3';'4';'5';'6';'7';'8'};
            t_row = cell2table(replicate);

            for r = 1:4 %loop r mra-datasets

                exportData = mra_all{r};
                le = size(exportData,2);
                t_exportData = array2table(exportData,'VariableNames',celllinenames_file(1:le));
                t_final = [t_row,t_exportData];

                % Export to Excel
                writetable(t_final,circadian_parameters_excel,'sheet',sheetnames{r});
                clear t_final
                clear exportData
            end

            clear mra_all

        end %if i == ii (mra)

    end %loop a reporter

    if i == 1
        detrended_all = dataSubset;
    elseif i == 2
        envelope_all = dataSubset;
    elseif i == 3
        norm_all = dataSubset;
    end

    varstoclear2 = {'dataSubset'};
    clear(varstoclear2{:})

end %loop i datasets

disp('f1_sort_lumicycle_and_mra_data.m is completed')

end %function