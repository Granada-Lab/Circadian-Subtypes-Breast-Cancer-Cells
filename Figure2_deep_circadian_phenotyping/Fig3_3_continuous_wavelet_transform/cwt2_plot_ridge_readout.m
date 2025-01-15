function cwt2_plot_ridge_readout(folderPath_global,folderPath_adaptive,AllFiles_global,AllFiles_adaptive,is_bmal,is_per,celllinenames_file,reporters)

%Carolin Ector, 02.11.2023

%Function plots the time-resolved wavelet spectrum parameters extracted from continuous ridges, averaged by cell line
%Function saves the mean, standard deviation, coefficient of variation, and median of each parameter across replicates by cell line to excel sheets

%input: stored in "cwt_analysis_pipeline.mat" & defined in cwt_analysis_pipeline.m (run in sequence before this script)
% folderPath_global: folder name that contains the ridge readout files extracted with a global ridge detection threshold.
% folderPath_adaptive: folder name that contains the ridge readout files extracted with an adaptive ridge detection threshold.
% AllFiles_global: all file names in 'folderPath_global'
% AllFiles_adaptive: all file names in 'folderPath_adaptive'
% is_bmal / is_per: files that contain variants of either "Bmal1" or "Per2" in their name.
% celllinenames_file: names of the cell lines being analysed, as written in the file names
% reporters: Bmal1-Luc or Per2-Luc

%define remaining parameters
cwtvalues = {'period';'amplitude'};
columns = [2,4]; %columns in ridge readout files for period, amplitude, and power, respectively.
replicate = {'1';'2';'3';'4';'5';'6';'7';'8'};

disp('cwt2_plot_ridge_readout.m is executed')

for v = 1:numel(cwtvalues) %loop v cwt-values

    col = columns(:,v);

    if v == 1
        folderPath_continuous = append(folderPath_adaptive,'_continuous/'); %normalized
        AllFileNames = AllFiles_adaptive;
    else
        folderPath_continuous = append(folderPath_global,'_continuous/'); %unnormalized
        AllFileNames = AllFiles_global;
    end

    %create empty arrays to store median or coefficient of variation of the time series, per replicate and cell model
    mediandata_bmal = nan(8,numel(celllinenames_file));
    mediandata_per = nan(8,numel(celllinenames_file));
    coeffvar_bmal = nan(8,numel(celllinenames_file));
    coeffvar_per = nan(8,numel(celllinenames_file));

    for c = 1:length(celllinenames_file) %loop c celllines

        %disp(celllinenames_file{c})

        for a = 1:numel(reporters) %loop a Luc-reporters

            % Filter and sort data for reporter group
            if a == 1
                filtered_files = AllFileNames(:,is_bmal);
            elseif a == 2
                filtered_files = AllFileNames(:,is_per);
            end

            % Use 'contains' to find which filenames include the cell line name
            % matches = logical array where 1 indicates a match
            matches = contains(filtered_files, celllinenames_file{c});

            % Extract only the filenames that match the cell line name
            matchingFiles = filtered_files(matches);

            filenames_all{a} = matchingFiles;

            %% extract and process data.

            for r = 1:numel(matchingFiles) %loop r replicates

                %load data for respective sample
                pathtomatchingfile = append(folderPath_continuous,'/',matchingFiles{r});
                [t_pyboatdata] = readtable(pathtomatchingfile);
                pyboatdata = table2array(t_pyboatdata(:,col));
                time = table2array(t_pyboatdata(:,1));

                %add NaN to missing time points
                mintime = min(time);

                %total recording time
                total_time = (0:0.16666667:137.7)';

                lengthdata = length(pyboatdata);
                totallength = length(total_time);

                if mintime ~= 0
                    coltoadd = numel(0:0.1666666667:mintime);
                    emptycols(1:coltoadd,:) = NaN;
                    pyboatdata = [emptycols;pyboatdata];
                    clear coltoadd
                    clear emptycols
                end

                if lengthdata < totallength
                    pyboatdata(end+1:totallength) = NaN;
                end

                if lengthdata < 288
                    pyboatdata(1:end,:) = NaN;
                end

                %save median period or amplitude across all timepoints per sample
                mediandata = median(pyboatdata,'all','omitnan');
                coeffvar = (std(pyboatdata,[],'all','omitnan'))/(mean(pyboatdata,'all','omitnan'));

                if a == 1
                    mediandata_bmal(r,c) = mediandata;
                    coeffvar_bmal(r,c) = coeffvar;
                else
                    mediandata_per(r,c) = mediandata;
                    coeffvar_per(r,c) = coeffvar;
                end

                data{a,r} = pyboatdata;

                vars = {'mediandata','coeffvar','pyboatdata'};
                clear(vars{:})

            end %loop r replicates

            vars = {'matchingFiles','filtered_files','matches'};
            clear(vars{:})

        end %loop a Luc-reporters

        sorted_data{c} = data;

        %save median and coeffvar values by cellline in cell array, replicate-by-replicate
        allfinalvaluesbyrep{1} = {mediandata_bmal;coeffvar_bmal};
        allfinalvaluesbyrep{2} = {mediandata_per;coeffvar_per};

        if c > 16 %no Per2-Luc data for U2OS-KO cell lines.
            max_aa = 1;
        else
            max_aa = numel(reporters);
        end

        clear data

    end %loop c celllines

    sorted_periods_amplitudes{v} = sorted_data;

    %save median values and coeffvar by cellline in cell array, replicate-by-replicate
    for aaa = 1:numel(reporters) %loop aa Luc-reportersv

        %define excel file & output sheet
        circadian_parameters_excel = append('extracted_circadian_parameters_by_replicate_',reporters{aaa},'.xlsx');
        metrics = {'_median'; '_coeffvar'};
        t_row = cell2table(replicate);

        %save median and coeffvar values
        for s = 1:2

            outputsheet = append(cwtvalues{v},metrics{s});
            valuestosave1 = allfinalvaluesbyrep{aaa}; %load values for respective reporter
            valuestosave2 = valuestosave1{s}; %load median or coeffvar

            le = size(valuestosave2,2);

            t_valuestosave = array2table(valuestosave2,'VariableNames',celllinenames_file(1:le));
            t_final2 = [t_row,t_valuestosave];
            writetable(t_final2,circadian_parameters_excel,'sheet',outputsheet);

            clear t_final2
            clear valuestosave

        end %loop s

    end %loop aaa

    vars = {'allfinalvaluesbyrep','mediandata_bmal','coeffvar_bmal','mediandata_per','coeffvar_per',};
    clear(vars{:})

end %loop v cwt-values

save('sorted_periods_amplitudes.mat', 'sorted_periods_amplitudes');
disp('cwt2_plot_ridge_readout.m is completed')

end %function


