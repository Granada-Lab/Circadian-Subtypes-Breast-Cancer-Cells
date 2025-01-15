function cwt3_extract_phase_difference(folderPath_adaptive,bmal,per,celllinenames_file)

%Carolin Ector, 31.10.2023

%Function calculates and saves the following Bmal1-Per2-phase difference parameters, sample-by-sample and per cell line
    % common ridge lengths between the Bmal1 and Per2 signal traces, later to be used for weighting 
    % median phase difference calculated by circular statistics
    % coefficient of variation of the phase difference calculated by circular statistics

% Function utilizes scripts of the Circular Statistics Toolbox Version 1.21.0.0 by Philipp Berens
    % circ_dist.m; circ_dist2.m; circ_mean.m; circ_median.m; circ_r.m; circ_var.m.

%Clock-TNBC Manuscript Fig. 2E, F

%input: stored in "cwt_analysis_pipeline.mat" & defined in cwt_analysis_pipeline.m (run in sequence before this script)
% folderPath_adaptive: folder name that contains the ridge readout files extracted with an adaptive ridge detection threshold.
% bmal: names for bmal1-luc reporter {'BLH','BMAL','Bmal','bmal'}
% per: names for per-luc reporter {'PER','Per','per'}
% celllinenames_file: names of the cell lines being analysed, as written in the file names

disp('cwt3_extract_phase_difference.m is executed')

%define remaining parameters
celllines = celllinenames_file(1:16,:); %exclude U2OS-KO cell lines (no Per2-Luc reporter)
folderPath_continuous = append(folderPath_adaptive,'_continuous/'); %normalized
fileList = dir(folderPath_continuous);

% Convert all file names to lowercase for case-insensitive comparison
AllFileNames = ({fileList.name});
is_bmal = false(1, length(AllFileNames));
is_per = false(1, length(AllFileNames));

% check for reporter names
for b = 1:numel(bmal)
    is_bmal = is_bmal | contains(AllFileNames, bmal{b});
end
for p = 1:numel(per)
    is_per = is_per | contains(AllFileNames, per{p});
end

% Filter and sort data for reporter group
filtered_files_bmal = AllFileNames(:,is_bmal);
filtered_files_per = AllFileNames(:,is_per);

%total recording time
total_time = (0:0.16666667:137.7)';

%create empty array to store values in, Bmal1-Per2-combination by Bmal1-Per2-combination
arraylength = 40; %maximal possible Bmal1-Per2-combinations
commonridge = NaN(arraylength, numel(celllines)); %will be used for weighting of the specific combination
c_median = NaN(arraylength, numel(celllines)); %median phase difference, sample-by-sample
c_var = NaN(arraylength, numel(celllines)); %coefficient of variation of the phase difference, sample-by-sample

for c = 1:length(celllines) %loop c celllines

    %disp(celllines{c})

    % Use 'contains' to find which filenames include the cell line name
    % matches = logical array where 1 indicates a match
    matches_bmal = contains(filtered_files_bmal,celllines{c});
    matches_per = contains(filtered_files_per, celllines{c});

    % Extract only the filenames that match the cell line name
    matchingfiles{1} = filtered_files_bmal(matches_bmal);
    matchingfiles{2} = filtered_files_per(matches_per);

    for a = 1:2 %loop a Luc-reporter

        matchingfile = matchingfiles{a};

        for r = 1:numel(matchingfile) %loop r replicates

            matchingfilename = matchingfile{r};

            pathtomatchingfile = append(folderPath_continuous,'/',matchingfilename);
            [t_pyboatdata] = readtable(pathtomatchingfile);
            phases = t_pyboatdata.phase;

            %add NaN to missing time points
            mintime = min(t_pyboatdata.time);

            if mintime ~= 0
                coltoadd = numel(0:0.1666666667:mintime);
                emptycols(1:coltoadd,:) = NaN;
                phases = [emptycols;phases];
                clear coltoadd
                clear emptycols
            end

            %add NaN to end of time series, if shorter than 137.7 hours
            if length(phases) < length(total_time)
                phases(end+1:length(total_time)) = NaN;
            end

            %save phases for each sample per reporter
            if a == 1
                phases_bmal(:,r) = phases;
                rr_b = numel(matchingfile); %number of Bmal1 samples
            elseif a == 2
                phases_per(:,r) = phases;
                rr_p = numel(matchingfile); %number of Per2 samples
            end

            varstoclear1 = {'t_pyboatdata','phases'};
            clear(varstoclear1{:})

        end %loop r replicates

    end %loop a Luc-reporter

    %% loop through different Bmal1-Per2 combinations

    plotID = (1:1:(rr_p*rr_b)); %possible Bmal1-Per2 phase difference combinations

    rep_b = (1:1:rr_b); %number of Bmal1 samples
    rep_p = (1:1:rr_p); %number of Per2 samples

    for rp = 1:rr_p %loop Per2

        nr_per = num2str(rep_p(:,rp)); %sample number Per2

        for rb = 1:rr_b %loop Bmal1

            nr_bmal = num2str(rep_b(:,rb)); %sample number Bmal1

            %find common rows between phases from the Bmal1 and phases from the Per2 reporter
            ph_per_bmal = [phases_per(:,rp),phases_bmal(:,rb),total_time];
            rowsToExclude = sum(~isnan(ph_per_bmal), 2) < 3;
            ph_per_bmal(rowsToExclude, :) = [];
            times{rp,rb} = ph_per_bmal(:,3);

            %calculate phase difference for the current combination
            diff_rel = (ph_per_bmal(:,1)-ph_per_bmal(:,2));
            phdiff1=atan2(sin(diff_rel),cos(diff_rel));

            commonridge1{rp,rb} = length(phdiff1)/6;

            if length(phdiff1) < 288 %disregard phase differences which are shorter than 288 timepoints (1 timepoint = 10 min -> 288 timepoints = 48 hours)
                phdiff2{rp,rb} = NaN;
                times{rp,rb} = NaN;
            else %save phase differences that are from common ridges longer than 48 hours
                phdiff2{rp,rb} = unwrap(phdiff1);
            end

            varstoclear2 = {'ph_per_bmal','rowsToExclude','diff_rel','nr_bmal'};
            clear(varstoclear2{:});
        end

        clear nr_per

    end

     %% calculate common ridge length, median phase difference, and coefficient of variation of the phase difference for each Bmal1-Per2 combination

    phdiff = cat(1, phdiff2(:));
    time = cat(1, times(:));

    % create empty array to store common ridge length for eachcombination
    commonridge2 = cell2mat(cat(1, commonridge1(:)));

    % assure equal array lengths across cell models
    if length(commonridge2) < 40 %40 is the maximum number of combinations
        commonridge2(end+1:40,:) = NaN;
    end

    % this will be
    commonridge(:,c) = commonridge2;

    varstoclear3 = {'phdiff2','subtitles','times','commonridge1','commonridge2'};
    clear(varstoclear3{:});

    for e = 1:numel(plotID) %loop e combination

        %load data for specific combination
        alpha = cell2mat(phdiff(e,:));
        TF = isnan(alpha);

        if TF == 1
            %parameter calculation not possible.
        else
            %calculate phase difference parameters:
            c_median(e,c) = circ_median(alpha);
            c_var(e,c) = circ_var(alpha);

            varstoclear4 = {'alpha','TF'};
            clear(varstoclear4{:});
        end

    end %loop e combination

    phdiff_all{c} = phdiff;
    time_all{c} = time;

    varstoclear5 = {'phdiff'};
    clear(varstoclear5{:})

end %loop c celllines

save('extracted_phase_differences.mat', 'phdiff_all', 'time_all');

%% save calculated parameters.
resultexcel_phdiff = append('extracted_circadian_parameters_by_replicate_Bmal1.xlsx');
allfinalvaluesbyrep = {c_median;c_var;commonridge};
outputsheets = {'phdiff_median';'phdiff_coeffvar';'phdiff_commonridge'};

for n = 1:numel(allfinalvaluesbyrep) %loop n final values ph. diff

    valuestosave = allfinalvaluesbyrep{n};
    outputsheet = outputsheets{n};

    nrofcombinations = (1:1:arraylength)';
    t_row = array2table(nrofcombinations);

    t_valuestosave = array2table(valuestosave,'VariableNames',celllines);
    t_final = [t_row,t_valuestosave];
    writetable(t_final,resultexcel_phdiff,'sheet',outputsheet);

end %loop n final values ph. diff

disp('cwt3_extract_phase_difference.m is completed')

end %function