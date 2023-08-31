%% Full mEPSC Analysis 

% READ!!!!!!!

%must adhere strictly to the file organization as given in the sample_data
%naming strucutre must also be consistent: default of clampex:
%YYYY_MM_DD_xxxx_filtered.abf
%YYYY = year of recording
%MM = month of recording
%DD = day of recording
%xxxx = number recotding that day (e.g. 3rd recording of a day is 0003,
%fifteenth recording is 0015, etc)
%filtered: ALL of my recordings I filter in clampfit, specifically I use a
%boxcar with 7 smoothing points. this is important because all the
%accurately marking events are based on this filtering. please ask
%questions if you cannot figure this out. in clampex, open file and then
%hit analyze->filter->select boxcar and select 7 smoothing points and
%ensure you filter the entire trace

%As presently constructed, this code will analyze the sample data for
%parvalbumin positive interneurons - also available is sample data for
%pyramidal neurons. To run this, you will need to switch the paths to grab
%the PYR files (organized the same as the PV-INs, see 'sample_data' folder)
%AND switch the analysis function to find_EPSCs_PYR.

%for each cell, I record 2-5 minutes of gap-free recordings, 1 minute each

%LAF, Boston University
%Work in progress, Last update 2023/07
%% Turn off warnings
warning('off','all')
%% Load in Experimemt Folder - can run line for PV or PYR for sample data depending on which cell type you are interested in
all_contents = dir('D:\Users\Luke_Fournier\Desktop\mEPSC_analysis_suite\sample_data\PV-IN\*Condition*');
%OR
all_contents = dir('D:\Users\Luke_Fournier\Desktop\mEPSC_analysis_suite\sample_data\PYR\*Condition*');
%all_contents should be a n x 1 struct where n = the number of conditions
%you have in a given experiment (e.g. control vs KO)
%% Add function path
addpath('D:\Users\Luke_Fournier\Desktop\mEPSC_analysis_suite\MATLAB\mEPSC_detection_functions');
%% Main Body of Code - run entire section at once
for i = 1:length(all_contents)
    working_folder = all_contents(i); %selects one of the n conditions
    working_folder_contents = dir(fullfile(strcat(working_folder.folder, '\', working_folder.name), '*cell*')); %this is why you must organize cells as they are in the sample folder
    %working_folder_contents should be an m x 1 struct where m = the number
    %of cells you have for a given condition
    
    %progressbar('Overall Progress') %init 3 bars
    for ii = 1:length(working_folder_contents) %iterate through all cells of the given condition
         working_cell = working_folder_contents(ii);
         cell_path = strcat(working_cell.folder, '\', working_cell.name);
         cell_recordings = dir(fullfile(cell_path, '*.abf')); %get all recordings, extension is .abf

         final_data_output = cell(7000,1); %simple preallocation, 7000 is arbitrary
         
         for iii = 1:length(cell_recordings) %iterate through all recordings of a given cell
             
          
            working_filename = strcat(cell_recordings(iii).folder, '\', cell_recordings(iii).name);
             
            %this small function helps correct for baseline drift, resets
            %baseline to 0 100 times over the course of your recording
            new_d_I = detrend_y_convert(working_filename);
            detrended_d_I = new_d_I'; %simple transpose of data

%% This Block runs the function that actually identifies the events

            raw_trace = detrended_d_I; %preserve raw trace
            trace = rmmissing(raw_trace); %removes any missing points from abf
        
            %THIS IS THE FUNCTION THAT IDs EVENTS - will need to change
            %depending on your desire to analyze mEPSCs from PV-INs or PYR
            
            %FOR PV-INs
            %[yi, yii,event_minima, event_times, event_amp, event_rise1090, event_decaytau, all_fits] = find_EPSCs_PVIN(trace);

            %FOR PYR
            [yi, yii,event_minima, event_times, event_amp, event_rise1090, event_decaytau, all_fits] = find_EPSCs_PYR(trace);

%% Plot Trace with all Events Marked

%THIS DISPLAYS:
%1. Raw detrended trace (blue)
%2. Markers on minima with number (red circle)
%3. Moving baseline (black)
%4. Fit for decay (pink transparent circles on each event decay)

close all
figure
plot(trace,'Color', [.4 .7 1])
hold on
plot(event_times, trace(event_times), 'o','MarkerFaceColor','r','MarkerEdgeColor','r')%, 'LineWidth', 2)

plot(yi, 'color', 'k')
%  hold on #un-comment this out if you want to see that smoothed trace
%  plot(yii)

combine_fits = cell2mat(all_fits);
scatter_fits = scatter(combine_fits(:,1),combine_fits(:,2),'MarkerFaceColor',[1, 0.4, 0.7],'MarkerEdgeColor',[1, 0.4, 0.7]); 
% Set property MarkerFaceAlpha and MarkerEdgeAlpha to <1.0
scatter_fits.MarkerFaceAlpha = .1;
scatter_fits.MarkerEdgeAlpha = .4;

for iv = 1:length(event_times)
    t = text(event_times(iv),trace(event_times(iv)),num2str(iv), 'fontweight', 'bold', 'Color', 'k', 'FontSize', 12);
end
xlim([1 20000])
%% Save fig where original .abf is so you can open it back up in matlab to check it out
save_fig_path = cell_recordings(iii).folder;
save_fig_title = strcat(cell_recordings(iii).name, '_Figure.fig');

if isfile(fullfile(save_fig_path, save_fig_title))
    delete(fullfile(save_fig_path, save_fig_title))
end

saveas(gcf, fullfile(save_fig_path, save_fig_title));

close all       
%% Make Cell Array for this recording's data
frequency = length(event_amp) / 60; %<------- CHANGE THIS TO LENGTH OF YOUR RECORDING (listed in seconds)
average_amp = mean(event_amp);
average_rise = nanmean(event_rise1090);
average_decay = nanmean(event_decaytau);
         
     working_data_cell_array = cell(8000,9);
     working_data_cell_array(1, 1:9) = {'Recording', 'Event Times', 'Event Amplitudes(pA)', 'Event Rise 10-90 (ms)', 'Event Decay Tau (ms)', 'Average Amplitude', 'Frequency (Hz)', 'Average Rise 10-90(ms)', 'Average Decay tau (ms)'};
     working_data_cell_array(2,1) = {cell_recordings(iii).name};
     working_data_cell_array(3:length(event_times)+1,1) = num2cell([2:1:length(event_times)]');
     working_data_cell_array(2:length(event_times)+1,2) = num2cell(event_times');
     working_data_cell_array(2:length(event_amp)+1,3) = num2cell(event_amp');
     working_data_cell_array(2:length(event_amp)+1,4) = num2cell(event_rise1090');
     working_data_cell_array(2:length(event_amp)+1,5) = num2cell(event_decaytau');
     working_data_cell_array(2,6) = {average_amp};
     working_data_cell_array(2,7) = {frequency};
     working_data_cell_array(2,8) = {average_rise};
     working_data_cell_array(2,9) = {average_decay};
     
     

xlsx_name = strcat(cell_recordings(iii).name, '_','Cell_output.xlsx');
savepath = strcat(working_cell.folder, '\', working_cell.name);
full_xlsxname = fullfile(savepath, xlsx_name);
if isfile(full_xlsxname)
    delete(full_xlsxname)
end
xlswrite(full_xlsxname, working_data_cell_array);
    %add this recording's data to your final data output for that cell
    %final_data_output = horzcat(final_data_output, working_data_cell_array)
    
clear yi
clear yii
clear event_minima
clear event_times
clear event_amp
clear event_rise1090
clear event_decaytau 
clear frequency
clear average_amp
clear average_rise
clear average_decay
clear all_fits

    
         end %third for loop, iii, each recording
         
    end %second for loop, ii
    
end %first for loop, i
 
%% Concatenate averaged results from each recording for each individual cell

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     STEP2      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%select entire Condition folder that has all the cell folders in it
%genotype_folder = uigetdir('Z:\Luke\Electrophysiology\PV-mC4_P40_all\4.PV-mC4_P40-60_mEPSCs_PV', 'Select a folder');
genotype_folder = uigetdir('Z:\Luke\Electrophysiology\PV-mC4_P40_all\5.PV-mC4_P40-60_mEPSCs_PYR', 'Select a folder');


all_cell_folders = dir([genotype_folder, '\*cell*']);

for i = 1:length(all_cell_folders)
    current_cell = strcat(all_cell_folders(i).folder, '\', all_cell_folders(i).name);
    remaining_xlsx = dir(fullfile(strcat(all_cell_folders(i).folder, '\', all_cell_folders(i).name), '*Cell_output.xlsx*'));
    remaining_xlsx_cell_array = cell(length(remaining_xlsx)+ 4,5);
    remaining_xlsx_cell_array(1,1:5) = {'Filename', 'Average Amplitude (mV)', 'Average Frequency (Hz)', 'Average Rise 10-90 (ms)', 'Average Decay tau (ms)'};
    remaining_xlsx_cell_array(end-1:end,1) = {'Cell Standard Deviation', 'Cell Average'};


    for ii = 1:length(remaining_xlsx)
        working_xlsx = readtable(strcat(remaining_xlsx(ii).folder, '\', remaining_xlsx(ii).name));
        working_xlsx = table2cell(working_xlsx);
        remaining_xlsx_cell_array(ii+1, 1) = working_xlsx(1,1);
        remaining_xlsx_cell_array(ii+1, 2) = working_xlsx(1,6);
        remaining_xlsx_cell_array(ii+1, 3) = working_xlsx(1,7);
        remaining_xlsx_cell_array(ii+1, 4) = working_xlsx(1,8);
        remaining_xlsx_cell_array(ii+1, 5) = working_xlsx(1,9);
    end

    %Amplitude
    remaining_xlsx_cell_array(end-1,2) = {std(cell2mat(remaining_xlsx_cell_array(2:end-3,2)))};
    remaining_xlsx_cell_array(end,2) = {mean(cell2mat(remaining_xlsx_cell_array(2:end-3,2)))};

    %Freq
    remaining_xlsx_cell_array(end-1,3) = {std(cell2mat(remaining_xlsx_cell_array(2:end-3,3)))};
    remaining_xlsx_cell_array(end,3) = {mean(cell2mat(remaining_xlsx_cell_array(2:end-3,3)))};
    
    %Rise
    remaining_xlsx_cell_array(end-1,4) = {std(cell2mat(remaining_xlsx_cell_array(2:end-3,4)))};
    remaining_xlsx_cell_array(end,4) = {mean(cell2mat(remaining_xlsx_cell_array(2:end-3,4)))};
    
   %Decay
    remaining_xlsx_cell_array(end-1,5) = {std(cell2mat(remaining_xlsx_cell_array(2:end-3,5)))};
    remaining_xlsx_cell_array(end,5) = {mean(cell2mat(remaining_xlsx_cell_array(2:end-3,5)))};
        

    Parts = strsplit(current_cell, '\');
    final_xlsx_name = cell2mat(strcat(Parts(end), '_FinalResults.xlsx'));
    full_xlsxname = fullfile(current_cell, final_xlsx_name);
    if isfile(full_xlsxname)
        delete(full_xlsxname)
    end
    xlswrite(full_xlsxname, remaining_xlsx_cell_array);
end

%% Get Amplitudes, IEIs, Rises, and Decays for all recordings that you kept - will need this to make cumulative frequnecy distributions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     STEP3      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

working_folder = uigetdir('Z:\Luke\Electrophysiology\PV-mC4_P40_all\4.PV-mC4_P40-60_mEPSCs_PV', 'Select a folder');
%working_folder = uigetdir('Z:\Luke\Electrophysiology\PV-mC4_P40_all\5.PV-mC4_P40-60_mEPSCs_PYR', 'Select a folder');
working_folder_contents = dir([working_folder, '\*cell*']);
    
    for ii = 1:length(working_folder_contents)
         working_cell = working_folder_contents(ii);
         cell_path = strcat(working_cell.folder, '\', working_cell.name);
         cell_path_parts = strsplit(cell_path, '\');
         recording_xlsxs = dir(fullfile(cell_path, '*_Cell_output.xlsx'));
         
        if length(recording_xlsxs) > 0
        working_cell_xlsx = readtable(strcat(recording_xlsxs(1).folder, '\', recording_xlsxs(1).name));
        working_cell_array = table2cell(working_cell_xlsx(:,1:5));
        %Find IEI
            for iv = 1:length(working_cell_array)-1
                working_cell_array(iv+1,6) = {(cell2mat(working_cell_array(iv+1,2)) - cell2mat(working_cell_array(iv,2)))*0.05};
            end 
            
         
            if length(recording_xlsxs) > 1
         
                for iii = 2:length(recording_xlsxs)
                    working_cell_xlsx = readtable(strcat(recording_xlsxs(iii).folder, '\', recording_xlsxs(iii).name));
                    temp_working_array = table2cell(working_cell_xlsx(:,1:5));
                        for iv = 1:length(temp_working_array)-1
                            temp_working_array(iv+1,6) = {(cell2mat(temp_working_array(iv+1,2)) - cell2mat(temp_working_array(iv,2)))*0.05};
                        end
                    working_cell_array = vertcat(working_cell_array, temp_working_array);
                end
            elseif length(recording_xlsxs) == 1
               working_cell_array = working_cell_array;
            end

            %AVG IEI
            working_cell_array(1,7) = {mean(cell2mat(working_cell_array(:,6)))};
         
            %AVG Freq
            working_cell_array(1,8) = {(1/cell2mat(working_cell_array(1,7)))*1000};
         
            %AVG Amplitude
            working_cell_array(1,9) = {abs(mean(cell2mat(working_cell_array(:,3))))};
            
            %AVG Rise
            working_cell_array(1,10) = {abs(mean(cell2mat(working_cell_array(:,4))))};
            
            %AVG Decay
            working_cell_array(1,11) = {abs(mean(cell2mat(working_cell_array(:,5))))};
           
         
            working_cell_array_header = {'Recording', 'Events', 'Amplitude (mV)', 'Rise (ms)', 'Decay (ms)' 'IEI (ms)', 'Average IEI (ms)', 'Average Frequency (Hz)', 'Average Amplitude (mV)' 'Average Rise (ms)', 'Average Decay (ms)'};
         
            working_cell_array = vertcat(working_cell_array_header, working_cell_array);
           

            xlsx_name = strcat(working_cell.name, '_','FinalResults_w_IEI.xlsx');
            full_xlsxname = fullfile(cell_path, xlsx_name);
            if isfile(full_xlsxname)
                delete(full_xlsxname)
            end
            xlswrite(full_xlsxname,working_cell_array);
            
        elseif length(recording_xlsxs) == 0 
            display(['There are no files to be concatenated in ', strcat(working_cell.folder, '\', working_cell.name)]);
        end
            
    end
    
%end

%% Compile all the averaged data for each condition

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     STEP4      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%working_folder = uigetdir('Z:\Luke\Electrophysiology\PV-mC4_P40_all\4.PV-mC4_P40-60_mEPSCs_PV', 'Select a folder');
working_folder = uigetdir('Z:\Luke\Electrophysiology\PV-mC4_P40_all\5.PV-mC4_P40-60_mEPSCs_PYR', 'Select a folder');
working_folder_contents = dir([working_folder, '\*cell*']);

%rearrange cells so they are in chronological order
all_cells = working_folder_contents;
[size_all_cells, ~] = size(all_cells);

for i = 1:size_all_cells
    temp_cell = strsplit(all_cells(i).name, '_');
    cell_num(i,1) =  string(temp_cell(2));
end
    
for i = 1:size_all_cells
    to_sort_all_cells(i,1) = str2num(cell_num(i,1));
end

sorted_all_cells = sortrows(to_sort_all_cells);

for i = 1:size_all_cells
    working_folder_contents(i).name = strcat('cell_', int2str(sorted_all_cells(i,1)));
end


%now actually make .xlsx

    Final_Results(1,1:5) = {'Cell', 'Average Amplitude', 'Average Frequency', 'Average Rise 10-90 (ms)', 'Average Decay tau (ms)'};
    
    for ii = 1:length(working_folder_contents)
         working_cell = working_folder_contents(ii);
         cell_path = strcat(working_cell.folder, '\', working_cell.name);
         cell_path_parts = strsplit(cell_path, '\');
         Final_Results(ii+1,1) = cell_path_parts(1,end);
         cell_final_xlsx = dir(fullfile(cell_path, '*_FinalResults.xlsx'));
         working_cell_xlsx = readtable(strcat(cell_final_xlsx.folder, '\', cell_final_xlsx.name));
         Final_Results(ii+1, 2) = table2cell(working_cell_xlsx(end,2));
         Final_Results(ii+1, 3) = table2cell(working_cell_xlsx(end,3));
         Final_Results(ii+1, 4) = table2cell(working_cell_xlsx(end,4));
         Final_Results(ii+1, 5) = table2cell(working_cell_xlsx(end,5));
    end
    
     xlsxname = strcat('_','FinalResults.xlsx');
     savepath = working_folder_contents.folder;
     full_xlsxname = strcat(savepath, xlsxname);
     if isfile(full_xlsxname)
         delete(full_xlsxname)
     end
     xlswrite(full_xlsxname,Final_Results);
         
     
%% Combine Amplitudes and IEIs for each genotype

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     STEP5      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%working_folder = uigetdir('D:\Users\Luke_Fournier\Desktop\mEPSC_analysis\sample_data\PV-IN', 'Select a folder');
working_folder = uigetdir('D:\Users\Luke_Fournier\Desktop\mEPSC_analysis\sample_data\PYR', 'Select a folder');
working_folder_contents = dir([working_folder, '\*cell*']);

%rearrange cells so they are in chronological order
all_cells = working_folder_contents;
[size_all_cells, ~] = size(all_cells);

for i = 1:size_all_cells
    temp_cell = strsplit(all_cells(i).name, '_');
    cell_num(i,1) =  string(temp_cell(2));
end
    
for i = 1:size_all_cells
    to_sort_all_cells(i,1) = str2num(cell_num(i,1));
end

sorted_all_cells = sortrows(to_sort_all_cells);

for i = 1:size_all_cells
    working_folder_contents(i).name = strcat('cell_', int2str(sorted_all_cells(i,1)));
end

%now make excel
 
 numcols = length(working_folder_contents)*2;
 place_amp_header = numcols/2 + 1;
 all_IEIs = cell(1,numcols);
 all_IEIs(1,1:numcols/2) = {'IEI (ms)'};
 all_IEIs(1,place_amp_header:end) = {'Amplitude (pA)'};
 


for i = 1:length(working_folder_contents)
    working_cell = working_folder_contents(i);
    working_IEI_file = dir(strcat(working_cell.folder, '\', working_cell.name, '\', '*IEI.xlsx')); 
    working_IEI_file = xlsread(strcat(working_IEI_file.folder, '\', working_IEI_file.name));
    working_IEI_file = num2cell(working_IEI_file(:,3:6));
    all_IEIs(2,i) = {working_cell.name};
    all_IEIs(2,i+(numcols/2)) = {working_cell.name};
    
    for ii = 1:length(working_IEI_file)
        all_IEIs(ii+2,i) = working_IEI_file(ii,4);
        all_IEIs(ii+2,i+(numcols/2)) = working_IEI_file(ii,1);
    end
end

     xlsxname = strcat('_','all_IEIs.xlsx');
     savepath = working_folder_contents.folder;
     full_xlsxname = strcat(savepath, xlsxname);
     if isfile(full_xlsxname)
         delete(full_xlsxname)
     end
     xlswrite(full_xlsxname,all_IEIs);

%% Combine Rise/Decays for each genotype

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     STEP6      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%working_folder = uigetdir('D:\Users\Luke_Fournier\Desktop\mEPSC_analysis\sample_data\PV-IN', 'Select a folder');
working_folder = uigetdir('D:\Users\Luke_Fournier\Desktop\mEPSC_analysis\sample_data\PYR', 'Select a folder');

working_folder_contents = dir([working_folder, '\*cell*']);

%rearrange cells so they are in chronological order
all_cells = working_folder_contents;
[size_all_cells, ~] = size(all_cells);

for i = 1:size_all_cells
    temp_cell = strsplit(all_cells(i).name, '_');
    cell_num(i,1) =  string(temp_cell(2));
end
    
for i = 1:size_all_cells
    to_sort_all_cells(i,1) = str2num(cell_num(i,1));
end

sorted_all_cells = sortrows(to_sort_all_cells);

for i = 1:size_all_cells
    working_folder_contents(i).name = strcat('cell_', int2str(sorted_all_cells(i,1)));
end

%now make excel
 
 numcols = length(working_folder_contents)*2;
 place_decay_header = numcols/2 + 1;
 all_kinetics = cell(1,numcols);
 all_kinetics(1,1:numcols/2) = {'Decay (ms)'};
 all_kinetics(1,place_decay_header:end) = {'Rise(ms)'};

for i = 1:length(working_folder_contents)
    working_cell = working_folder_contents(i);
    working_kin_file = dir(strcat(working_cell.folder, '\', working_cell.name, '\', '*IEI.xlsx')); 
    working_kin_file = xlsread(strcat(working_kin_file.folder, '\', working_kin_file.name));
    working_kin_file = num2cell(working_kin_file(:,4:5));
    all_kinetics(2,i) = {working_cell.name};
    all_kinetics(2,i+(numcols/2)) = {working_cell.name};
    
    for ii = 1:length(working_kin_file)
        all_kinetics(ii+2,i) = working_kin_file(ii,2);
        all_kinetics(ii+2,i+(numcols/2)) = working_kin_file(ii,1);
    end
end

     xlsxname = strcat('_','all_Kinetics.xlsx');
     savepath = working_folder_contents.folder;
     full_xlsxname = strcat(savepath, xlsxname);
     if isfile(full_xlsxname)
         delete(full_xlsxname)
     end
     xlswrite(full_xlsxname,all_kinetics);
    
%%
%EOC