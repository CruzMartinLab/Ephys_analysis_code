function [yi, yii,event_minima, event_times, event_amp, event_rise1090, event_decaytau, all_fits] = find_EPSCs_PVIN(trace)

y = trace'; %<---------- RAW TRACE

yi = smoothdata(y, 'sgolay', 4000); %<--------MOVING BASELINE - needed to determine amplitude of each event; should cut through the middle of the noise
yii = smoothdata(y, 'sgolay', 50); %<---------smoothed trace, used later for a couple parameters
yi = yi+1;%brings basleine up by 1 pA, found this to be more accurate
TF = islocalmin(yii, 'MinProminence', 5, 'MinSeparation', 100); %<--------editable dependning on what you want your min event size to be
%TF variable holes all local minima of trace (i.e. events); over-marks a ton obvioulsy,
%so from here we pare down to discard all the false-positives


%This makes array that gives the mV value and the x-coordinate for each
%event
peaks = [sum(TF),2];
peak_fill = 1; 
for i = 1:length(y)
    TF_binary = TF(i);
    if TF_binary == 1 %i.e. if data point is one of the identified minima
        [~, min_idx] = min(y(i-20:i+20)); %find min of points around it, helps control for noise
        min_idx = min_idx + i - 20; %gets actual idx of min relative to entitre trace
        min_value = mean(y(min_idx-5:min_idx + 5)); %we take average to account for noise around the min
        peaks(peak_fill,1) = min_value; %fill array
        peaks(peak_fill,2) = min_idx;%fill array
        peak_fill = peak_fill+1;
    end
end


%So now you have a m x 2 array where m = the number of minima detected
%However, this array gives the the mV value, not the amplitude of each
%minima relative to the baseline we set
%'deltas' is the array that gives you the actual amplitudes, as calculated
%by subtracting the minimum from yi, which is the line that cuts through
%the noise of the data (i.e. is the baseline)

deltas = [length(peaks)];
for i = 1:length(peaks)
    deltas(i,1) = yi(peaks(i,2)) - peaks(i,1);
end

all_data = horzcat(peaks,deltas);%concat. the data, just for organization-sake


%Now you have all the amplitudes in deltas, but the minimum function
%overmarks due to noise. in this step, we are finding the local noise on a
%second-by-second basis and eliminating events based on this noise. in this
%way, the code constantly adapts to the noise (i.e. how wide your baseline
%is)

%Anything in this section marked with a ** means you can edit it to play
%with the stringency of your event marking. I have found that these
%parameters work well for PV-INs
trace_parts = round(length(trace)/60); %splits trace into 60 even sections (~1s/section)
start_subtrace = 1;
for p = 1:59 %iterating through trace, second-by-second
    end_subtrace = p*trace_parts;
    working_subtrace = trace(start_subtrace:end_subtrace);    
    trace_rm = working_subtrace;
    indices = find(trace_rm < -10); %find noise by looking at all data points >-10mv to capture (roughly) the baseline)
    trace_rm(indices) = [];
    noise = std(trace_rm); %noise = standard deviation of the trace, a somewhat arbitrary value
    if noise < 1.25 %** this corresponds to a very flat baseline, clean recording
        for i = 1:length(all_data)
            if deltas(i,1) < 5 && all_data(i,2) > start_subtrace && all_data(i,2) < end_subtrace %** 5 is your minimum event size in this case
                all_data(i,1) = NaN; %if it is less than 5pA, remove this event from the all_data array
                all_data(i,2) = NaN;
                all_data(i,3) = NaN;
            end
        end
    elseif noise > 3.125 %** this corresponds to a really noisy of a recording, shakiness in baseline
        for i = 1:length(all_data)
            if deltas(i,1) < 10  && all_data(i,2) > start_subtrace && all_data(i,2) < end_subtrace %** now your minumum is 10pA
                all_data(i,1) = NaN;
                all_data(i,2) = NaN;
                all_data(i,3) = NaN;
            end
        end
    elseif noise > 1.25 && noise < 3.125 %** this corresponds to pretty clean recording with decently flat baseline
        for i = 1:length(all_data)
            if deltas(i,1) < noise*3.2 && all_data(i,2) > start_subtrace && all_data(i,2) < end_subtrace %** min corresponds to the noise 
                all_data(i,1) = NaN;
                all_data(i,2) = NaN;
                all_data(i,3) = NaN;
            end
        end
    end  

    start_subtrace = end_subtrace;
end

%for last bit of trace
working_subtrace = trace(end_subtrace:end);
trace_rm = working_subtrace;
indices = find(trace_rm < -10);
trace_rm(indices) = [];
noise = std(trace_rm);
if noise < 1.25
    for i = 1:length(all_data)
        if deltas(i,1) < 5 && all_data(i,2) > end_subtrace
            all_data(i,1) = NaN;
            all_data(i,2) = NaN;
            all_data(i,3) = NaN;
        end
    end
elseif noise > 3.125
    for i = 1:length(all_data) 
        if deltas(i,1) < 10 && all_data(i,2) > end_subtrace
            all_data(i,1) = NaN;
            all_data(i,2) = NaN;
            all_data(i,3) = NaN;
        end
    end
elseif noise > 1.25 && noise < 3.125  
    for i = 1:length(all_data)
        if deltas(i,1) < noise*3.2 && all_data(i,2) > end_subtrace
            all_data(i,1) = NaN;
            all_data(i,2) = NaN;
            all_data(i,3) = NaN;
        end
    end
end 

%Now we have replaced all rows with amplitudes < given min in deltas with NaN. In
%this step, we are deleting all the rows with NaN, because we no longer
%need them.

to_plot = [length(all_data),3];
to_plot_filler = 1;
for i = 1:length(all_data)
    if isnan(all_data(i,1))
        continue
    else
        to_plot(to_plot_filler,1) = all_data(i,1);
        to_plot(to_plot_filler,2) = all_data(i,2);
        to_plot(to_plot_filler,3) = all_data(i,3);
        to_plot_filler = to_plot_filler + 1;
    end
end

%%

%Finally, sometimes the noise exceeds the min, so an event may be marked despite it just being noise. However, the noise has a super
%fast decay, far shorter than a real mEPSC. Here, we eliminate those by taking only the minima where the 
%+ 10 points ahead of the peak is not within 2mV of that smoothed baseline
%we originally drew (yi)

for i = 1:length(to_plot)-1
    if isnan(to_plot(i,1))
          continue
    else   
        working_point = to_plot(i,:);
        start_window = working_point(2);
        end_window = start_window+10;
        points = start_window:end_window;

        for ii = 1:length(points)
            if yi(points(ii)) - yii(points(ii)) < 4  %** can play with this number
                current_point = points(ii);
                %if you quickly return to baseline (<4pA of baseline 10
                %points after minima was detected, i.e is decayed super
                %quick), look 5 more points ahead and see if it stays the
                %same
                add5 = current_point+5; %OG: 5
                if yi(add5) - yii(add5) < 5 
                    to_plot(i,1) = NaN;
                    to_plot(i,2) = NaN;
                    to_plot(i,3) = NaN;
                end
            end
        end    
    end
end


% Here we are doing something similar, but now looking behdind the detected
% minima instead of after it (i.e. we are looking at the rise to determine
% real event vs noise).
% This will also help eliminate noise, goes 65 points back from peak and
% looks to see if it is above the baseline, if it is it is probably just
% noise (real events typically won't rise that qucikly
for ii = 2:length(to_plot)
   if isnan(to_plot(ii,1))
       continue
   else
   working_peak = to_plot(ii,1);
   working_peak_idx = to_plot(ii,2);
   check_behind = working_peak_idx - 65; %**could play with this
       if  working_peak - yii(check_behind) > 0
            to_plot(ii,1) = NaN;
            to_plot(ii,2) = NaN;
            to_plot(ii,3) = NaN;
       end
   end
end


%Now making a final_results array since all the false-positives have been
%removed

final_results = [length(to_plot),5];
final_results_filler = 1;
for i = 1:length(to_plot)
    if isnan(to_plot(i,1))
        continue
    else
        final_results(final_results_filler,1) = to_plot(i,1);
        final_results(final_results_filler,2) = to_plot(i,2);
        final_results(final_results_filler,3) = to_plot(i,3);
        final_results_filler = final_results_filler + 1;
    end
end

%Find rise time 10-90 (time taken to rise from 10% of total amplitude to
% 90% of total amplitude; don for each event
for i = 3:length(final_results) %starting at 3 to account for events really close to start of recording
    working_delta = final_results(i,3);
    rise90_value = final_results(i,1) + (working_delta*0.1); 
    rise10_value = final_results(i,1) + (working_delta*0.9);
    rise90_idx = final_results(i,2);
    rise10_idx = final_results(i,2);
    while y(rise90_idx) < rise90_value
        rise90_idx = rise90_idx - 1;
    end
    
    while y(rise10_idx) < rise10_value
        if sum(y(rise10_idx-500:rise10_idx) < rise10_value) == 501
            rise10_idx = 0;
            break
        end
        rise10_idx = rise10_idx - 1;
    end
    
    if rise10_idx == 0 | rise90_idx == 0
        final_results(i,4) = NaN;
    else
    final_results(i,4) = (rise90_idx - rise10_idx)*0.05;
    end
end


%Find decay tau (63.2% decay back to baseline) - done for each event
all_fits = cell(1,1);

for i = 1:length(final_results)-2 %ending on second penultimate event to account for events really close to end of recording
    working_delta = final_results(i,3); %working on the current one
    decaytaumin_idx = final_results(i,2); %this is the min
    decaytaubl_idx = final_results(i,2); %this is what will change each loop until we get back to baseline
    %need to find where it reaches back to baseline
    while mean(yii(decaytaubl_idx-10:decaytaubl_idx+10)) < yi(decaytaubl_idx) %while smoothed trace < baseline
            decaytaubl_idx = decaytaubl_idx+1; %if smoothed is still less than baseline, move to next point
    end
   
    current_tofit_idc = (decaytaumin_idx:decaytaubl_idx)'; %indices to fit (peak to baseline)
    current_tofit = y(decaytaumin_idx:decaytaubl_idx); %actual pA values 
    
    if length(current_tofit_idc) < 2
        event_tau = NaN;
        all_fits{i,1} = NaN;
        all_fits{i,2} = NaN;
    else
        f = fit(current_tofit_idc,current_tofit, 'exp1', 'Normalize', 'on', 'Upper', [0,0]); %fit, need to normalize to avoid errors
        
        all_fits{i,1} = current_tofit_idc;
        all_fits{i,2} = f(current_tofit_idc);
        %save fits for each event in a separate cell array, can output and plot later in main code

        %now need to actually find tau - 63% of return to baseline
        fitted_yvalues = f(current_tofit_idc);

        pA_delta = abs(fitted_yvalues(1) - fitted_yvalues(end));
        perc63 = pA_delta*.628;
        tau_endpt = fitted_yvalues(1) + perc63;

        working_fitted_yvalues_idx = 1;

            while fitted_yvalues(working_fitted_yvalues_idx) < tau_endpt
                working_fitted_yvalues_idx = working_fitted_yvalues_idx + 1;
            end
            event_tau = working_fitted_yvalues_idx*0.05;
    end
    
    final_results(i,5) = event_tau;
   
end   

event_minima = final_results(:,1);
event_times = final_results(:,2);
event_amp = final_results(:,3);
event_rise1090 = final_results(:,4);
event_decaytau = final_results(:,5);
    
end
