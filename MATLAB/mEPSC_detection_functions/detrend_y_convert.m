function new_d_I = detrend_y_convert(working_file)

[d,si,h]=abfload_fcollman(working_file); %function for reading in .abf files to MATLAB, open source (I did not write abdload_fcollman)
    
d_I = squeeze(d(:,1,:));

new_d_I = nan(1200000,1); %preallocate, 60s, 20kHz sampling = 1200000 data points/recording
start = 1;
ending = 12000;

for i = 1:99                                  %
    working_d_I = d_I(start:ending);
    out = detrend(working_d_I,1);
    new_d_I(start:ending) = out;
    start = ending+1;
    ending = ending+12000;
end


%helps eliminate any positive noise (uncontrollable oscillations from
%computer, for example) - simply cuts off any data points above 6pA
noise_values = new_d_I > 6;            %
new_d_I(noise_values) = [];            %
end