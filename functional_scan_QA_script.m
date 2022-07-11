%% a script to run QA on functional scans individually and as part of a group 
% to flag up outliers with large fluctuations in signal values - an
% indicator of physics related artefacts from the scanner or surrounds that
% may effect analysis of scans
% this script is run on workstation with multiple nodes - so can use the
% parfor function to speed up processing, 
% if used on a convential PC change all instances of parfor to 'for'

% cd to directory that contains functional scans to run QA on and use 'dir'
% function to generate list of files
clear
cd ''
subs = struct2cell(dir("sub-CC*")).';
subs = subs([1:108,110:end],1);

% number of subjects/rows in the list
num_sub = 651;

for z = 1:num_sub
    
    % cd to directory that contains functional scans to run QA on
    cd '';
    % loop through each of the scan files (in generated 'sub' variable) -
    % isolate sections of the filepath inc. subject ID
    file_path = strcat('C:...', subs{z}, '\swu', subs{z}, '_epi_rest.nii');
    file_char = char(file_path);
    full_dir = file_char(1:64);
    ID2 = file_char(75:80);
    ID = str2double(ID2)

% go into the paths one by one
cd(full_dir)

% get volumes from that particular subject, using the epi file name
vols = cellstr(file_path);%sub;
M=nifti(vols);

% calculates the mean of each slice for each volume
parfor v = 1:261
    for s = 1:32
        slicemean(v,s) = mean(mean(M.dat(:,:,s,v)));
    end
end

% transposes the matrix
slicemean_t = slicemean(:,:)';

% mean centre over timeseries within each slice (subtracts the mean from
% each value, to normalise the score, resetting the mean to zero)
parfor n = 1:32
    slicemean_norm(n,:) = slicemean_t(n,:) - mean(slicemean_t(n,:));   
end

% calculate mean for each volume, and other relevant stats
slicemean_mean = mean(slicemean(:,:));
slicemean_norm_mean = mean(slicemean_norm);
slicemean_std = std(slicemean(:,:));
slicemean_norm_std = std(slicemean_norm);
parfor www = 1:length(slicemean(:,:))
   slicemean_mean2(www) = mean(slicemean(www,:));
   slicemean_std2(www) = std(slicemean(www,:));
end
% fast fourier transform
slicemean_fft = (abs(fft(slicemean_norm')))';

% calculates the range of scores accross volumes, within each slice (in the
% temporal domain)
parfor xnx = 1:length(slicemean(1,:))
    all_val_t(xnx) = max(slicemean(:,xnx)) - min(slicemean(:,xnx));
end
% calculates the maximum range accross all slices, telling us within this
% subject, the greatest variability over time.
max_range_temporal = max(all_val_t);

% does the same as above, but over slices, within each volume, telling us
% the greatest range of activation, within a volume, accross the slices
% (thus in the spatial domain)
parfor nxn = 1:length(slicemean(:,1))
    all_val_s(nxn) = max(slicemean(nxn,:)) - min(slicemean(nxn,:));
end
% calculates the maximum
max_range_spatial = max(all_val_s);

% create a results file with key info for automatic outlier
% identification
    results(z,1) = ID;

% add the max values for each participant
    results(z,2) = max_range_temporal;
    results(z,3) = max_range_spatial;
        
% calculate the coefficient of variation for each participant, looking at
% the average of each slice. This is run at each slice accross time (cov)
% and at each scan accross slices (cov2)
cov = slicemean_std ./ slicemean_mean;
cov2 = slicemean_std2 ./ slicemean_mean2;
% calculate the mean average, max average and median average accross all
% slices, and accross all volumes
covmax = max(cov);
cov2max = max(cov2);
covmean = mean(cov);
cov2mean = mean(cov2);
covmed = median(cov);
cov2med = median(cov2);

% create a matrix containing the key values. cov is temporal, cov2 is
% spatial
cov_results(z,1) = ID;
cov_results(z,2) = covmax;
cov_results(z,3) = covmean;
cov_results(z,4) = covmed;
cov_results(z,5) = cov2max;
cov_results(z,6) = cov2mean;
cov_results(z,7) = cov2med;

figure(z);

% specify the file into a variable to later save the plots_prepro in a loop
h = figure(z);

% hides the figures, so they don't dominate the screen
set(h, 'Visible', 'off');

% plots_prepro the standard deviation
    subplot(5,1,1), plot(slicemean_std);
    subplot(5,1,1), title('standard deviation of each volume');
    subplot(5,1,1), xlabel('slice');
    subplot(5,1,1), ylabel('std');
    %subplot(5,1,2), axis([1 length(B) min(B) max(B)])    
%     colorbar;
    
% plots_prepro the normalised mean
    subplot(5,1,2), imagesc(slicemean_norm);
    subplot(5,1,2), title('signal intensity (slice mean corrected)');
    subplot(5,1,2), xlabel('volume');
    subplot(5,1,2), ylabel('slice');
    colorbar;
    %subplot(3,1,2), axis([1 length(A_norm) min(min(A_norm)) max(max(A_norm))])
    
    % plots_prepro the FFT
    subplot(5,1,3), imagesc(slicemean_fft(:,2:(length(slicemean_fft)-1)));
    subplot(5,1,3), axis([1 50 1 size(slicemean_fft,1)]);
    subplot(5,1,3), title('Fast Fourier Transform of above');
    subplot(5,1,3), xlabel('Number of cycles in timecourse');
    subplot(5,1,3), ylabel('slice');
    colorbar;

    % plots_prepro the cov
    subplot(5,1,4), plot(cov);
    subplot(5,1,4), title('coefficient of variance temporal');
    subplot(5,1,4), xlabel('slice');
    subplot(5,1,4), ylabel('cov');
    subplot(5,1,4), axis([1 length(cov) min(cov) max(cov)]);
%     colorbar;
    

% plots_prepro the cov2
    subplot(5,1,5), plot(cov2);
    subplot(5,1,5), title('coefficient of variance spatial');
    subplot(5,1,5), xlabel('volume');
    subplot(5,1,5), ylabel('cov');
    subplot(5,1,5), axis([1 length(cov2) min(cov2) max(cov2)]);
%     colorbar;

    % cd into the directory to save plots_prepro, then save each plot
    % in a loop
    cd 'C:...\plots_prepro';
    saveas(h,sprintf('FIG%d.jpg',z));

    clear M max* slice* all* cov cov2* covmax covmean covmed vols
end

%%%%%% NEEDS UPDATING %%%%%%%
% saves key files
cd '';
save results results;
save cov_results cov_results;

%% From here use the stats calculated for each subjects functional scans
% make group comparisons to identify outliers

num_sub = 651;

% loads key files
cd '';
load results;
load cov_results;

%%%%%% NEEDS UPDATING %%%%%%%
% set threshold for outlier identification (SD)
thresh_val = 2.5;

% generates a threshold of mean plus thresh_val stds for results table
temp_mean = mean(results(:,2));
temp_std = std(results(:,2));
thresh_temp = temp_mean + (thresh_val.*temp_std);

spat_mean = mean(results(:,3));
spat_std = std(results(:,3));
thresh_spat = spat_mean + (thresh_val.*spat_std);

covmax_mean = mean(cov_results(:,2));
covmax_std = std(cov_results(:,2));
cmax_thresh = covmax_mean + (thresh_val.*covmax_std);

covmean_mean = mean(cov_results(:,3));
covmean_std = std(cov_results(:,3));
cmean_thresh = covmean_mean + (thresh_val.*covmean_std);

covmed_mean = mean(cov_results(:,4));
covmed_std = std(cov_results(:,4));
cmed_thresh = covmed_mean + (thresh_val.*covmed_std);

cov2max_mean = mean(cov_results(:,5));
cov2max_std = std(cov_results(:,5));
c2max_thresh = cov2max_mean + (thresh_val.*cov2max_std);

cov2mean_mean = mean(cov_results(:,6));
cov2mean_std = std(cov_results(:,6));
c2mean_thresh = cov2mean_mean + (thresh_val.*cov2mean_std);

cov2med_mean = mean(cov_results(:,7));
cov2med_std = std(cov_results(:,7));
c2med_thresh = cov2med_mean + (thresh_val.*cov2med_std);

% saves histograms of all relevant measures
figure(1000); histogram(cov_results(:,2));
title('Maximum Cov (temporal)')
xlabel('Coefficient of variation (maximum)');
ylabel('Number of participants');

figure(1001); histogram(cov_results(:,3));
title('Mean Cov (temporal)')
xlabel('Coefficient of variation (mean)');
ylabel('Number of participants');

figure(1002); histogram(cov_results(:,4));
title('Median Cov (temporal)')
xlabel('Coefficient of variation (median)');
ylabel('Number of participants');

figure(1003); histogram(cov_results(:,5));
title('Maximum Cov (spatial)')
xlabel('Coefficient of variation (maximum)');
ylabel('Number of participants');

figure(1004); histogram(cov_results(:,6));
title('Mean Cov (spatial)')
xlabel('Coefficient of variation (mean)');
ylabel('Number of participants');

figure(1005); histogram(cov_results(:,7));
title('Median Cov (spatial)')
xlabel('Coefficient of variation (median)');
ylabel('Number of participants');

figure(1006); histogram(results(:,2));
title('Max Range (temporal)')
xlabel('Max range');
ylabel('Number of participants');

figure(1007); histogram(results(:,3));
title('Max Range (spatial)')
xlabel('Max range');
ylabel('Number of participants');
j = figure(1000);
% 
saveas(j,'maxab_norm_cov.jpg');
jj = figure(1001);
saveas(jj,'mean_cov.jpg');
jjj = figure(1002);
saveas(jjj,'median_cov.jpg');
l = figure(1003);
saveas(l,'maxab_cov2.jpg');
ll = figure(1004);
saveas(ll,'mean_cov2.jpg');
lll = figure(1005);
saveas(lll,'median_cov2.jpg');
rangt = figure(1006);
saveas(rangt,'range_temporal.jpg');
rangs = figure(1007);
saveas(rangs,'range_spatial.jpg');

% runs through all subjects, if their max range (temporal) is above the threshold it
% is saved into an outliers file, if not it goes into an included file

outliers_temp = zeros(10,2);
for loop = 1:num_sub
   % to regenerate ID for each subject - to add them to the outlier/
   % not-outlier lists - again isolate the ID from the file path, this is
   % pointlessly long-winded, generate a single list of IDs and filepaths
   % at start of script and call on that throughout the script
   
    clear file_char
    clear ID3
    clear ID4
    
    file_path = strcat('C:...\epi_rest\', subs{loop}, '\', subs{loop}, '_epi_rest.nii');
    file_char = char(file_path);
    full_dir = file_char(1:64);
    ID3 = file_char(72:77);
    ID4 = str2double(ID3)

    if results(loop,2) > thresh_temp
        outliers_temp(loop,1) = ID4;
        outliers_temp(loop,2) = results(loop,2);
    else included_temp(loop,1) = ID4;
        included_temp(loop,2) = results(loop,2);
    end
end
clear loop

% removes zero values 
outliernew_temp = outliers_temp(outliers_temp(:,1)>0,:);
includednew_temp = included_temp(included_temp(:,1)>0,:);
% orders files in ascending order
outlier_files_temp = sortrows(outliernew_temp,1);
included_files_temp = sortrows(includednew_temp,1);

outliers_spat = zeros(10,2);
% same for spatial range
for loop = 1:num_sub
    clear file_char
    clear ID3
    clear ID4

    file_path = strcat('C:...\epi_rest\', subs{loop}, '\', subs{loop}, '_epi_rest.nii');
    file_char = char(file_path);
    full_dir = file_char(1:64);
    ID3 = file_char(72:77);
    ID4 = str2double(ID3)

    if results(loop,3) > thresh_spat
        outliers_spat(loop,1) = ID4;
        outliers_spat(loop,2) = results(loop,3);
    else included_spat(loop,1) = ID4;
        included_spat(loop,2) = results(loop,3);
    end
end
clear loop

% removes zero values
outliernew_spat = zeros(10,2);
outliernew_spat = outliers_spat(outliers_spat(:,1)>0,:);
includednew_spat = included_spat(included_spat(:,1)>0,:);
% orders files in ascending order
outlier_files_spat = sortrows(outliernew_spat,1);
included_files_spat = sortrows(includednew_spat,1);

% covmax
% runs through all subjects, if their max value is above the threshold it
% is saved into an outliers file, if not it goes into an included file. We
% do not care about min values as low cov means low variability.

outliers_cmax = zeros(10,2);
for loop = 1:num_sub
    clear file_char
    clear ID3
    clear ID4
    
    file_path = strcat('C:...\epi_rest\', subs{loop}, '\', subs{loop}, '_epi_rest.nii');
    file_char = char(file_path);
    full_dir = file_char(1:64);
    ID3 = file_char(72:77);
    ID4 = str2double(ID3)

    if cov_results(loop,2) > cmax_thresh
        outliers_cmax(loop,1) = ID4;
        outliers_cmax(loop,2) = cov_results(loop,2);
    else included_cmax(loop,1) = ID4;
        included_cmax(loop,2) = cov_results(loop,2);
    end
end
clear loop

% removes zero values
outliernew_cmax = outliers_cmax(outliers_cmax(:,1)>0,:);
includednew_cmax = included_cmax(included_cmax(:,1)>0,:);
% orders files in ascending order
outlier_files_cmax = sortrows(outliernew_cmax,1);
included_files_cmax = sortrows(includednew_cmax,1);

outliers_cmean = zeros(10,2);
% cov mean
for loop = 1:num_sub
    clear file_char
    clear ID3
    clear ID4
    
    file_path = strcat('C:...\epi_rest\', subs{loop}, '\', subs{loop}, '_epi_rest.nii');
    file_char = char(file_path);
    full_dir = file_char(1:64);
    ID3 = file_char(72:77);
    ID4 = str2double(ID3)

    if cov_results(loop,3) > cmean_thresh
        outliers_cmean(loop,1) = ID4;
        outliers_cmean(loop,2) = cov_results(loop,3);
    else included_cmean(loop,1) = ID4;
        included_cmean(loop,2) = cov_results(loop,3);
    end
end
clear loop

% removes zero values
outliernew_cmean = outliers_cmean(outliers_cmean(:,1)>0,:);
includednew_cmean = included_cmean(included_cmean(:,1)>0,:);
% orders files in ascending order
outlier_files_cmean = sortrows(outliernew_cmean,1);
included_files_cmean = sortrows(includednew_cmean,1);

outliers_cmed = zeros(10,2);
% covmed
for loop = 1:num_sub
    clear file_char
    clear ID3
    clear ID4

    file_path = strcat('C:...\epi_rest\', subs{loop}, '\', subs{loop}, '_epi_rest.nii');
    file_char = char(file_path);
    full_dir = file_char(1:64);
    ID3 = file_char(72:77);
    ID4 = str2double(ID3)

    if cov_results(loop,4) > cmed_thresh
        outliers_cmed(loop,1) = ID4;
        outliers_cmed(loop,2) = cov_results(loop,4);
    else included_cmed(loop,1) = ID4;
        included_cmed(loop,2) = cov_results(loop,4);
    end
end
clear loop

% removes zero values
outliernew_cmed = outliers_cmed(outliers_cmed(:,1)>0,:);
includednew_cmed = included_cmed(included_cmed(:,1)>0,:);
% orders files in ascending order
outlier_files_cmed = sortrows(outliernew_cmed,1);
included_files_cmed = sortrows(includednew_cmed,1);

outliers_c2max = zeros(10,2);
%c2max
for loop = 1:num_sub
    clear file_char
    clear ID3
    clear ID4

    file_path = strcat('C:...\epi_rest\', subs{loop}, '\', subs{loop}, '_epi_rest.nii');
    file_char = char(file_path);
    full_dir = file_char(1:64);
    ID3 = file_char(72:77);
    ID4 = str2double(ID3)

    if cov_results(loop,5) > c2max_thresh
        outliers_c2max(loop,1) = ID4;
        outliers_c2max(loop,2) = cov_results(loop,5);
    else included_c2max(loop,1) = ID4;
        included_c2max(loop,2) = cov_results(loop,5);
    end
end
clear loop

% removes zero values
outliernew_c2max = outliers_c2max(outliers_c2max(:,1)>0,:);
includednew_c2max = included_c2max(included_c2max(:,1)>0,:);
% orders files in ascending order
outlier_files_c2max = sortrows(outliernew_c2max,1);
included_files_c2max = sortrows(includednew_c2max,1);

outliers_c2mean = zeros(10,2);
% cov2 mean
for loop = 1:num_sub
    clear file_char
    clear ID3
    clear ID4
   
    file_path = strcat('C:...\epi_rest\', subs{loop}, '\', subs{loop}, '_epi_rest.nii');
    file_char = char(file_path);
    full_dir = file_char(1:64);
    ID3 = file_char(72:77);
    ID4 = str2double(ID3)

    if cov_results(loop,6) > c2mean_thresh
        outliers_c2mean(loop,1) = ID4;
        outliers_c2mean(loop,2) = cov_results(loop,6);
    else included_c2mean(loop,1) = ID4;
        included_c2mean(loop,2) = cov_results(loop,6);
    end
end
clear loop

% removes zero values
outliernew_c2mean = outliers_c2mean(outliers_c2mean(:,1)>0,:);
includednew_c2mean = included_c2mean(included_c2mean(:,1)>0,:);
% orders files in ascending order
outlier_files_c2mean = sortrows(outliernew_c2mean,1);
included_files_c2mean = sortrows(includednew_c2mean,1);

outliers_c2med = zeros(10,2);
% cov2med
for loop = 1:num_sub
    clear file_char
    clear ID3
    clear ID4
    
    cd 'C:...\epi_rest';

    file_path = strcat('C:...\epi_rest\', subs{loop}, '\', subs{loop}, '_epi_rest.nii');
    file_char = char(file_path);
    full_dir = file_char(1:64);
    ID3 = file_char(72:77);
    ID4 = str2double(ID3)

    if cov_results(loop,7) > c2med_thresh
        outliers_c2med(loop,1) = ID4;
        outliers_c2med(loop,2) = cov_results(loop,7);
    else included_c2med(loop,1) = ID4;
        included_c2med(loop,2) = cov_results(loop,7);
    end
end
clear loop

% removes zero values
outliernew_c2med = outliers_c2med(outliers_c2med(:,1)>0,:);
includednew_c2med = included_c2med(included_c2med(:,1)>0,:);
% orders files in ascending order
outlier_files_c2med = sortrows(outliernew_c2med,1);
included_files_c2med = sortrows(includednew_c2med,1);

% creates a grand file of all outliers, checking for duplicates and only
% adding if not already present in the file
grand_out = outlier_files_temp;

for aaa = 1:length(outlier_files_spat(:,1))
    if outlier_files_spat(aaa,1) ~= grand_out(:,1)
        nnn = length(grand_out(:,1)) + 1;
        grand_out(nnn,1) = outlier_files_spat(aaa,1);
        grand_out(nnn,2) = outlier_files_spat(aaa,2);
    end
end

for aaa = 1:length(outlier_files_cmax(:,1))
    if outlier_files_cmax(aaa,1) ~= grand_out(:,1)
        % includes the below file_char to consistenly update the grand_out fie
        % to write to the next available file_char, starting at 1 beause its
        % original length is zero, so 0+1 = 1st row.
        nnn = length(grand_out(:,1)) + 1;
        grand_out(nnn,1) = outlier_files_cmax(aaa,1);
        grand_out(nnn,2) = outlier_files_cmax(aaa,2);
    end
end

for aaa = 1:length(outlier_files_cmean(:,1))
    if outlier_files_cmean(aaa,1) ~= grand_out(:,1)
        nnn = length(grand_out(:,1)) + 1;
        grand_out(nnn,1) = outlier_files_cmean(aaa,1);
        grand_out(nnn,2) = outlier_files_cmean(aaa,2);
    end
end

for aaa = 1:length(outlier_files_cmed(:,1))
    if outlier_files_cmed(aaa,1) ~= grand_out(:,1)
        nnn = length(grand_out(:,1)) + 1;
        grand_out(nnn,1) = outlier_files_cmed(aaa,1);
        grand_out(nnn,2) = outlier_files_cmed(aaa,2);
    end
end

for aaa = 1:length(outlier_files_c2max(:,1))
    if outlier_files_c2max(aaa,1) ~= grand_out(:,1)
        nnn = length(grand_out(:,1)) + 1;
        grand_out(nnn,1) = outlier_files_c2max(aaa,1);
        grand_out(nnn,2) = outlier_files_c2max(aaa,2);
    end
end

for aaa = 1:length(outlier_files_c2mean(:,1))
    if outlier_files_c2mean(aaa,1) ~= grand_out(:,1)
        nnn = length(grand_out(:,1)) + 1;
        grand_out(nnn,1) = outlier_files_c2mean(aaa,1);
        grand_out(nnn,2) = outlier_files_c2mean(aaa,2);
    end
end

for aaa = 1:length(outlier_files_c2med(:,1))
    if outlier_files_c2med(aaa,1) ~= grand_out(:,1)
        nnn = length(grand_out(:,1)) + 1;
        grand_out(nnn,1) = outlier_files_c2med(aaa,1);
        grand_out(nnn,2) = outlier_files_c2med(aaa,2);
    end
end

%identifies included files (non-outliers) for double checking
for aaa = 1:length(cov_results)
    if cov_results(aaa,1) ~= grand_out(:,1)
        grand_in3(aaa,1) = cov_results(aaa,1);
        grand_in3(aaa,2) = 99999;
    elseif isempty(grand_out)
        grand_in3(aaa,1) = cov_results(aaa,1);
        grand_in3(aaa,2) = 99999;
    end
end

% removes zero values
grand_in2 = grand_in3(grand_in3(:,1)>0,:);
% orders ascending
grand_in = sortrows(grand_in2,1);

% generates random numbers ranging from 1 to the number of included files,
% then uses these to randomly select control participants for the blind
% manual check
% for check = 1:length(grand_out(:,1))
%     xxx = length(grand_in(:,1));
%     random = randi(xxx,1);
%     control_list(check,1) = grand_in(random,1);
%     control_list(check,2) = grand_in(random,2);
% end

% puts the files together 
% stack = [grand_out; control_list];
% % randomises rows
% stack(:,3) = rand(length(stack(:,1)),1);
% stack_rand = sortrows(stack,3);

% saves the relevant files
cd 'C:...\checks_prepro';
% save stack_rand stack_rand;
save grand_out grand_out
save grand_in grand_in
