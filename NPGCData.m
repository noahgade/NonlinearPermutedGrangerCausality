load("~/Desktop/Research/Rodu/NPGC/Application/Data/crcns-ac1/wehr/Stimuli/fragments/13.mat")
res = zeros(64000, 20);
for idx = 25:44
    filename1 = "~/Desktop/Research/Rodu/NPGC/Application/Data/crcns-ac1/wehr/Results/20020608-mw-002/20020608-mw-002-";
    filename = strcat(strcat(filename1, sprintf("%03d", idx)), ".mat");
    load(filename);
    res(:, idx - 24) = double(response.trace) * response.scale;
end

time_res = transpose((1:length(response.trace)) / response.sf);
time_stim = transpose((1:length(stimulus.samples)) / triggers.param.sf);
stim_new = resample(stimulus.samples, response.sf, triggers.param.sf);

plot(time_stim, stimulus.samples);
hold on
plot(time_res(1:60001), stim_new);

data_out = horzcat(time_res(1000:59000), stim_new(1000:59000), res(1000:59000, :));
save("~/Desktop/Research/Rodu/NPGC/Application/dat.mat", "data_out");


