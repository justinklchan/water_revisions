clear
% close all

Ns=1920;
N_pilot = Ns*7;
fs=48e3;
Gi=0;
inc=fs/Ns;
f_range = [1000, 4000]; %6000
exp_num = 4;

if Ns==480
    Ncs=34;
elseif Ns==960
    Ncs=67;
elseif Ns==1920
    Ncs=134;
elseif Ns==3840
    Ncs=269;
end

folder_name = '../Air_data3/1920';
preamble=dlmread('sending_signal/naiser.txt')/30000;

sounding_file = strcat(folder_name, '/Alice-Sounding-',int2str(exp_num),'.txt');
sending_signal=dlmread(sounding_file)/30000;
% sending_signal = sending_signal(2*fs + 2400 + 1:end);   

rx_file = strcat(folder_name, '/Bob-Sounding-',int2str(exp_num),'-0-bottom.txt');
recv_dat=dlmread(rx_file)/30000;
rx_snr = strcat(folder_name, '/Bob-SNRs-',int2str(exp_num),'.txt');
recv_SNR=dlmread(rx_snr);

% figure
% plot(preamble)
% figure
% plot(sending_signal)

nbin1=round(f_range(1)/inc) + 1;
nbin2=round(f_range(2)/inc) ;
subcarrier_number = nbin2 - nbin1+ 1;
valid_carrier = [];
for i = nbin1:nbin2
   valid_carrier = [valid_carrier, i];
end
f_seq = linspace(0, fs, Ns);

[f_begin, f_end, data_rate] = fre_bin_select(recv_SNR, 10, f_seq(valid_carrier), 1, fs, 0.6);

%% the cross correlation of preamble
dat = recv_dat;
[acor,lag]=xcorr(dat,preamble);
[pks,locs,w,p]=findpeaks(acor,'MinPeakHeight',1,'MinPeakDistance',2*fs);
locs=lag(locs);


figure
hold on
plot(lag, acor)
scatter(locs, pks, 'rx')

%% symbol parameter
first_gap = 960;
CP = 0; %Ncs*5;
N_pre = length(preamble);

%% begin process the fft data
pilot_idx = 1 + length(preamble) +  first_gap +CP; 
idx=locs(1) + length(preamble) + first_gap + CP;
bias = 10;
pilot_symbol = dat(idx+1-bias:idx+N_pilot-bias);

figure
plot(recv_dat)
figure
plot(pilot_symbol)

pilot_spectrums = [];
pilot_gts = [];

%% sactter plot method for SNR
for i =1:7
    each_pilot = pilot_symbol(1 + (i-1)*Ns : i*Ns);
    pilot_gt = sending_signal(pilot_idx-1: pilot_idx + Ns - 2);  
    each_fft = fft(each_pilot);
%     figure
%     plot(each_pilot)
%     xlim([1000,5000])
    each_gt = fft(pilot_gt);
    pilot_spectrums = [pilot_spectrums, each_fft];
    pilot_gts = [pilot_gts, each_gt];
end

mean_signal_level = mean(abs(pilot_spectrums), 2);   
mean_gt_level = mean(abs(pilot_gts), 2);  
snr_bins = snr_calculate(pilot_spectrums, pilot_gts, valid_carrier, f_seq)';
[f_begin, f_end, data_rate] = fre_bin_select(snr_bins, 10, f_seq(valid_carrier), 1,fs,0.6);
[f_begin, f_end, data_rate]

figure
hold on
plot(f_seq(valid_carrier), snr_bins)
xlim([1000,5000])
plot(f_seq(valid_carrier), recv_SNR)
xlim([1000,5000])
legend('freq SNR', 'SNR in phone process')