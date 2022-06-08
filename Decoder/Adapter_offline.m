clear
close all

Ns=960;

fs=48e3;
Gi=0;

f_range = [1000, 4000]; %6000


folder_name = strcat('../preamble_exp');
preamble=dlmread('sending_signal/naiser_240.txt');

sounding_file = strcat(folder_name, '/sending_signal');
sending_signal=dlmread(sounding_file)'/30000;

rx_file = strcat(folder_name, '/dist/30-bottom.txt');
recv_dat=dlmread(rx_file)/30000;
recv_dat = recv_dat(8*fs:end);

filter_order = 128;
wn = [(1000-300)/(fs/2), (4000+200)/(fs/2)];    
b = fir1(filter_order, wn, 'bandpass');   
y_after_fir=filter(b,1,recv_dat);
delay_fir = filter_order/2;
recv_dat = y_after_fir(delay_fir+1:end);

figure
plot(recv_dat)
figure
plot(sending_signal)

inc=fs/Ns;
nbin1=round(f_range(1)/inc) + 1;
nbin2=round(f_range(2)/inc) ;
subcarrier_number = nbin2 - nbin1+ 1;
valid_carrier = [];
for i = nbin1:nbin2
   valid_carrier = [valid_carrier, i];
end
f_seq = linspace(0, fs, Ns);

%% the cross correlation of preamble
dat = recv_dat;
[acor,lag]=xcorr(dat,preamble);
[pks,locs,w,p]=findpeaks(acor,'MinPeakHeight',0.5,'MinPeakDistance',24600);
locs=lag(locs);


figure
hold on
plot(lag, acor)
scatter(locs, pks, 'rx')

%% symbol parameter
first_gap = 960;
CP = 67; %Ncs*5;
N_pre = length(preamble);
N_pilot = (Ns+CP)*7;

PN_seq = [1, -1, -1, -1, -1, -1, 1, -1];
L = 8;
Nu = 960;
N0 = 240;

for i = 1:3
    %% begin process the fft data
    now_locs = locs(i);

    pilot_idx = 1 + length(preamble) +  first_gap +CP; 
    idx=now_locs + length(preamble) + first_gap +CP;
    bias = 80;
    pilot_symbol = dat(idx+1-bias:idx+N_pilot-bias);
    
    

    preamble_recv = dat(now_locs - 1200 + 1  :  now_locs + 1200 +length(preamble));
    [begin_idx, max_idx, peak, Mn] = naiser_corr3(preamble_recv, Nu, N0, L, PN_seq);
    new_loc = now_locs - 1200 + begin_idx +240;
    idx2=new_loc+ length(preamble) + first_gap +CP;
    pilot_symbol2 = dat(idx2+1-bias:idx2+N_pilot-bias);
    
    figure
    hold on
    plot(pilot_symbol)
    plot(pilot_symbol2)

    %% sactter plot method for SNR
    pilot_spectrums = [];
    pilot_spectrums2 = [];
    pilot_gts = [];
    
    index1 = 1;
    for i = 1:7
        each_pilot = pilot_symbol(index1 : index1 + Ns - 1);
        each_pilot2 = pilot_symbol2(index1 : index1 + Ns - 1);
        pilot_gt = sending_signal(pilot_idx: pilot_idx + Ns - 1);  
        
        pilot_idx = pilot_idx + Ns+CP;
        index1 = index1 + Ns + CP;
        
        each_fft = fft(each_pilot);
        each_fft2 = fft(each_pilot2);
    %     figure(34)
    %     hold on
    %     plot(f_seq, abs(each_fft))
    %     figure(35)
    %     hold on
    %     plot(each_pilot)
    %     figure(36)
    %     hold on
    %     plot(each_pilot)
    
        each_gt = fft(pilot_gt);
        pilot_spectrums = [pilot_spectrums, each_fft];
        pilot_spectrums2 = [pilot_spectrums2, each_fft2];
        pilot_gts = [pilot_gts, each_gt];
    end
    
    mean_signal_level = mean(abs(pilot_spectrums), 2);   
    mean_gt_level = mean(abs(pilot_gts), 2);  
    
    
    % figure
    % subplot(211)
    % plot(mean_signal_level)
    % subplot(212)
    % plot(mean_gt_level)
    
    snr_bins = snr_calculate(pilot_spectrums, pilot_gts, valid_carrier, f_seq, 0)';
    snr_bins2 = snr_calculate(pilot_spectrums2, pilot_gts, valid_carrier, f_seq, 0)';
    
    [f_begin, f_end, data_rate] = fre_bin_select(snr_bins, 12, f_seq(valid_carrier), 1,fs,0.7);
    [f_begin, f_end, data_rate] = fre_bin_select(snr_bins2, 12, f_seq(valid_carrier), 1,fs,0.7);
    
%     figure
%     hold on
%     % plot(f_seq0(valid_carrier0), snr_bins_naiser)
%     % xlim([1000,5000])
%     plot(f_seq(valid_carrier), snr_bins2)
%     xlim([1000,5000])
%     plot(f_seq(valid_carrier), snr_bins)
%     xlim([1000,5000])
%     plot(f_seq(valid_carrier), recv_SNR)
%     xlim([1000,5000])
%     legend('freq SNR shifted', 'freq SNR')

end