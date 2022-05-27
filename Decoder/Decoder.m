clear
close all
% clc

%% configuration of OFDM symbols
Ns=9600; %length of OFDM symbol
fs=48e3; % sampling rate
inc=fs/Ns; % frequency spacing
% length of cyclic prefix
if Ns==960
    Ncs=67;
elseif Ns==1920
    Ncs=135;
elseif Ns==4800
    Ncs=336;
elseif Ns==9600
    Ncs=672;
end
CP = Ncs;
sym_len=Ns+Ncs; %length of symbol
first_gap = 960; % zero interval between preamble and OFDM symbol
tap_num = 1200; % length of equalizer

offset = 40 ; % 5m -120 10m -180
fre_offset = 20;
    
% read the raw data
save_name = strcat('../Air_data6/', int2str(Ns));
folder_name = strcat(save_name,'/sync_file/');  % exp_lake

% read the preamble
preamble=dlmread('sending_signal/naiser.txt')/30000;
preamble = preamble';
N_pre = length(preamble);

% configuration of coding 
code_rate = [2 3]; % code rate of convolution code
code_in_differential  = 1; % whether apply the 
mic = 'bottom'; % which mic to use for decodinig
Repeat = 40; % expriment time

bandwisth_pick = [];
Error_packets =0;
Valid_Packets = 0;

Error_bits_time_uncoded = 0;
Error_bits_fre_uncoded = 0;
Total_bits_uncoded = 0;

Error_bits_diff = 0;
Total_bits_diff = 0;

Error_bits_coded = 0;
Total_bits_coded = 0;

% sending data after encoding
encode_data_gt0  =dlmread(strcat('sending_signal/encode_data_', int2str(code_rate(1)), '_', int2str(code_rate(2)), '.txt'));
% sending data before encoding
uncode_data_gt0  =dlmread(strcat('sending_signal/uncode_data_', int2str(code_rate(1)), '_', int2str(code_rate(2)), '.txt'));

visual_debug = 1;

for r = 0
    if(~exist(strcat(folder_name, 'Alice-DataAdapt-',int2str(r),'.txt'), 'file'))
        continue
    end
    %% read the recv raw packets 
    recv_name =strcat('Bob-DataRx-',int2str(r),'-', mic,'.txt');
    rx_file = strcat(strcat(folder_name, recv_name));
    recv_dat=dlmread(rx_file)/30000;
    dat = recv_dat;
    
    %% load the sending configuartion and ground truth for transmitting bits
    sending_file = strcat('Alice-DataAdapt-',int2str(r),'.txt');
    sending_bits_file = strcat('Alice-BitsAdapt-',int2str(r),'.txt');
    padding_bits_file = strcat('Alice-BitsAdapt_Padding-',int2str(r),'.txt');
    fill_order = dlmread(strcat(folder_name, 'Alice-Bit_Fill_Adapt-',int2str(r),'.txt'));
    sending_signal=dlmread(strcat(folder_name, sending_file))/30000;
    sending_bit = dlmread(strcat(folder_name, sending_bits_file));
    padding_bits = dlmread(strcat(folder_name, padding_bits_file));

    if(size(sending_bit, 1) == 1)
        sending_bit = sending_bit';
    end
        

    %% frequency selection choice
    f_seq = linspace(0, fs, Ns);

    freq_name = strcat('Bob-FreqEsts-',int2str(r),'.txt');
    freq_name2 = strcat('Alice-FeedbackFreqs-',int2str(r),'.txt');
    freq_name3 = strcat('Alice-ExactFeedbackFreqs-',int2str(r),'.txt');
    freq_est = dlmread(strcat(folder_name, freq_name));
    freq_send = dlmread(strcat(folder_name, freq_name2));
    freq_exact = dlmread(strcat(folder_name, freq_name3));
    SNR_phone = dlmread(strcat(folder_name, 'Bob-SNRs-',int2str(r),'.txt'));

    if(freq_send(1) == -1 || freq_send(2) == -1)
        disp(strcat('error frequenct and FSK in Packet:', int2str(r)))
        continue;
    end
    
    %% select the valid frequency bin
    f_range = [freq_exact(1), freq_exact(end)]; %600
    bandwisth_pick = [bandwisth_pick f_range(2) - f_range(1) + inc];

    nbin1=round(f_range(1)/inc) + 1;
    nbin2=round(f_range(2)/inc) + 1;
    subcarrier_number = nbin2 - nbin1+ 1;
    valid_carrier = [];
    for i = nbin1:nbin2
        valid_carrier = [valid_carrier, i];
    end
    

    %% load the coding data from the disk
    code_number = length(padding_bits); % ????????
    assert(mod(code_number, subcarrier_number) == 0)
    Nsyms=code_number/subcarrier_number + 1; 
    assert(Nsyms  == length(fill_order) + 1)
    
    uncode_number = length(sending_bit); %floor(length(sending_bit)*code_rate(1)/code_rate(2)); ???????
    assert(sum(fill_order) == uncode_number)
    sending_uncode = uncode_data_gt0(1 : uncode_number*code_rate(1)/code_rate(2));
    sending_encode = sending_bit; %encode_data_gt0(1 : uncode_number);


    %%  ------------------------- begin demodulation ------------------------------

        
    blocklen=length(preamble)+ first_gap + sym_len*(Nsyms); 
    pilots = 1;
    locs = find_chirp(dat,preamble,fs,r,visual_debug);
    
        figure
        hold on
        plot(dat)
        plot([locs, locs], [-1, 1])
    
    idx=locs+1;
    filter_order = 96;
    if(blocklen +filter_order + tap_num + idx > length(dat))
        block=dat(idx:end);
    else
        block=dat(idx:idx+blocklen +filter_order + tap_num);
    end

    %% band-pass filter 
    filter_order = 96;
%     wn = [(f_range(1)-200)/(fs/2), (f_range(2)+200)/(fs/2)];    
    wn = [(1000-200)/(fs/2), (4000+200)/(fs/2)];    
    b = fir1(filter_order, wn, 'bandpass');   
    y_after_fir=filter(b,1,block);
    delay_fir = filter_order/2;
    block_fir = y_after_fir(delay_fir+1:end);
    block = block_fir;
    
    pkg_idx = length(preamble) + first_gap + CP + 1;  % the index for the received symbols
    gt_idx = length(preamble) + first_gap + CP+1;  % the index for the received symbols

    
    if(visual_debug )
        figure
        hold on
        plot(block_fir(1: pkg_idx + 1000))
        plot(sending_signal(1+offset: gt_idx + offset+1000 )/5, 'r--')
    end

    % time equalizer coding
    data_to_decode2 = [];
    data_to_gt2 =[];
    % fre euqalizer coding
    data_to_decode1= [];
    data_to_gt1 = [];
    % differential coding
    differential_bits = [];
    differential_bits_gt = [];

    % to save the FFT demod results of  ground-truth symbol and rx symbol in method2
    points_gt =[];
    points_pred = [];
    points_gt1 =[];
    points_pred1 = [];

    training_tx1 = [];
    training_rx1  = [];
    training_tx = [];
    training_rx = [];
    
    spect_rx =[];
    spect_gt = [];


    %% demodulation for each symbol
    for sym=1:Nsyms
        rx = block(pkg_idx-fre_offset: pkg_idx-fre_offset + Ns - 1);
        gt = sending_signal(gt_idx: gt_idx + Ns - 1);  
        

        rx_fft = fft(rx);
        gt_fft = fft(gt);

        spect_rx = [spect_rx, rx_fft];
        spect_gt = [spect_gt, gt_fft];
       
        %% training process
        if(sym <= pilots)
            %% frequency-domain training data
            training_tx1 = [training_tx1, gt];
            training_rx1 = [training_rx1, rx];
           
            %% time-domain training data
            len_rx = Ns+tap_num-1;
            symbol_tx = sending_signal(gt_idx : gt_idx + Ns - 1);
            symbol_rx=block(pkg_idx-offset : pkg_idx+len_rx-offset-1);
            training_tx = [training_tx, symbol_tx];
            training_rx = [training_rx, symbol_rx];

            if(sym == pilots) 
                %% training of frequency domain
                Hs = [];
                for t = 1:size(training_rx1, 2)
                    [H0, G0] = cir_estimation_fre(training_tx1(:, t), training_rx1(:, t), nbin1, nbin2); %use the training symbol to calculate the phase shifting
                    Hs = [Hs, H0];
                end
                H= mean(Hs, 2);
                G = 1./H;
                
                for t = 1:1:size(training_rx1, 2)
                    % decode
                    [dat_gt1, dat_pred1, sig_gt1, sig_pred1, acc1] = symbol_decode_fft(training_tx1(:, t), training_rx1(:, t), nbin1, nbin2,valid_carrier,G);
                    if(sym == 1)
                        differential_bits = [differential_bits; dat_gt1];
                        differential_bits_gt = [differential_bits_gt; dat_gt1];
                    end
                    % record decode result
                    data_to_decode1 = [data_to_decode1; dat_pred1];
                    data_to_gt1 = [data_to_gt1; dat_gt1];
                    points_gt1 = [points_gt1, sig_gt1];
                    points_pred1 = [points_pred1, sig_pred1];
                end
                training_tx1 = [];
                training_rx1 = [];   

                %% training of time 
                g = equalizer_estimation_multiple(training_tx, training_rx, tap_num);

                if(visual_debug )
                    figure
                    subplot(211)
                    plot(training_tx)
                    subplot(212)
                    plot(training_rx)
                    figure
                    plot(g)
                end
                
                for t = 1:size(training_rx, 2)
                    % decode
                    symbol_pred = equalizer_recover(training_rx(:, t), g);
                    err = immse(training_tx(:, t),symbol_pred);
                    [dat_gt, dat_pred, sig_gt, sig_pred, acc] = symbol_decode(training_tx(:, t), symbol_pred, nbin1, nbin2, valid_carrier);
                    
                    % record decode result
                    data_to_decode2 = [data_to_decode2; dat_pred];
                    data_to_gt2 = [data_to_gt2; dat_gt];
                    points_gt = [points_gt, sig_gt];
                    points_pred = [points_pred, sig_pred];
                end
                % reset the training data cache
                training_tx = [];
                training_rx = [];
            end
        %% equalizer process
        else
            %% equalizer of fre
            [dat_gt1, dat_pred1, sig_gt1, sig_pred1, acc1] = symbol_decode_fft(gt, rx, nbin1, nbin2,valid_carrier,G);
            points_gt1 = [points_gt1, sig_gt1];
            points_pred1 = [points_pred1, sig_pred1];    
            data_to_decode1 = [data_to_decode1; dat_pred1];
            data_to_gt1 = [data_to_gt1; dat_gt1];

            %% equalizer of time
            len_rx = Ns+tap_num-1;
            symbol_gt = sending_signal(gt_idx : gt_idx + Ns - 1); 
            symbol_rx=block(pkg_idx-offset : pkg_idx+len_rx-offset-1);
            symbol_pred = equalizer_recover(symbol_rx, g);
            
            [dat_gt, dat_pred, sig_gt, sig_pred, acc] = symbol_decode(symbol_gt, symbol_pred, nbin1, nbin2, valid_carrier);
            points_gt = [points_gt, sig_gt];
            points_pred = [points_pred, sig_pred];
            data_to_decode2 = [data_to_decode2; dat_pred];
            data_to_gt2 = [data_to_gt2; dat_gt];        
        end


        pkg_idx = pkg_idx + sym_len;
        gt_idx = gt_idx + sym_len;
    end
    
%     visual_scatter( points_gt1, points_pred1, f_seq, valid_carrier)
        %% debug with the scatter plot image
%         visual_scatter( points_gt, points_pred, f_seq, valid_carrier, 100+50*r);
%         visual_scatter( points_gt1, points_pred1, f_seq, valid_carrier, 50*r+125);
%         signal_level = mean(abs(spect_rx), 2);
%         figure
%         plot(f_seq, signal_level)

    if(visual_debug)
        signal_level = mean(abs(spect_rx), 2);
        figure(r+100)
        hold on
        plot(f_seq(1000/inc+1:4000/inc), SNR_phone)
        %xline(f_range(1))
        %xline(f_range(2))
        plot(f_seq, mag2db(signal_level))
        xlim([1000 5000])
        legend('SNR in phones', 'recv')
        title(int2str(r))
        hold off
    end

    differential_based_time = 1;
    data_to_decode = [];
    data_to_gt = [];

    for i = 1:size(points_gt, 2)-1
        if(differential_based_time)
            delta_phase = wrapToPi(phase(points_pred(:, i+1)./points_pred(:, i)));
            delta_phase_gt = wrapToPi(phase(points_gt(:, i+1)./points_gt(:, i)));
        else
            delta_phase = wrapToPi(phase(points_pred1(:, i+1)./points_pred1(:, i)));
            delta_phase_gt = wrapToPi(phase(points_gt1(:, i+1)./points_gt1(:, i)));
        end
       
        decode_bit2 = ((delta_phase >= pi/2) | (delta_phase <= -pi/2));
        truth_bit2 = ((delta_phase_gt >= pi/2) | (delta_phase_gt <= -pi/2));
        differential_bits = [differential_bits; decode_bit2];
        data_to_decode = [data_to_decode; decode_bit2(1:fill_order(i))];
        data_to_gt =[data_to_gt; truth_bit2(1:fill_order(i))];
        differential_bits_gt = [differential_bits_gt; truth_bit2];
    end
    assert(sum(data_to_gt ~= sending_encode) == 0);

    
    %% differential BER
    interleaving = 0;
    error_bit_diff = sum(differential_bits~=differential_bits_gt);
    total_bit_diff = length(differential_bits_gt);
    Error_bits_diff = Error_bits_diff + error_bit_diff;
    Total_bits_diff = Total_bits_diff + total_bit_diff;

    if(visual_debug) 
        title_name = strcat(int2str(r), '---');
        diffs=hmap(reshape(differential_bits, subcarrier_number, Nsyms), reshape(differential_bits_gt, subcarrier_number, Nsyms), title_name, r);
    end
    
    error_bit_coded = 0;
    total_bit_coded = 0;

    if code_in_differential
        
        print_BE = sum(data_to_gt ~= data_to_decode);

        [decode_err2, decode_total2] = viterbi_decode(data_to_decode, code_rate, sending_uncode);

        
        Error_bits_coded = Error_bits_coded + decode_err2;
        Total_bits_coded = Total_bits_coded + decode_total2;

        error_bit_coded = decode_err2;
        total_bit_coded = decode_total2;
    end
    
    %% fre ber
    error_bit_fre = sum(data_to_decode1~=data_to_gt1);
    total_bit_fre = length(data_to_gt1);
    Error_bits_fre_uncoded = Error_bits_fre_uncoded + error_bit_fre;
    Total_bits_uncoded = Total_bits_uncoded + total_bit_fre;
    %% time ber
    error_bit_time = sum(data_to_decode2~=data_to_gt2);
    total_bit_time = length(data_to_gt2);
    Error_bits_time_uncoded = Error_bits_time_uncoded + error_bit_time;

    Valid_Packets = Valid_Packets + 1;
    if(error_bit_coded > 0 )
        Error_packets = Error_packets  +1;
        disp(strcat('Error packets-', int2str(r)))
    end

    if(visual_debug)
        print_info = strcat('Packet ', int2str(r), '; Uncoded Fre error bit: ', int2str(error_bit_fre),  '; Time error bit: ', int2str(error_bit_time), '; Total bit: ', int2str(total_bit_fre));
        disp(print_info)
        print_info = strcat('Uncoded Diff error bit: ', int2str(error_bit_diff), '; Total diff bit: ', int2str(total_bit_diff), '; Coded error bit: ',int2str(error_bit_coded), '; Coded total bit: ', int2str(total_bit_coded) );
        disp(print_info)
    end

end
Valid_Packets
disp('PER')
Error_packets/Valid_Packets
disp('Code')
Error_bits_coded/Total_bits_coded
disp('diff')
Error_bits_diff/Total_bits_diff
disp('time')
Error_bits_time_uncoded/Total_bits_uncoded
disp('Freq')
Error_bits_fre_uncoded/Total_bits_uncoded

bandwisth_pick = bandwisth_pick*code_rate(1)/code_rate(2);
figure
cdfplot(bandwisth_pick)
title('CDF of code rate with conv code')
%xline(500*code_rate(1)/code_rate(2))
%xline(1500*code_rate(1)/code_rate(2))
%xline(2000*code_rate(1)/code_rate(2))
xlim([0 3200])
ylim([0 1])


% save_name2 = strcat('../new_result/lake20/L', int2str(tap_num));
% if(~visual_debug)
% %     dlmwrite(strcat(save_name2,'Bitrate-',exp_num,'.txt'), bandwisth_pick)
%     dlmwrite(strcat(save_name2,'_Packet_err-',exp_num,'.txt'), [Error_packets(1) Valid_Packets])
%     dlmwrite(strcat(save_name2,'_Bit_err_time-',exp_num,'.txt'), [Error_bits_time(1) Total_bits(1)])
%     dlmwrite(strcat(save_name2,'_Bit_err_diff-',exp_num,'.txt'), [Error_bits_diff(1) Total_bits(1)])
%     dlmwrite(strcat(save_name2,'_Bit_err_code-',exp_num,'.txt'), [Error_bits_code(1) Total_bits_code(1)])
% end


