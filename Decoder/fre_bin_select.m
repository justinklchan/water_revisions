function [f_begin, f_end, data_rate, new_snr] = fre_bin_select(SNR, threshold, f_seq, neighbor, fs, ratio)
    total_bins = length(SNR);
    group_fft = [];
    f_middle = [];
    for h =0:neighbor:total_bins-1
        group_fft = [group_fft mean(SNR(h+1:h + neighbor))];
        f_middle = [f_middle; f_seq(h+1), f_seq(h+neighbor)];
    end
    
    
    new_bins = length(group_fft);
    group_fft_bk = group_fft;
    total_bins = new_bins;
    
    range = [0, 0];
    for L = new_bins : -1: 1
        incre = mag2db(new_bins/L)/2*ratio;
        best_idx = -1;
        best_vally = -1000;
        for i = 1:new_bins - L + 1
            seg = group_fft(i : i + L -1);
            valley = min(seg) + incre;
            if(valley < threshold )
                continue
            else
                if(valley > best_vally)
                    best_idx = i;
                    best_vally = valley;
                end
            end
        end
        if(best_idx > 0)
            range = [best_idx, best_idx+L-1];
            break
        end
    end

    if(range(1) == 0)
        f_begin = -1;
        f_end = -1;
        disp('channel is too weak to support commnucation')
        return
    end
%     for i = 1:new_bins
%         new_power = group_fft + mag2db(new_bins/total_bins);
%         [min_value, min_idx] = min(new_power);
%         if(min_value > threshold)
%             break;
%         end
%         if(total_bins <= 1)
%             f_begin = -1;
%             f_end = -1;
%             disp('channel is too weak to support commnucation')
%             return
%         else
%             f_middle(min_idx, :) = [];
%             group_fft(min_idx) = [];
%             total_bins = total_bins - 1;
%         end
%     end
    L = range(2) - range(1) + 1;
    new_power = group_fft_bk + mag2db(new_bins/L)/2;
    new_snr = new_power(range(1):range(2));
    figure
    hold on
    plot(f_seq, new_power)
    plot([1000, 4000], [threshold threshold])
    plot([f_seq(range(1)), f_seq(range(1))], [10 30])
    plot([f_seq(range(2)), f_seq(range(2))], [10 30])
    
    
    f_begin = f_middle(range(1), 1);
    f_end = f_middle(range(2), end);
    BW = round(f_end - f_begin + 50);
    title(strcat('new power ', int2str(BW)))
    data_rate = (neighbor*L*fs)/(1920*(1+0.07));
    

end