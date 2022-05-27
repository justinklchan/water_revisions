function [snr_list] = snr_calculate(fft_rx, fft_tx, valid_carrier, if_norm)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    snr_list = [];

    for i = valid_carrier
        send_sym = fft_tx(i, :);
        send_sym = send_sym./abs(send_sym);
        
        recv_sym = fft_rx(i, :);
        if(if_norm) 
            recv_sym = recv_sym./abs(recv_sym);
        end
        H_i = ( recv_sym*(send_sym'))/(send_sym*(send_sym'));

        recv_valid = H_i*send_sym;
       

        noise_level = recv_sym - recv_valid;
        SNR_i = sumsqr(abs(recv_valid))/sumsqr(abs(noise_level));
        snr_list = [snr_list mag2db(SNR_i)/2];

        %% visualize the scatter plot
%         data = recv_sym./send_sym;
%         x = real(data);
%         y = imag(data);
% 
%         figure(i+10)
%         scatter(x, y)
%         xlim([-2, 2])
%         ylim([-2, 2])


%         bit0_idx = find(real(send_sym)>0);
%         bit1_idx = find(real(send_sym)<0);
%         
% 
%         bit0_x = real(recv_sym(bit0_idx));
%         bit0_y = imag(recv_sym(bit0_idx));
%         bit1_x = real(recv_sym(bit1_idx));
%         bit1_y = imag(recv_sym(bit1_idx));
% 
%         valid0_x = real(1*H_i);
%         valid0_y = imag(1*H_i);
%         valid1_x = real(-1*H_i);
%         valid1_y = imag(-1*H_i);
        

%         figure(i)
%         hold on
%         scatter(bit0_x, bit0_y, 'bx');
%         scatter(bit1_x, bit1_y, 'rx');
% 
%         scatter(valid0_x, valid0_y, 100, 'bo');
%         scatter(valid1_x, valid1_y, 100, 'ro');
%         title(strcat('my-decode, freq = ', int2str(f_seq(i))));

    end

end