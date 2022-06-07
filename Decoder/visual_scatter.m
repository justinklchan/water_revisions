function bers = visual_scatter( points_gt, points_pred, f_seq, valid_carrier)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    bers = [];
    for i = 1:size(points_gt, 1)
            ps_gt =  points_gt(i, :);
            ps_pred =  points_pred(i, :);

            dat_gt=pskdemod(ps_gt,2);
            dat_pred=pskdemod(ps_pred,2);
            ber = 1- sum(dat_gt==dat_pred)/length(dat_gt);
            bers = [bers ber];
            disp(strcat(int2str(valid_carrier(i)), '-BER for frequency-',   int2str(f_seq(valid_carrier(i))), '-Hz is_', num2str(ber)))
            range = mean(abs(ps_gt));
            ps_gt = ps_gt./range;

            bit0_idx = find(real(ps_gt)>0);
            bit1_idx = find(real(ps_gt)<0);

%             scatter(real(ps_gt(bit0_idx)), imag(ps_gt(bit0_idx)), 'bx');
%             scatter(real(ps_gt(bit1_idx)), imag(ps_gt(bit1_idx)), 'ro');
%             xlim([-1.5 1.5])
%             ylim([-1.5 1.5])
%             title(strcat('ground-truth, freq = ', int2str(f_seq(valid_carrier(i)))));
%             subplot(122)
%             hold on;
            figure(100 + i)
            hold on
            ps_pred = ps_pred./range;
            scatter(real(ps_pred(1)), imag(ps_pred(1)), 'k^')
            scatter(real(ps_pred(bit0_idx)), imag(ps_pred(bit0_idx)), 'bx');
            scatter(real(ps_pred(bit1_idx)), imag(ps_pred(bit1_idx)), 'ro');
            xlim([-1.5 1.5])
            ylim([-1.5 1.5])
            title(strcat('my-decode, freq = ', int2str(f_seq(valid_carrier(i)))));
    end
end

