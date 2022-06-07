function [begin_idx,max_idx, peak, Mn] = naiser_corr3(signal, Nu, N0, L, PN_seq)
    total_L = length(signal);
    preamble_L = (Nu + N0) *L;
    N_both = Nu + N0;
    Mn = [];
    num = 0;
    if(total_L-preamble_L < 0)
        peak = -1;
        begin_idx = -1;
        return
    end
    for i =  1 : total_L-preamble_L
        seg = signal(i: i+preamble_L - 1);
        Pd = 0;
        Rd = 0;
        
        valid_seg = seg( 1+N0: Nu+N0);
        for k = 1: L -1
            bk = PN_seq(k)*PN_seq(k+1);
            seg1 = seg((k-1)*N_both + N0+1: (k)*N_both);
            seg2 = seg((k)*N_both + N0+1: (k+1)*N_both); 

            temp_P = bk*sum(seg1.*seg2);
            
            Pd = Pd + temp_P;
            if(k == 1)
                Rd = Rd + sum(seg1.*seg1);
                Rd = Rd + sum(seg2.*seg2);
            else
                Rd = Rd + sum(seg2.*seg2);
            end
           
        end
%         Rd = sumsqr(seg)*Nu/N_both;
%         Rd = sumsqr(valid_seg);
        corr = Pd/Rd;
        Mn = [Mn, corr];

        num = num + 1;
    end
    
    [peak, begin_idx] = max(Mn);
    max_idx= begin_idx;
    shoulder= peak*0.8;
    right = -1; 
    left = -1;
    for i = begin_idx:length(Mn)-1
        if(Mn(i) >= shoulder && Mn(i+1) <= shoulder)
            right = i;
            break;
        end
    end

    for i = begin_idx:-1:2
        if(Mn(i) >= shoulder && Mn(i-1) <= shoulder)
            left = i;
            break;
        end
    end


    figure(234)
    hold off
    plot(Mn);
%     hold on
%     scatter(1201, Mn(1201), 'rx')
%     ylim([-1, 1])
    
%     xline(right)
%     xline(left)
%     xline((right+left)/2)
    if(right~=-1 && left~=-1 && peak > 0.45)
        disp(strcat('peak width (80%): ', int2str(right - left)))
        begin_idx = round(0.5*right+0.5*left) ;
    else
        begin_idx = -1;
    end
    

end