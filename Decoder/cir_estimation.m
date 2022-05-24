function [h] = cir_estimation(tx,rx,tap_num)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    if((length(rx)+ tap_num - 1) ~= length(tx))
        disp('warinig tx and rx different size')
        h= [];
        return;
    end
    if(tap_num > length(rx))
        disp('tap number is too large')
        h= [];
        return;
    end

    P = length(rx);
    L= tap_num;
    M = zeros(P, L);

    for i=1:P
        M(i, :) = tx(L+i-1:-1:i);
    end

    h = pinv((M'* M) + 1e-3*eye(L))*M'*rx;

end