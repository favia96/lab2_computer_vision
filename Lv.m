function pixels = Lv(inpic, smooth_flag, var, shape)
    % some difference operators, feel free to use
    sdx = [-1 0 1];
    sdy = sdx';
    cdx = [-0.5 0 0.5];
    cdy = cdx';
    rob_pos_dig = [-1 0; 0 1];
    rob_neg_dig = [0 -1; 1 0];

    dxmask = cdx; % or sdx, or rob_pos_dig
    dymask = cdy; % or sdy, or rob_neg_dig

    if (nargin < 5)
        shape = 'same';
    end

    if (smooth_flag == true) % smoothing in case flag true
        inpic = discgaussfft(inpic, var);
    end

    Lx = filter2(dxmask, inpic, shape);
    Ly = filter2(dymask, inpic, shape);
    pixels = sqrt(Lx.^2 + Ly.^2);

end