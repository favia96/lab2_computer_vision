function curves = extractedge(inpic, var, threshold, shape)
    if (nargin < 4) 
        shape = 'same';
    end

    lv = Lv(inpic, true, var, shape); % 1st partial derivatives and gradient
    Lvv = Lvvtilde(discgaussfft(inpic, var), shape); % 2nd derivative
    Lvvv = Lvvvtilde(discgaussfft(inpic, var), shape); % 3rd derivative
    Lv_mask = (lv > threshold) - 0.4; % threshold on gradient (1st derivatives)
    Lvvv_mask = (Lvvv < 0) - 0.4; % check negativity of 3rd derivative

    % sort zerocrossings with negative Lvvv and Lv above threshold
    curves = zerocrosscurves(Lvv, Lvvv_mask); % when 2nd derivative == 0 and 3rd < 0
    curves = thresholdcurves(curves, Lv_mask); % thresholding with 1st
    overlaycurves(inpic, curves);

end