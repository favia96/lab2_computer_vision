function [linepar, acc] = houghedgeline(pic, scale, gradmagnthreshold, nrho, ntheta, nlines, verbose)

    curves = extractedge(pic, scale, gradmagnthreshold, 'same'); %first plot
    magnitude = Lv(pic, true, scale, 'same');

    [linepar, acc] = houghline(curves, magnitude, nrho, ntheta, gradmagnthreshold, nlines);
    visual_output(pic, linepar); % shows edge lines --> fundamental

    if verbose == 1
        figure('name', 'gradient magnitude')
        showgrey(magnitude);
    else    
        figure('name', 'gradient magnitude')
        showgrey(magnitude);
        
        %figure('name', 'parameter space')
        %showgrey(acc);
        
        figure('name', 'parameter space with smoothing')
        showgrey(binsepsmoothiter(acc, 0.5, 1));
    end

end