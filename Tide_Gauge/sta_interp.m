function [a_int, p_int] = sta_interp(x,y,a,p,xq,yq)
    %Linear interpolation of z from scatter points
    %check for outliers and remove

    if size(x,1) ~= size(a,1)
       x = x'; y = y';
    end
    % convert to complex number
    c = a.*exp(1i*deg2rad(p));
    
    % find median, standard deviations and absolute differences
    Cmed    = median(c,2,'omitnan');
    Cstd    = std(c,0,2,'omitnan');
    Cdiff   = abs(c - repmat(Cmed,1,size(c,2)));
    
    % Flag where all surrounding nodes return bad value for the station
    bad = zeros(size(x,1),1);
    for ii = 1:size(x,1)
        I = (Cdiff(ii,:) > 2.5*Cstd(ii) | isnan(c(ii,:)));
        if all(I)
            bad(ii) = 1;
        end
    end
    
    % Find values outside of 2.5 standard deviations
    x(Cdiff > 2.5*repmat(Cstd,1,size(c,2)) | isnan(c)) = [];
    y(Cdiff > 2.5*repmat(Cstd,1,size(c,2)) | isnan(c)) = [];
    c(Cdiff > 2.5*repmat(Cstd,1,size(c,2)) | isnan(c)) = [];

    % create scattered interpolant
    F = scatteredInterpolant(x(:),y(:),c(:));
    
    % do interpolation
    c_int = F(xq,yq);
    
    % make c_int of bad nans
    c_int(bad == 1) = NaN;
    
    % get back amp and phase
    a_int = abs(c_int);
    p_int = rad2deg(angle(c_int));
    p_int(p_int < 0) = p_int(p_int < 0) + 360;
end

