function [s_segment] = smooth_coastline(segment,ispoly,window,plot_on)
% smooth polygons and coastline by applying window pt moving average 
% kjr apr. 2017
[rows, ~] = find(isnan(segment(:,1))); % segments and polys are padded with NaNs
if rows(end) == length(segment)
    rows(end) = [];
end
for i = 1 : length(rows)
    
    if i == length(rows)
        iseg(:,1) = segment(rows(i)+1 : end,1);
        iseg(:,2) = segment(rows(i)+1 : end,2);
    else
        iseg(:,1) = segment(rows(i)+1 : rows(i+1) -1,1);
        iseg(:,2) = segment(rows(i)+1 : rows(i+1) -1,2);
    end
    % pad with NaNs
    s_segment(rows(i),1:2) = NaN;
    if i ~= length(rows)
        s_segment(rows(i+1),1:2)= NaN;
        
        s_segment(rows(i)+1:rows(i+1)-1,1)=smooth(iseg(:,1),window); %apply moving average
        s_segment(rows(i)+1:rows(i+1)-1,2)=smooth(iseg(:,2),window);
        
        if(ispoly==1) 
          s_segment(rows(i+1)-1,1) = s_segment(rows(i)+1,1); % first point equals last
          s_segment(rows(i+1)-1,2) = s_segment(rows(i)+1,2); 
        end
    else % do not pad end with NaNs..that will occur in General_DistMesh_FP2.m
        
        s_segment(rows(i)+1:length(segment),1)=smooth(iseg(:,1),window); 
        s_segment(rows(i)+1:length(segment),2)=smooth(iseg(:,2),window);
        
        if(ispoly==1) 
           s_segment(length(segment),1) = s_segment(rows(i)+1,1); %first point equals last 
           s_segment(length(segment),2) = s_segment(rows(i)+1,2);
        end
    end

    if plot_on == 1
        %plot(iseg(:,1),iseg(:,2),'r-x'); hold on;
        plot(s_segment(:,1),s_segment(:,2),'b-');
    end
    iseg = [] ;
end
end