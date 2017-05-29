function segment = smooth_coastline(segment,window,plot_on) %ispoly,
% smooth polygons and coastline by applying window pt moving average 
% kjr apr. 2017
% modified by WJP June 23 2017
Isnan = find(isnan(segment(:,1)));
Isnan = vertcat(0,Isnan); 

if length(Isnan) - 1 == 0
    iseg = segment;
    if iseg(1,1) == iseg(end,1)
        % Is polygon so need pad to make cyclic
        iseg = [iseg(end-floor(window/2):end,:); iseg; ...
                iseg(1:floor(window/2),:)];
        iseg(:,1) = smooth(iseg(:,1),window);
        iseg(:,2) = smooth(iseg(:,2),window);
        segment = iseg(floor(window/2)+1:end-floor(window/2)-1,:);
        segment(end,:) = segment(1,:);
    else
        segment(:,1) = smooth(segment(:,1),window);
        segment(:,2) = smooth(segment(:,2),window);
    end
else
    for i = 1 : length(Isnan) - 1
        iseg = segment(Isnan(i)+1:Isnan(i+1)-1,:);
        if ~isempty(iseg)
            if iseg(1,1) == iseg(end,1)
                % Is polygon so need pad to make cyclic
                iseg = [iseg(end-floor(window/2):end,:); iseg; ...
                        iseg(1:floor(window/2),:)];
                iseg(:,1) = smooth(iseg(:,1),window);
                iseg(:,2) = smooth(iseg(:,2),window);
                segment(Isnan(i)+1:Isnan(i+1)-1,:) = ...
                       iseg(floor(window/2)+1:end-floor(window/2)-1,:);
                segment(Isnan(i+1)-1,:) = segment(Isnan(i)+1,:);
            else
                segment(Isnan(i)+1:Isnan(i+1)-1,1) = smooth(iseg(:,1),window);
                segment(Isnan(i)+1:Isnan(i+1)-1,2) = smooth(iseg(:,2),window); 
            end
        end
    end
end
% [rows, ~] = find(isnan(segment(:,1))); % segments and polys are padded with NaNs
% if rows(end) == length(segment)
%     rows(end) = [];
% end
% for i = 1 : length(rows)
%     
%     if i == length(rows)
%         iseg(:,1) = segment(rows(i)+1 : end,1);
%         iseg(:,2) = segment(rows(i)+1 : end,2);
%     else
%         iseg(:,1) = segment(rows(i)+1 : rows(i+1) -1,1);
%         iseg(:,2) = segment(rows(i)+1 : rows(i+1) -1,2);
%     end
%     % pad with NaNs
%     s_segment(rows(i),1:2) = NaN;
%     if i ~= length(rows)
%         s_segment(rows(i+1),1:2)= NaN;
%         
%         s_segment(rows(i)+1:rows(i+1)-1,1)=smooth(iseg(:,1),window); %apply moving average
%         s_segment(rows(i)+1:rows(i+1)-1,2)=smooth(iseg(:,2),window);
%         
%         if(ispoly==1) 
%           s_segment(rows(i+1)-1,1) = s_segment(rows(i)+1,1); % first point equals last
%           s_segment(rows(i+1)-1,2) = s_segment(rows(i)+1,2); 
%         end
%     else % do not pad end with NaNs..that will occur in General_DistMesh_FP2.m
%         
%         s_segment(rows(i)+1:length(segment),1)=smooth(iseg(:,1),window); 
%         s_segment(rows(i)+1:length(segment),2)=smooth(iseg(:,2),window);
%         
%         if(ispoly==1) 
%            s_segment(length(segment),1) = s_segment(rows(i)+1,1); %first point equals last 
%            s_segment(length(segment),2) = s_segment(rows(i)+1,2);
%         end
%     end

if plot_on == 1
    plot(segment(:,1),segment(:,2));
end
%    iseg = [] ;
%end
end