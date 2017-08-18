    function [cons_pos] = impoly_snap(new_pos,mainland) 
        [idx,dis]= knnsearch(mainland,new_pos); 

        if(dis <= 0.25) 
            disp('snapping point')
            cons_pos = mainland(idx,:); 
        else
            cons_pos = new_pos; 
        end
    return
    end