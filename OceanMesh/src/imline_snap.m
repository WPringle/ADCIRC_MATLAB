    function [cons_pos] = imline_snap(new_pos,mainland) 
        [idx,dis]= knnsearch(mainland,new_pos); 
        disp('blah here'); 
        dis
        if(dis < 0.01) 
            disp('snapping point')
            dis
            cons_pos = mainland(idx,:); 
        else
            disp('did not snap'); 
            dis
            cons_pos = new_pos; 
        end
    return
    end