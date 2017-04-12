function [vs,ids,ide] = extdom_polygon_fix( etbv, ev, pv, iedbeg, ipsbeg, ipsend )
% Still very lousy
% 
figure(1);
hold off
simpplot(pv,ev);
hold(gca); 

ter = 1 ;

ipb  = ipsend ;
ipst = ipsbeg ;
idcur = iedbeg ;

sk = 1 ;
vs(sk,:) = pv(etbv(ipb,idcur),:) ;
ids(sk) = etbv(ipb,idcur) ;
ide = [];
while ( ter )
    for i = 1:length(ipst)
        ipst_n = ipst(i); idcur_n = idcur(i); 
        ikk_n = etbv(ipst_n,idcur_n); 
        sk_n = sk; vs_n = vs; ids_n = ids; ide_n = ide;
        while length(ipst_n) == 1
            ide_n(sk_n) = idcur_n; 
            [idcur_n,sk_n,vs_n,ids_n,ipst_n] = ...
                         get_next_node(ikk_n,idcur_n,sk_n,vs_n,ids_n);
            ikk_n = etbv(ipst_n,idcur_n) ;
            if ( idcur_n == iedbeg ) 
                if length(ipst) == 1
                    break;
                end
            end
        end
        if length(ipst) == 1
            sk = sk_n; vs = vs_n; ids = ids_n; ide = ide_n;
            idcur = idcur_n; ipst = ipst_n;
        else
            skc{i} = sk_n; vsc{i} = vs_n; idsc{i} = ids_n; 
            idec{i} = ide_n; idcurc{i} = idcur_n; ipstc{i} = ipst_n;
            if i == length(ipst)
                if ids(end) == 183088
                end
                % We need to check which answer will "not" lead 
                % us back to the last node
                keep_j = get_possible_routes;
                if ~isempty(keep_j)
                    if length(keep_j) == 1
                        j = keep_j;
                    else
                        % User select ml
                        plot_current_area
                        ml = input_validate('Choose an option: ', 1:length(keep_j));
                        set(h1,'Visible','off')
                        set(h2,'Visible','off')
                        set(h3,'Visible','off')
                        j = keep_j(ml);
                    end
                else
                    broken = 0;
                    for j = 1:length(ipstc)
                        if isempty(intersect(ids,idsc{j}(sk+1:end)))
                            broken = 1;
                            break;
                        end
                    end
                    if broken == 0 
                        % we wanna just skip out and delete this edge
                        ter = 0; idsc{j} = ids(1);
                        idec{j} = iedbeg; vsc{j} = [NaN NaN];
                    end
                end
                sk = skc{j}; vs = vsc{j}; ids = idsc{j}; ide = idec{j};
                idcur = idcurc{j}; ipst = ipstc{j};
                clear skc vsc idsc idcurc ipstc idec            
            end
        end
        if ( idcur == iedbeg ) 
            ter = 0 ; break;
        end
    end
    
end

    function [idcur,sk,vs,ids,ipst] = get_next_node(ikk,idcur,sk,vs,ids)       
        sk = sk + 1 ;
        vs(sk,:) = pv(ikk,:) ;
        ids(sk) = ikk ;
        itemp = find(~(etbv - ikk)) ; 
        iednext = ceil(itemp/2)  ;
        
        %
        if ~isempty(iednext(find((iednext - idcur))))
            idcur = iednext(find((iednext - idcur))) ; 
            [order,ipst] = find(etbv(:,idcur)' - ikk) ;
            [~,order] = sort(order);
            ipst = ipst(order);
            if length(ipst) > 1
            end
        else
            ipst = 1;
            idcur = iedbeg;
        end
    end

    function keep_j = get_possible_routes
        keep_j = [];
        for j = 1:length(ipstc)
            sk_n = skc{j}; vs_n = vsc{j}; ids_n = idsc{j}; 
            for ii = 1:length(ipstc{j})
                idcur_n = idcurc{j}(ii); ipst_n = ipstc{j}(ii);
                ikk_n = etbv(ipst_n,idcur_n); 
                [~,sk_n,vs_n,ids_n,~] = ...
                     get_next_node(ikk_n,idcur_n,sk_n,vs_n,ids_n);
            end
            if isempty(intersect(ids,ids_n(sk+1:end)))
                keep_j(end+1,:) = j;
            elseif ~isempty(find(idsc{j}(sk+1:end) == ids(1),1))
                keep_j(end+1,:) = j;
                idx = find(idsc{j}(2:end) == ids(1));
                idsc{j} = idsc{j}(1:idx+1);
                skc{j} = idx+1;
                vsc{j} = vsc{j}(1:idx+1,:);
                idcurc{j} = iedbeg;
            end
        end
    end

    function plot_current_area              
%         figure(1);
        idx = knnsearch(pv,pv(ids(end),:),'k',1000);
%         enow = ev(idx,:); 
%         pnow = pv(unique(enow),:);
%         tri = delaunay(pnow(:,1),pnow(:,2));
%         simpplot(pnow,tri)
%         hold on
        % change view to center in on on current decision node 
        axis([min(pv(idx,1)) max(pv(idx,1)) ...
              min(pv(idx,2)) max(pv(idx,2))])
        % plotting recent edge
        h3 = plot(vs(max(1,sk-20):end,1),vs(max(1,sk-20):end,2),'-r');
        for jj = 1:size(keep_j,1)
            h1(jj) = plot(vsc{keep_j(jj)}(sk:end,1),...
                 vsc{keep_j(jj)}(sk:end,2),'-o');
            if length(vsc{keep_j(jj)}) > sk
                h2(jj) = text(vsc{keep_j(jj)}(sk+1,1),...
                         vsc{keep_j(jj)}(sk+1,2),num2str(jj));
            else
                h2(jj) = text(vsc{keep_j(jj)}(sk,1),...
                         vsc{keep_j(jj)}(sk,2),['must choose ' num2str(jj)]);
            end
        end  
    end
end

function x = input_validate(prompt,valid)
backspace = sprintf('\b');
prompt_len = length(prompt);
while 1
    x_str = input(prompt,'s');
    x = str2num(x_str);
    if ismember(x,valid)
        break
    end
    num_bs = prompt_len + length(x_str) + 1;
    fprintf(repmat(backspace,[1 num_bs]))
end
end