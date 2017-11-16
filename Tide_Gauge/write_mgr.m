function write_mgr(filename,M,coef_cell,SNRmin,PEmin,Amin)
% write_mgr(filename,coef,SNRmin,PEmin,Amin)
% writes out .mgr file from M information on location and
% coef structure output from ut_solv (appends to end of filename.mgr)
% will only output constituents that satisfy the min parameters specified
% will reconstruct and output seasonal variation for you
% 
fid = fopen([filename '.mgr'],'a');
% write out station info
fprintf(fid,'_______________________________________________________________________________________\n');
fprintf(fid,'Station ID = %s \t: %s\n',M.ID,M.Name);
fprintf(fid,'Station Location \t\t: %fN \t %fE\n',M.Lat,M.Lon);
fprintf(fid,'Start time \t\t\t\t: %s\n',M.Start_Time);
fprintf(fid,'End time \t\t\t\t: %s\n',M.End_Time);
fprintf(fid,'Gapless length (years) \t: %.2f\n',M.GaplessLength);
fprintf(fid,'Number of Constituents \t: %d\n',M.CNUM);
fprintf(fid,'Percent Tidal Variance \t: %.2f\n',M.PTV);
fprintf(fid,'Constituent Name \t\t: Amplitude [m]\t\t Phase [GMT deg/date] \t Percent Energy [%%]\n');
fprintf(fid,'---------------------------------------------------------------------------------------\n');

% Get Total TV
Total_TV = 0;
for cc = coef_cell
   Total_TV = cc{1}.diagn.TVallc + Total_TV; 
end
jj = 0;
% Evaluate Seasonal effects
for cc = coef_cell
    coef = cc{1};
    i = strfind(coef.name,'H1'); i = find(~cellfun('isempty',i));
    j = strfind(coef.name,'H2'); j = find(~cellfun('isempty',j));
    if ~isempty(i) && ~isempty(j)
        S_PE = (coef.diagn.PE(i) + coef.diagn.PE(j))*...
                coef.diagn.TVallc/Total_TV;
        if coef.diagn.SNR(i) ~= Inf && coef.diagn.SNR(j) ~= Inf && ...
           coef.diagn.SNR(i) + coef.diagn.SNR(j) > SNRmin && ...
           round(S_PE,2) >= round(PEmin,2)
            % Reconstructing seasonal signal
            t_eval = datenum('01-01-2017'):datenum('12-31-2017');
            sl_fit = ut_reconstr_plus(t_eval',coef,'cnstit',{'S'});
            S_amp = (max(sl_fit) - min(sl_fit))/2;
            if S_amp > Amin 
                jj = jj + 1;
                [~,st] = max(sl_fit);
                S_phs    = datestr(t_eval(st),'dd-mmm');
                name{jj}  = 'M2 Seasonal Variation';
                C(jj,1)   = S_amp;
                C(jj,2)   = coef.A_ci(i) + coef.A_ci(j);
                C(jj,3)   = 0;
                C(jj,4)   = int32((coef.g_ci(i) + coef.g_ci(j))*365.25/360);
                C(jj,5)   = S_PE;
            end
        end
        break
    end
end
    
% Get the constituents from all coefs
for cc = coef_cell
    coef = cc{1};
    for ii = 1:length(coef.name)
         PEn = coef.diagn.PE(ii)*coef.diagn.TVallc/Total_TV;
         if coef.diagn.SNR(ii) > SNRmin && ...
            round(PEn,2) >= round(PEmin,2) && coef.A(ii) > Amin 
            jj = jj + 1;
            name{jj} = coef.name{ii};
            C(jj,1)  = coef.A(ii);
            C(jj,2)  = coef.A_ci(ii);
            if strcmp(name{jj},'SA')
                % Utide gives SA phase compared to calender year rather
                % than the Vernal equinox
                C(jj,3)  = coef.g(ii);
            else
                C(jj,3)  = coef.g(ii);
            end
            C(jj,4)  = coef.g_ci(ii);
            C(jj,5)  = PEn;
         end
    end
end
% sort in order of PE
[C, s] = sortrows(C,5,'descend'); name = name(s);
% make sure any unique values are deleted
[name, s] = unique(name,'stable'); C = C(s,:);

% now write out
for ii = 1:length(name)
    if length(name{ii}) == 2
        fprintf(fid,['%s \t\t\t\t\t\t %.4f ' ,char(177), ...
                     ' %.4f \t %7.3f ', char(177), ...
                     ' %5.3f \t\t %5.2f\n'],name{ii},C(ii,:));
    elseif length(name{ii}) < 5
        fprintf(fid,['%s \t\t\t\t\t %.4f ' ,char(177), ...
                     ' %.4f \t %7.3f ', char(177), ...
                     ' %5.3f \t\t %5.2f\n'],name{ii},C(ii,:));
                 
    else
        fprintf(fid,['%s \t %.4f ' ,char(177), ...
                     ' %.4f \t  %s ', char(177), ...
                     ' %i \t\t\t %5.2f\n'],name{ii},C(ii,1:2),...
                     S_phs,C(ii,4:5));
    end
end
fclose(fid);
end