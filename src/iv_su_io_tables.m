%% this is the Matlab codes to construct SUTand IOT in FABIO-CHN

% Created: Quanliang Ye
% Date: 24-October-2024
% Email: yequanliang1993@gmail.com
% Version: 2.3.0

% Note:
%      You are more than welcome to comments on the model and codes.

%% Configure default information
clear

% root path
path_root = pwd+"\";

% current version
path_parts = split(path_root,'\');
version ="_"+path_parts{end-1};
clear path_parts

r = 31;  % 31 provinces
ts = 2013-1990+1; % time series, 1990-2013
p_c = 84; % number of products in FABIO-CHN
s_c = 75; % number of processes in FABIO-CHN
%% data balance check

% -----------------------------------------------------------------
% important to make sure which trade data are using, optimized only or
% optimized plus benchmarked
benchmarked = 1;
% ------------------------------------------------------------------
if benchmarked == 0
    load(path_root+"opti_results"+version+".mat");
elseif benchmarked == 1
    load(path_root+"opti_results_benchmarked"+version+".mat");
end

check_bal_sup_use = zeros(ts,31,p_c);
for t = 1:ts % year
    for m = 1:r % province
        for j_c = 1:p_c % product code
            data_chn = squeeze(commod_bal_chn(t,j_c,m,:));
            tr_vol = squeeze(Trade_max(t,j_c,:,:));

            temp_supply = data_chn(1)+... % production
                data_chn(2)+... % international import
                data_chn(3)+... % stock addition
                sum(tr_vol(1:r,m)); % inter-provincial import

            temp_use = data_chn(5)+... % provincial use
                data_chn(4)+... % internatioanl export
                sum(tr_vol(m,1:r)); % 

            check_bal_sup_use(m,j_c) = temp_supply-temp_use;

            clear data_chn tr_vol temp_supply temp_use
        end
    end
end

%% supply tables
load(path_root+"fabio_chn"+version+".mat",...
    'meta','sup_tab','sup_tab_chn','prod_c','use_tab','proc_est','feed_ani_31');

Trade_max(find(Trade_max<0)) = 0;

prod_c(isnan(prod_c)) = 130;

sup_tab_chn_v = zeros(ts,r,p_c,s_c+r+1+1);  
% supple table with certain values
% in the last dimension:    1:s_c representing local production
%                           s_c+1:s_c+r representing inter-provincial import
%                           s_c+r+1 representing international import
%                           s_c+r+1+1 representing stock variation
for t = 1:ts % year
    for m = 1:r
        data_chn = squeeze(commod_bal_chn(t,:,m,:));
        data_chn(isnan(data_chn)) = 0;
        
        temp = squeeze(commod_bal_chn(t,:,m,1));
        temp(isnan(temp)) = 0;
        
        temp_sup = squeeze(sup_tab_chn(t,m,:,:));
        temp_sup(isnan(temp_sup)) = 0;           
        sup_tab_chn_v(t,m,:,1:s_c)= diag(temp)*temp_sup;
        temp_check(m,t) = sum(sum(temp,'omitnan'))-sum(sum(diag(temp)*temp_sup,'omitnan')); % value check
        
        sup_tab_chn_v(t,m,:,s_c+1:s_c+r) = squeeze(Trade_max(t,:,1:r,m));
        sup_tab_chn_v(t,m,:,s_c+1+r) = squeeze(Trade_max(t,:,end,m));
        sup_tab_chn_v(t,m,:,s_c+2+r) = data_chn(:,3);

        clear temp temp_sup data_chn 
    end
end

%% Use tables
use_tab_chn_v = zeros(ts,r,p_c,s_c+3+r+1);  % use table with certain values; 
% dimension 4: processes + food use + other use + balance+trade
ani_c = [97 99 101:106];  % feeds distributed by animals

for t = 1:ts % year
    for m = 1:r % province
        temp = squeeze(commod_bal_chn(t,:,m,:));
        temp(isnan(temp)) = 0;
        
        temp_use = zeros(p_c,s_c+3);
        % dimension 4: processes + food use + other use + balance
        for j_c = 1:p_c
            temp_p2s_1 = find(use_tab(j_c,:) == 1);   % 1 products also used for production (including seed and loss)
            if isempty(temp_p2s_1) == 0
                if j_c == 70 || j_c == 80
                    temp_out = find(sup_tab(j_c,:) == 1);
                    temp_ratio = squeeze(sup_tab_chn(t,m,j_c,temp_out));
                    temp_use(j_c,temp_p2s_1) = sum(temp(j_c,7:8))*temp_ratio;
                    
                    clear temp_out temp_ratio
                else
                    temp_use(j_c,temp_p2s_1) = sum(temp(j_c,7:8));
                end
            end
            
            temp_p2s_2 = find(use_tab(j_c,:) == 2); % 2 for processing
            if isempty(temp_p2s_2) == 0
                if j_c == 1 || j_c == 4
                    temp_out_1 = proc_est(t,j_c,m);
                    temp_out_1(isnan(temp_out_1)) = 0;
                    temp_out_2 = temp(j_c,9) - temp_out_1;
                    temp_out_2(find(temp_out_2 < 0)) = 0;
                    temp_use(j_c,temp_p2s_2) = horzcat(temp_out_2,temp_out_1);
                    
                    clear temp_out_1 temp_out_2
                else
                    temp_use(j_c,temp_p2s_2) = temp(j_c,9);
                end
            elseif isempty(temp_p2s_2) == 1 && temp(j_c,9) ~=0
                temp_use(j_c,s_c+2) = temp(j_c,9);
            end
            
            temp_p2s_3 = find(use_tab(j_c,:) == 3); % 3 for feed
            if isempty(temp_p2s_3) == 0
                temp_feed = squeeze(sum(feed_ani_31(t,ani_c,prod_c(j_c,:),m),3));
                temp_feed(1) = temp_feed(1)+sum(feed_ani_31(t,98,prod_c(j_c,:),m)); % plus feed for buffaloes
                temp_feed(2) = temp_feed(2)+sum(feed_ani_31(t,100,prod_c(j_c,:),m)); % plus feed for goats
                temp_use(j_c,temp_p2s_3) = temp_feed*temp(j_c,6)/sum(temp_feed);
                clear temp_feed
            end
            
            clear temp_p2s_1 temp_p2s_2 temp_p2s_3
        end
        
        temp_use(isnan(temp_use)) = 0;
        temp_bal = temp(:,5)-sum(temp_use,2)-temp(:,10)-temp(:,11);
        temp_use(:,s_c+1:end) = temp_use(:,s_c+1:end)+horzcat(temp(:,10),temp(:,11),temp_bal);
        temp_check_a(m,t) = sum(temp(:,5),'omitnan')-sum(sum(temp_use,'omitnan')); % value check
        
        % combine with intra- and inter-national trade
        temp_use = horzcat(temp_use,squeeze(Trade_max(t,:,m,1:r)),squeeze(Trade_max(t,:,m,end))');
        
        use_tab_chn_v(t,m,:,:)= temp_use;
        clear temp temp_use
    end
end


%% link trade into supply and use table
clear temp_check temp_check_a
sup_tab_chn_tr = zeros(ts,p_c*r,s_c*r+r+1+1);
use_tab_chn_tr = zeros(ts,p_c*r+p_c,s_c*r); % dimension 2: the last row for imported products
io_tab_chn_Z = zeros(ts,p_c*r+p_c,p_c*r); % dimension 2: the last row for imported products
io_tab_chn_Y = zeros(ts,p_c*r+p_c,3*r+1); % dimension 2: the last row for imported products; dimension 3: the last column for exported prodcuts

for t = 1:ts % year
    temp_sup_tr = zeros(p_c*r,s_c*r+r+1+1);
    temp_use_tr = zeros(p_c*r+p_c,s_c*r);
    temp_y_tr = zeros(p_c*r+p_c,3*r+1);
    
    for m = 1:r
        temp_sup_m = squeeze(sup_tab_chn_v(t,m,:,:));
        temp_use_m = squeeze(use_tab_chn_v(t,m,:,:));
        temp_commod_bal_m = squeeze(commod_bal_chn(t,:,m,:));
        temp_commod_bal_m(isnan(temp_commod_bal_m)) = 0;
        
        temp_sup_tr((m-1)*p_c+(1:p_c),(m-1)*s_c+(1:s_c)) = temp_sup_m(:,1:s_c);
        temp_sup_tr((m-1)*p_c+(1:p_c),s_c*r+1:end) = temp_sup_m(:,s_c+1:end);
        
        for j = 1:p_c
            temp_trade_j = squeeze(Trade_max(t,j,:,:));
            temp_sup_j = temp_trade_j([1:r end],m);
            temp_sup_j(m) = temp_commod_bal_m(j,1)-sum(temp_trade_j(m,[1:r end])); %************* how to confirm the local supply?
            if temp_sup_j(m) <= 0 
                temp_sup_j(m) = temp_commod_bal_m(j,1);
            end
            
            for i = 1:s_c+3
                if temp_use_m(j,i) ~=0
                    if sum(temp_sup_j) == 0
                        temp_sup_j(m) = 1;
                    end
                    
                    temp_use_j = temp_use_m(j,i)*temp_sup_j/sum(temp_sup_j);
                    if i <= s_c
                        temp_use_tr(([1:r]-1)*p_c+j,(m-1)*s_c+i) = temp_use_j(1:end-1);
                        temp_use_tr(r*p_c+j,(m-1)*s_c+i) = temp_use_j(end);
                    else
                        temp_y_tr(([1:r]-1)*p_c+j,(m-1)*3+i-s_c) = temp_use_j(1:end-1);
                        temp_y_tr(r*p_c+j,(m-1)*3+i-s_c) =temp_use_j(end);
                    end
                    
                    temp_check_a((t-1)*p_c+j,i) = sum(temp_use_j)-temp_use_m(j,i);

                    clear temp_use_j
                end
            end

            clear temp_trade_j temp_sup_j
        end
        temp_y_tr((m-1)*p_c+(1:p_c),end) = squeeze(Trade_max(t,:,m,end));

        temp_check((t-1)*r+m,:) = [sum(sum(temp_use_m(:,1:s_c)))-sum(sum(temp_use_tr(:,(m-1)*s_c+(1:s_c)))) ...
            sum(sum(temp_use_m(:,s_c+1:s_c+3)))-sum(sum(temp_y_tr(:,(m-1)*3+(1:3))))];
        
        clear temp_sup_m temp_use_m temp_commod_bal_m
    end
    
    transf = diag(ones(1,r*s_c)./sum(temp_sup_tr(:,1:r*s_c)));
    transf(isinf(transf)) = 0;
    temp_io_Z = temp_use_tr*transf*temp_sup_tr(:,1:r*s_c)';
    
    
    sup_tab_chn_tr(t,:,:) = temp_sup_tr;
    use_tab_chn_tr(t,:,:) = temp_use_tr;
    io_tab_chn_Z(t,:,:) = temp_io_Z;
    io_tab_chn_Y(t,:,:) = temp_y_tr;
    clear temp_sup_tr temp_use_tr temp_io_Z transf 
end

%% save data
if benchmarked ==0
    output_name = path_root+"sut_iot_hybrid"+version+".mat";
elseif benchmarked == 1
    output_name = path_root+"sut_iot_hybrid_benchmarked"+version+".mat";
end
save(output_name,'sup_tab_chn_tr','sup_tab_chn_v','use_tab_chn_tr','use_tab_chn_v','io_tab_chn_Z','io_tab_chn_Y','-v7.3')


