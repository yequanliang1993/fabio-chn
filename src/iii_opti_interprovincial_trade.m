%% this is the Matlab codes to estimate trade within China for FABIO-CHN

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
p0 = 130; % products in FABIO
s0 = 121; % processes in FABIO

% code for harbour provinces in China: Tianjin(2), Liaoning (6), Shandong (15), Guangdong (19)
harb_pro = [2 6 15 19]; % harbour provinces for international shipping

%% 
% assuming the export only happens in provinces with surplus products, 
% whereas import only happens in provinces with insufficient production

load(path_root+"fabio_chn"+version+".mat")

p_c = size(sup_tab,1); % products in Chinese FABIO
prod_c(isnan(prod_c)) = 130;

commod_bal_chn = zeros(ts,size(sup_tab,1),r,length(element));
Trade_max = zeros(ts,size(sup_tab,1),r+5,r+5);  
% 1 dimension: time serious; 
% 2 dimension: products in China;
% 3 dimension: outsourcing province (incl. export);
% 4 dimension: receiving provinces (incl. import)

% textile products
text_prod = xlsread(path_root+"Data\Production of textile.xlsx",'B3:B33');
text_prod = text_prod/sum(text_prod);

food_demand = xlsread(path_root+"Data\per-capital_consuption_of_products.xlsx");
food_demand(isnan(food_demand)) = 0;

% we assume all the stock variation is regarded as other use
commod_bal(:,:,end) = commod_bal(:,:,end)-commod_bal(:,:,3);
commod_bal(:,:,5) = commod_bal(:,:,5)-commod_bal(:,:,3);
commod_bal(:,:,3) = 0;

%% Inter-provincial/international trade of crops
prob_prod = []; % year; product code; trade type 
index_p_vs_d = zeros(ts,p_c,r);
for t = 1:ts  % year
    for j_c = 1:p_c
        if any([1 2 3] == proc_type(prod_c(j_c,1),1)) && j_c ~= 80 && j_c ~= p_c   % processing type and product with no optimal solutions
            data_fao = squeeze(sum(commod_bal(t,prod_c(j_c,:),:),2,'omitnan'));
            data_fao(find(data_fao<0)) = 0;
            
            data_chn = zeros(length(element),r);
            data_chn(1,:) = squeeze(prod_stat(t,prod_c(j_c,1),:));  % provincial production
            data_chn(2,:) = squeeze(imp_stat(t,prod_c(j_c,1),:));   % provincial import
            data_chn(4,:) = squeeze(exp_stat(t,prod_c(j_c,1),:));   % provincial export
            
            data_chn(6,:) = squeeze(sum(feed_ani_31(t,:,prod_c(j_c,1),:),'omitnan'));    % provincial feeds
            data_chn(7,:) = squeeze(seed_est(t,prod_c(j_c,1),:));   % provincial seed
            r_loss = data_fao(8)/sum(data_fao(1:3));   % loss of product assumed by a fixed fraction of the availability (i.e., production, import,stock variation);
            
            if proc_type(prod_c(j_c,1),1) == 1 || proc_type(prod_c(j_c,1),1) == 3
                r_proc = data_fao(9)/data_fao(5); % we assume processing use as a fix fraction of domestic supply
            elseif proc_type(prod_c(j_c,1),1) == 2
                data_chn(9,:) = proc_est(t,prod_c(j_c,1),:); % processed products for the outputs
            end
            
            if sum(food_demand(j_c,:)) ~= 0
                data_chn(10,:) = food_demand(j_c,:).*pop(t,:);
            else
                data_chn(10,:) = pop(t,:)'.*data_fao(10)/sum(pop(t,:),'omitnan'); % food demand
            end
            r_oth = data_fao(11)/data_fao(5); % we assume other use as a fix fraction of domestic supply
            
            if proc_type(prod_c(j_c,1),1) == 1 || proc_type(prod_c(j_c,1),1) == 3
                r_f_s_f = sum(data_fao([6 7 10]))/data_fao(5); % the total fraction of feed, seed, and food use in domestic supply
            elseif proc_type(prod_c(j_c,1),1) == 2
                r_f_s_f = sum(data_fao([6 7 9 10]))/data_fao(5); % the total fraction of feed, seed, processing and food use in domestic supply
            end
            
            % balancing national total value with values in FAO
            sf_data_chn = data_fao./sum(data_chn,2,'omitnan');
            sf_data_chn(isnan(sf_data_chn)) = 0;
            sf_data_chn(isinf(sf_data_chn)) = 0;
            data_chn = data_chn.*sf_data_chn;
            data_chn(isnan(data_chn)) = 0;
            clear sf_data_chn
            
            if r_f_s_f ~= 0
                if proc_type(prod_c(j_c,1),1) == 1 || proc_type(prod_c(j_c,1),1) == 3
                    data_chn(5,:) = sum(data_chn([6 7 10],:))/r_f_s_f;  % domestic supply
                    data_chn(9,:) = data_chn(5,:)*r_proc; % product for processing
                elseif proc_type(prod_c(j_c,1),1) == 2
                    data_chn(5,:) = sum(data_chn([6 7 9 10],:))/r_f_s_f;  % domestic supply
                end
                
                data_chn(11,:) = data_chn(5,:)*r_oth; % product for other use
                data_chn(8,:) = data_chn(5,:)-sum(data_chn(6:11,:),'omitnan'); % product loss
                
                % balancing national total value with values in FAO
                sf_data_chn = data_fao./sum(data_chn,2,'omitnan');
                sf_data_chn(isnan(sf_data_chn)) = 0;
                sf_data_chn(isinf(sf_data_chn)) = 0;
                data_chn = data_chn.*sf_data_chn;
                data_chn(isnan(data_chn)) = 0;
            elseif r_f_s_f == 0 % product used for neither of feed, seed, processing and food, usually for loss and other use
                temp_p_i_s_e = sum(data_chn(1:3,:))-data_chn(4,:);
                data_chn(5,:) = data_fao(5)*temp_p_i_s_e/sum(temp_p_i_s_e); % we distribute the domestic supply by the production+import+stock-export
                
                if t == 23 && j_c == 61
                    data_chn(5,:) = data_fao(5)*text_prod;
                    clear text_prod
                end
                
                data_chn(9,:) = data_chn(5,:)*r_proc; % product for processing
                data_chn(11,:) = data_chn(5,:)*r_oth; % product for other use
                data_chn(8,:) = data_chn(5,:)-sum(data_chn(6:11,:),'omitnan'); % product loss
                
                % balancing national total value with values in FAO
                sf_data_chn = data_fao./sum(data_chn,2,'omitnan');
                sf_data_chn(isnan(sf_data_chn)) = 0;
                sf_data_chn(isinf(sf_data_chn)) = 0;
                data_chn = data_chn.*sf_data_chn;
                data_chn(isnan(data_chn)) = 0;
                clear sf_data_chn
                clear temp_p_i_s_e
            end
                
            % data available of import and export
            if sum(imp_stat(t,prod_c(j_c,1),:),3,'omitnan') ~= 0 && sum(exp_stat(t,prod_c(j_c,1),:),3,'omitnan') ~= 0
                v_n = r*r;  % number of variables
                tr_type = 1; % both import and export data are available
            elseif sum(imp_stat(t,prod_c(j_c,1),:),3,'omitnan') == 0 && sum(exp_stat(t,prod_c(j_c,1),:),3,'omitnan') ~= 0
                v_n = r*r+r*4;  % number of variables
                tr_type = 2; % export data are available
            elseif sum(imp_stat(t,prod_c(j_c,1),:),3,'omitnan') ~= 0 && sum(exp_stat(t,prod_c(j_c,1),:),3,'omitnan') == 0
                v_n = r*r+r*4;  % number of variables
                tr_type = 3; % import data are available
            elseif sum(imp_stat(t,prod_c(j_c,1),:),3,'omitnan') == 0 && sum(exp_stat(t,prod_c(j_c,1),:),3,'omitnan') == 0
                v_n = r*r+4*r+4*r;  % number of variables
                tr_type = 4; % neither import nor export data are available
            end
                
            % optimization objective: mimimal cost of transportation
            f = zeros(1,v_n);
            for m = 1:r
                f((m-1)*r+(1:r)) = cost(m,:);
                
                if tr_type == 3 || tr_type == 4
                    f((m-1)*4+r*r+(1:4)) = cost(m,harb_pro);  % cost for export
                end
            end
            
            if tr_type == 4
                for m = 1:4
                    f((m-1)*r+r*r+4*r+(1:r)) = cost(harb_pro(m),:);
                end
            end
            
            %constraint 1: positive values
            A = -eye(v_n);
            b = 0.001*ones(v_n,1);
            
            % constrain 2: provincial supply and use equial
            Aeq = zeros(r,v_n);
            beq = zeros(r,1);
            
            temp_rest = data_chn(5,:)-sum(data_chn(1:2,:))+data_chn(4,:); % residual between use and supply
            index_p_d = find(temp_rest<0); % in provinces, product j's production more than demand
            index_d_p = find(temp_rest>0); % in provinces, product j's production less than demand
            index_p_vs_d(t,j_c,index_p_d) = 9999;
            index_p_vs_d(t,j_c,index_d_p) = -9999;
            
            if tr_type == 1 % with export and import
                for m = 1:r
                    if any(index_p_d == m)
                        Aeq(m,(m-1)*r+(1:r)) = 1;       % intra-national export
                        Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                        beq(m) = -temp_rest(m);
                    elseif any(index_d_p == m)
                        Aeq(m,(index_p_d-1)*r+m) = 1;   % intra-national import only
                        Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                        beq(m) = temp_rest(m);
                    end
                end
                
                % constrain 3: national import and export
                Aeq_exp = zeros(1,v_n);
                beq_exp = 0;
                Aeq_imp = zeros(1,v_n);
                beq_imp = 0;
            elseif tr_type == 2    % with export no import
                for m = 1:r
                    if any(index_p_d == m)
                        Aeq(m,(m-1)*r+(1:r)) = 1;       % intra-national export
                        Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                        %                     Aeq(m,(m-1)*4+(1:4)+r*r) = 1;   % international export only, no international import
                        beq(m) = -temp_rest(m);
                    elseif any(index_d_p == m)
                        Aeq(m,(index_p_d-1)*r+m) = 1;   % intra-national import only
                        Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                        Aeq(m,([1:4]-1)*r+m+r*r) = 1;  % international import only, no international export
                        beq(m) = temp_rest(m);
                    end
                end
                
                % constrain 3: national import and export
                Aeq_exp = zeros(1,v_n);
                beq_exp = 0;
                Aeq_imp = zeros(1,v_n);
                Aeq_imp(r*r+1:end) = 1;
                beq_imp = data_fao(2);
                
                if r_loss == 0
                    prob_prod = vertcat(prob_prod,[t,j_c,tr_type]);
                end
            elseif tr_type == 3   % no export with import
                for m = 1:r
                    if any(index_p_d == m)
                        Aeq(m,(m-1)*r+(1:r)) = 1;       % intra-national export
                        Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                        Aeq(m,(m-1)*4+(1:4)+r*r) = 1;   % international export only, no international import
                        beq(m) = -temp_rest(m);
                    elseif any(index_d_p == m)
                        Aeq(m,(index_p_d-1)*r+m) = 1;   % intra-national import only
                        Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                        %                         Aeq(m,([1:4]-1)*r+m+r*r) = 1;  % international import only, no international export
                        beq(m) = temp_rest(m);
                    end
                end
                
                % constrain 3: national import and export
                Aeq_exp = zeros(1,v_n);
                Aeq_exp(r*r+1:end) = 1;
                beq_exp = data_fao(4);
                Aeq_imp = zeros(1,v_n);
                %                 Aeq_imp(r*r+1:end) = 1;
                beq_imp = 0;
                
                if r_loss == 0
                    prob_prod = vertcat(prob_prod,[t,j_c,tr_type]);
                end
            elseif tr_type == 4
                for m = 1:r
                    if any(index_p_d == m)
                        Aeq(m,(m-1)*r+(1:r)) = 1;       % intra-national export
                        Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                        Aeq(m,(m-1)*4+(1:4)+r*r) = 1;   % international export only, no international import
                        beq(m) = -temp_rest(m);
                    elseif any(index_d_p == m)
                        Aeq(m,(index_p_d-1)*r+m) = 1;   % intra-national import only
                        Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                        Aeq(m,([1:4]-1)*r+m+r*r+4*r) = 1;  % international import only, no international export
                        beq(m) = temp_rest(m);
                    end
                end
                
                % constrain 3: national import and export
                Aeq_exp = zeros(1,v_n);
                Aeq_exp(r*r+1:r*r+4*r) = 1;
                beq_exp = data_fao(4);
                Aeq_imp = zeros(1,v_n);
                Aeq_imp(r*r+4*r+1:end) = 1;
                beq_imp = data_fao(2);
            end
            j_c
            temp_sol = linprog(f,A,b,vertcat(Aeq,Aeq_exp,Aeq_imp),vertcat(beq,beq_exp,beq_imp));
            temp_sol(find(temp_sol<0)) = 0;
                
            Trade_max(t,j_c,1:r,1:r) = reshape(temp_sol(1:r*r),r,r)';
            if tr_type == 1
                Trade_max(t,j_c,1:r,end) = data_chn(4,:);
                Trade_max(t,j_c,end,1:r) = data_chn(2,:);

                commod_bal_chn(t,j_c,:,:) = data_chn';
            elseif  tr_type == 2 % with export no import
                Trade_max(t,j_c,1:r,end) = data_chn(4,:);
                Trade_max(t,j_c,r+1:r+4,1:r) = reshape(temp_sol(r*r+1:end),r,4)';
                Trade_max(t,j_c,end,1:r) = sum(reshape(temp_sol(r*r+1:end),r,4)');

                data_chn(2,:) = sum(reshape(temp_sol(r*r+1:end),r,4)');
                commod_bal_chn(t,j_c,:,:) = data_chn';
            elseif tr_type == 3
                Trade_max(t,j_c,1:r,r+1:r+4) = reshape(temp_sol(r*r+1:end),4,r)';
                Trade_max(t,j_c,1:r,end) = sum(reshape(temp_sol(r*r+1:end),4,r)',2);
                Trade_max(t,j_c,end,1:r) = data_chn(2,:);

                data_chn(4,:) = sum(reshape(temp_sol(r*r+1:end),4,r)',2);
                commod_bal_chn(t,j_c,:,:) = data_chn';
            elseif tr_type == 4
                Trade_max(t,j_c,1:r,r+1:r+4) = reshape(temp_sol(r*r+1:r*r+4*r),4,r)';
                Trade_max(t,j_c,1:r,end) = sum(reshape(temp_sol(r*r+1:r*r+4*r),4,r)',2);
                Trade_max(t,j_c,r+1:r+4,1:r) = reshape(temp_sol(r*r+4*r+1:end),r,4)';
                Trade_max(t,j_c,end,1:r) = sum(reshape(temp_sol(r*r+4*r+1:end),r,4)');

                data_chn(2,:) = sum(reshape(temp_sol(r*r+4*r+1:end),r,4)');
                data_chn(4,:) = sum(reshape(temp_sol(r*r+1:r*r+4*r),4,r)',2);
                commod_bal_chn(t,j_c,:,:) = data_chn';
            end
        end
    end
    clear j_c v_n tr_type f A b data_fao r_f_s_f r_loss r_oth r_proc...
        data_chn Aeq beq index_p_d index_d_p Aeq_imp Aeq_exp beq_imp ...
        beq_exp temp_sol temp stock temp_rest
end


% check negative values
for t = 23
    for j_c = 1:p_c
        check_neg_val = squeeze(commod_bal_chn(t,j_c,:,[1,2,4:11]));
        check_neg_val(check_neg_val<0) = 0;
        commod_bal_chn(t,j_c,:,[1,2,4:11]) = check_neg_val;
        clear check_neg_val
    end
end

%% address the problemtic products
unique(prob_prod(:,end)) % only for the trade type 2 i.e., with export data available 
unique(prob_prod(:,2)) % only for product 73 (beef) 74 (mutton) 81 (honey)

for k = 1:size(prob_prod,1)
    t = prob_prod(k,1); % year
    j_c = prob_prod(k,2); % product
        
    data_fao = squeeze(sum(commod_bal(t,prod_c(j_c,:),:),2,'omitnan'));
    
    data_chn = zeros(length(element),r);
    data_chn(1,:) = squeeze(prod_stat(t,prod_c(j_c,1),:));  % provincial production
    data_chn(2,:) = squeeze(imp_stat(t,prod_c(j_c,1),:));   % provincial import
    data_chn(4,:) = squeeze(exp_stat(t,prod_c(j_c,1),:));   % provincial export
    
    data_chn(6,:) = squeeze(sum(feed_ani_31(t,:,prod_c(j_c,1),:),'omitnan'));    % provincial feeds
    data_chn(7,:) = squeeze(seed_est(t,prod_c(j_c,1),:));   % provincial seed
    r_loss = data_fao(8)/sum(data_fao(1:3));   % loss of product assumed by a fixed fraction of the availability (i.e., production, import,stock variation);
    
    if proc_type(prod_c(j_c,1),1) == 1 || proc_type(prod_c(j_c,1),1) == 3
        r_proc = data_fao(9)/data_fao(5); % we assume processing use as a fix fraction of domestic supply
    elseif proc_type(prod_c(j_c,1),1) == 2
        data_chn(9,:) = proc_est(t,prod_c(j_c,1),:); % processed products for the outputs
    end
    
    if sum(food_demand(j_c,:)) ~= 0
        data_chn(10,:) = food_demand(j_c,:).*pop(t,:);
    else
        data_chn(10,:) = pop(t,:)'.*data_fao(10)/sum(pop(t,:),'omitnan'); % food demand
    end
    r_oth = data_fao(11)/data_fao(5); % we assume other use as a fix fraction of domestic supply
    
    if proc_type(prod_c(j_c,1),1) == 1 || proc_type(prod_c(j_c,1),1) == 3
        r_f_s_f = sum(data_fao([6 7 10]))/data_fao(5); % the total fraction of feed, seed, and food use in domestic supply
    elseif proc_type(prod_c(j_c,1),1) == 2
        r_f_s_f = sum(data_fao([6 7 9 10]))/data_fao(5); % the total fraction of feed, seed, processing and food use in domestic supply
    end
    
    % balancing national total value with values in FAO
    sf_data_chn = data_fao./sum(data_chn,2,'omitnan');
    sf_data_chn(isnan(sf_data_chn)) = 0;
    sf_data_chn(isinf(sf_data_chn)) = 0;
    data_chn = data_chn.*sf_data_chn;
    data_chn(isnan(data_chn)) = 0;
    clear sf_data_chn
    
    if proc_type(prod_c(j_c,1),1) == 1 || proc_type(prod_c(j_c,1),1) == 3
        data_chn(5,:) = sum(data_chn([6 7 10],:))/r_f_s_f;  % domestic supply
        data_chn(9,:) = data_chn(5,:)*r_proc; % product for processing
    elseif proc_type(prod_c(j_c,1),1) == 2
        data_chn(5,:) = sum(data_chn([6 7 9 10],:))/r_f_s_f;  % domestic supply
    end
    
    data_chn(11,:) = data_chn(5,:)*r_oth; % product for other use
    data_chn(8,:) = data_chn(5,:)-sum(data_chn(6:11,:),'omitnan'); % product loss
    
    % balancing national total value with values in FAO
    sf_data_chn = data_fao./sum(data_chn,2,'omitnan');
    sf_data_chn(isnan(sf_data_chn)) = 0;
    sf_data_chn(isinf(sf_data_chn)) = 0;
    data_chn = data_chn.*sf_data_chn;
    data_chn(isnan(data_chn)) = 0;
    clear sf_data_chn
    
    % data available of import and export
    v_n = r*r+4*r;  % number of variables
    tr_type = 2; % export data are available
       
    % optimization objective: mimimal cost of transportation
    f = zeros(1,v_n);
    for m = 1:r
        f((m-1)*r+(1:r)) = cost(m,:);
    end
    for m = 1:4 % harbour provinces to 31 provinces
        f((m-1)*r+r*r+(1:r)) = cost(harb_pro(m),:);
    end
    
    %constraint 1: positive values
    A = -eye(v_n);
    b = 0.001*ones(v_n,1);
    
    % constrain 2: provincial supply and use equial
    Aeq = zeros(r,v_n);
    beq = zeros(r,1);
    
    % address stock variation
    if data_fao(3) == 0
        data_chn(3,:) = zeros(1,r);
    else
        data_chn(1,:) = data_chn(1,:)*(data_fao(1)+data_fao(3))/data_fao(1); % assume stock is part of production
    end
        
    index_p_d = []; % in provinces, product j's production more than demand
    index_d_p = []; % in provinces, product j's production less than demand
    for m = 1:r
        if sum(data_chn(1:3,m))-data_chn(4,m) > data_chn(5,m)
            index_p_d = [index_p_d,m];
        elseif sum(data_chn(1:3,m))-data_chn(4,m) < data_chn(5,m)
            index_d_p = [index_d_p,m];
        end
    end
    
    index_p_vs_d(t,j_c,index_p_d) = 9999;
    index_p_vs_d(t,j_c,index_d_p) = -9999;
    
    for m = 1:r
        if any(index_p_d == m)
            Aeq(m,(m-1)*r+(1:r)) = 1;       % intra-national export
            Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
            %                     Aeq(m,(m-1)*4+(1:4)+r*r) = 1;   % international export only, no international import
            beq(m) = sum(data_chn(1:3,m))-data_chn(4,m)-data_chn(5,m);
        elseif any(index_d_p == m)
            Aeq(m,(index_p_d-1)*r+m) = 1;   % intra-national import only
            Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
            Aeq(m,([1:4]-1)*r+m+r*r) = 1;  % international import only, no international export
            beq(m) = data_chn(5,m)-(sum(data_chn(1:3,m))-data_chn(4,m)); 
        end
    end
    
    % constrain 3: national import and export
    Aeq_exp = zeros(1,v_n);
%     Aeq_exp(r*r+1:r*r+4*r) = 1;
    beq_exp = 0;
    Aeq_imp = zeros(1,v_n);
    Aeq_imp(r*r+1:end) = 1;
    beq_imp = data_fao(2);
    
    temp_sol = linprog(f,A,b,vertcat(Aeq,Aeq_exp,Aeq_imp),vertcat(beq,beq_exp,beq_imp));
    temp_sol(find(temp_sol<0)) = 0;
    
    Trade_max(t,j_c,1:r,1:r) = reshape(temp_sol(1:r*r),r,r)';
    Trade_max(t,j_c,1:r,end) = data_chn(4,:);
    Trade_max(t,j_c,r+1:r+4,1:r) = reshape(temp_sol(r*r+1:end),r,4)';
    Trade_max(t,j_c,end,1:r) = sum(reshape(temp_sol(r*r+1:end),r,4)');
        
    data_chn(2,:) = sum(reshape(temp_sol(r*r+1:end),r,4)');
    commod_bal_chn(t,j_c,:,:) = data_chn';
    
    clear j_c v_n tr_type f A b data_fao r_f_s_f r_loss r_oth r_proc...
        data_chn Aeq beq index_p_d index_d_p Aeq_imp Aeq_exp beq_imp beq_exp temp_sol
end

% check negative values
for t = 23
    for j_c = 1:p_c
        check_neg_val = squeeze(commod_bal_chn(t,j_c,:,[1,2,4:11]));
        check_neg_val(find(check_neg_val<0)) = 0;
        commod_bal_chn(t,j_c,:,[1,2,4:11]) = check_neg_val;
        clear check_neg_val
    end
end

%% for 'Hides and skins' (80) and 'Wood fuel' (84)
for t = 1:ts
    for j_c = [80, 84]
        data_fao = squeeze(sum(commod_bal(t,prod_c(j_c,:),:),2,'omitnan'));
        data_fao(find(data_fao<0)) = 0;

        data_chn = zeros(length(element),r);
        data_chn(1,:) = squeeze(prod_stat(t,prod_c(j_c,1),:));  % provincial production
        data_chn(2,:) = squeeze(imp_stat(t,prod_c(j_c,1),:));   % provincial import
        data_chn(4,:) = squeeze(exp_stat(t,prod_c(j_c,1),:));   % provincial export

        data_chn(6,:) = squeeze(sum(feed_ani_31(t,:,prod_c(j_c,1),:),'omitnan'));    % provincial feeds
        data_chn(7,:) = squeeze(seed_est(t,prod_c(j_c,1),:));   % provincial seed

        if proc_type(prod_c(j_c,1),1) == 1 || proc_type(prod_c(j_c,1),1) == 3
            r_proc = data_fao(9)/data_fao(5); % we assume processing use as a fix fraction of domestic supply
        elseif proc_type(prod_c(j_c,1),1) == 2
            data_chn(9,:) = proc_est(t,prod_c(j_c,1),:); % processed products for the outputs
        end

        if sum(food_demand(j_c,:)) ~= 0
            data_chn(10,:) = food_demand(j_c,:).*pop(t,:);
        else
            data_chn(10,:) = pop(t,:)'.*data_fao(10)/sum(pop(t,:),'omitnan'); % food demand
        end
        r_oth = data_fao(11)/data_fao(5); % we assume other use as a fix fraction of domestic supply

        % balancing national total value with values in FAO
        sf_data_chn = data_fao./sum(data_chn,2,'omitnan');
        sf_data_chn(isnan(sf_data_chn)) = 0;
        sf_data_chn(isinf(sf_data_chn)) = 0;
        data_chn = data_chn.*sf_data_chn;
        data_chn(isnan(data_chn)) = 0;
        clear sf_data_chn

        temp_p_i_s_e = sum(data_chn(1:3,:))-data_chn(4,:);
        data_chn(5,:) = data_fao(5)*temp_p_i_s_e/sum(temp_p_i_s_e); % we distribute the domestic supply by the production+import+stock-export

        data_chn(9,:) = data_chn(5,:)*r_proc; % product for processing
        data_chn(11,:) = data_chn(5,:)*r_oth; % product for other use
        data_chn(8,:) = data_chn(5,:)-sum(data_chn(6:11,:),'omitnan'); % product loss

        % balancing national total value with values in FAO
        sf_data_chn = data_fao./sum(data_chn,2,'omitnan');
        sf_data_chn(isnan(sf_data_chn)) = 0;
        sf_data_chn(isinf(sf_data_chn)) = 0;
        data_chn = data_chn.*sf_data_chn;
        data_chn(isnan(data_chn)) = 0;
        clear sf_data_chn
        clear temp_p_i_s_e

        temp_rest = data_chn(5,:)-sum(data_chn(1:2,:))+data_chn(4,:); % residual between use and supply
        data_chn(2,:) = data_fao(2)/(data_fao(2)-data_fao(4))*temp_rest;
        data_chn(4,:) = data_chn(2,:)-temp_rest;

        Trade_max(t,j_c,1:r,end) = data_chn(4,:);
        Trade_max(t,j_c,end,1:r) = data_chn(2,:);
        commod_bal_chn(t,j_c,:,:) = data_chn';

        clear temp_rest  data_fao data_chn
    end
    clear t
end

%% oil and oil cake production
o_c_prod(:,1) = 39:56;  
o_c_prod(:,2) = [13:16 37 17:18 1 4 19 13:16 37 17:19];  
% dimension 1: oil or cake products in Chinese FABIO; 
% dimension 2: original products in Chinese FABIO

for t = 1:ts  % year
    for i = 1:length(o_c_prod) 
        j_c = o_c_prod(i,1);
        v_n = r*r+4*r+r*4;  % number of variables
        
        % optimization objective: mimimal cost of transportation
        f = zeros(1,v_n);
        for m = 1:r
            f((m-1)*r+(1:r)) = cost(m,:);
            f((m-1)*4+r*r+(1:4)) = cost(m,harb_pro);  % cost for export
        end

        for m = 1:4
            f((m-1)*r+r*r+4*r+(1:r)) = cost(harb_pro(m),:);
        end
        
        %constraint 1: positive values
        A = -eye(v_n);
        b = zeros(v_n,1);
            
        
        % constrain 2: provincial supply and use equal
        data_fao = squeeze(sum(commod_bal(t,prod_c(j_c,:),:),2,'omitnan'));
        
        data_chn = zeros(length(element),r);
        if proc_type(prod_c(o_c_prod(i,2),1),1) == 3
            temp_proc_ori = squeeze(commod_bal_chn(t,o_c_prod(i,2),:,9));
            temp_proc_ori_oth = squeeze(proc_est(t,o_c_prod(i,2),:));  % processing requirement for other output, e.g., rice for alcoholic beverages
            temp_proc_ori(isnan(temp_proc_ori)) = 0;
            temp_proc_ori_oth(isnan(temp_proc_ori_oth)) = 0;
            
            data_chn(1,:) = (temp_proc_ori-temp_proc_ori_oth)*tcf_fao(prod_c(j_c,1),2)/100;
        else
            temp_proc_ori = squeeze(commod_bal_chn(t,o_c_prod(i,2),:,9));
            temp_proc_ori(isnan(temp_proc_ori)) = 0;
            
            data_chn(1,:) = temp_proc_ori*tcf_fao(prod_c(j_c,1),2)/100;
        end
        clear temp_proc_ori temp_proc_ori_oth
            
        data_chn(6,:) = squeeze(sum(feed_ani_31(t,:,prod_c(j_c,1),:),'omitnan'));    % provincial feeds
        data_chn(7,:) = squeeze(seed_est(t,prod_c(j_c,1),:));   % provincial seed
        r_loss = data_fao(8)/sum(data_fao(1:3));   % loss of product assumed by a fixed fraction of the availability (i.e., production, import,stock variation);
        r_proc = data_fao(9)/data_fao(5); % we assume processing use as a fix fraction of domestic supply
        if sum(food_demand(j_c,:)) ~= 0
            data_chn(10,:) = food_demand(j_c,:).*pop(t,:);
        else
            data_chn(10,:) = pop(t,:)'.*data_fao(10)/sum(pop(t,:),'omitnan'); % food demand
        end
        r_oth = data_fao(11)/data_fao(5); % we assume other use as a fix fraction of domestic supply
        
        r_f_s_f = sum(data_fao([6 7 10]))/data_fao(5); % the total fraction of feed, seed, and food use in domestic supply

        
        % balancing national total value with values in FAO
        sf_data_chn = data_fao./sum(data_chn,2,'omitnan');
        sf_data_chn(isnan(sf_data_chn)) = 0;
        sf_data_chn(isinf(sf_data_chn)) = 0;
        data_chn = data_chn.*sf_data_chn;
        data_chn(isnan(data_chn)) = 0;
        clear sf_data_chn

        if r_f_s_f ~= 0
            data_chn(5,:) = sum(data_chn([6 7 10],:))/r_f_s_f;  % domestic supply
        else
            data_chn(5,:) = data_fao(5)/31;
        end
        data_chn(9,:) = data_chn(5,:)*r_proc; % product for processing
        data_chn(11,:) = data_chn(5,:)*r_oth; % product for other use
        data_chn(8,:) = data_chn(5,:)-sum(data_chn(6:11,:),'omitnan'); % product loss
        
        if sum(data_chn(1,:)) == 0 && sum(data_chn(5,:)) ~= 0
            data_chn(1,:) = data_chn(5,:)*data_fao(1)/data_fao(5); % assume stock is part of production
        end
        
        % balancing national total value with values in FAO
        sf_data_chn = data_fao./sum(data_chn,2,'omitnan');
        sf_data_chn(isnan(sf_data_chn)) = 0;
        sf_data_chn(isinf(sf_data_chn)) = 0;
        data_chn = data_chn.*sf_data_chn;
        data_chn(isnan(data_chn)) = 0;
        clear sf_data_chn
        
        
        Aeq = zeros(r,v_n);
        beq = zeros(r,1);
        
        temp_rest = data_chn(5,:)-sum(data_chn(1:2,:))+data_chn(4,:); % residual between use and supply
        index_p_d = find(temp_rest<0); % in provinces, product j's production more than demand
        index_d_p = find(temp_rest>0); % in provinces, product j's production less than demand
        index_p_vs_d(t,j_c,index_p_d) = 9999;
        index_p_vs_d(t,j_c,index_d_p) = -9999;
        
        if isempty(index_p_d) == 1 && data_fao(4) ~= 0
            for m = 1:r
                Aeq(m,(m-1)*r+(1:r)) = -1;       % intra-national export
                Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                Aeq(m,(m-1)*4+(1:4)+r*r) = -1;   % international export
                
                Aeq(m,([1:r]-1)*r+m) = 1;   % intra-national import
                Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                Aeq(m,([1:4]-1)*r+m+r*r+4*r) = 1;  
                beq(m) = temp_rest(m);
            end
        elseif isempty(index_d_p) == 1 && data_fao(2) ~= 0
            for m = 1:r
                Aeq(m,(m-1)*r+(1:r)) = -1;       % intra-national export
                Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                Aeq(m,(m-1)*4+(1:4)+r*r) = -1;   % international export
                
                Aeq(m,([1:r]-1)*r+m) = 1;   % intra-national import
                Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                Aeq(m,([1:4]-1)*r+m+r*r+4*r) = 1;  
                beq(m) = temp_rest(m);
            end
        else
            for m = 1:r
                if any(index_p_d == m)
                    Aeq(m,(m-1)*r+(1:r)) = 1;       % intra-national export
                    Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                    Aeq(m,(m-1)*4+(1:4)+r*r) = 1;   % international export only, no international import
                    beq(m) = -temp_rest(m);
                elseif any(index_d_p == m)
                    Aeq(m,(index_p_d-1)*r+m) = 1;   % intra-national import only
                    Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                    Aeq(m,([1:4]-1)*r+m+r*r+4*r) = 1;  % international import only, no international export
                    beq(m) = temp_rest(m);
                end
            end
        end
        
        % constrain 3: national import and export
        Aeq_exp = zeros(1,v_n);
        Aeq_exp(r*r+1:r*r+4*r) = 1;
        beq_exp = data_fao(4);
        Aeq_imp = zeros(1,v_n);
        Aeq_imp(r*r+4*r+1:end) = 1;
        beq_imp = data_fao(2);
        
        temp_sol = linprog(f,A,b,vertcat(Aeq,Aeq_exp,Aeq_imp),vertcat(beq,beq_exp,beq_imp));
        temp_sol(find(temp_sol<0)) = 0;
        
        Trade_max(t,j_c,1:r,1:r) = reshape(temp_sol(1:r*r),r,r)';
        Trade_max(t,j_c,1:r,r+1:r+4) = reshape(temp_sol(r*r+1:r*r+4*r),4,r)';
        Trade_max(t,j_c,1:r,end) = sum(reshape(temp_sol(r*r+1:r*r+4*r),4,r)',2);
        Trade_max(t,j_c,r+1:r+4,1:r) = reshape(temp_sol(r*r+4*r+1:end),r,4)';
        Trade_max(t,j_c,end,1:r) = sum(reshape(temp_sol(r*r+4*r+1:end),r,4)');
        
        data_chn(2,:) = sum(reshape(temp_sol(r*r+4*r+1:end),r,4)');
        data_chn(4,:) = sum(reshape(temp_sol(r*r+1:r*r+4*r),4,r)',2);
        commod_bal_chn(t,j_c,:,:) = data_chn';
  
        clear j_c v_n f A b data_fao r_f_s_f r_loss r_oth r_proc...
            data_chn Aeq beq index_p_d index_d_p Aeq_imp Aeq_exp beq_imp beq_exp temp_sol
    end
end

clear o_c_prod

% check negative values
for t = 23
    for j_c = 1:p_c
        check_neg_val = squeeze(commod_bal_chn(t,j_c,:,[1,2,4:11]));
        check_neg_val(find(check_neg_val<0)) = 0;
        commod_bal_chn(t,j_c,:,[1,2,4:11]) = check_neg_val;
        clear check_neg_val
    end
end


%% Inter-provincial/international trade of live animals to meet the production of livestock
ls_prod(:,1) = 62:69;
ls_prod(:,2) = [73:77 77 77 77];
% dimension 1: live animal product;
% dimension 2: output livestock

j_lso = 77; % code of livestock others
for t =  1:ts  % year
    for i = 1:size(ls_prod,1)  % product codes of live animals in Chinese FABIO
        j_c = ls_prod(i,1);
        if sum(imp_stat(t,prod_c(j_c,1),:),3,'omitnan') ~= 0 && sum(exp_stat(t,prod_c(j_c,1),:),3,'omitnan') ~= 0
            v_n = r*r;  % number of variables
            tr_type = 1; % both import and export data are available
        elseif sum(imp_stat(t,prod_c(j_c,1),:),3,'omitnan') == 0 && sum(exp_stat(t,prod_c(j_c,1),:),3,'omitnan') ~= 0
            v_n = r*r+4*r;  % number of variables
            tr_type = 2; % export data are available
        elseif sum(imp_stat(t,prod_c(j_c,1),:),3,'omitnan') ~= 0 && sum(exp_stat(t,prod_c(j_c,1),:),3,'omitnan') == 0
            v_n = r*r+4*r;  % number of variables
            tr_type = 3; % import data are available
        elseif sum(imp_stat(t,prod_c(j_c,1),:),3,'omitnan') == 0 && sum(exp_stat(t,prod_c(j_c,1),:),3,'omitnan') == 0
            v_n = r*r+4*r+4*r;  % number of variables
            tr_type = 4; % neither import nor export data are available
        end
            
        % optimization objective: mimimal cost of transportation 
        f = zeros(1,v_n);
        for m = 1:r
            f((m-1)*r+(1:r)) = cost(m,:);
            
            if tr_type == 3 || tr_type == 4
                f((m-1)*4+r*r+(1:4)) = cost(m,harb_pro);  % cost for export
            end
        end
        
        if tr_type == 2
            for m = 1:4
                f((m-1)*r+r*r+(1:r)) = cost(harb_pro(m),:); % cost for import
            end
        elseif tr_type == 4
            for m = 1:4
                f((m-1)*r+r*r+4*r+(1:r)) = cost(harb_pro(m),:);
            end
        end
        
        
        %constraint 1: positive values
        A = -eye(v_n);
        b = zeros(v_n,1);
        
        
        data_fao = squeeze(sum(commod_bal(t,prod_c(j_c,:),:),2,'omitnan'));
        data_chn = zeros(length(element),r);
        data_chn(1,:) = squeeze(prod_stat(t,prod_c(j_c,1),:));
        data_chn(2,:) = squeeze(imp_stat(t,prod_c(j_c,1),:));
        data_chn(4,:) = squeeze(exp_stat(t,prod_c(j_c,1),:));
        
        % balancing national total value with values in FAO
        sf_data_chn = data_fao./sum(data_chn,2,'omitnan');
        sf_data_chn(isnan(sf_data_chn)) = 0;
        sf_data_chn(isinf(sf_data_chn)) = 0;
        data_chn = data_chn.*sf_data_chn;
        data_chn(isnan(data_chn)) = 0;
        clear sf_data_chn
       
        % live animal requirement in each province
        meat_ls = find(sup_tab(j_lso,:) == 1);
        if any([5:8] == i)
            temp_prod = squeeze(prod_stat(t,prod_c(ls_prod(i,2),1),:));
            temp_ratio = squeeze(sup_tab_chn(t,:,j_lso,meat_ls(i-4)))';
            temp_prod(isnan(temp_prod)) = 0;
            temp_ratio(isnan(temp_ratio)) = 0;
            
            data_chn(5,:) = temp_prod.*temp_ratio*1000/tcf_fao(prod_c(j_c,1),2);
            clear temp_prod temp_ratio
        else
            data_chn(5,:) = squeeze(prod_stat(t,prod_c(ls_prod(i,2),1),:))*1000/tcf_fao(prod_c(j_c,1),2);
        end
        data_chn(5,:) = data_chn(5,:)*(data_fao(1)+data_fao(2)-data_fao(4))/sum(data_chn(5,:),'omitnan');
        data_chn(9,:) = data_chn(5,:);
        
        % constrain 2: provincial supply and use equial
        Aeq = zeros(r,v_n);
        beq = zeros(r,1);
        
        temp_rest = data_chn(5,:)-sum(data_chn(1:2,:))+data_chn(4,:); % residual between use and supply
        index_p_d = find(temp_rest<0); % in provinces, product j's production more than demand
        index_d_p = find(temp_rest>0); % in provinces, product j's production less than demand
        index_p_vs_d(t,j_c,index_p_d) = 9999;
        index_p_vs_d(t,j_c,index_d_p) = -9999;
        
        if tr_type == 1
            for m = 1:r
                if any(index_p_d == m)
                    Aeq(m,(m-1)*r+(1:r)) = 1;       % intra-national export
                    Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                    %                     Aeq(m,(m-1)*4+(1:4)+r*r) = 1;   % international export only, no international import
                    beq(m) = -temp_rest(m);
                elseif any(index_d_p == m)
                    Aeq(m,(index_p_d-1)*r+m) = 1;   % intra-national import only
                    Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                    %                     Aeq(m,([1:4]-1)*r+m+r*r+4*r) = 1;  % international import only, no international export
                    beq(m) = temp_rest(m);
                end
            end
            
            % constrain 3: national import and export
            Aeq_exp = zeros(1,v_n);
            beq_exp = 0;
            Aeq_imp = zeros(1,v_n);
            beq_imp = 0;
        elseif tr_type == 2
            for m = 1:r
                if any(index_p_d == m)
                    Aeq(m,(m-1)*r+(1:r)) = 1;       % intra-national export
                    Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                    %                     Aeq(m,(m-1)*4+(1:4)+r*r) = 1;   % international export only, no international import
                    beq(m) = -temp_rest(m);
                elseif any(index_d_p == m)
                    Aeq(m,(index_p_d-1)*r+m) = 1;   % intra-national import only
                    Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                    Aeq(m,([1:4]-1)*r+m+r*r) = 1;  % international import only, no international export
                    beq(m) = temp_rest(m);
                end
            end
            
            % constrain 3: national import and export
            Aeq_exp = zeros(1,v_n);
            beq_exp = 0;
            Aeq_imp = zeros(1,v_n);
            Aeq_imp(r*r+1:end) = 1;
            beq_imp = data_fao(2);
        elseif tr_type == 3
            for m = 1:r
                if any(index_p_d == m)
                    Aeq(m,(m-1)*r+(1:r)) = 1;       % intra-national export
                    Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                    Aeq(m,(m-1)*4+(1:4)+r*r) = 1;   % international export only, no international import
                    beq(m) = -temp_rest(m);
                elseif any(index_d_p == m)
                    Aeq(m,(index_p_d-1)*r+m) = 1;   % intra-national import only
                    Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                    %                         Aeq(m,([1:4]-1)*r+m+r*r) = 1;  % international import only, no international export
                    beq(m) = temp_rest(m);
                end
            end
            
            % constrain 3: national import and export
            Aeq_exp = zeros(1,v_n);
            Aeq_exp(r*r+1:end) = 1;
            beq_exp = data_fao(4);
            Aeq_imp = zeros(1,v_n);
            beq_imp = 0;
        elseif tr_type == 4
            for m = 1:r
                if any(index_p_d == m)
                    Aeq(m,(m-1)*r+(1:r)) = 1;       % intra-national export
                    Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                    Aeq(m,(m-1)*4+(1:4)+r*r) = 1;   % international export only, no international import
                    beq(m) = -temp_rest(m);
                elseif any(index_d_p == m)
                    Aeq(m,(index_p_d-1)*r+m) = 1;   % intra-national import only
                    Aeq(m,(m-1)*r+m) = 0; % no trade from region m to region m
                    Aeq(m,([1:4]-1)*r+m+r*r+4*r) = 1;  % international import only, no international export
                    beq(m) = temp_rest(m);
                end
            end
            
            % constrain 3: national import and export
            Aeq_exp = zeros(1,v_n);
            Aeq_exp(r*r+1:r*r+4*r) = 1;
            beq_exp = data_fao(4);
            Aeq_imp = zeros(1,v_n);
            Aeq_imp(r*r+4*r+1:end) = 1;
            beq_imp = data_fao(2);
        end
        
        temp_sol = linprog(f,A,b,vertcat(Aeq,Aeq_exp,Aeq_imp),vertcat(beq,beq_exp,beq_imp));
        temp_sol(find(temp_sol<0)) = 0;
              
        Trade_max(t,j_c,1:r,1:r) = reshape(temp_sol(1:r*r),r,r)';
        if tr_type == 1 
            Trade_max(t,j_c,1:r,end) = data_chn(4,:);
            Trade_max(t,j_c,end,1:r) = data_chn(2,:);
            
            commod_bal_chn(t,j_c,:,:) = data_chn';
        elseif tr_type == 2
            Trade_max(t,j_c,1:r,end) = data_chn(4,:);
            Trade_max(t,j_c,r+1:r+4,1:r) = reshape(temp_sol(r*r+1:end),r,4)';
            Trade_max(t,j_c,end,1:r) = sum(reshape(temp_sol(r*r+1:end),r,4)');
            
            data_chn(2,:) = sum(reshape(temp_sol(r*r+1:end),r,4)');
            commod_bal_chn(t,j_c,:,:) = data_chn';
        elseif tr_type == 3
            Trade_max(t,j_c,1:r,r+1:r+4) = reshape(temp_sol(r*r+1:end),4,r)';
            Trade_max(t,j_c,1:r,end) = sum(reshape(temp_sol(r*r+1:end),4,r)',2);
            Trade_max(t,j_c,end,1:r) = data_chn(2,:);
            
            data_chn(4,:) = sum(reshape(temp_sol(r*r+1:end),4,r)',2);
            commod_bal_chn(t,j_c,:,:) = data_chn';
        elseif tr_type == 4
            Trade_max(t,j_c,1:r,r+1:r+4) = reshape(temp_sol(r*r+1:r*r+4*r),4,r)';
            Trade_max(t,j_c,1:r,end) = sum(reshape(temp_sol(r*r+1:r*r+4*r),4,r)',2);
            Trade_max(t,j_c,r+1:r+4,1:r) = reshape(temp_sol(r*r+4*r+1:end),r,4)';
            Trade_max(t,j_c,end,1:r) = sum(reshape(temp_sol(r*r+4*r+1:end),r,4)');
            
            data_chn(2,:) = sum(reshape(temp_sol(r*r+4*r+1:end),r,4)');
            data_chn(4,:) = sum(reshape(temp_sol(r*r+1:r*r+4*r),4,r)',2);
            commod_bal_chn(t,j_c,:,:) = data_chn';
        end
                
        clear j_c v_n tr_type f A b data_fao data_chn meat_ls Aeq beq index_p_d index_d_p Aeq_imp Aeq_exp beq_imp beq_exp temp_sol
    end
end


tr_mn_all = {};
for j_c = 1:84
    tr_mn_all{j_c,1} = squeeze(Trade_max(23,j_c,:,:));
end

% check negative values
for t = 23
    for j_c = 1:p_c
        check_neg_val = squeeze(commod_bal_chn(t,j_c,:,[1,2,4:11]));
        check_neg_val(find(check_neg_val<0)) = 0;
        commod_bal_chn(t,j_c,:,[1,2,4:11]) = check_neg_val;
        clear check_neg_val
    end
end

%% save the optimized data, with high uncertainty
output_name = path_root+"opti_results"+version+".mat";
save(output_name,'commod_bal_chn','Trade_max')

output_name = path_root+"province_exp_imp"+version+".mat";
save(output_name,'index_p_vs_d');

%% adjust inter-provincial trade of main agri-food products
% Note:
% Dimensions of Trade_max
%      1st dimension: years 1990-2013
%      2nd dimension: FABIO-CHN products
%      3nd dimension: 31 provinces + 4 harbor provinces + total import
%      4th dimension: 31 provinces + 4 harbor provinces + total export

load(path_root+"trade_structure_benchmarking"+version+".mat")

% load main crops for China
main_crops = xlsread(path_root+"Data\main_crops.xlsx",'d2:D85');
main_crops(isnan(main_crops)) = 0;

year0 = 1989;
for t =1:ts
    if t+year0 <= 2008
        tr_str = tr_str_2007;
    elseif t+year0 > 2008 && t+year0<=2010
        tr_str = tr_str_2010;
    elseif t+year0>2010
        tr_str = tr_str_2012;
    end
    
    exp_str = tr_str.exp_str;
    imp_str = tr_str.imp_str;
    tr_str = tr_str.tr_str;

    for p = 1:p_c
        if main_crops(p) == 1
            data_chn = squeeze(commod_bal_chn(t,p,:,:));
            data_chn(isnan(data_chn)) = 0;
            data_fao = sum(data_chn);
            
            % data available of import and export
            if sum(imp_stat(t,prod_c(p,1),:),3,'omitnan') ~= 0 && sum(exp_stat(t,prod_c(p,1),:),3,'omitnan') ~= 0
                tr_type = 1; % both import and export data are available
            elseif sum(imp_stat(t,prod_c(p,1),:),3,'omitnan') == 0 && sum(exp_stat(t,prod_c(p,1),:),3,'omitnan') ~= 0
                tr_type = 2; % export data are available
            elseif sum(imp_stat(t,prod_c(p,1),:),3,'omitnan') ~= 0 && sum(exp_stat(t,prod_c(p,1),:),3,'omitnan') == 0
                tr_type = 3; % import data are available
            elseif sum(imp_stat(t,prod_c(p,1),:),3,'omitnan') == 0 && sum(exp_stat(t,prod_c(p,1),:),3,'omitnan') == 0
                tr_type = 4; % neither import nor export data are available
            end
            
            
            tr_vol = squeeze(Trade_max(t,p,...
                [1:r size(Trade_max,3)],...
                [1:r size(Trade_max,4)]));
            tr_vol(isnan(tr_vol)) = 0;
            
            tr_ratio = sum(sum(tr_vol))/...
                sum(commod_bal_chn(t,p,:,1));
            
            tr_vol_adj = tr_vol;
            
            if tr_type == 1 || tr_type == 2
                temp_rest = sum(data_chn(:,4:5),2)-sum(data_chn(:,1:3),2);
                tr_vol_adj(1:r,1:r) = sum(abs(temp_rest)).*tr_str(:,1:r);
            else
                temp_rest = sum(data_chn(:,4:5),2)-sum(data_chn(:,1:3),2);
                
                temp_exp = (sum(abs(temp_rest))+sum(data_chn(:,4))).*tr_str(:,end);
                sf_exp = data_fao(4)/sum(temp_exp);
                tr_vol_adj(1:r,end) = temp_exp*sf_exp;
                
                tr_vol_adj(1:r,1:r) = sum(abs(temp_rest))*tr_str(:,1:r);
                clear sf_exp temp_exp
            end
                
            Trade_max(t,p,1:r,1:r) = tr_vol_adj(1:r,1:r);
            Trade_max(t,p,1:r,end) = tr_vol_adj(1:r,end);
            
            for m = 1:r 
                temp_commod_bal = data_chn(m,:);
%                 temp_commod_bal_a = temp_commod_bal;
                
                % import
                temp_commod_bal(2) = tr_vol_adj(end,m);
                
                % export
                temp_commod_bal(4) = tr_vol_adj(m,end);
                
                % stock addition
                temp_commod_bal(3) = temp_commod_bal(5)+...
                    sum(tr_vol_adj(m,:))-temp_commod_bal(1)-...
                    sum(tr_vol_adj(:,m));

                commod_bal_chn(t,p,m,:) = temp_commod_bal;
                
                check_bal_sup_use(m,p) = sum(temp_commod_bal(1:3))+sum(tr_vol_adj(1:r,m))-...
                    sum(temp_commod_bal(4:5))-sum(tr_vol_adj(m,1:r));
                
                clear temp_commod_bal temp_commod_bal_a 
            end
            
            clear tr_u tr_v tr_vol_adj tr_vol_adj1 tr_vol_adj2
            clear index_p_d index_d_p tr_vol
            
        end
    end
    clear imp_str exp_str
end

% check negative values
for t = 23
    for j_c = 1:p_c
        check_neg_val = squeeze(commod_bal_chn(t,j_c,:,[1,2,4:11]));
        check_neg_val(find(check_neg_val<0)) = 0;
        commod_bal_chn(t,j_c,:,[1,2,4:11]) = check_neg_val;
        clear check_neg_val
    end
end   


path_output = path_root+"opti_results_benchmarked"+version+".mat";
save(path_output,'commod_bal_chn','Trade_max')
            

        

