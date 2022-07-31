%% supply tables
clear
cd('C:\Users\YeQ\Documents\MATLAB\Chinese FABIO')
load('FABIO_CHN', 'meta','sup_tab','sup_tab_chn','prod_c','use_tab','proc_est','feed_ani_31');
load opti_results

r = 31;
ts = 29;
p_c = size(sup_tab,1);
s_c = size(sup_tab,2);
prod_c(isnan(prod_c)) = 130;

sup_tab_chn_v = zeros(ts,r,p_c,s_c);  % supple table with certain values
for t = 1:24 % data available for 1990-2013
    for m = 1:r
        temp = squeeze(commod_bal_chn(t,:,m,1));
        temp(isnan(temp)) = 0;
        
        temp_sup = squeeze(sup_tab_chn(t,m,:,:));
        temp_sup(isnan(temp_sup)) = 0;           
        sup_tab_chn_v(t,m,:,:)= diag(temp)*temp_sup;
        temp_check(m,t) = sum(sum(temp,'omitnan'))-sum(sum(diag(temp)*temp_sup,'omitnan')); % value check

        clear temp temp_sup
    end
end

%% Use tables
use_tab_chn_v = zeros(ts,r,p_c,s_c+4);  % use table with certain values; 
% dimension 4: processes + food use + stock + other use+balance
ani_c = [97 99 101:106];  % feeds distributed by animals

for t = 1:24% 1:ts
    for m = 1:r
        temp = squeeze(commod_bal_chn(t,:,m,:));
        temp(isnan(temp)) = 0;
        
        temp_use = zeros(p_c,s_c+4);
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
                temp_use(j_c,s_c+3) = temp(j_c,9);
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
        temp_use(:,s_c+1:end) = temp_use(:,s_c+1:end)+horzcat(temp(:,10),-temp(:,3),temp(:,11),temp_bal);
        use_tab_chn_v(t,m,:,:)= temp_use;
        temp_check_a(m,t) = sum(temp(:,5),'omitnan')-sum(temp(:,3),'omitnan')-sum(sum(temp_use,'omitnan')); % value check

%         horzcat(sum(temp_use,2),sum(temp(:,[3 6:11]),2));
%         ans(:,3) = ans(:,1)-ans(:,2);
        clear temp temp_use
    end
end


%% link trade into supply and use table
clear temp_check temp_check_a
sup_tab_chn_tr = zeros(ts,p_c*r,s_c*r);
use_tab_chn_tr = zeros(ts,p_c*r+p_c,s_c*r); % dimension 2: the last row for imported products
io_tab_chn_Z = zeros(ts,p_c*r+p_c,p_c*r); % dimension 2: the last row for imported products
io_tab_chn_Y = zeros(ts,p_c*r+p_c,4*r+1); % dimension 2: the last row for imported products; dimension 3: the last column for exported prodcuts

for t = 1:24
    temp_sup_tr = zeros(p_c*r,s_c*r);
    temp_use_tr = zeros(p_c*r+p_c,s_c*r);
    temp_y_tr = zeros(p_c*r+p_c,4*r+1);
    
    for m = 1:r
        temp_sup_m = squeeze(sup_tab_chn_v(t,m,:,:));
        temp_use_m = squeeze(use_tab_chn_v(t,m,:,:));
        temp_commod_bal_m = squeeze(commod_bal_chn(t,:,m,:));
        temp_commod_bal_m(isnan(temp_commod_bal_m)) = 0;
        
        temp_sup_tr((m-1)*p_c+(1:p_c),(m-1)*s_c+(1:s_c)) = temp_sup_m;
        
        for j = 1:p_c
            temp_trade_j = squeeze(Trade_max(t,j,:,:));
            temp_sup_j = temp_trade_j([1:r end],m);
            temp_sup_j(m) = temp_commod_bal_m(j,1)-sum(temp_trade_j(m,[1:r end])); %************* how to confirm the local supply?
            if temp_sup_j(m) <= 0 
                temp_sup_j(m) = temp_commod_bal_m(j,1);
            end
            
            
            for i = 1:s_c+4
                if temp_use_m(j,i) ~=0
                    if sum(temp_sup_j) == 0
                        temp_sup_j(m) = 1;
                    end
                    
                    temp_use_j = temp_use_m(j,i)*temp_sup_j/sum(temp_sup_j);
                    if i <= s_c
                        temp_use_tr(([1:r]-1)*p_c+j,(m-1)*s_c+i) = temp_use_j(1:end-1);
                        temp_use_tr(r*p_c+j,(m-1)*s_c+i) = temp_use_j(end);
                    else
                        temp_y_tr(([1:r]-1)*p_c+j,(m-1)*4+i-s_c) = temp_use_j(1:end-1);
                        temp_y_tr(r*p_c+j,(m-1)*4+i-s_c) =temp_use_j(end);
                    end
                    
                    temp_check_a((t-1)*p_c+j,i) = sum(temp_use_j)-temp_use_m(j,i);

                    clear temp_use_j
                end
            end

            clear temp_trade_j temp_sup_j
        end
        temp_y_tr((m-1)*p_c+(1:p_c),end) = squeeze(Trade_max(t,:,m,end));

        temp_check((t-1)*r+m,:) = [sum(sum(temp_use_m(:,1:s_c)))-sum(sum(temp_use_tr(:,(m-1)*s_c+(1:s_c)))) sum(sum(temp_use_m(:,s_c+1:end)))-sum(sum(temp_y_tr(:,(m-1)*4+(1:4))))];
        
        clear temp_sup_m temp_use_m temp_commod_bal_m
    end
    
    transf = diag(ones(1,r*s_c)./sum(temp_sup_tr));
    transf(isinf(transf)) = 0;
    temp_io_Z = temp_use_tr*transf*temp_sup_tr';
    
    
    sup_tab_chn_tr(t,:,:) = temp_sup_tr;
    use_tab_chn_tr(t,:,:) = temp_use_tr;
    io_tab_chn_Z(t,:,:) = temp_io_Z;
    io_tab_chn_Y(t,:,:) = temp_y_tr;
    clear temp_sup_tr temp_use_tr temp_io_Z transf 
end

cd('C:\Users\YeQ\Documents\MATLAB\Chinese FABIO')
save('IOT_physical','sup_tab_chn_tr','sup_tab_chn_v','use_tab_chn_tr','use_tab_chn_v','io_tab_chn_Z','io_tab_chn_Y')

%% check the balance for any commodity

i = 46; % product
m = 23; % region
t = 2012-1989; % year 2012

temp_cbs = squeeze(commod_bal_chn(t,i,m,:));


temp_prod = sum(sup_tab_chn_tr(t,(m-1)*p_c+i,:));

temp_use_do = sum(use_tab_chn_tr(t,(m-1)*p_c+i,(m-1)*s_c+(1:s_c)));
temp_imp_do = sum(sum(use_tab_chn_tr(t,([1:r]-1)*p_c+i,(m-1)*s_c+(1:s_c))))-temp_use_do;
temp_imp_int = sum(use_tab_chn_tr(t,r*p_c+i,(m-1)*s_c+(1:s_c)));

temp_exp_do = sum(use_tab_chn_tr(t,(m-1)*p_c+i,:))-temp_use_do;

temp_fd_do = sum(io_tab_chn_Y(t,(m-1)*p_c+i,(m-1)*4+(1:4)));
temp_fd_exp = sum(io_tab_chn_Y(t,(m-1)*p_c+i,1:end-1))-temp_fd_do;
temp_fd_imp_do = sum(sum(io_tab_chn_Y(t,([1:r]-1)*p_c+i,(m-1)*4+(1:4))))-temp_fd_do;
temp_fd_imp_int = sum(sum(io_tab_chn_Y(t,r*p_c+i,(m-1)*4+(1:4))));

temp_exp_int = sum(io_tab_chn_Y(t,(m-1)*p_c+i,end));

temp_prod+temp_imp_do+temp_imp_int-temp_exp_do-temp_exp_int+temp_fd_imp_do+temp_fd_imp_int-temp_fd_exp

temp_use_do+temp_imp_do+temp_imp_int+temp_fd_do+temp_fd_imp_int+temp_fd_imp_do

squeeze(sum(commod_bal_chn(2012-1989,:,:,:),3));


%% heatmap figure
clear 
cd('C:\Users\YeQ\Documents\MATLAB\Chinese FABIO')
load('IOT_physical','io_tab_chn_Z','io_tab_chn_Y')

r = 31;
s = 84;
t = 23;

for m = 1:r
    reg{m} = sprintf('r%d',m);
end

for j = 1:s
    prod{j} = sprintf('c%d',j);
end

temp_z = squeeze(io_tab_chn_Z(t,:,:));
temp_z(:,r*s+1) = squeeze(io_tab_chn_Y(t,:,end));

temp_rr = zeros(r+1,r+1);
temp_ss = zeros(s,s);
for i = 1:r
    for m = 1:r
        temp_rr(i,m) = sum(sum(temp_z((i-1)*s+(1:s),(m-1)*s+(1:s))));
        temp_ss = temp_ss+temp_z((i-1)*s+(1:s),(m-1)*s+(1:s));
    end
    
    temp_rr(i,end) = sum(temp_z((i-1)*s+(1:s),end));
    temp_rr(end,i) = sum(temp_z(end,(i-1)*s+(1:s)));
end

fs = 20;
lw = 1;
figure
reg{r+1} = 'Export';
x_v = reg;
reg{r+1} = 'Import';
y_v = reg;
h = heatmap(x_v,y_v,log(temp_rr));
h.Title = 'log(Z^r^r)';
h.XLabel = 'Province';
h.YLabel = 'Province';
h.FontSize=fs;

cd('C:\Users\YeQ\Documents\MATLAB\Chinese FABIO')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 11 10];
print('Heatmap_Z_rr','-dpng','-r300')



figure
h = heatmap(prod,prod,log(temp_ss));
h.Title = 'log(Z_c_c)';
h.XLabel = 'Commodity';
h.YLabel = 'Commodity';
h.FontSize=fs;

cd('C:\Users\YeQ\Documents\MATLAB\Chinese FABIO')
fig = gcf;
fig.PaperUnits = 'inches';
fig.PaperPosition = [0 0 11 10];
print('Heatmap_Z_cc','-dpng','-r300')

%% 
load MRIOT_2012_42sectors

temp_rr_m = zeros(r,r);
for m = 1:r
    for n = 1:r
        temp_rr_m(m,n) = sum(IO.Z((m-1)*42+1,(n-1)*42+(1:42)));
        
        if temp_rr(m,n) ~= 0
            tr_mn(m,n) = 1;
        end
    end
end

temp_rr_m_tr = temp_rr_m -diag(diag(temp_rr_m));
temp_rr_tr = temp_rr-diag(diag(temp_rr));
