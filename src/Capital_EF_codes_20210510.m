%% data preparation
clear
r = 31;
s = 42;
ts = 2017-1994;

cd('C:\Users\YeQ\Documents\MATLAB\Capital_EF_codes_20210510')
c_axp = xlsread('Conordance_with_42_sectors_in_Wang','3 assets to 42 MRIO sectors','C3:Ar5');
c_axp(isnan(c_axp)) = 0;
save('c_axp','c_axp')
c_ixs = xlsread('Conordance_with_42_sectors_in_Wang','37 KLEMS to 42 MRIO sectors','C3:AR39');
c_ixs(isnan(c_ixs)) = 0;
save('c_ixs','c_ixs')

c_c2mrio = xlsread('Conordance_with_42_sectors_in_Wang','45 carbon sectors to 42 MRIO','b3:AQ47');
c_c2mrio(isnan(c_c2mrio)) = 0;
save('c_c2mrio','c_c2mrio')


% value added by sectors
cd('C:\Users\YeQ\Documents\Ph.D UT\Phd in UT\Data China\Annual value added by province\Value_Added_By_Sectors')
sec_n = {'Agri.xls','Industry.xls','Construction.xls','Transportation.xls','Wholesale.xls',...
    'Accommodation.xls','Financial.xls','Real_Estate.xls','Other_Services.xls'};
vd_s = zeros(ts+2,r,length(sec_n));
for i = 1:length(sec_n)
    temp = xlsread(sec_n{i},'C5:AA35')'*10000;
    vd_s(:,:,i) = temp([ts+2:-1:1],:);
    clear temp
end
clear sec_n

% value added by categories
cd('C:\Users\YeQ\Documents\Ph.D UT\Phd in UT\Data China\Annual value added by province\Value_Added_By_Categories')
vd_c = zeros(ts+2,r,4);
for i = 1:4
    temp = xlsread(sprintf('C%d.xls',i),'C5:AA35')'*10000;
    vd_c(:,:,i) = temp([ts+2:-1:1],:);
    clear temp
end

cd('C:\Users\YeQ\Documents\MATLAB\Capital_EF_codes_20210510')
c_mrio2vd =  xlsread('Conordance_with_42_sectors_in_Wang','42 sectors to vd sectors','C2:K43'); % concordance to link value added sectors into 30 MRIO sectors 
c_mrio2vd(isnan(c_mrio2vd)) = 0;

save('vd_info','vd_c','vd_s','c_mrio2vd')


% per-capita expenditure concordance
cd('C:\Users\YeQ\Documents\MATLAB\Capital_EF_codes_20210510')
c_exp = xlsread('Conordance_with_42_sectors_in_Wang','42 sectors to exp types','C2:j43');
c_exp(isnan(c_exp)) = 0;

% per-capita expenditure of rural population, in unit: 10000 capita
cd('C:\Users\YeQ\Documents\Ph.D UT\Phd in UT\Data China\Annual expenditures of urban and rural population')
temp_exp_r = xlsread('Annual expenditures of population','Rural','A3:J777');
exp_r = zeros(ts+2,r,size(c_exp,2)); % per-capita expenditure of rural population
for n = 1:ts+2
    index = find(temp_exp_r(:,1) == n+1994);
    exp_r(n,:,:) = temp_exp_r(index,3:end)/10000; % unit in 10000 Yuan per capita
    clear index
end
clear temp_exp_r

% per-capita expenditure of urban population
temp_exp_u = xlsread('Annual expenditures of population','Urban','A3:J777');
exp_u = zeros(ts+2,r,size(c_exp,2)); % per-capita expenditure of urban population
for n = 1:ts+2
    index = find(temp_exp_u(:,1) == n+1994);
    exp_u(n,:,:) = temp_exp_u(index,3:end)/10000; % unit in 10000 Yuan per capita
    clear index
end
clear temp_exp_u

% final expenditure and gross capital formation by province
cd('C:\Users\YeQ\Documents\Ph.D UT\Phd in UT\Data China\Annual value added by province\Value_Added_By_Final_Demand')
fd_c = {'Expenditure_rural.xls','Expenditure_urban.xls','Expenditure_gov.xls','GFCF.xls','Stock_Changes.xls','Net_Outflow.xls'};
fd_chn = zeros(ts+2,r,6);
for i = 1:6
    temp = xlsread(fd_c{i},'C5:AA35')'*10000;
    fd_chn(:,:,i) = temp([ts+2:-1:1],:);
    clear temp
end

for n = 1:ts+2
    sf = squeeze(sum(vd_c(n,:,:),3))./squeeze(sum(fd_chn(n,:,:),3));
    fd_chn(n,:,:) = squeeze(fd_chn(n,:,:)).*sf';
end
fd_chn(:,:,6) = [];

cd('C:\Users\YeQ\Documents\MATLAB\Capital_EF_codes_20210510')
save('final_expenditure','exp_r','exp_u','c_exp','fd_chn')


% annual export/import by product categories
cd('C:\Users\YeQ\Documents\Ph.D UT\Phd in UT\Data China\Annual import and export by products')
temp_exp = xlsread('Annual import and export by products','National Export by product','A3:AA102'); % export in 10000 Yuan
temp_exp(find(temp_exp(:,1) == 0),:) = [];
temp_imp = xlsread('Annual import and export by products','National Import by products','A3:AA102'); % import in 10000 Yuan
temp_imp(find(temp_imp(:,1) == 0),:) = [];

cd('C:\Users\YeQ\Documents\MATLAB\Capital_EF_codes_20210510')
c_expimp = xlsread('Conordance_with_42_sectors_in_Wang','42 sectors to expt_impt types','C50:o68'); % concordance to change the categories of products
c_expimp(isnan(c_expimp)) = 0;
export = temp_exp(:,3:end)'*c_expimp;
import = temp_imp(:,3:end)'*c_expimp;
clear temp_exp temp_imp

c_mrio2expimp =  xlsread('Conordance_with_42_sectors_in_Wang','42 sectors to expt_impt types','C2:O43'); % concordance to link export/import products into 30 MRIO sectors 
c_mrio2expimp(isnan(c_mrio2expimp)) = 0;

cd('C:\Users\YeQ\Documents\MATLAB\Capital_EF_codes_20210510')
save('expt_impt','export','import','c_mrio2expimp')


% population 1995-2019
cd('C:\Users\YeQ\Documents\Ph.D UT\Phd in UT\Data China\Population')
pop_r = xlsread('Population rural 1990-2018','M5:AK35'); % rural population, in 10000 capita
pop_u = xlsread('Population urban  1990-2018','Q5:AO35'); % urban population, in 10000 capita
pop_to = xlsread('Population total 1986-2018','L5:Aj35'); % total population, in 10000 capita

for n = 1:ts
    if sum(pop_u(:,n),'omitnan') == 0
        pop_u(:,n) = pop_to(:,n)-pop_r(:,n);
    end
end

% population 2020-2040
cd('C:\Users\YeQ\Documents\MATLAB\Capital_model_IOT\Data');
temp = xlsread('Pop_TOTAL.csv');
temp_r = temp(6:end,11:31);
temp_r = temp_r(([1:31]-1)*5+2,:)/10000;

u_r = sum(pop_u(:,2017-1994))/sum(pop_to(:,2017-1994)); % urbanizaton rate in 2019
u_ra = 0.814; % urbanizaton rate in 2040
r_u_r = nthroot(u_ra/u_r,23);

r_r_reg = pop_r(:,2019-1994)./sum(pop_r(:,2019-1994));

for n = 1:21
    r_u_r_a(1,n) = u_r*(r_u_r^n);
end

temp_pop_r = (1-r_u_r_a).*sum(temp_r).*r_r_reg;
temp_pop_u = temp_r-temp_pop_r;

pop_u(:,2020-1994:2040-1994) = temp_pop_u;
pop_r(:,2020-1994:2040-1994) = temp_pop_r;
pop_to(:,2020-1994:2040-1994) = temp_r;


plot(pop_to');

meta_pop.Years = 1995:2040;
meta_pop.notes{1,1} = {'Population between 1995-2019 from NBSC; Population between 2020-2040 from Yidan Chen 2020'};
meta_pop.notes{2,1} = {'Urbanization rates increase from 61.06% in 2019 to 81.4% 2040, with annual growth rate 1.38% per year'};
cd('C:\Users\YeQ\Documents\MATLAB\Capital_model_IOT')
save('pop_chn','pop_r','pop_to','pop_u','meta_pop')

cd('C:\Users\YeQ\Documents\MATLAB\Capital_EF_codes_20210510')
cd('C:\Users\YeQ\Documents\MATLAB\Capital_model_IOT')
save('pop_chn','pop_r','pop_u','pop_to','meta_pop')
%% adjust MRIO tables
clear

load MRIO_1978_2017
load cpi
load meta
load final_expenditure
load expt_impt
load vd_info
load pop_chn
load ki_chn

r = 31;
s = 42;
ts = 2017-1994;

for t = 1:ts+2
    if t<23 % period of 1995-2016
        t0 = t+17;
    else
        t0 = ts+17; % for year 2017, 2018 and 2019
    end
    temp_z = IO{t0,2};
    temp_y = IO{t0,4};
    temp_expt = IO{t0,6};
    temp_expt = temp_expt(1:r*s);
    temp_imp_z = IO{t0,3};
    temp_imp_y = IO{t0,5};
    temp_vd = IO{t0,7};
    temp_x = IO{t0,8};
        
    temp_z(find(temp_z<10^-6)) = 0;
    temp_y(find(temp_y<10^-6)) = 0;
    temp_expt(find(temp_expt<10^-6)) = 0;
    temp_imp_z(find(temp_imp_z<10^-6)) = 0;
    temp_imp_y(find(temp_imp_y<10^-6)) = 0;
    temp_vd(find(temp_vd<10^-6)) = 0;
    temp_x(find(temp_x<10^-6)) = 0;
    
    % adjust processes
    y_t = zeros(size(temp_y)); % estimated domestic final demand in year t+1994
    impt_y_t = zeros(size(temp_imp_y)); % estimated imported final demand in year t+1994
    vd_n = zeros(size(temp_vd)); % estimated value added in year t+1994
    for m = 1:r
        % rural expenditure
        yr_t = reshape(temp_y(:,(m-1)*5+1),s,r); % domestic rural expenditure of province m in year t
        impt_yr_t = temp_imp_y(:,(m-1)*5+1); % imported rural expenditure of province m in year t
        exp_t = squeeze(exp_r(t,m,:))*pop_r(m,t)*10000; % total rural expenditure in year t
        sf = ones(s,1); % sacling factor of rural expenditure
        for i = 1:size(c_exp,2)
            index_i = find(c_exp(:,i) == 1);
            sf(index_i) = exp_t(i)/(sum(sum(yr_t(index_i,:)))+sum(impt_yr_t(index_i)));
            clear index_i
        end
        yr_t_a = yr_t.*sf;
        impt_yr_t_a = impt_yr_t.*sf;
        sf_a = fd_chn(t,m,1)/(sum(sum(yr_t_a))+sum(impt_yr_t_a)); % balanced into year n's statistical final demand
        yr_t_a = yr_t_a*sf_a;
        impt_yr_t_a = impt_yr_t_a*sf_a;
        clear exp_t sf sf_a
        check_y(t,m,1) = sum(sum(yr_t_a))/sum(sum(yr_t));
        
        % urban expenditure
        yu_t = reshape(temp_y(:,(m-1)*5+2),s,r); % domestic urban expenditure of province m in year t
        impt_yu_t = temp_imp_y(:,(m-1)*5+2); % imported urban expenditure of province m in year t
        exp_t = squeeze(exp_u(t,m,:))*pop_u(m,t)*10000; % total urban expenditure in year t
        sf = ones(s,1); % sacling factor of urban expenditure
        for i = 1:size(c_exp,2)
            index_i = find(c_exp(:,i) == 1);
            sf(index_i) = exp_t(i)/(sum(sum(yu_t(index_i,:)))+sum(impt_yu_t(index_i)));
            clear index_i
        end
        yu_t_a = yu_t.*sf;
        impt_yu_t_a = impt_yu_t.*sf;
        sf_a = fd_chn(t,m,2)/(sum(sum(yu_t_a))+sum(impt_yu_t_a)); % balanced into year n's statistical final demand
        yu_t_a = yu_t_a*sf_a;
        impt_yu_t_a = impt_yu_t_a*sf_a;
        clear exp_t sf sf_a
        check_y(t,m,2) = sum(sum(yu_t_a))/sum(sum(yu_t));
        
        % expenditure of government
        yg_t = reshape(temp_y(:,(m-1)*5+3),s,r); % domestic final expenditure of government of province m in year t
        impt_yg_t = temp_imp_y(:,(m-1)*5+3); % imported final expenditure of government of province m in year t
        sf = (yr_t_a+yu_t_a)./(yr_t+yu_t); % scaling factor according to the changes in rural/urban expenditures
        sf(isnan(sf)) = 0;
        sf(isinf(sf)) = 0;
        yg_t_a = yg_t.*sf;
        impt_yg_t_a = impt_yg_t.*(impt_yr_t_a+impt_yu_t_a)./(impt_yr_t+impt_yu_t);
        impt_yg_t_a(isnan(impt_yg_t_a)) = 0;
        sf_a = fd_chn(t,m,3)/(sum(sum(yg_t_a))+sum(impt_yg_t_a)); % balanced into year t's statistical final demand
        yg_t_a = yg_t_a*sf_a;
        impt_yg_t_a = impt_yg_t_a*sf_a;
        clear sf sf_a
        check_y(t,m,3) = sum(sum(yg_t_a))/sum(sum(yg_t));
        
        % gross fixed capital formation
        gcf_t = reshape(temp_y(:,(m-1)*5+4),s,r); % domestic gross capital formation of province m in year t
        impt_gcf_t = temp_imp_y(:,(m-1)*5+4); % imported gross capital formation of province m in year t
        sf_a = fd_chn(t,m,4)/(sum(sum(gcf_t))+sum(impt_gcf_t)); % balanced into year n's statistical final demand
        gcf_t_a = gcf_t*sf_a;
        impt_gcf_t_a = impt_gcf_t*sf_a;
        clear sf sf_a
        check_y(t,m,4) = sum(sum(gcf_t_a))/sum(sum(gcf_t));
        
        % Stock changes
        ys_t = reshape(temp_y(:,(m-1)*5+5),s,r); % domestic stock changes of province m in year t
        impt_ys_t = temp_imp_y(:,(m-1)*5+5); % imported stock changes of province m in year t
        sf_a = fd_chn(t,m,5)/(sum(sum(ys_t))+sum(impt_ys_t)); % balanced into year n's statistical final demand
%         if sf_a < 0
%             sf_a = 1;
%         end
        ys_t_a = ys_t*sf_a;
        impt_ys_t_a = impt_ys_t*sf_a;
        clear sf sf_a
        check_y(t,m,5) = sum(sum(ys_t_a))/sum(sum(ys_t));
        
         
        % value added
        vd_m_t = temp_vd(:,(m-1)*s+(1:s)); % value added of year t
        vd_s_m_t = squeeze(vd_s(t,m,:)); % value added of year t by sectors 
        sf = ones(s,1); % sacling factor of value added
        for i = 1:size(c_mrio2vd,2)
            index_i = find(c_mrio2vd(:,i) == 1);
            sf(index_i) = vd_s_m_t(i)/sum(sum(vd_m_t(:,index_i)));
            clear index_i
        end
        vd_m_t_a = vd_m_t.*sf';
        vd_m_t_a = vd_m_t_a.*(squeeze(vd_c(t,m,:))./sum(vd_m_t_a,2));   

        
        y_t(:,(m-1)*5+(1:5)) = horzcat(reshape(yr_t_a,s*r,1),reshape(yu_t_a,s*r,1),...
            reshape(yg_t_a,s*r,1),reshape(gcf_t_a,s*r,1),reshape(ys_t_a,s*r,1));
        impt_y_t(:,(m-1)*5+(1:5)) = horzcat(impt_yr_t_a,impt_yu_t_a,impt_yg_t_a,impt_gcf_t_a,impt_ys_t_a);
        vd_t(:,(m-1)*s+(1:s)) = vd_m_t_a;
        
        clear yr_t_a yu_t_a yg_t_a gcf_t_a ys_t_a yr_t yu_t yg_t gcf_t ys_t  ...
            impt_yr_t_a impt_yu_t_a impt_yg_t_a impt_gcf_t_a impt_ys_t_a ...
            impt_yr_t impt_yu_t impt_yg_t impt_gcf_t impt_ys_t ...
            vd_m_t_a vd_m_t
        
    end
    
    
    % export
    expt_t = reshape(temp_expt(1:r*s),s,r); % exports of year t
    sf = c_mrio2expimp*(export(t,:)'./sum(c_mrio2expimp'*expt_t,2));
    expt_t = reshape(expt_t.*sf,s*r,1);
    clear sf
    
    A = temp_z./temp_x';
    A(isnan(A)) = 0;
    A(isinf(A)) = 0;
    L = inv(eye(r*s)-A);
    z_t = A*diag(L*(sum(y_t,2)+expt_t));
    x_t = max((sum(z_t)+sum(vd_t))',L*(sum(y_t,2)+expt_t));
    
    impt_z_t = x_t'-sum(z_t)-sum(vd_t);
    impt_z_t(find(impt_z_t<0)) = 0;
%     impt_z_a = reshape(impt_z',s,r)';
    sf = (sum(sum(y_t))+sum(expt_t)-sum(sum(vd_t)))/sum(impt_z_t);
    if sf>0
        impt_z_t = impt_z_t.*sf;
    end
    sf_impt(t) = sf;
    x_t = max((sum(z_t)+sum(vd_t)+impt_z_t)',L*(sum(y_t,2)+expt_t));
    clear sf
    
    u_proxy = x_t-sum(y_t,2)-expt_t; % total outputs (i.e., row sum of Z)
    v_proxy = x_t'-sum(vd_t)-impt_z_t; % total intermediate inputs (i.e., column sum of Z) and final demand (i.e., column sum of Y)
    z_t = gras(z_t,u_proxy,v_proxy',0.1e-5); % balance the MRIO tables by GRAS method
    
    check(t,1) = max(sum(z_t)+sum(vd_t)+impt_z_t-x_t')-min(sum(z_t)+sum(vd_t)+impt_z_t-x_t');
    check(t,2) = max(sum(z_t,2)+sum(y_t,2)+expt_t-x_t)-min(sum(z_t,2)+sum(y_t,2)+expt_t-x_t);
    
    
    IOT(t).Z = z_t;
    IOT(t).Y = y_t;
    IOT(t).Exp = expt_t;
    IOT(t).Imp_Z = impt_z_t;
    IOT(t).Imp_Y = impt_y_t;
    IOT(t).vd = vd_t;
    IOT(t).x = x_t;

    clear t0 temp_z temp_y temp_expt temp_imp_z temp_imp_y temp_vd temp_x
    
    t
end

% check negative values
for n = 1:ts+2
    index(n).z = find(IOT(n).Z<0);
    index(n).y = find(IOT(n).Y(:,([1:r]-1)*5+(1:4)')<0);
    index(n).expt = find(IOT(n).Exp<0);
    index(n).x = find(IOT(n).x<0);
    index(n).impt_z = find(IOT(n).Imp_Z<0);
    index(n).impt_y = find(IOT(n).Imp_Y(:,([1:r]-1)*5+(1:4)')<0);
    index(n).vd = find(IOT(n).vd<0);
end


y_ts= zeros(ts+2,1);
x_ts = zeros(ts+2,1);
z_ts = zeros(ts+2,1);
expt_ts = zeros(ts+2,1);
for n = 1:ts+2
    y_ts(n) = sum(sum(IOT(n).Y));
    x_ts(n) = sum(IOT(n).x);
    z_ts(n) = sum(sum(IOT(n).Z));
    expt_ts(n) = sum(IOT(n).Exp);
end

figure
plot(horzcat(z_ts,y_ts,expt_ts,x_ts));
legend('z','y','expt','x')

save('MRIOT_CHN_42sectors','IOT','cpi','meta');

%% Carbon intensity
clear
load MRIOT_CHN_42sectors
load pop_chn % population
load c_c2mrio % concordance
load F_chn % direct co2 emissions

ts = 2017-1994;
r = 31;
s = 42;
s_c = 45;

f_co2 = zeros(ts,r,s);
s_tpi = find(c_c2mrio(end-2,:) == 1);
s_rra = find(c_c2mrio(end-1,:) == 1);
s_sev = find(c_c2mrio(end,:) == 1);
for n = 1:ts
    if n > 2
        x_n = reshape(IOT(n).x,s,r)'; % province by sectors
        for m = setdiff([1:r],[26])
            F_m = squeeze(F_co2(n,m,:));
            
            if sum(F_m)~=0
                temp_f = F_m(1:end-3)'*c_c2mrio(1:end-3,:)./x_n(m,:);
                temp_f(s_tpi) = F_m(end-2)/sum(x_n(m,s_tpi));
                temp_f(s_rra) = F_m(end-1)/sum(x_n(m,s_rra));
                temp_f(s_sev) = F_m(end)/sum(x_n(m,s_sev));
                temp_f(isnan(temp_f)) = 0;
                temp_f(isinf(temp_f)) = 0;
                f_co2(n,m,:) = temp_f;
                clear temp_f
            else
                f_co2(n,m,:) = squeeze(f_co2(n-1,m,:));
            end
            
            clear F_m
        end
        f_co2(n,26,:) = sum(squeeze(f_co2(n,:,:)).*x_n)./sum(x_n([1:25 27:r],:)); % national average for Tibet
        
        if n == 3
            f_co2(1,:,:) =  f_co2(n,:,:);
            f_co2(2,:,:) =  f_co2(n,:,:);
        end
    end
end

save('MRIO_co2_42sector','f_co2','Fhh_co2','F_co2')

%% calculation annual capital depreciation
clear
load MRIOT_CHN_42sectors
load ki_chn; % capital investment
load depr % depreciation rates
load c_axp
load c_ixs

ts = 2017-1994;
s = 42;
r = 31;

cd_ss = zeros(r,ts,ts,s,s);
for m  = 1:r
    for n = 1:ts % each capital consuming year
        for t = 1:n % each capital investment year
            Yk = IOT(t).Y(:,(m-1)*5+4);
%             impt_y = IOT(t).Imp_Y(:,(m-1)*5+4);
            
            Ki = squeeze(ki_chn(t,m,:,:))';
            gfcf_a = Ki;
            gfcf_a(isnan(gfcf_a))=0;
            gfcf_p = sum(reshape(Yk,[s,r]),2);
            gfcf_p(gfcf_p<0) = 0;
            cfc_v = IOT(n).vd(3,(m-1)*s+1:m*s);
            
            % converting newly capital investment data based on IOT: from assets (3) to the capital goods producing sectors (42)
            c_axp_v1 = c_axp; % products producing assets: adjusted concordance - considering the distribution of sectors producing the same assets based on the production structure described in Yk 
            for a = 1:size(Ki,2)
                temp = find(c_axp(a,:)==1);
                c_axp_v1(a,temp) = c_axp(a,temp).*gfcf_p(temp)'./sum(gfcf_p(temp));
            end
            clear temp;
            c_axp_v1(isnan(c_axp_v1)) = 0;
            gfcf_ap = gfcf_a*c_axp_v1;
            gfcf_ap = gfcf_ap.*gfcf_p'./sum(gfcf_ap,'omitnan');       % scale to Yk
            gfcf_ap(isnan(gfcf_ap))=0;
            clear c_axp_v1
            
            % adjust when Ki = 0 whereas gfcf of sector is not 0
            for j=1:s
                if gfcf_p(j)>0 && sum(gfcf_ap(:,j))==0
                    gfcf_ap(:,j) = gfcf_p(j).*(sum(gfcf_a,2)./sum(sum(gfcf_a)))';
                    j
                end
            end
            
            
            cd=zeros(size(Ki,1),s);
            for j = 1:s
                if sum(gfcf_ap(:,j))>0
                    temp = find(c_axp(:,j)>0);
                    temp1 = zeros(size(Ki,2),1); % address: one sector producing more than one assets
                    if any(temp)
                        temp1(temp)=sum(gfcf_a(:,temp))./sum(sum(gfcf_a(:,temp)));
                        if sum(temp1,'omitnan')==0
                            temp1=repmat(1/size(Ki,2),[size(Ki,2),1]);
                        end
                    else
                        temp1= sum(gfcf_a,1)./sum(sum(gfcf_a,1));
                    end
                    cdr = ((1-depr).^(n-t)).*depr;
                    cd(:,j) = sum(gfcf_ap(:,j)*c_axp(:,j)'*diag(temp1).*cdr,2);
                    
                    clear cdr temp temp1
                end
            end
            
            % distributing across capital-consumning sectors
            c_ixs_v1=c_ixs;
            for j = 1:size(Ki,1) % adjusting on the consumption side - going through each capital consuming sector
                temp = find(c_ixs(j,:)==1); % one KLEMS sector is often matched to multiple IOT sectors
                c_ixs_v1(j,temp) = c_ixs(j,temp).*(cfc_v(temp)./sum(cfc_v(temp))); %assuming sectors with higher capital compensations also consume more capital goods
                clear temp
            end
            
            cd_n_t = cd'*c_ixs_v1; % rows: capital producting sectors; columns: capital consuming (i.e. investing) sectors
            cd_n_t(isnan(cd_n_t))=0;
            cd_n_t = cd_n_t.*(sum(cd)'./sum(cd_n_t,2)); % adjusts for the case where one IOT sector is matched to multiple KLEMS sectors, by constraining on the amount of capital goods produced
                            
            cd_ss(m,n,t,:,:) = cd_n_t; 
            clear Yk impt_y Ki gfcf_a gfcf_p cfc_v c_ixs_v1 gfcf_ap cd cd_n_t
        end
    end
end

figure
temp = squeeze(sum(sum(sum(cd_ss,1,'omitnan'),4),5)).*cpi(end,end)./cpi(end,:);
bar(temp,'stacked');

save('MRIO_capital_depreciation','cd_ss')

%% supply chain-wide CO2 emissions embodied in capital depreciation
clear
load MRIOT_CHN_42sectors
load MRIO_capital_depreciation
load MRIO_co2_42sector

r = 31;
s = 42;
ts = 23;

for n = 1:ts
    temp_z = IOT(n).Z;
    temp_y = IOT(n).Y;
    temp_expt = IOT(n).Exp;
    temp_imp_y = IOT(n).Imp_Y;
    temp_x = IOT(n).x;
    vd(n).vd = IOT(n).vd;
    Imp_Z(n).Imp_Z = IOT(n).Imp_Z;
    
    Z(n).Z = temp_z;
    for m = 1:r
        Yc(n).Yc(:,m) = sum(temp_y(:,(m-1)*5+(1:3)),2);
        Yk(n).Yk(:,m) = temp_y(:,(m-1)*5+4);
        
        Imp_Yc(n).Imp_Yc(:,m) = sum(temp_imp_y(:,(m-1)*5+(1:3)),2);
        Imp_Yk(n).Imp_Yk(:,m) = sum(temp_imp_y(:,(m-1)*5+4),2);
    end
    Y(n).Y = temp_y;
    Imp_Y(n).Imp_Y = temp_imp_y;
    Expt(n).Expt = temp_expt;
    x(n).x = temp_x;
    
    clear temp_z temp_y temp_expt temp_x temp_imp_y
end
clear IOT
    
Fk = zeros(ts,ts,r,s*r); % supply chain-wide carbon emissions of capital consumption, 1st for capital depreciation year; 2nd for capital investment year
Fka = zeros(ts,ts,r,s*r); % supply chain-wide carbon emissions of capital consumption, estimated by current technology assumption 
PBEk_p_s_cd  = zeros(ts,ts,r*s,r*s); % production based carbon emissions of capital depreciation
%3rd dimension: basic production sectors
%4rd dimension:regarding the capital depreciation by capital using sectors
for t = 1:ts % production year of the depreciated capital goods
    A = Z(t).Z./x(t).x';
    A(isnan(A)) = 0;
    A(isinf(A)) = 0;
    L = inv(eye(r*s)-A);
    
    f0 = reshape(squeeze(f_co2(t,:,:))',r*s,1)';
    for n = t:ts % year of final consumption
        f0n = reshape(squeeze(f_co2(n,:,:))',r*s,1)';
        for i = 1:r
            gfcf_p = Yk(t).Yk(:,i); %i: investing province
%             gfcf_impt = Imp_Yk(t).Imp_Yk(:,i); % imported capital 
            if sum(sum(cd_ss(i,n,t,:,:),'omitnan'))>0 % capital goods produced in year t and depreciated in year n, invested by province i
                temp = squeeze(cd_ss(i,n,t,:,:)); % rows (dimension 4): capital goods producing sectors; columns: capital goods consuming sectors
                temp(isnan(temp))=0;
                temp0 = reshape(gfcf_p,[s,r]); % capital goods purchased in year t by province i; columns: province where the capital goods were purchased from
                temp0 = temp0./sum(temp0,2);
                temp0(isnan(temp0)) = 0;
                temp0(isinf(temp0)) = 0;
                temp1 = zeros(r*s,s);
                for m = 1:r % 'finished' capital goods producer
                    for j = 1:s % sectors that produced the 'finished' capital goods
                        temp1((m-1)*s+j,:) = temp0(j,m).*temp(j,:); % for the 'same' capital goods produced in the same year, assuming the cross-province distributions of depreciation and investment are idential
                    end
                end
                clear temp temp0
                %                     test_k1 =  test_k1+sum(sum(temp1,'omitnan'));
                
                temp2 = f0'.*L*temp1;
                for m = 1:r % m: provinces producing for the capital goods; i: provinces using the capital goods; columns: capital using sectors in province i
                    Fk(n,t,m,(i-1)*s+1:i*s) = sum(temp2((m-1)*s+1:m*s,:));
                end
                PBEk_p_s_cd(n,t,:,(i-1)*s+(1:s)) = temp2;
                
                temp2a = f0n'.*L*temp1;
                for m = 1:r % m: countries producing for the capital goods; i: countries using the capital goods; columns: capital using sectors in country i
                    Fka(n,t,m,(i-1)*s+1:i*s) = sum(temp2a((m-1)*s+1:m*s,:));
                end
                clear temp1 temp2 temp2a
            end
            clear gfcf_p
        end
        clear f0n
    end
    
    clear L A f0
    t
end

figure
temp = squeeze(sum(sum(Fk,3),4));
bar(temp,'stacked');

squeeze(sum(sum(Fk,2),3))-squeeze(sum(sum(PBEk_p_s_cd,2),3));

save('MRIO_CF_depr','Fk','Fka','PBEk_p_s_cd','Z','x','Y','Expt','Imp_Y','vd','-v7.3')

%% carbon footprint reassessment
clear
load MRIO_co2_42sector
load MRIO_CF_depr
load MRIOT_CHN_42sectors

r=31;
s=42;
ts=23;

for n = 1:ts
    for m = 1:r
        Yc(n).Yc(:,m) = sum(Y(n).Y(:,(m-1)*5+(1:3)),2);
        Yk(n).Yk(:,m) = Y(n).Y(:,(m-1)*5+4);
    end
end


EF_c = zeros(ts,r,s,r,s); % conventional consumption-based carbon emissions of final consumption; 2nd-3rd: province of final consumption and sectors; 4th-5th: where impacts occurred
EF_gfcf = zeros(ts,r,s,r,s); % conventional consumption-based carbon emissions of gross fixed capital formation; 2nd-3rd: province of GFCF and sectors; 4th-5th: where impacts occurred
EF_expt = zeros(ts,r,s,s);
EFk_c_s_p = zeros(ts,ts,r,s,r,s,2); % assuming capitals are used for producing noncapital goods only
    % 3rd: province of final consumption; 
    % 5th: where impacts occurred
    % 7th: for local consumption or for other places' consumption
EFk_c_s_p_a = zeros(ts,ts,r,s,r,s,2); % using current technology assumption
fk = zeros(ts,ts,r,r*s);
fk_a = zeros(ts,ts,r,r*s);
for n = 1:ts % year of capital investment
    tic
    A = Z(n).Z./x(n).x';
    A(isnan(A)) = 0;
    A(isinf(A)) = 0;
    L = inv(eye(r*s)-A);
    
    yc = zeros(r*s,r*s);
    yk = zeros(r*s,r*s);
    c_m = repmat(eye(s),r,1);
    for m = 1:r
        yc(:,(m-1)*s+(1:s)) = Yc(n).Yc(:,m).*c_m;
        yk(:,(m-1)*s+(1:s)) = Yk(n).Yk(:,m).*c_m;
    end
    expt = Expt(n).Expt.*c_m;
    
    f0 = reshape(squeeze(f_co2(n,:,:))',s*r,1)';
    temp_ef_c = f0'.*L*yc;
    temp_ef_gfcf = f0'.*L*yk;
    temp_ef_expt = f0'.*L*expt;
    for m = 1:r % region where EP occur
        for i = 1:r % region where final goods consumed
            EF_c(n,i,:,m,:) = temp_ef_c((m-1)*s+(1:s),(i-1)*s+(1:s))';
            EF_gfcf(n,i,:,m,:) = temp_ef_gfcf((m-1)*s+(1:s),(i-1)*s+(1:s))';
        end
        EF_expt(n,m,:,:) = temp_ef_expt((m-1)*s+(1:s),:);
    end
    clear temp_ef_c temp_ef_gfcf f0 temp_ef_expt

    x_fc = L*(sum(Yc(n).Yc,2)+Expt(n).Expt); % total outputs of final consumption
    for t = 1:n        
        fk_t = squeeze(Fk(n,t,:,:))./sum(x_fc,2)';
%         fk_t_a = squeeze(Fka(n,t,:,:))./sum(x_fc,2)';
        % carbon intensities for capital goods depreciated in year n (rows: province where ...
        % environemntal impacts occurred directly; columns: province-sector -where capital goods are consumed)
        % Note: depreciated capitals are used to produce both noncapital goods consumed in year n and capital goods to be consumed...
        % in year n, n+1 ... 
        fk_t(isnan(fk_t))=0;
        fk_t(isinf(fk_t))=0;
%         fk_t_a(isnan(fk_t_a))=0;
%         fk_t_a(isinf(fk_t_a))=0;
        fk(n,t,:,:) = fk_t;
%         fk_a(n,t,:,:)=fk_t_a;
        
        
        % by actual technology
        for m = 1:r
            fk1 = zeros(r,r*s);
            fk1(:,(m-1)*s+1:m*s) = fk_t(:,(m-1)*s+1:m*s);
            temp_do = fk1(m,:).*c_m'*L*yc;
            temp_fgn = sum(fk1(setdiff([1:r],m),:)).*c_m'*L*yc;
            for i = 1:r 
                EFk_c_s_p(n,t,i,:,m,:,1) = temp_do(:,(i-1)*s+1:i*s)'; %dimension3: province of final consumption; %dimension 5: where capital is consumed; 6: impacts at capital consuming province
                EFk_c_s_p(n,t,i,:,m,:,2) = temp_fgn(:,(i-1)*s+1:i*s)'; %dimension3: province of final consumption; %dimension 5: where capital is consumed; 6: impacts outside of capital consuming province
            end
            clear temp_do fk1 temp_fgn
        end
        
        % by current technology assumption
%         for m = 1:r
%             fk1 = zeros(r,r*s);
%             fk1(:,(m-1)*s+1:m*s) = fk_t_a(:,(m-1)*s+1:m*s);
%             temp_do = fk1(m,:).*c_m'*L*yc;
%             temp_fgn = sum(fk1(setdiff([1:r],m),:)).*c_m'*L*yc;
%             for i = 1:r
%                 EFk_c_s_p_a(n,t,i,:,m,:,1) = temp_do(:,(i-1)*s+1:i*s)'; %dimension3: province of final consumption; %dimension 5: where capital is consumed; 6: impacts at capital consuming province
%                 EFk_c_s_p_a(n,t,i,:,m,:,2) = temp_fgn(:,(i-1)*s+1:i*s)';
%             end
%             clear temp_do fk1 temp_fgn
%         end
%         
        clear fk_t fk_t_a
    end
    n
    toc
    clear A L y x_fc
end


temp = squeeze(sum(sum(sum(sum(sum(EFk_c_s_p,3),4),5),6),7));
figure
bar(temp,'stacked');

for n = 1:ts
    yr{n,1} = num2str(n+1994);
end

ts_a = ts;
m = 1:r;
ef_c = squeeze(sum(sum(sum(sum(EF_c(1:ts_a,m,:,:,:),2),3),4),5))+...
    squeeze(sum(sum(Fhh_co2(1:ts_a,m,:),2),3));  % unit in 10000 t
ef_p = squeeze(sum(sum(F_co2(1:ts_a,m,:),2),3))+squeeze(sum(sum(Fhh_co2(1:ts_a,m,:),2),3));
ef_gfcf = squeeze(sum(sum(sum(sum(EF_gfcf(1:ts_a,m,:,:,:),2),3),4),5));
efk_c = squeeze(sum(sum(sum(sum(sum(sum(EFk_c_s_p(1:ts_a,:,m,:,:,:,:),2),3),4),5),6),7));
figure
plot([ef_p ef_c ef_c+ef_gfcf])
hold on
plot(ef_c+efk_c,'--k')
ax=gca;
ax.XLim = [0.25 ts_a+0.75];
ax.XTick = [1:ts_a];
ax.XTickLabel = yr(1:ts_a);
ax.XTickLabelRotation = 90;
ax.XColor = 'black';
ax.YColor = 'black';
legend('PBA','CBA (fc)','CBA (all)','CBA+KCBA','Location','northwest');

% save('MRIO_CF_reassess','EF_c','EF_gfcf','EF_expt','EFk_c_s_p','EFk_c_s_p_a','fk','fk_a','-v7.3')
save('MRIO_CF_reassess','EF_c','EF_gfcf','EF_expt','EFk_c_s_p','fk','-v7.3')

%% carbon footprint reassessment: capital used for producing all the outputs
clear
load MRIO_co2_42sector
load MRIO_CF_depr
load MRIOT_CHN_42sectors

r=31;
s=42;
ts=23;

for n = 1:ts
    for m = 1:r
        Yc(n).Yc(:,m) = sum(Y(n).Y(:,(m-1)*5+(1:3)),2);
        Yk(n).Yk(:,m) = Y(n).Y(:,(m-1)*5+4);
    end
end


EF_c = zeros(ts,r,s,r,s); % conventional consumption-based carbon emissions of final consumption; 2nd-3rd: province of final consumption and sectors; 4th-5th: where impacts occurred
EF_gfcf = zeros(ts,r,s,r,s); % conventional consumption-based carbon emissions of gross fixed capital formation; 2nd-3rd: province of GFCF and sectors; 4th-5th: where impacts occurred
EF_expt = zeros(ts,r,s,s);

% capitals are used for producing total outputs
EFk_c_s_fc = zeros(ts,ts,r,s,r,s); 
EFk_c_s_gfcf = zeros(ts,ts,r,s,r,s);
    % 3rd: province of final consumption / GFCF; 
    % 5th: where impacts occurred
EFk_c_s_expt = zeros(ts,ts,s,r,s);
    % 3rd: sectors that are exported; 
    % 4th: where impacts occurred

fk_x = zeros(ts,ts,r,r*s); % co2 intensity of total outputs
for n = 1:ts % year of capital investment
    tic
    A = Z(n).Z./x(n).x';
    A(isnan(A)) = 0;
    A(isinf(A)) = 0;
    L = inv(eye(r*s)-A);
    
    yc = zeros(r*s,r*s);
    yk = zeros(r*s,r*s);
    c_m = repmat(eye(s),r,1);
    for m = 1:r
        yc(:,(m-1)*s+(1:s)) = Yc(n).Yc(:,m).*c_m;
        yk(:,(m-1)*s+(1:s)) = Yk(n).Yk(:,m).*c_m;
    end
    expt = Expt(n).Expt.*c_m;
    
    f0 = reshape(squeeze(f_co2(n,:,:))',s*r,1)';
    temp_ef_c = f0'.*L*yc;
    temp_ef_gfcf = f0'.*L*yk;
    temp_ef_expt = f0'.*L*expt;
    for m = 1:r % region where EP occur
        for i = 1:r % region where final goods consumed
            EF_c(n,i,:,m,:) = temp_ef_c((m-1)*s+(1:s),(i-1)*s+(1:s))';
            EF_gfcf(n,i,:,m,:) = temp_ef_gfcf((m-1)*s+(1:s),(i-1)*s+(1:s))';
        end
        EF_expt(n,m,:,:) = temp_ef_expt((m-1)*s+(1:s),:);
    end
    clear temp_ef_c temp_ef_gfcf f0 temp_ef_expt

    x_n = x(n).x; % total outputs of final consumption
    for t = 1:n        
        fk_t = squeeze(Fk(n,t,:,:))./sum(x_n,2)';
        % carbon intensities for capital goods depreciated in year n (rows: province where ...
        % environemntal impacts occurred directly; columns: province-sector -where capital goods are consumed)
        % Note: depreciated capitals are used to produce both noncapital goods consumed in year n and capital goods to be consumed...
        % in year n, n+1 ... 
        fk_t(isnan(fk_t))=0;
        fk_t(isinf(fk_t))=0;
        fk_x(n,t,:,:) = fk_t;        
        
        % by actual technology
        for m = 1:r
            fk1 = zeros(1,r*s);
            fk1(:,(m-1)*s+1:m*s) = sum(fk_t(:,(m-1)*s+1:m*s));
            temp_fc = fk1.*c_m'*L*yc;
            temp_gfcf = fk1.*c_m'*L*yk;
            temp_expt = fk1.*c_m'*L*expt;
            for i = 1:r 
                EFk_c_s_fc(n,t,i,:,m,:) = temp_fc(:,(i-1)*s+1:i*s)'; %dimension3: province of final consumption; %dimension 5: where capital is consumed; 6: impacts at capital consuming province
                EFk_c_s_gfcf(n,t,i,:,m,:) = temp_gfcf(:,(i-1)*s+1:i*s)'; %dimension3: province of final consumption; %dimension 5: where capital is consumed; 6: impacts outside of capital consuming province
            end
            EFk_c_s_expt(n,t,:,m,:) = temp_expt';
            clear fk1 temp_fc temp_gfcf temp_expt
        end

        clear fk_t fk_t_a
    end
    n
    toc
    clear A L y x_n
end


temp = squeeze(sum(sum(sum(sum(EFk_c_s_fc,3),4),5),6));
figure
bar(temp,'stacked');

for n = 1:ts
    yr{n,1} = num2str(n+1994);
end

ts_a = ts;
m = 1:r;
ef_c = squeeze(sum(sum(sum(sum(EF_c(1:ts_a,m,:,:,:),2),3),4),5))+...
    squeeze(sum(sum(Fhh_co2(1:ts_a,m,:),2),3));  % unit in 10000 t
ef_p = squeeze(sum(sum(F_co2(1:ts_a,m,:),2),3))+squeeze(sum(sum(Fhh_co2(1:ts_a,m,:),2),3));
ef_gfcf = squeeze(sum(sum(sum(sum(EF_gfcf(1:ts_a,m,:,:,:),2),3),4),5));
ef_expt = squeeze(sum(sum(sum(EF_expt,2),3),4));
efk_fc = squeeze(sum(sum(sum(sum(sum(EFk_c_s_fc(1:ts_a,:,m,:,:,:),2),3),4),5),6));
efk_gfcf = squeeze(sum(sum(sum(sum(sum(EFk_c_s_gfcf(1:ts_a,:,m,:,:,:),2),3),4),5),6));
efk_expt = squeeze(sum(sum(sum(sum(EFk_c_s_expt(1:ts_a,:,:,m,:),2),3),4),5));
figure
plot([ef_p ef_c ef_c+ef_gfcf])
hold on
plot(ef_c+efk_fc+efk_gfcf+ef_expt+efk_expt,'--k')  % PBE considering capital
plot(ef_c+efk_fc+efk_gfcf,'--r') % CBE considering capital
ax=gca;
ax.XLim = [0.25 ts_a+0.75];
ax.XTick = [1:ts_a];
ax.XTickLabel = yr(1:ts_a);
ax.XTickLabelRotation = 90;
ax.XColor = 'black';
ax.YColor = 'black';
legend('PBA','CBA (fc)','CBA (all)','CBA+KCBA','Location','northwest');

% cd('C:\Users\YeQ\Documents\MATLAB\Capital_model_IOT')
% save('MRIO_CF_reassess','EF_c','EF_gfcf','EF_expt','EFk_c_s_p','EFk_c_s_p_a','fk','fk_a','-v7.3')
save('MRIO_CF_reassess_cd_for_output','EF_c','EF_gfcf','EF_expt','EFk_c_s_fc','EFk_c_s_gfcf','EFk_c_s_expt','fk_x','-v7.3')




%%
load MRIO_capital_depreciation
load CF_2012
n = 2010-1994;
squeeze(sum(sum(cd_ss(:,n,:,:,:),1,'omitnan'),3));

n = 2007-1994;
cba = squeeze(sum(sum(EF_c(n,:,:,:,:),3),5))+squeeze(sum(sum(EF_gfcf(n,:,:,:,:),3),5));
cba_to = sum(cba,2)/100;

n = 2012-1994;
cba = squeeze(sum(EF_c(n,:,:,:),4))+squeeze(sum(EF_gfcf(n,:,:,:),4));
cba_to = sum(cba,2)/100;

cba_to_12 = squeeze(sum(sum(EF_c_2012(:,:,:),2),3))/100+squeeze(sum(sum(EF_gfcf_2012(:,:,:),2),3))/100;

%% carbon footprint in 2012
clear
cd('C:\Users\YeQ\Documents\Ph.D UT\Phd in UT\Data China\MRIOTs')
load MRIO_CHN_2012

cd('C:\Users\YeQ\Documents\MATLAB\Capital_model_IOT')
load ('MRIO_co2_42sector','F_co2','Fhh_co2')

cd('C:\Users\YeQ\Documents\MATLAB\Capital_model_IOT\Data')
c_c2mrio = xlsread('Concordances','45 carbon sectors to 42 (2012)','b3:AQ47');
c_c2mrio(isnan(c_c2mrio)) = 0;
sum(c_c2mrio,2);

r = 31;
s = 42;
n = 2012-1994;

f_2012 = zeros(r,s);
s_rra = find(c_c2mrio(end-1,:) == 1);
s_sev = find(c_c2mrio(end,:) == 1);
x_n = reshape(IO.x,s,r)'; % province by sectors
for m = setdiff([1:r],[26])
    F_m = squeeze(F_co2(n,m,:));
    
    if sum(F_m)~=0
        temp_f = F_m(1:end-2)'*c_c2mrio(1:end-2,:)./x_n(m,:);
        temp_f(s_rra) = F_m(end-1)/sum(x_n(m,s_rra));
        temp_f(s_sev) = F_m(end)/sum(x_n(m,s_sev));
        temp_f(isnan(temp_f)) = 0;
        temp_f(isinf(temp_f)) = 0;
        f_2012(m,:) = temp_f;
        clear temp_f
    end
    
    clear F_m
end
f_2012(26,:) = sum(f_2012(:,:).*x_n)./sum(x_n([1:25 27:r],:)); % national average for Tibet


IO.Y(find(IO.Y<0)) = 0;
IO.Exp(find(IO.Exp<0)) = 0;
for m = 1:r
    Yc(:,m) = sum(IO.Y(:,(m-1)*5+(1:3)),2);
    Yk(:,m) = IO.Y(:,(m-1)*5+4);
end

A = IO.Z./IO.x';
A(isnan(A)) = 0;
A(isinf(A)) = 0;
L = inv(eye(r*s)-A);

y = zeros(r*s,r*s);
c_m = repmat(eye(s),r,1);
for m = 1:r
    y(:,(m-1)*s+(1:s)) = Yc(:,m).*c_m;
end
clear c_m

f0 = reshape(f_2012(:,:)',s*r,1)';
temp_ef_c = f0'.*L*Yc;
temp_ef_gfcf = f0'.*L*Yk;
temp_ef_expt = f0'.*L*IO.Exp;
for m = 1:r
    EF_c_2012(:,m,:) = temp_ef_c((m-1)*s+(1:s),:)';
    EF_gfcf_2012(:,m,:) = temp_ef_gfcf((m-1)*s+(1:s),:)';
    EF_expt_2012(m,:) = temp_ef_expt((m-1)*s+(1:s));
end
clear temp_ef_c temp_ef_gfcf f0

cd('C:\Users\YeQ\Documents\MATLAB\Capital_model_IOT')
save('CF_2012','EF_c_2012','EF_gfcf_2012','EF_expt_2012')

%% carbon footprints of capital investment, by asset and sectors
clear
load MRIOT_CHN_42sectors
load ki_chn; % capital investment
load c_axp
load MRIO_co2_42sector

ts = 2017-1994;
s = 42;
s_k = 37;
r = 31;

EF_ki = zeros(r,ts,3,s_k); % carbon footprint of capital investment sectors
% dimension 3rd: capital asset type
% dimension 4th: capitla investing sectors
Ki_gfcf = zeros(r,ts,s,s_k);
for t = 1:ts % each capital consuming year
    A = IOT(t).Z./IOT(t).x';
    A(isnan(A)) = 0;
    A(isinf(A)) = 0;
    L = inv(eye(r*s)-A);
    
    f0 = reshape(squeeze(f_co2(t,:,:))',r*s,1)';
    for m  = 1:r
        Yk = IOT(t).Y(:,(m-1)*5+4);
        
        Ki = squeeze(ki_chn(t,m,:,:))';
        gfcf_a = Ki;
        gfcf_a(isnan(gfcf_a))=0;
        gfcf_p = sum(reshape(Yk,[s,r]),2);
        gfcf_p(gfcf_p<0) = 0;
        
        % converting newly capital investment data based on IOT: from assets (3) to the capital goods producing sectors (42)
        c_axp_v1 = c_axp; % products producing assets: adjusted concordance - considering the distribution of sectors producing the same assets based on the production structure described in Yk
        for a = 1:size(Ki,2)
            temp = find(c_axp(a,:)==1);
            c_axp_v1(a,temp) = c_axp(a,temp).*gfcf_p(temp)'./sum(gfcf_p(temp));
        end
        clear temp;
        c_axp_v1(isnan(c_axp_v1)) = 0;
        gfcf_ap = gfcf_a*c_axp_v1;
        gfcf_ap = gfcf_ap.*gfcf_p'./sum(gfcf_ap,'omitnan');       % scale to Yk
        gfcf_ap(isnan(gfcf_ap))=0;
        
        % adjust when Ki = 0 whereas gfcf of sector is not 0
        for j=1:s
            if gfcf_p(j)>0 && sum(gfcf_ap(:,j))==0
                gfcf_ap(:,j) = gfcf_p(j).*(sum(gfcf_a,2)./sum(sum(gfcf_a)))';
                j
            end
        end
        gfcf_ap = gfcf_ap';
        Ki_gfcf(m,t,:,:) = gfcf_ap;
        
        % allocate across provinces
        temp0 = reshape(Yk,[s,r]); % capital goods purchased in year t by province i; columns: province where the capital goods were purchased from
        temp0 = temp0./sum(temp0,2);
        temp0(isnan(temp0)) = 0;
        temp0(isinf(temp0)) = 0;
        
        temp1 = zeros(r*s,s_k);
        for i = 1:r % 'finished' capital goods producer
            for j = 1:s % sectors that produced the 'finished' capital goods
                temp1((i-1)*s+j,:) = temp0(j,i).*gfcf_ap(j,:); % for the 'same' capital goods produced in the same year, assuming the cross-province distributions of depreciation and investment are idential
            end
        end
        clear temp temp0 gfcf_ap Ki gfcf_av c_axp_v1
        
        % CF estimation
        temp2 = f0'.*L*temp1;
        for i = 1:s_k
            temp_ef = sum(reshape(temp2(:,i),s,r),2);
            EF_ki(m,t,:,i) = c_axp*temp_ef;
            clear temp_ef
        end
       
        clear Yk gfcf_p temp2
    end
    clear L A f0
end


cd('C:\Users\YeQ\Documents\MATLAB\Capital_model_IOT')
save('MRIO_CF_by_3assets_37sector','EF_ki','Ki_gfcf')
