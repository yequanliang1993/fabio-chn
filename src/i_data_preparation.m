%% this is the Matlab codes to prepare data input for FABIO-CHN

% Created: Quanliang Ye
% Date: 22-October-2024
% Email: yequanliang1993@gmail.com
% Version: 2.3.0

% Note:
%   You are more than welcome to comments on the model and codes.

%% Configure default information
clear

% root path
path_root = pwd+"\";

% current version
path_parts = split(path_root,'\');
version ="_"+path_parts{end-1};
clear path_parts

% number of products in FABIO, including 130 agri-food products
p0 = 130;

% number of provinces in china, 31 provinces
r = 31;

% number of years, for raw data, between 1990 and 2018
ts = 2018-1990+1;


%% match product codes between FABIO and FABIO-CHN 

% import data
path_faostat = path_root+"Data\FAOSTAT\";
filename = "Supply-use table (refered to FABIO supporting infor).xlsx";
data = importdata(path_faostat+filename);
meta.fabio_crop = data.textdata.DataAvailableInChina(2:end,[1 3:end]);

% matching FABIO product code with FABIO_CHN product code
prod_c = data.data.DataAvailableInChina(find(data.data.SupplyTableCHN(2:131,5) == 1),8:11); % product codes of Chinese FABIO

% supply table concordance, no real data
sup_tab = data.data.SupplyTableCHN(2:131,6:126);  % format of the supple table (concordance between products and processes)
sup_tab(find(data.data.SupplyTableCHN(2:131,5) == 0),:) = []; % remove the rows of products without data
sup_tab(:,find(data.data.SupplyTableCHN(1,6:126) == 0)) = []; % remove the columns of processes withourt data
sup_tab(isnan(sup_tab)) = 0;
sup_tab_chn = zeros(ts,r,size(sup_tab,1),size(sup_tab,2));
% sup_tab_chn dimension 1: year? 1990-2018
% sup_tab_chn dimension 2: region
% sup_tab_chn dimension 3: products
% sup_tab_chn dimension 4: processes
for t = 1:ts
    for m = 1:r
        sup_tab_chn(t,m,:,:) = sup_tab;
    end
end

% use table concordance, no real data
use_tab = data.data.UseTable(2:131,6:126);  % format of the supple table (concordance between products and processes)
use_tab(find(data.data.UseTable(2:131,5) == 0),:) = []; % remove the rows of products without data
use_tab(:,find(data.data.UseTable(1,6:126) == 0)) = []; % remove the columns of processes withourt data
use_tab(isnan(use_tab)) = 0;

% technical conversion factors of products
tcf_fao = data.data.TechConvFactor(:,1:3);

% type of processing (single, multiple or others)
proc_type = data.data.TechConvFactor(:,4:5); 

%% Commodity balance sheets (CBS) from FAOSTAT --- Crops
filename = "Commodity_balance_sheet_crop.xlsx";
data = importdata(path_faostat+filename);
element = {'Production';'Import Quantity';'Stock Variation';'Export Quantity';...
    'Domestic supply quantity';'Feed';'Seed';'Losses';'Processing';...
    'Food supply quantity (tonnes)';'Other uses'}; % elements in FAOSTAT
item_crop = unique(data.textdata(2:end,8));
ts_a = max(data.data(:,8))-min(data.data(:,8))+1; % only avaiable between 1990 and 2013
commod_bal = zeros(ts,p0,length(element));   % commodity balance sheets, 1st dimension: time series; 2nd dimension: number of crops; 3rd dimension: elements

data_avi_fao = zeros(p0,length(element)); % check the data available for each product
crop_a = [];       % for crop in a short name included in crops in a longer name, e.g., 'Cottonseed' also in 'Cottonseed oil'
for i = 1:length(element)
    index1_c = find(contains(data.textdata(2:end,6),element(i)));
    for j = 1:p0
        index2_c = find(contains(data.textdata(2:end,8),meta.fabio_crop(j,2)));
        if isempty(index2_c) == 0
            index_c = intersect(index1_c,index2_c);
            if isempty(index_c) == 0
                if size(index_c,1) == ts_a
                    commod_bal(1:ts_a,j,i) = data.data(index_c,end)/10000;   % unit: 10000 tons
                    data_avi_fao(j,i) = 1;
                else
                    k = length(crop_a);
                    if k == 0
                        crop_a(k+1) = j;
                    else
                        if any(crop_a == j)
                            crop_a;
                        else
                            crop_a(k+1) = j;
                        end
                    end
                end
            end
        end
    end
end

data_avi_fao(crop_a,:) = zeros(length(crop_a),length(element)); % omit data for further adjustment

% futher adjust crop_a
for i = 1:length(element)
    index1_c = find(contains(data.textdata(2:end,6),element(i)));
        
    temp = find(contains(meta.fabio_crop(:,2),meta.fabio_crop(crop_a(1),2)));
    index2_c1 = find(contains(data.textdata(2:end,8),meta.fabio_crop(temp(1),2)));
    index2_c2 = find(contains(data.textdata(2:end,8),meta.fabio_crop(temp(2),2)));
    index2_c = setdiff(index2_c1,index2_c2);
    index_c = intersect(index1_c,index2_c);
    if isempty(index_c) == 0
        commod_bal(1:ts_a,temp(1),i) = data.data(index_c,end)/10000;
        data_avi_fao(temp(1),i) = 1;
    end
    clear temp index2_c1 index2_c2 index2_c index_c
    
    temp = find(contains(meta.fabio_crop(:,2),meta.fabio_crop(crop_a(2),2)));
    index2_c1 = find(contains(data.textdata(2:end,8),meta.fabio_crop(temp(1),2)));
    index2_c2 = find(contains(data.textdata(2:end,8),meta.fabio_crop(temp(2),2)));
    index2_c3 = find(contains(data.textdata(2:end,8),meta.fabio_crop(temp(3),2)));
    index2_c = setdiff(index2_c1,vertcat(index2_c2,index2_c3));
    index_c = intersect(index1_c,index2_c);
    if isempty(index_c) == 0
        commod_bal(1:ts_a,temp(1),i) = data.data(index_c,end)/10000;
        data_avi_fao(temp(1),i) = 1;
    end
    clear temp index2_c1 index2_c2 index2_c3 index2_c index_c
    
    temp = find(contains(meta.fabio_crop(:,2),meta.fabio_crop(crop_a(3),2)));
    index2_c1 = find(contains(data.textdata(2:end,8),meta.fabio_crop(temp(1),2)));
    index2_c2 = find(contains(data.textdata(2:end,8),meta.fabio_crop(temp(2),2)));
    index2_c3 = find(contains(data.textdata(2:end,8),meta.fabio_crop(temp(3),2)));
    index2_c = setdiff(index2_c3,vertcat(index2_c1,index2_c2));
    index_c = intersect(index1_c,index2_c);
    if isempty(index_c) == 0
        commod_bal(1:ts_a,temp(3),i) = data.data(index_c,end)/10000;
        data_avi_fao(temp(3),i) = 1;
    end
    clear temp index2_c1 index2_c2 index2_c3 index2_c index_c
        
    temp = find(contains(meta.fabio_crop(:,2),meta.fabio_crop(crop_a(4),2)));
    index2_c1 = find(contains(data.textdata(2:end,8),meta.fabio_crop(temp(1),2)));
    index2_c2 = find(contains(data.textdata(2:end,8),meta.fabio_crop(temp(2),2)));
    index2_c = setdiff(index2_c1,index2_c2);
    index_c = intersect(index1_c,index2_c);
    if isempty(index_c) == 0
        commod_bal(1:ts_a,temp(1),i) = data.data(index_c,end)/10000;
        data_avi_fao(temp(1),i) = 1;
    end
    clear temp index2_c1 index2_c2 index2_c index_c index1_c
end

clear crop_a

% balancing data of coffee
j = 45; 
t = [7:19];
commod_bal(t,j,5) = commod_bal(t,j,10);
commod_bal(t,j,1) = commod_bal(t,j,10)+commod_bal(t,j,4)-commod_bal(t,j,2)-commod_bal(t,j,3);

%% Commodity balance sheets (CBS) from FAOSTAT --- Livestock
filename = "Commodity_balance_sheet_livestock.xlsx";
data = importdata(path_faostat+filename);
item_livest = unique(data.textdata(2:end,8));
ts_a = max(data.data(:,8))-min(data.data(:,8))+1; % years covered in FAOSTAT

crop_a = [];
for i = 1:length(element)
    index1_c = find(contains(data.textdata(2:end,6),element(i)));
    for j = 1:p0
        index2_c = find(contains(data.textdata(2:end,8),meta.fabio_crop(j,2)));
        if isempty(index2_c) == 0
            index_c = intersect(index1_c,index2_c);
            if isempty(index_c) == 0
                if size(index_c,1) == ts_a
                    commod_bal(1:ts_a,j,i) = data.data(index_c,end)/10000;   % unit: 10000 tons
                    data_avi_fao(j,i) = 1;
                else
                    k = length(crop_a);
                    if k == 0
                        crop_a(k+1) = j;
                    else
                        if any(crop_a == j)
                            crop_a;
                        else
                            crop_a(k+1) = j;
                        end
                    end
                end
            end
        end
    end
end

clear crop_a index1_c index2_c index_c

%% Commodity balance sheets (CBS) from FAOSTAT --- Live animals
filename = "Live_animal_slaughtered.xlsx";    % regarded as element "Production"
data = importdata(path_faostat+filename);            % only contains the production of live animals
item_animal = unique(data.textdata(2:end,8));
ts_a = max(data.data(:,8))-min(data.data(:,8))+1;  % time series of data from FAO (1991-2018) match with statistics data during (1990-2017) 

crop_a = [];
for j = 1:p0 
    if j == 102  % 102 is 'Poultry Birds' as a total of chicken, duck, geese and guinea fowls
        index_c = vertcat(find(contains(data.textdata(2:end,8),'chicken')),find(contains(data.textdata(2:end,8),'duck')),...
            find(contains(data.textdata(2:end,8),'goose and guinea fowl')));
        commod_bal(1:ts_a,j,1) = sum(reshape(data.data(index_c,end),ts_a,3),2,'omitnan')/10;   % unit: 10000 stocks
        
        data_avi_fao(j,1) = 1;
    else
        index_c = find(contains(data.textdata(2:end,8),meta.fabio_crop(j,2)));
        if isempty(index_c) == 0
            if size(index_c,1) == ts_a
                if j == 108   % 108 'Rabbits and hares' in unit of 1000 heads
                    commod_bal(1:ts_a,j,1) = data.data(index_c,end)/10;   % unit: 10000 heads
                    data_avi_fao(j,1) = 1;
                else
                    commod_bal(1:ts_a,j,1) = data.data(index_c,end)/10000; % unit: 10000 heads
                    data_avi_fao(j,1) = 1;
                end
            else
                k = length(crop_a);
                if k == 0
                    crop_a(k+1) = j;
                else
                    if any(crop_a == j)
                        crop_a;
                    else
                        crop_a(k+1) = j;
                    end
                end
            end
        end
    end
end

clear crop_a index_c 

%% Commodity balance sheets (CBS) from FAOSTAT --- Live animals trade
filename = "Live_animal_trade from FAO.xlsx";
data = importdata(path_faostat+filename);            % only contains the production of live animals
item_animal = unique(data.textdata(2:end,8));
ts_a = max(data.data(:,8))-min(data.data(:,8))+1;  % time series of data from FAO (1991-2018) match with statistics data during (1990-2017) 
meta.fabio_crop(97:110,2) = {'Cattle','Buffaloes','Sheep','Goats','Pigs','Poultry Birds',...
    'Horses','Asses','Mules','Camels','Camelids, other','Rabbits and hares','Rodents, other','Live animals, other'};

for i = [2 4]
    index1_c = find(contains(data.textdata(2:end,6),element(i)));
    for j = 1:p0
        if j == 102  % 102 is 'Poultry Birds' as a total of chicken, duck, geese and guinea fowls
            index2_c = vertcat(find(contains(data.textdata(2:end,8),'Chickens')),find(contains(data.textdata(2:end,8),'Ducks')),...
                find(contains(data.textdata(2:end,8),'Geese and guinea fowls')),find(contains(data.textdata(2:end,8),'Turkeys')));
            if isempty(index2_c) == 0
                index_c = intersect(index1_c,index2_c);
                if isempty(index_c) == 0
                    for t = 1:max(data.data(index_c,7))-min(data.data(index_c,7))+1  % the time-series of data are not identical 
                        temp = index_c(find(data.data(index_c,7) == t+1989));
                        commod_bal(t,j,i) = sum(data.data(temp,end),'omitnan')/10;   % unit: 10000 tons
                        data_avi_fao(j,i) = 1;
                    end
                end
            end
        else
            index2_c = find(contains(data.textdata(2:end,8),meta.fabio_crop(j,2)));
            if isempty(index2_c) == 0
                index_c = intersect(index1_c,index2_c);
                if isempty(index_c) == 0
                    if j == 108   % 108 'Rabbits and hares' in unit of 1000 heads
                        commod_bal(1:length(index_c),j,i) = data.data(index_c,end)/10;   % unit: 10000 heads
                        data_avi_fao(j,i) = 1;
                    else
                        commod_bal(1:length(index_c),j,i) = data.data(index_c,end)/10000; % unit: 10000 heads
                        data_avi_fao(j,i) = 1;
                    end
                end
            end
        end
    end
end

clear index_c index1_c index2_c

%% Commodity balance sheets (CBS) from FAOSTAT --- Fishery
filename = "Fishery products.xlsx";
data = importdata(path_faostat+filename);
item_fish = unique(data.textdata(2:end,8));
ts_a = max(data.data(:,8))-min(data.data(:,8))+1;

j = 127;
for i = 1:length(element)
    index_c = find(contains(data.textdata(2:end,6),element(i)));
    if isempty(index_c) == 0
        temp_item = unique(data.textdata(index_c+1,8));
        commod_bal(1:ts_a,j,i) = sum(reshape(data.data(index_c,end),ts_a,length(temp_item)),2)/10000;
        data_avi_fao(j,i) = 1;
    end
    
    clear index_c temp_item
end

%% Commodity balance sheets (CBS) from FAOSTAT --- Forestery
filename = "Forestery products.xlsx";
data = importdata(path_faostat+filename);
item_forest = unique(data.textdata(2:end,8));
ts_a = max(data.data(:,8))-min(data.data(:,8))+1;
j = 128;
for i = 1:length(element)
    index_c = find(contains(data.textdata(2:end,6),element(i)));
    if isempty(index_c) == 0      
        for t = 1:max(data.data(index_c,7))-min(data.data(index_c,7))+1  % the time-series of data are not identical
            temp = index_c(find(data.data(index_c,7) == t+1989));
            commod_bal(t,j,i) = sum(data.data(temp,end),'omitnan')/10000;   % unit: 10000 m^3
            data_avi_fao(j,i) = 1;
            
            clear temp
        end
    end
    
    clear index_c
end
commod_bal(1:ts_a,j,5) = squeeze(sum(commod_bal(1:ts_a,j,1:2),3))-squeeze(commod_bal(1:ts_a,j,4));
commod_bal(1:ts_a,j,11) = commod_bal(1:ts_a,j,5);

save(path_root+"data_avi_fao"+version+".mat",'data_avi_fao')

%% China statistics data for production --- Crops
path_chn_crop = path_root+"Data\China_Data\Crop production in provinces\";

prod_stat = zeros(ts,p0,r); % dimension 1: time period 1990-2018; dimension 2: 130 agricultural products, same with FABIO; dimension 3: 31 procinves
for i = 1:r
    fname = sprintf("%d.xls",i);
    temp = importdata(path_chn_crop+fname);
    meta.regions_chn{i,1} = temp.textdata{2,1}(8:end);

    for j = 1:p0
        if isempty(meta.fabio_crop{j,4}(1:end)) == 0 
            index_c = find(contains(temp.textdata(5:44,1),meta.fabio_crop(j,4)));
            if isempty(index_c) == 0
                for t = 1:ts
                    if j == 13
                        temp1 = temp.data(index_c(1),ts-t+1); % production of all tubers
                        temp2 = temp.data(index_c(1)+1,ts-t+1); % production of potatoes
                        temp1(isnan(temp1)) = 0;
                        temp2(isnan(temp2)) = 0;
                        temp3 = temp1-temp2; % production of tubers excluding potatoes as 'Roots, other'
                        temp3(find(temp3<0)) = 0;
                        
                        prod_stat(t,j,i) = temp3;   % unit: 10000 tons
                        clear temp1 temp2 temp3
                    elseif j == 17
                        temp1 = temp.data(index_c(1),ts-t+1); % production of all beans
                        temp2 = temp.data(index_c(1)+3,ts-t+1); % production of all soybeans
                        temp1(isnan(temp1)) = 0;
                        temp2(isnan(temp2)) = 0;
                        temp3 = temp1-temp2; % production of beans excluding soybeans as 'beans other'
                        temp3(find(temp3<0)) = 0;
                        
                        prod_stat(t,j,i) = temp3;   % unit: 10000 tons
                        clear temp1 temp2 temp3
                    elseif j == 30
                        temp1 = temp.data(index_c(1),ts-t+1); % production of all oil crops
                        temp2 = sum(temp.data(index_c(1)+[1:4],ts-t+1),'omitnan'); % production of peanuts, sunflowers, rapeseeds, sesame
                        temp1(isnan(temp1)) = 0;
                        temp2(isnan(temp2)) = 0;
                        temp3 = temp1-temp2; % production of oil crops excluding peanuts, sunflowers, rapeseeds, sesame as 'oil crops other'
                        temp3(find(temp3<0)) = 0;
                        
                        prod_stat(t,j,i) = temp3;   % unit: 10000 tons
                        clear temp1 temp2 temp3
                    else
                        prod_stat(t,j,i) = temp.data(index_c(1),ts-t+1);   % unit: 10000 tons
                    end
                end
            end
        end
    end
end

% missing data of barley
j = 3; 
temp_r = squeeze(prod_stat(13,j,:))./sum(prod_stat(13,j,:),'omitnan');
prod_stat(1:12,j,:) = squeeze(commod_bal(1:12,j,1)).*temp_r';
clear temp_r j 

% missing data of beans
j = 17; 
temp_r = squeeze(prod_stat(2,j,:))./sum(prod_stat(2,j,:),'omitnan');
prod_stat(1,j,:) = squeeze(commod_bal(1,j,1)).*temp_r';
clear temp_r j

% cotton seed;
j_1 = 63; 
% cotton lint
j_2 = 96; 
prod_stat(:,j_1,:) = squeeze(prod_stat(:,j_2,:))/tcf_fao(j_2,2)*tcf_fao(j_1,2);
clear j_1 j_2

%% China statistics data for production --- Animal slaughtered and products
path_chn_ani_prod = path_root+"Data\China_Data\Animal Slaughtered and products\";

filename = "Cattle and Buffaloes.xls";
data = importdata(path_chn_ani_prod+filename);
j = 97; % product code 97 is Cattle.
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
end
clear j data 


filename = "Sheeps and goats.xls";
data = importdata(path_chn_ani_prod+filename);
j = 99; % product code 99 is Sheep
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
end
clear data j


filename = "pigs for pork.xls";
data = importdata(path_chn_ani_prod+filename);
j = 101; % product code 101 is Pigs
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
end
clear data j


filename = "Poultry.xls";
data = importdata(path_chn_ani_prod+filename);
j = 102; % product code 102 is Poultary
for t = 1:ts
    if t == 1  % for year without producation data, we distribute the FAO data by the provincial sharing of the nesrest year
        prod_stat(t,j,:) = commod_bal(t,j,1)*(data.data(1:r,28)/sum(data.data(1:r,28),'omitnan'));
    else
        prod_stat(t,j,:) = data.data(1:r,ts-t+1);
    end
end
clear data j


filename = "Horses.xls";
data = importdata(path_chn_ani_prod+filename);
j = 103; % product code 103 is Horses
for t = 1:ts
    if t == 8 || t == 9
        temp = sum(data.data(1:r,[7 10]),2,'omitnan')./sum(sum(data.data(1:r,[7 10]),'omitnan'));
        prod_stat(t,j,:) = commod_bal(t,j,1)*temp;
        
        clear temp
    else
        prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
    end
end
clear data j


filename = "Donkey.xls";
data = importdata(path_chn_ani_prod+filename);
j = 104; % product code 104 is Horses
for t = 1:ts
    if t == 8 || t == 9
        temp = sum(data.data(1:r,[7 10]),2,'omitnan')./sum(sum(data.data(1:r,[7 10]),'omitnan'));
        prod_stat(t,j,:) = commod_bal(t,j,1)*temp;
        
        clear temp
    else
        prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
    end
end
clear data j


filename = "Mules.xls";
data = importdata(path_chn_ani_prod+filename);
j = 105; % product code 105 is Mules
for t = 1:ts
    if t == 8 || t == 9
        temp = sum(data.data(1:r,[7 10]),2,'omitnan')./sum(sum(data.data(1:r,[7 10]),'omitnan'));
        prod_stat(t,j,:) = commod_bal(t,j,1)*temp;
        
        clear temp
    else
        prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
    end
end
clear data j


filename = "All milk.xls"; % milk including cow milk
data1 = importdata(path_chn_ani_prod+filename);
filename = "Milk.xls"; % cow milk only
data2 = importdata(path_chn_ani_prod+filename);
j = 111; % product code 111 is Milk
j_c = 70; % product code of Chinese FABIO
for t = 1:ts
    if t == 10 % missing year: distributing by the provincial share in the nearest year
        temp = sum(data1.data(1:r,[9 11]),2,'omitnan')./sum(sum(data1.data(1:r,[9 11]),'omitnan'));
        temp_a = squeeze(commod_bal(t,j,1))*temp;
        temp_milk = max(temp_a,data2.data(1:r,ts-t+1));
        prod_stat(t,j,:) = temp_milk;
        
        clear temp temp_a
    else  
        prod_stat(t,j,:) = max(data1.data(1:r,ts-t+1),data2.data(1:r,ts-t+1));   % unit: 10000 tonnes.
    end
    
    
    % adjust in supply tables of China
    temp1 = data2.data(1:r,ts-t+1);
    temp2 = squeeze(prod_stat(t,j,:));
    temp1(isnan(temp1)) = 0;
    temp2(isnan(temp2)) = 0;
    temp_r = temp1./temp2; % fraction of cow milk in total milk
    temp_r(find(temp_r > 1)) = 1;
    
    % distributing total milk to different diary animals according to the
    % live animal numbers
    milk_ls = find(sup_tab(j_c,:) == 1);
    
    sup_tab_chn(t,:,j_c,milk_ls(1)) = temp_r;
    sup_tab_chn(t,:,j_c,milk_ls(2)) = 1-temp_r;
    
    clear temp1 temp2 temp_r milk_ls
end
clear data1 data2 j j_c



filename = "Poultry eggs.xls";
data = importdata(path_chn_ani_prod+filename);
j = 113; % product code 113 is Eggs
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end


filename = "Wool of sheep.xlsx";
data1 = importdata(path_chn_ani_prod+filename);
filename = "Wool of goat.xlsx";
data2 = importdata(path_chn_ani_prod+filename);
j = 114; % product code 114 is Wool
for t = 1:ts
    temp1 = data1.data(1:r,ts-t+1)/10000;   % unit: 10000 tonnes.
    temp2 = data2.data(1:r,ts-t+1)/10000;   % unit: 10000 tonnes.
    temp1(isnan(temp1)) = 0;
    temp2(isnan(temp2)) = 0;
  
    prod_stat(t,j,:) = temp1+temp2;   % should be separated in use table.
    
    clear temp1 temp2 
end
clear data1 data2 j 


filename = "Beef.xls";
data = importdata(path_chn_ani_prod+filename);
j = 115; % product code 115 is Bovine meat
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end


filename = "Mutton.xlsx";
data = importdata(path_chn_ani_prod+filename);
j = 116; % product code 116 is Mutton & Goat Meat
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end



filename = "Pigmeat.xls";
data = importdata(path_chn_ani_prod+filename);
j = 117; % product code 117 is Pigmeat
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end


filename = "Poultry meat.xls";
data = importdata(path_chn_ani_prod+filename);
j = 118; % product code 118 is Poultry meat
for t = 1:ts
    if any([8 9] == t)
        temp = sum(data.data(1:r,[7 10]),2,'omitnan')/sum(sum(data.data(1:r,[7 10]),2,'omitnan'));
        prod_stat(t,j,:) = squeeze(commod_bal(t,j,1))*temp;
        clear temp
    elseif t == 19
        temp = sum(data.data(1:r,[18 20]),2,'omitnan')/sum(sum(data.data(1:r,[18 20]),2,'omitnan'));
        prod_stat(t,j,:) = squeeze(commod_bal(t,j,1))*temp;
        clear temp
    else
        prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
    end
end


filename = "Livestock meat total.xls";
data1 = importdata(path_chn_ani_prod+filename);
filename = "Meat of horses, donkey, mules, camel and rabbit.xlsx";
data2 = importdata(path_chn_ani_prod+filename);
j = 119; % product code 119 is Meat, other
j_c = 77; % product code 119 is Meat, other

% we distribute provicial meat of others into meat of horses, donkey, mules, 
% camel and rabbit by 1995 data
r_ani = data2.data(find(data2.data(:,1) == 1995),2:5)./sum(data2.data(find(data2.data(:,1) == 1995),2:5),2,'omitnan');
temp = find(isnan(r_ani(:,1)));
r_ani(temp,1:3) = repmat([1/3],length(temp),3);
clear temp

ani_ls = find(sup_tab(j_c,:) == 1);
for t = 1:ts
    temp = data1.data(1:r,ts-t+1);   % unit: 10000 tonnes.
    temp(isnan(temp)) = 0;
    temp_1 = squeeze(sum(prod_stat(t,j-(1:4),:),2,'omitnan'));
    temp_2 = temp-temp_1;
    temp_2(find(temp_2<0)) = 0;

    prod_stat(t,j,:) = temp_2;   % unit: 10000 tonnes.
     
    sup_tab_chn(t,:,j_c,ani_ls) = r_ani;
    
    clear temp temp_1 temp_2
end
clear data1 data2 j j_c r_ani ani_ls



filename = "Honey.xlsx";
data = importdata(path_chn_ani_prod+filename);
j = 125; % product code 125 is Honey
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end


filename = "Silk.xlsx";
data = importdata(path_chn_ani_prod+filename);
j = 126; % product code 126 is Silk
temp1 = commod_bal(:,j,1);
for t = 1:ts
    temp2(t) = data.data(end,ts-t+1);
end
temp = temp1./temp2'; % FAO report silk production whereas National Bureau of Statistics reported production of cocoons
ratio = sum(temp([1 4:7 10:24]),'omitnan')/length([1 4:7 10:24]);   % average conversion factor from cocoon to silk
for t = 1:ts
    if t == 2 | t == 3  %  for years without producation data, we distribute the FAO data by the provincial sharing of the nesrest year
        prod_stat(t,j,:) = data.data(1:r,26)./sum(data.data(1:r,26),'omitnan')*commod_bal(t,j,1);
    elseif t == 8 | t == 9
        prod_stat(t,j,:) = sum(data.data(1:r,[20 23]),2,'omitnan')./sum(sum(data.data(1:r,[20 23]),'omitnan'))*commod_bal(t,j,1);
    else
        prod_stat(t,j,:) = data.data(1:r,ts-t+1)*ratio;   % unit: 10000 tonnes.
    end
end
clear temp1 temp2     

        
filename = "Fishery product.xlsx";
data = importdata(path_chn_ani_prod+filename);
j = 127; % product code 127 is Fish, Seafood
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end
clear data


%% China statistics data for production --- Animal in stock
path_chn_ani_stock = path_root+"Data\China_Data\Animal in-stock\";
filename = "Camels.xls";
data = importdata(path_chn_ani_stock+filename);
j = 106; % product code 106 is Camels
for t = 1:ts
    temp = data.data(1:r,ts-t+1)./sum(data.data(1:r,ts-t+1),'omitnan');
    prod_stat(t,j,:) = commod_bal(t,j,1)*temp;;   % distribute the camels slaughted by the in-stock number
    
    clear temp
end
clear data j

%% Estimate animal products from FAOSTAT
filename = "Offal_of_animals.xlsx";
data_offal = importdata(path_faostat+filename);
filename = "Fats_of_animals.xlsx";
data_fat = importdata(path_faostat+filename);
ani_item = {'cattle';'buffaloes';'sheep';'goats';'pigs';'horses';'camels'};

ani_prod = zeros(ts,length(ani_item),4); % animal products; dimension 4: offal(1), fat(2),hides and skin(3),meat production(4)
for i = 1:length(ani_item)
    index1_c = find(contains(data_offal.textdata(2:end,8),ani_item{i}));
    ani_prod(2:end,i,1) = data_offal.data(index1_c,end)/10000;
    
    if i ~= 6
        index2_c = find(contains(data_fat.textdata(2:end,8),ani_item{i}));
        ani_prod(2:end,i,2) = data_fat.data(index2_c,end)/10000;
    end
    clear index1_c index2_c
end

filename = "Hides_skins_of_animals.xlsx";
data_hid_skin = importdata(path_faostat+filename);
filename = "Freshmeat_animal.xlsx";
data_freshm = importdata(path_faostat+filename);
ani_item = {'cattle';'buffalo';'sheep';'goat';'pig';'horse';'camel'};
for i = 1:length(ani_item)
    index1_c = find(contains(data_hid_skin.textdata(2:end,8),ani_item{i}));
    if isempty(index1_c) == 0
        ani_prod(:,i,3) = data_hid_skin.data(index1_c,end)/10000;
    end
    
    index2_c = find(contains(data_freshm.textdata(2:end,8),ani_item{i}));
    ani_prod(:,i,4) = data_freshm.data(index2_c,end)/10000;
    
    clear index1_c index2_c
end
clear data_offal data_fat data_hid_skin data_freshm ani_item

unit_offal = squeeze(ani_prod(:,:,1))./squeeze(ani_prod(:,:,4)); % unit-based offal production of one unit meat prodcution
unit_offal(1,:) = unit_offal(2,:);
unit_fat = squeeze(ani_prod(:,[1:5 7],2))./squeeze(ani_prod(:,[1:5 7],4)); % unit-based fat production of one unit meat prodcution
unit_fat(1,:) = unit_fat(2,:);
unit_hid_skin = squeeze(ani_prod(:,1:4,3))./squeeze(ani_prod(:,1:4,4)); % unit-based hides and skins production of one unit meat prodcution

j_offal = 78; % codes in FABIO_CHN
j_fat = 79; % codes in FABIO_CHN
j_hid_skin = 80; % codes in FABIO_CHN
offal_ls = find(sup_tab(j_offal,:) == 1);
fat_ls = find(sup_tab(j_fat,:) == 1);
hid_skin_ls = find(sup_tab(j_hid_skin,:) == 1);
for t = 1:ts
    temp_beef = squeeze(prod_stat(t,115,:));
    temp_beef(isnan(temp_beef)) = 0;
    
    temp_mutton = squeeze(prod_stat(t,116,:));
    temp_mutton(isnan(temp_mutton)) = 0;
    
    temp_pork = squeeze(prod_stat(t,117,:));
    
    temp_meat_oth = squeeze(prod_stat(t,119,:));
    temp_meat_oth(isnan(temp_meat_oth)) = 0;
    temp_horse = temp_meat_oth.*sup_tab_chn(t,:,77,68)'; % meat from horse
    temp_cam = temp_meat_oth.*sup_tab_chn(t,:,77,71)'; % meat from camel

    ani_meat = horzcat(temp_beef,temp_mutton,temp_pork,temp_horse,temp_cam);
    ani_meat(isnan(ani_meat)) = 0;
    
    clear temp_beef temp_cat temp_buf temp_sheep temp_goat temp_pork temp_horse temp_cam temp_mutton temp_meat_oth
    
    % total offal, fat, hide and skin production in each provinces
    ani_offal = ani_meat.*unit_offal(t,[1 3 5:7]);
    ani_fat = ani_meat(:,[1:3 5]).*unit_fat(t,[1 3 5 6]);
    ani_hid_skin = ani_meat(:,1:2).*unit_hid_skin(t,[1 3]);
    
    % total offal, fat, hide and skin production in each provinces
    prod_stat(t,120,:) = sum(ani_offal,2)*commod_bal(t,120,1)/sum(sum(ani_offal,'omitnan'));
    prod_stat(t,121,:) = sum(ani_fat,2)*commod_bal(t,121,1)/sum(sum(ani_fat,'omitnan'));
    prod_stat(t,122,:) = sum(ani_hid_skin,2)*commod_bal(t,122,1)/sum(sum(ani_hid_skin,'omitnan'));
    
    sup_tab_chn(t,:,j_offal,offal_ls) = ani_offal./sum(ani_offal,2,'omitnan');
    sup_tab_chn(t,:,j_fat,fat_ls) = ani_fat./sum(ani_fat,2);
    sup_tab_chn(t,:,j_hid_skin,hid_skin_ls) = ani_hid_skin./sum(ani_hid_skin,2);
    
    clear ani_meat ani_offal ani_fat ani_hid_skin
end 
clear j_offal j_fat j_hid_skin unit_offal unit_fat unit_hid_skin ani_prod offal_ls fat_ls hid_skin_ls
 
%% China statistics data for production --- Other products
path_chn_oths = path_root+"Data\China_Data\Other products\";

filename = "Oats_Tomato_Chilli(Pimento)_1995.xlsx";
data = importdata(path_chn_oths+filename);
temp_tomato = data.data(1:r,2)./sum(data.data(1:r,2),'omitnan'); % tomato
temp_chilli = data.data(1:r,3)./sum(data.data(1:r,3),'omitnan'); % chilli

filename = "All types of vegetable.xlsx";
data1 = importdata(path_chn_oths+filename);

for t = 1:ts
    temp_1 = squeeze(commod_bal(t,31,1))*temp_tomato;  % distribute to provinces based on 1995 data
    temp_2 = data1.data(1:r,ts-t+1);  % all vegetables
    temp_1(isnan(temp_1)) = 0;
    temp_2(isnan(temp_2)) = 0;
    temp_3 = temp_2-temp_1;  % vegetables other
    temp_3(find(temp_3<0)) = 0;
    
    prod_stat(t,31,:) = temp_1;
    prod_stat(t,33,:) = temp_3;
    prod_stat(t,50,:) = squeeze(commod_bal(t,50,1))*temp_chilli;
    
    clear temp_1 temp_2 temp_3
end
clear data data1


filename = "Coconut.xls";
data = importdata(path_chn_oths+filename);
j = 26; % product code 26 is Coconuts
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1)/10000;   % unit: 10000 tonnes.
end


filename = "All types of cirtus.xls";
data = importdata(path_chn_oths+filename);
j1 = 34; % oranges and mandarines
j2 = 36; % grapefruit
j3 = 37; % citrus, others
for t = 1:ts
    temp_1 = data.data(1:r,ts-t+1)/10000; % mandarines 1, unit: 10000 tonnes
    temp_2 = data.data(41:71,ts-t+1)/10000; % mandarines 2, unit: 10000 tonnes
    temp_3 = data.data(80:110,ts-t+1)/10000; % oranges, unit: 10000 tonnes
    temp_1(isnan(temp_1)) = 0;
    temp_2(isnan(temp_2)) = 0;
    temp_3(isnan(temp_3)) = 0; 
    prod_stat(t,j1,:) = temp_1+temp_2+temp_3;   % unit: 10000 tonnes.
    
    temp_4 = data.data(119:149,ts-t+1)/10000; % grapefruit unit: 10000 tonnes.
    temp_4(isnan(temp_4)) = 0 ;
    prod_stat(t,j2,:) = temp_4;
    
    temp_5 = data.data(158:188,ts-t+1)/10000; % all types of cirtus
    temp_5(isnan(temp_5)) = 0;
    temp_6 = temp_5-temp_4-temp_3-temp_2-temp_1;
    temp_6(find(temp_6 < 0)) = 0;
    prod_stat(t,j3,:) = temp_6;
    
    clear temp_1 temp_2 temp_3 temp_4 temp_5 temp_6 
end


filename = "Bananas.xls";
data = importdata(path_chn_oths+filename);
j = 38; % product code 38 is Bananas
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end


filename = "Apples.xls";
data = importdata(path_chn_oths+filename);
j = 40; % product code 40 is Apples
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end


filename = "Pineapples.xlsx";
data = importdata(path_chn_oths+filename);
j = 41; % product code 41 is Pineapples
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end


filename = "Dates.xlsx";
data = importdata(path_chn_oths+filename);
j = 42; % product code 42 is Red dates
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end


filename = "Grapes.xls";
data = importdata(path_chn_oths+filename);
j = 43; % product code 43 is Grapes
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end


filename = "All types of fruits.xlsx";
data = importdata(path_chn_oths+filename);
j = 44; % product code 44 is fruits others
for t = 1:ts
    temp1 = data.data(1:r,ts-t+1);
    temp2 = squeeze(sum(prod_stat(t,34:43,:),2,'omitnan')); % fruits in the dataset already
    temp1(isnan(temp1)) = 0;
    temp2(isnan(temp2)) = 0;
    temp3 = temp1 - temp2;
    temp3(find(temp3<0)) = 0;
    
    prod_stat(t,j,:) = temp3;   % unit: 10000 tonnes.
    clear temp1 temp2 temp3
end


filename = "Coffee.xls";
data = importdata(path_chn_oths+filename);
j = 45; % product code 45 is Coffee beans
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1)/10000;   % unit: 10000 tonnes.
end


filename = "Tea.xlsx";
data = importdata(path_chn_oths+filename);
j = 47; % product code 47 is Tea
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 tonnes.
end


filename = "Pepper.xlsx";
data = importdata(path_chn_oths+filename);
j = 49; % product code 49 is Pepper
for t = 1:ts
    if t<4
        prod_stat(t,j,:) = data.data(1:r,ts-t+1)/10000;   % unit: 10000 tonnes.
    else
        temp = data.data(1:r,27)./sum(data.data(1:r,27),'omitnan'); % distributing the producation from FAO to each province
        prod_stat(t,j,:) = squeeze(commod_bal(t,j,1))*temp;
    end
end


filename = "Sisal.xls";
data = importdata(path_chn_oths+filename);
j = 56; % product code 56 is Sisal
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1)/10000;   % unit: 10000 tonnes.
end


filename ="Rubber.xls";
data = importdata(path_chn_oths+filename);
j = 60; % product code 60 is Rubber
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1)/10000;   % unit: 10000 tonnes.
end


filename = "Sugar.xls";
data = importdata(path_chn_oths+filename);
j = 67; % product code 6 is 'Sugar, Refined Equiv'
for t = 1:ts
    if any([13:15] == t)
        temp = sum(data.data(1:r,[12 16]),2,'omitnan')./sum(sum(data.data(1:r,[12 16]),'omitnan'));
        prod_stat(t,j,:) = commod_bal(t,j,1)*temp;
    else
        prod_stat(t,j,:) = data.data(1:r,ts-t+1);
    end
    clear temp
end


filename = "Wine.xlsx";
data1 = importdata(path_chn_oths+filename);
data1.data(isnan(data1.data)) = 0;
filename = "Beer.xls";
data2 = importdata(path_chn_oths+filename);
data2.data(isnan(data2.data)) = 0;
filename = "Alcoholic Beverages Total (including beer and wine).xls";
data3 = importdata(path_chn_oths+filename);
data3.data(isnan(data3.data)) = 0;

j1 = 91; % product code 91 is Wine
j2 = 92; % product code 92 is Beer
j3 = 93; % product code 94 is Beverages, Fermented
for t = 1:ts
    prod_stat(t,j2,:) = data2.data(1:r,ts-t+1);   % unit: 10000 m^3.
end

temp3 = data3.data;
temp2 = data2.data;
temp = temp3-temp2;   % other alcoholic beverages including wine
temp(find(temp<0)) = 0;
clear temp3 temp2
for t = 1:ts
    if any([1:10] == t)
        temp0 = data1.data(:,19)./temp(:,19);
        temp0(isnan(temp0)) = 0 ;
        prod_stat(t,j1,:) = temp(1:r,ts-t+1).*temp0(1:r);   % unit: 10000 m^3.
        prod_stat(t,j3,:) = temp(1:r,ts-t+1).*(1-temp0(1:r));  % unit: 10000 m^3.
    elseif any([24:29] == t)
        temp0 = data1.data(:,7)./temp(:,7);
        temp0(isnan(temp0)) = 0 ;
        prod_stat(t,j1,:) = temp(1:r,ts-t+1).*temp0(1:r);
        prod_stat(t,j3,:) = temp(1:r,ts-t+1).*(1-temp0(1:r));
    else
        temp0 = data1.data(:,ts-t+1);
        temp0(isnan(temp0)) = 0 ;
        prod_stat(t,j1,:) = temp0(1:r);
        temp1 = temp(1:r,ts-t+1)-temp0(1:r);
        temp1(find(temp1<0)) = 0;
        prod_stat(t,j3,:) = temp1;
    end
end
clear data1 data2 data3 j1 j2 j3 temp temp0


filename = "Alcohol (96%).xlsx";
data = importdata(path_chn_oths+filename);
j = 95;  % product code 95 is Alcohol, non-food
for t = 1:ts
    if any([1:10] == t)
        temp = data.data(1:r,19)./sum(data.data(1:r,19),'omitnan');
        temp(isnan(temp)) = 0;
        prod_stat(t,j,:) = commod_bal(t,j,1)*temp; % unit: 10000 m^3.
    elseif any([17:19] == t)
        temp = sum(data.data(1:r,[14,10]),2,'omitnan')/sum(sum(data.data(1:r,[14,10]),2,'omitnan'));
        prod_stat(t,j,:) = commod_bal(t,j,1)*temp; % unit: 10000 m^3.
    else
        prod_stat(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 m^3.
    end
end


filename = "Forestry value added.xls";
data = importdata(path_chn_oths+filename);
j = 128; % product code 128 is Wood fuels
for t = 1:ts
    prod_stat(t,j,:) = data.data(1:r,ts-t+1)*10000;   % unit: 10000 Yuan.
end

%% Chinese statistics data of import and export
path_chn_imp_exp = path_root+"Data\China_Data\Import and Export of products\";
fname = "Import.xlsx";
data = importdata(path_chn_imp_exp+fname);

imp_stat = zeros(ts,p0,r); % dimension 1: time period 1990-2018; dimension 2: 130 agricultural products, same with FABIO; dimension 3: 31 procinves
for j = 1:p0
    if isempty(meta.fabio_crop{j,5}(1:end)) == 0
        index_c = find(contains(data.textdata(3:end,2),meta.fabio_crop{j,5}(1:end)));
        if isempty(index_c) == 0
            for t = 1:ts
                if j == 17
                    temp1 = data.data(index_c,ts-t+1)/10000;   % all beans, unit: 10000 tons
                    temp2 = data.data(index_c+2*r,ts-t+1)/10000; % soybeans, unit: 10000 tons
                    temp1(isnan(temp1)) = 0;
                    temp2(isnan(temp2)) = 0;
                    temp3 = temp1-temp2;
                    temp3(find(temp3<0)) = 0;
                    imp_stat(t,j,:) = temp3;   % unit: 10000 tons
                    clear temp1 temp2 temp3
                elseif j == 117
                    imp_stat(t,j,:) = data.data(index_c(1:r),ts-t+1)/10000;   % unit: 10000 tons
                else
                    imp_stat(t,j,:) = data.data(index_c,ts-t+1)/10000;   % unit: 10000 tons
                end
            end
        end
    end
end
clear data


fname = "Export.xlsx";
data = importdata(path_chn_imp_exp+fname);

exp_stat = zeros(ts,p0,r); % dimension 1: time period 1990-2018; dimension 2: 130 agricultural products, same with FABIO; dimension 3: 31 procinves
for j = 1:p0
    if isempty(meta.fabio_crop{j,6}(1:end)) == 0
        index_c = find(contains(data.textdata(3:end,2),meta.fabio_crop{j,6}(1:end)));
        if isempty(index_c) == 0
            for t = 1:ts
                if j == 17
                    temp1 = data.data(index_c,ts-t+1)/10000;   % all beans, unit: 10000 tons
                    temp2 = data.data(index_c+2*r,ts-t+1)/10000; % soybeans, unit: 10000 tons
                    temp1(isnan(temp1)) = 0;
                    temp2(isnan(temp2)) = 0;
                    temp3 = temp1-temp2;
                    temp3(find(temp3<0)) = 0;
                    exp_stat(t,j,:) = temp3;   % unit: 10000 tons
                    clear temp1 temp2 temp3
                elseif j == 117
                    exp_stat(t,j,:) = data.data(index_c(1:r),ts-t+1)/10000;   % unit: 10000 tons
                else
                    exp_stat(t,j,:) = data.data(index_c,ts-t+1)/10000;   % unit: 10000 tons
                end
            end
        end
    end
end
clear data


%% feeds for animals
% animal in-stock at the end of each year
ani_stock = zeros(ts,p0,r);  % animal in-stock

filename = "Cattle and Buffaloes.xlsx";
data1 = importdata(path_chn_ani_stock+filename);
filename = "Buffaloes.xlsx";
data2 = importdata(path_chn_ani_stock+filename);
filename = "Cattle_and_Buffaloes_FAO.xlsx";
data3 = importdata(path_chn_ani_stock+filename);
temp_fao = reshape(data3.data(:,end),ts,2)/10000; % column 1 for buffaloes, column 2 for cattles
j_1 = 97; % product code 97 is Cattle in FABIO.
j_2 = 98; % product code 98 is Buffaloes in FABIO.
for t = 1:ts
    if sum(data2.data(1:r,ts-t+1),'omitnan') == 0
        r_buf = temp_fao(t,1)/sum(temp_fao(t,:)); % national share of buffaloes
        
        temp_prov = data2.data(1:r,12)/sum(data2.data(1:r,12),'omitnan');  % distribute by the data in 2007
        
        temp_buf = sum(data1.data(1:r,ts-t+1),'omitnan')*r_buf.*temp_prov;
        
        clear r_buf temp_prov
    else
        temp_buf = data2.data(1:r,ts-t+1);
    end
    
    temp_buf(isnan(temp_buf)) = 0;
    temp_cat = data1.data(1:r,ts-t+1)-temp_buf;
    temp_cat(find(temp_cat<0)) = 0;
    
    ani_stock(t,j_1,:) = temp_cat;   % unit: 10000 heads.
    ani_stock(t,j_2,:) = temp_buf;
    
    clear temp_buf temp_cat
end
clear j_1 j_2 data1 data2 data3


filename = "sheeps.xls";
data1 = importdata(path_chn_ani_stock+filename);
filename = "Goats.xls";
data2 = importdata(path_chn_ani_stock+filename);
j_1 = 99;  % product code 99 is Sheeps in FABIO.
j_2 = 100; % product code 100 is Goats in FABIO.
for t = 1:ts
    ani_stock(t,j_1,:) = data1.data(1:r,ts-t+1);   % unit: 10000 heads.
    ani_stock(t,j_2,:) = data2.data(1:r,ts-t+1);   % unit: 10000 heads.
end
clear j_1 j_2 data1 data2


filename = "Pigs.xls";
data = importdata(path_chn_ani_stock+filename);
j = 101; % product code 101 is pigs in FABIO.
for t = 1:ts
    ani_stock(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
end
clear j data


filename = "Poultary birds.xls";
data = importdata(path_chn_ani_stock+filename);
j = 102; % product code 102 is Poultary birds in FABIO.
for t = 1:ts
    ani_stock(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
end
clear j data


filename = "Horse.xls";
data = importdata(path_chn_ani_stock+filename);
j = 103; % product code 103 is Horse in FABIO.
for t = 1:ts
    ani_stock(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
end
clear j data


filename = "Donkeys.xls";
data = importdata(path_chn_ani_stock+filename);
j = 104; % product code 104 is Donkey in FABIO.
for t = 1:ts
    ani_stock(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
end
clear j data


filename = "Mules.xls";
data = importdata(path_chn_ani_stock+filename);
j = 105; % product code 105 is Mules in FABIO.
for t = 1:ts
    ani_stock(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
end
clear j data


filename = "Camels.xls";
data = importdata(path_chn_ani_stock+filename);
j = 106; % product code 106 is Camels in FABIO.
for t = 1:ts
    ani_stock(t,j,:) = data.data(1:r,ts-t+1);   % unit: 10000 heads.
end
clear j data


filename = "Animal feeds.xls";
data = importdata(path_root+"Data\"+filename);

prov_graz = [5 26 28:31];  % we assume the main system is grazing in Inner Mongolia (5), Tibet(26), Gansu(28), Qinghai(29), Ningxia(30), Xinjiang(31), 

prod_grai = [1:14];  % product codes for grain crops
prod_oth = [15:p0];  % product codes for other products
feed_c = data.data.IndustrialSystems(1,:); % codes of specific feed crops
ani_c = data.data.IndustrialSystems(:,1);  % codes of specific live animals

% we distribute the natioanl feeds data of crop products to 31 provinces
% according to provincial animal in-stock numbers and per-yr feed
% requirement quantities; we distinguish grazing and industrial system.
feed_ani = zeros(p0,p0,r);  % 1st dimension: animals; 2nd dimension: feed crops; 3rd dimension: provinces
for i = 1:r
    for j = 1:length(ani_c)
        if isnan(ani_c(j)) == 0
            for k = 1:length(feed_c)
                if isnan(feed_c(k)) == 0
                    if any(prov_graz == i)
                        feed_ani(ani_c(j),feed_c(k),i) = data.data.GrazingSystems(j,k); % unit: ton/(head*yr)
                    else
                        feed_ani(ani_c(j),feed_c(k),i) = data.data.IndustrialSystems(j,k);
                    end
                end
            end
            
            for k = 1:p0
                if any(setdiff(prod_grai,feed_c) == k)
                    if any(commod_bal(:,k,6) ~= 0)
                        if any(prov_graz == i)
                            feed_ani(ani_c(j),k,i) = data.data.GrazingSystems(j,7);
                        else
                            feed_ani(ani_c(j),k,i) = data.data.IndustrialSystems(j,7);
                        end
                    end
                elseif any(setdiff(prod_oth,feed_c) == k)
                    if any(commod_bal(:,k,6) ~= 0)
                        if any(prov_graz == i)
                            feed_ani(ani_c(j),k,i) = data.data.GrazingSystems(j,end);
                        else
                            feed_ani(ani_c(j),k,i) = data.data.IndustrialSystems(j,end);
                        end
                    end
                end
            end         
        end
    end
end
feed_ani(98,:,:) = feed_ani(97,:,:); % we assume 98 {'Buffaloes'} share the same feed requirement with cattle.

clear prod_grai prod_oth

%% %% Chinese statistics, product for animal feed
ani_c_a = vertcat(ani_c, [98]);
ani_c_a(isnan(ani_c_a)) = [];

feed_ani_31 = zeros(ts,p0,p0,r);  % 2nd dimension: live animals; 3rd dimension: feed crops;
for t = 1:ts
    for j = 1:p0  % for each crop as feeds
        if commod_bal(t,j,6) ~= 0
            temp_feed_fao = commod_bal(t,j,6);
            for i = 1:r
                for k = 1:length(ani_c_a)
                    temp_ani_1 = prod_stat(t,ani_c_a(k),i); % animal slaughted
                    temp_ani_2 = ani_stock(t,ani_c_a(k),i); % animal in-stock
                    temp_ani_1(isnan(temp_ani_1)) = 0;
                    temp_ani_2(isnan(temp_ani_2)) = 0;
                    temp_feed = feed_ani(ani_c_a(k),j,i);
                    temp_feed(isnan(temp_feed)) = 0;
                    temp_feed_v(k,i) = temp_ani_1*temp_feed*0.5+temp_ani_2*temp_feed; % we assume animals slaughted contributing 50% of feeds
                    clear temp_ani_1 temp_ani_2 temp_feed
                end
            end

            feed_ani_31(t,ani_c_a,j,:) = temp_feed_v*temp_feed_fao/sum(sum(temp_feed_v,'omitnan'));
            
            clear temp_feed_fao temp_feed_v
        end                    
    end
end

clear ani_c ani_c_a feed_c



%% Chinese statisticas, product for seeds

path_chn_seed = path_root+"Data\China_Data\Crop sown areas in provinces\";
sown_area = zeros(ts,p0,r); % dimension 1: time period 1990-2018; dimension 2: 130 agricultural products, same with FABIO; dimension 3: 31 procinves
seed_est = zeros(ts,p0,r);

for i = 1:r
    fname = sprintf("%d.xls",i);
    temp = importdata(path_chn_seed+fname);

    for j = 1:p0
        if isempty(meta.fabio_crop{j,4}(1:end)) == 0 
            index_c = find(contains(temp.textdata(5:44,1),meta.fabio_crop(j,4)));
            if isempty(index_c) == 0
                for t = 1:ts
                    if j == 13
                        temp1 = temp.data(index_c(1),ts-t+1); % production of all tubers
                        temp2 = temp.data(index_c(1)+1,ts-t+1); % production of potatoes
                        temp1(isnan(temp1)) = 0;
                        temp2(isnan(temp2)) = 0;
                        temp3 = temp1-temp2; % production of tubers excluding potatoes as 'Roots, other'
                        temp3(find(temp3<0)) = 0;
                        
                        sown_area(t,j,i) = temp3;   % unit: 1000 hectares
                        clear temp1 temp2 temp3
                    elseif j == 17
                        temp1 = temp.data(index_c(1),ts-t+1); % production of all beans
                        temp2 = temp.data(index_c(1)+3,ts-t+1); % production of all soybeans
                        temp1(isnan(temp1)) = 0;
                        temp2(isnan(temp2)) = 0;
                        temp3 = temp1-temp2; % production of beans excluding soybeans as 'beans other'
                        temp3(find(temp3<0)) = 0;
                        
                        sown_area(t,j,i) = temp3;   % unit: 1000 hectares
                        clear temp1 temp2 temp3
                    elseif j == 30
                        temp1 = temp.data(index_c(1),ts-t+1); % production of all oil crops
                        temp2 = sum(temp.data(index_c(1)+[1:4],ts-t+1),'omitnan'); % production of peanuts, sunflowers, rapeseeds, sesame
                        temp1(isnan(temp1)) = 0;
                        temp2(isnan(temp2)) = 0;
                        temp3 = temp1-temp2; % Oil crop, others
                        temp3(find(temp3<0)) = 0;
                        
                        sown_area(t,j,i) = temp3;   % unit: 1000 hectares
                        clear temp1 temp2 temp3
                    else
                        sown_area(t,j,i) = temp.data(index_c(1),ts-t+1);   % unit: 10000 tons
                    end
                end
            end
        end
    end
    
    clear temp
end

filename = "Harvested areas of barley.xlsx";
data = importdata(path_chn_seed+filename);
j = 3;
temp_r = squeeze(sown_area(13,j,:))./sum(sown_area(13,j,:),'omitnan');
sown_area(1:12,j,:) = data.data(1:12,10)/1000.*temp_r';
clear data j temp_r

j_1 = 63; % sown area of cottonseed same with that of cotton lint
j_2 = 96; 
sown_area(:,j_1,:) = sown_area(:,j_2,:);
clear j_1 j_2

% estimation for seed use
for j = 1:p0
    temp = squeeze(sown_area(:,j,:)); % sown areas from NBSC
    temp_a = vertcat(temp(2:ts,:),zeros(1,r)); % seed in year t for sown area in year t+1
    temp_r = temp_a./sum(temp_a,2,'omitnan');
    
    temp_1 = squeeze(commod_bal(:,j,7)).*temp_r;
    
    seed_est(:,j,:) = temp_1;
end

% seeds for eggs
seed_est(:,113,:) = squeeze(prod_stat(:,113,:))*0.024;

seed_est(isnan(seed_est)) = 0;


% change the format of sown_area from FABIO to FABIO_CHN
prod_c(isnan(prod_c)) = 130;
sown_area_a = zeros(ts,size(prod_c,1),r);
for j_c = 1:size(prod_c,1)
    sown_area_a(:,j_c,:) = squeeze(sum(sown_area(:,prod_c(j_c,:),:),2));
end

data = importdata(path_chn_seed+'Sown areas of other products.xlsx');
for j_c = 1:size(prod_c,1)
    index_j_c = find(data.data.Sheet1(:,1)==j_c);
    if isempty(index_j_c) == 0
        sown_area_a(23,j_c,:) = data.data.Sheet1(index_j_c,4);
    end
    clear index_j_c
end
sown_area_a(isnan(sown_area_a)) = 0;

%% processing of products with multiple procedures,estimate by least-squares optimization
proc_est = zeros(t,p0,r); % dimension 1: time period 1990-2018; dimension 2: 130 agricultural products, same with FABIO; dimension 3: 31 procinves

% 1. sugarcan and sugar beets for sugar
temp_prod = [15 16]; % 15 for sugarcane; 16 for sugarbeets
temp_tcf = [0.1104 0.1012];  % technical conversion fators of sugarcane and sugarbeets to refined sugar
temp_out = [67];  % output product 

v_n = length(temp_prod)*r; % number of variables

C = zeros(r,v_n);
for m = 1:r
    C(m,([1:length(temp_prod)]-1)*r+m) = temp_tcf;
end

for t = 1:ts
    % profuction of sugar balanced by FAO data
    d = squeeze(prod_stat(t,temp_out,:))*commod_bal(t,temp_out,1)/squeeze(sum(prod_stat(t,temp_out,:),'omitnan'));
    d(isnan(d)) = 0;
    
    % constraint 1
    A = -eye(v_n);
    b = zeros(v_n,1);
    
    % constrint 2
    Aeq = zeros(length(temp_prod),v_n);
    for i = 1:length(temp_prod)
        Aeq(i,(i-1)*r+(1:r)) = 1;
    end
    beq = squeeze(commod_bal(t,temp_prod,9));
    
    temp_x = lsqlin(C,d,A,b,Aeq,beq);
    
    for i = 1:length(temp_prod)
        proc_est(t,temp_prod(i),:) = temp_x((i-1)*r+(1:r));
    end
end
clear temp_prod temp_tcf temp_out v_n C d A b Aeq beq temp_x

% 2. barley, maize, millet, sorghum for beer
temp_prod = [3 4 7 8]; 
temp_tcf = [4.48 5.5 4.5 4.85];  % technical conversion fators of sugarcane and sugarbeets to refined sugar
temp_out = [92];  % output product 

v_n = length(temp_prod)*r; % number of variables

C = zeros(r,v_n);
for m = 1:r
    C(m,([1:length(temp_prod)]-1)*r+m) = temp_tcf;
end

for t = 1:ts
    % profuction of beer balanced by FAO data
    d = squeeze(prod_stat(t,temp_out,:))*commod_bal(t,temp_out,1)/squeeze(sum(prod_stat(t,temp_out,:),'omitnan'));
    d(isnan(d)) = 0;
    
    % constraint 1
    A = -eye(v_n);
    b = zeros(v_n,1);
    
    % constrint 2
    Aeq = zeros(length(temp_prod),v_n);
    for i = 1:length(temp_prod)
        Aeq(i,(i-1)*r+(1:r)) = 1;
    end
    beq = squeeze(commod_bal(t,temp_prod,9));
    
    beq(2) = beq(2) - commod_bal(t,80,1)/0.027;  % 'Maize Germ Oil' from maize
    beq(find(beq<0)) = 0;
    
    temp_x = lsqlin(C,d,A,b,Aeq,beq);
    
    for i = 1:length(temp_prod)
        proc_est(t,temp_prod(i),:) = temp_x((i-1)*r+(1:r));
    end
end
clear temp_prod temp_tcf temp_out v_n C d A b Aeq beq temp_x

% 3.rice, wheat, apple for beverage alcohol
temp_prod = [1 2 40]; 
temp_tcf = [0.5025 0.68 0.65];  % technical conversion fators of sugarcane and sugarbeets to refined sugar
temp_out = [93];  % output product 

v_n = length(temp_prod)*r; % number of variables

C = zeros(r,v_n);
for m = 1:r
    C(m,([1:length(temp_prod)]-1)*r+m) = temp_tcf;
end

for t = 1:ts
    % profuction of beer balanced by FAO data
    d = squeeze(prod_stat(t,temp_out,:))*sum(commod_bal(t,temp_out:temp_out+1,1))/squeeze(sum(prod_stat(t,temp_out,:),'omitnan'));
    d(isnan(d)) = 0;
    
    % constraint 1
    A = -eye(v_n);
    b = zeros(v_n,1);
    
    % constrint 2
    Aeq = zeros(length(temp_prod),v_n);
    for i = 1:length(temp_prod)
        Aeq(i,(i-1)*r+(1:r)) = 1;
    end
    beq = squeeze(commod_bal(t,temp_prod,9));
    
    beq(1) = beq(1) - commod_bal(t,79,1)/0.012;
    beq(find(beq<0)) = 0;
    
    temp_x = lsqlin(C,d,A,b,Aeq,beq);
    
    for i = 1:length(temp_prod)
        proc_est(t,temp_prod(i),:) = temp_x((i-1)*r+(1:r));
    end
end
clear temp_prod temp_tcf temp_out v_n C d A b Aeq beq temp_x

% 4. grape for wine   ***************************** question marks
c_grape = 43;
tcf_g_w = 0.74;
c_wine = 91;
temp_prod_w = squeeze(commod_bal(:,c_wine,1))./squeeze(sum(prod_stat(:,91,:),3,'omitnan'));
temp_proc_g_w = squeeze(prod_stat(:,c_wine,:)).*temp_prod_w/tcf_g_w;
temp_proc_g_w = temp_proc_g_w.*(squeeze(commod_bal(:,c_grape,9))./sum(temp_proc_g_w,2,'omitnan'));
proc_est(:,c_grape,:) = temp_proc_g_w;

clear c_grape tcf_g_w c_wine temp_proc_g_w tmep_prod_w

%% Chinese statistics population
path_chn_pop = path_root+"Data\China_Data\Population\";

data = importdata(path_chn_pop+"Population total 1986-2018.xls");

pop = zeros(ts,r);
for t = 1:ts
    pop(t,:) = data.data(:,ts-t+1);   % unit: 10000 capita
end


%% Chinese data, trade-cost between capital cities

filename = 'Inter-provincial trade cost.xlsx';
data = importdata(path_root+"Data\China_Data\Inter-provincial trade cost\"+filename);
cost = data.data;   % unit: Yuan per ton

%% save data

save(path_root+"fabio_chn"+version+".mat",'commod_bal','cost','data_avi_fao','element','exp_stat',...
    'feed_ani_31','imp_stat','meta','pop','prod_c','prod_stat','proc_type','proc_est','seed_est',...
    'sown_area','sown_area_a','sup_tab','sup_tab_chn','tcf_fao','use_tab')

