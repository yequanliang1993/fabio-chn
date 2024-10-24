%% this is the Matlab codes to calculate trade structure within China for FABIO-CHN

% Created: Quanliang Ye
% Date: 22-October-2024
% Email: yequanliang1993@gmail.com
% Version: 2.3.0

% Note:
%
% This code is used to generate trade structure of agriculture products in
% benchmarking years (2007, 2010, 2012, 2015, 2017). 
% We use trade structures in these benchmarking years (of which MRIO tables
% are available) as proxy to estimate inter-provincial trade of key
% agri-food products during the period of 1990-2013

% This approach maybe cause high uncertainty. The reason is data of 
% international import is integrated in intermediate inputs. But, it is 
% better than applying no constrains at all.

%% Configure default information
clear

% root path
path_root = pwd+"\";

% current version
path_parts = split(path_root,'\');
version ="_"+path_parts{end-1};
clear path_parts

% path to MRIOTs of China
path_chn_mrio = path_root+"Data\China_Data\MRIOTs\";

%% benchmarking year 2007
load(path_chn_mrio+"MRIO_CHN_2007.mat")

IO.Z(find(IO.Z<0)) = 0;
IO.Y(find(IO.Y<0)) = 0;
IO.Exp(find(IO.Exp<0)) = 0;

r = length(meta.regions);
s = length(meta.sectors);

sum_mat_z = zeros(r*s,r);
sum_mat_y = zeros(r*5,r);
for m = 1:r
  sum_mat_z((m-1)*s+(1:s),m) = 1;
  sum_mat_y((m-1)*5+(1:5),m) = 1;
end

% for different trade patterns
trade_z = zeros(r,r);
trade_y = zeros(r,r);
trade_exp = zeros(r,1);
for m = 1:r
    trade_z(m,:) = IO.Z((m-1)*s+1,:)*sum_mat_z;
    trade_y(m,:) = IO.Y((m-1)*s+1,:)*sum_mat_y;
    trade_exp(m) = IO.Exp((m-1)*s+1);
end
trade_z = trade_z-diag(diag(trade_z));
trade_y = trade_y-diag(diag(trade_y));

tibet = 26;
exp_str = horzcat(trade_z+trade_y,trade_exp)./sum(horzcat(trade_z+trade_y,trade_exp),2);
exp_str(tibet+1:r+1,:) = exp_str(tibet:r,:);
exp_str(tibet,:) = 0;
exp_str(:,tibet+1:r+2) = exp_str(:,tibet:end);
exp_str(:,tibet) = 0;

imp_str = (trade_z+trade_y)./sum(trade_z+trade_y);
imp_str(tibet+1:r+1,:) = imp_str(tibet:end,:);
imp_str(tibet,:) = 0;
imp_str(:,tibet+1:r+1) = imp_str(:,tibet:end);
imp_str(:,tibet) = 0;

tr_str = horzcat(trade_z+trade_y,trade_exp)./...
    sum(sum(horzcat(trade_z+trade_y,trade_exp)));
tr_str(tibet+1:r+1,:) = tr_str(tibet:r,:);
tr_str(tibet,:) = 0;
tr_str(:,tibet+1:r+2) = tr_str(:,tibet:end);
tr_str(:,tibet) = 0;


tr_str_2007.exp_str = exp_str;
tr_str_2007.imp_str = imp_str;
tr_str_2007.tr_str = tr_str;

clear tr_str exp_str i imp_str IO m meta r s sum_mat_y sum_mat_z trade_exp trade_y year yr trade_z

%% benchmarking year 2010
load(path_chn_mrio+"MRIO_CHN_2010.mat")

r = length(meta.regions);
s = length(meta.sectors);

IO.Z(find(IO.Z<0)) = 0;
IO.Y(find(IO.Y<0)) = 0;
IO.Exp(find(IO.Exp<0)) = 0;

sum_mat_z = zeros(r*s,r);
sum_mat_y = zeros(r*2,r);
for m = 1:r
  sum_mat_z((m-1)*s+(1:s),m) = 1;
  sum_mat_y((m-1)*2+(1:2),m) = 1;
end

% for different trade patterns
trade_z = zeros(r,r);
trade_y = zeros(r,r);
trade_exp = zeros(r,1);
for m = 1:r
    trade_z(m,:) = IO.Z((m-1)*s+1,:)*sum_mat_z;
    trade_y(m,:) = IO.Y((m-1)*s+1,:)*sum_mat_y;
    trade_exp(m) = IO.Exp((m-1)*s+1);
end
trade_z = trade_z-diag(diag(trade_z));
trade_y = trade_y-diag(diag(trade_y));

tibet = 26;
exp_str = horzcat(trade_z+trade_y,trade_exp)./sum(horzcat(trade_z+trade_y,trade_exp),2);
exp_str(tibet+1:r+1,:) = exp_str(tibet:r,:);
exp_str(tibet,:) = 0;
exp_str(:,tibet+1:r+2) = exp_str(:,tibet:end);
exp_str(:,tibet) = 0;

imp_str = (trade_z+trade_y)./sum(trade_z+trade_y);
imp_str(tibet+1:r+1,:) = imp_str(tibet:end,:);
imp_str(tibet,:) = 0;
imp_str(:,tibet+1:r+1) = imp_str(:,tibet:end);
imp_str(:,tibet) = 0;

tr_str = horzcat(trade_z+trade_y,trade_exp)./...
    sum(sum(horzcat(trade_z+trade_y,trade_exp)));
tr_str(tibet+1:r+1,:) = tr_str(tibet:r,:);
tr_str(tibet,:) = 0;
tr_str(:,tibet+1:r+2) = tr_str(:,tibet:end);
tr_str(:,tibet) = 0;


tr_str_2010.exp_str = exp_str;
tr_str_2010.imp_str = imp_str;
tr_str_2010.tr_str = tr_str;

clear exp_str i imp_str IO m meta r s sum_mat_y sum_mat_z trade_exp trade_y year yr trade_z

%% benchmarking year 2012
load(path_chn_mrio+"MRIO_CHN_2012.mat")

IO.Z(find(IO.Z<0)) = 0;
IO.Y(find(IO.Y<0)) = 0;
IO.Exp(find(IO.Exp<0)) = 0;

r = length(meta.regions);
s = length(meta.sectors);

sum_mat_z = zeros(r*s,r);
sum_mat_y = zeros(r*5,r);
for m = 1:r
  sum_mat_z((m-1)*s+(1:s),m) = 1;
  sum_mat_y((m-1)*5+(1:5),m) = 1;
end

% for different trade patterns
trade_z = zeros(r,r);
trade_y = zeros(r,r);
trade_exp = zeros(r,1);
for m = 1:r
    trade_z(m,:) = IO.Z((m-1)*s+1,:)*sum_mat_z;
    trade_y(m,:) = IO.Y((m-1)*s+1,:)*sum_mat_y;
    trade_exp(m) = IO.Exp((m-1)*s+1);
end
trade_z = trade_z-diag(diag(trade_z));
trade_y = trade_y-diag(diag(trade_y));

tibet = 26;
exp_str = horzcat(trade_z+trade_y,trade_exp)./sum(horzcat(trade_z+trade_y,trade_exp),2);
% exp_str(tibet+1:r+1,:) = exp_str(tibet:r,:);
% exp_str(tibet,:) = 0;
% exp_str(:,tibet+1:r+2) = exp_str(:,tibet:end);
% exp_str(:,tibet) = 0;

imp_str = (trade_z+trade_y)./sum(trade_z+trade_y);
% imp_str(tibet+1:r+1,:) = imp_str(tibet:end,:);
% imp_str(tibet,:) = 0;
% imp_str(:,tibet+1:r+1) = imp_str(:,tibet:end);
% imp_str(:,tibet) = 0;

tr_str = horzcat(trade_z+trade_y,trade_exp)./...
    sum(sum(horzcat(trade_z+trade_y,trade_exp)));
% tr_str(tibet+1:r+1,:) = tr_str(tibet:r,:);
% tr_str(tibet,:) = 0;
% tr_str(:,tibet+1:r+2) = tr_str(:,tibet:end);
% tr_str(:,tibet) = 0;

tr_str_2012.exp_str = exp_str;
tr_str_2012.imp_str = imp_str;
tr_str_2012.tr_str = tr_str;

% assume trade structure in Tibet in 2007 and 2012 is same as the structure
% in 2012
tr_str_2007.exp_str(tibet,:) = tr_str_2012.exp_str(tibet,:);
tr_str_2010.exp_str(tibet,:) = tr_str_2012.exp_str(tibet,:);

tr_str_2007.imp_str(:,tibet) = tr_str_2012.imp_str(:,tibet);
tr_str_2010.imp_str(:,tibet) = tr_str_2012.imp_str(:,tibet);

tr_str_2007.tr_str(tibet,:) = tr_str_2012.tr_str(tibet,:);
tr_str_2007.tr_str(:,tibet) = tr_str_2012.tr_str(:,tibet);
tr_str_2007.tr_str = tr_str_2007.tr_str*1/sum(sum(tr_str_2007.tr_str));
tr_str_2010.tr_str(tibet,:) = tr_str_2012.tr_str(tibet,:);
tr_str_2010.tr_str(:,tibet) = tr_str_2012.tr_str(:,tibet);
tr_str_2010.tr_str = tr_str_2010.tr_str*1/sum(sum(tr_str_2010.tr_str));

clear exp_str i imp_str IO m meta r s sum_mat_y sum_mat_z trade_exp trade_y year yr trade_z

%% save data
path_output = path_root+"trade_structure_benchmarking"+version+".mat";
save(path_output, 'tr_str_2007','tr_str_2010','tr_str_2012')
