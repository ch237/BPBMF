%%
%%small data
clear all;close all;
load kinship_idxi.mat
% load umls_idxi_new
% load nation_idxi
xi_tensor=double(xi_tensor);
xi_tensor(isnan(xi_tensor))=0;
para.R=max(id{3});% # of ralations
para.N=max(max(id{1}),max(id{2}));% # of entities
para.K=10;%rank
para.TrainFrac=0.9;
para.batchfrac=0.05;
para.itermax=300;
para.burnin=floor(para.itermax*3/4);
RelationOnlineMean1
%%
%%large data
clear all; close all;
load NELL50K.mat
para.R=max(max(idtest{3}),max(idtrain{3}));% # of ralations
para.N=max(max(max(idtrain{1}),max(idtrain{2})),max(max(idtest{1}),max(idtest{2})));% # of entities
para.K=10;%rank
para.TrainFrac=0.9;
para.batchfrac=1;
para.itermax=100;
para.burnin=floor(para.itermax*3/4);
Large_onlineEM1