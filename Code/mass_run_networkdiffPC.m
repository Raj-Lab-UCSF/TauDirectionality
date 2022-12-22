close all
clearvars -except Xcomp
%clearvars -except norms_struc corrs_struc norms_gene corrs_gene
%generate theoretical seeding by putting 1s in each area of the brain
%alterately

%Load connectivity matrix
    %load 212 x 212 genetic connectivity matrix
    load ('/home/chx2002/naomicode/Copythis/connectivity212.mat');
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              %load 212 x 212 structural connectivity matrix
    %load ('/Users/xia/Google Drive/Code/diseaseSpread/ipsi_contraAv.mat');
    %Coexpress = ipsi_contraAv;
    
    %load proximity connectivity matrix
    %load ('/home/chx2002/naomicode/Copythis/proxconMat.mat');
    %Coexpress = proxconMat;
    
for i=1:size(Coexpress,1)
    Coexpress(i,i)=0;
end
Lgen=genLplcns(Coexpress);
L=Lgen;

%specify your parameters
t=0:0.1:10; 
b=1;
    
%%%%Load seeding resolution and disease end states

%Load EC-APP data
[num,text] = xlsread('/home/chx2002/naomicode/Copythis/loadDF.xlsx','ECAPPLoadDF');
range = 1:10; %you can choose whether to include 10, the PC
APP_6mo = num(range,4);
APP_13mo = num(1:7,8);
arealabel = num(range,1);

%Load J20 data
%[num,text] = xlsread('/Users/xia/Google Drive/Code/loadDF.xlsx','J20LoadDF');
%range = 1:9; %J20 data does not have PC
%APP_6mo = num(range,4); %seed
%APP_13mo = num(1:6,8); %
%arealabel = num(range,1);
    
%seeding vector
%seed=zeros(size(L,1),1);
%seed(arealabel)=APP_6mo;
   

    %if you want to seed 100 different permutations with a certain location
    %restriction
    
    seed = zeros(212,1000);
    var = size(seed,2);
    seed_ind = [];
    for i=1:var
        seed_ind = seed_phys(seed_ind);
        seed(seed_ind(:,1),i) = seed_ind(:,2);
    end
    
    
    
%Structure size adjustments
%weights of constituent areas to area you compare with

%matVoxWeights = num(range,3)./num(range,6);
matVoxWeights = ones(10,1);

%size of final-partition brain structure relative to the whole brain
matArea10 = num(range,6);

%choose filter
mat_filter6 = [1 3 4 5 6 7];
mat_filter7 = [1 3 4 5 6 7 10];
if length(range)==9
    choosefilter = mat_filter6;
else
    choosefilter = mat_filter7;
end

matArea7 = matArea10(choosefilter);

for i=1:length(t)
    matArea7_T(:,i)=matArea7;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%
%run network diffusion model for each seeding vector variation
for i=1:size(seed,2)
    Xo=seed(:,i);
    Xt=spreadModeldelta(Xo,t,L,b);
    Xtcorr2=Xt(arealabel,:); %you get Xt over time with desired resolution of area labels
     
%       %Unionize Xtcorr2 from 10 regions -> 7 regions with appropriate
%       weights
        Xtcorr(1,:) = matVoxWeights(1)*Xtcorr2(1,:) + matVoxWeights(2)*Xtcorr2(2,:); %EC
        Xtcorr(2,:) = Xtcorr2(3,:); %Dentate Gyrus
        Xtcorr(3,:) = Xtcorr2(4,:);  %CA1
        Xtcorr(4,:) = Xtcorr2(5,:);  %CA2
        Xtcorr(5,:) = Xtcorr2(6,:);  %CA3
        Xtcorr(6,:) = matVoxWeights(7)*Xtcorr2(7,:) + matVoxWeights(8)*Xtcorr2(8,:) + matVoxWeights(9)*Xtcorr2(9,:);%RC
        
        if length(range)~=9
        Xtcorr(7,:) = Xtcorr2(10,:); %PC
        end
        
        %Normalize it to the area of the structure
        if t~=0
        Xtcorr=Xtcorr./(matArea7_T)*100;
        end
        
    [norms,corrs]=min_norm(Xtcorr,APP_13mo);%you get 1 value per Xt 
    
    %select either norm or correlation
    timeindex=find(norms==min(norms));
    norms(find(isnan(norms)))=0;
        
    %timeindex=find(corrs==max(corrs));
    %corrs(find(isnan(corrs)))=0;
    
    max_R_Xt(:,i)=Xtcorr(:,timeindex);    
    max_T_time(i)=timeindex;
    
        maxnorm(i)=norms(timeindex); %store norm
            
    
        %maxcorrs(i)=corrs(timeindex); %store corrs
        
    figure(1);
    corrs(find(isnan(corrs)))=0;
    plot(t,norms);
    norma_all(:,i)=norms;
    Xt_all(:,:,i)=Xt;
    Xtcorr_all(:,i)=corrs;
    hold on %plot all graphs on a single figure
    
end



    

    