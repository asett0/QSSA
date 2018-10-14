%Obtains data for chiral QSSA Algorithm for varying alpha

%Output/s: 

%data - Matrix containing success probabilities and run times for chosen
%graph type for different values of alpha. 

%Input/s: 

%gtype - Graph type input as string ('comp','hyper','bi')

function data=runalpha(gtype)

%Random values of alpha between 0 and pi for all graphs
alphavec=zeros(100,1);
for ii=2:length(alphavec)
    alphavec(ii)=pi*rand;
end

%Chosen graph size for each type of graph
switch gtype
    
    case 'comp'
        n=201;

    case 'hyper'
        n=8;
    
    case 'bi'
        n=100;

    otherwise
        
        error('Graph type not correctly chosen')

end

%search function is run for all values of alpha

data=zeros(length(alphavec),3);

for kk=1:length(alphavec)
    [data(kk,1),data(kk,2),data(kk,3)]=search(n,gtype,'chiral',alphavec(kk));
end

datatemp = data(any(data,2),:);
alphavec=alphavec(any(data,2),:);
data=[datatemp,alphavec];

[~,idx] = sort(data(:,4)); % sort alpha in ascending order
data = data(idx,:);   % sort the whole matrix using the sort indices

end
