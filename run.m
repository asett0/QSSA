%Obtains data for standard and chiral QSSA Algorithm

%Output/s: 

%data - Matrix containing success probabilities and run times for chosen
%graph type for different values of n. 

%Input/s: 

%gtype - Graph type input as string ('comp','hyper','bi')
%atype - Algorithm type input as string ('normal','chiral')


function data=run(gtype,atype)

%Values of n for each graph type
switch gtype
    
    case 'comp'
        n=301;
        index=3:2:n;

    case 'hyper'
        n=12;
        index=12:2:n; 
        
    case 'bi'
        n=250;
        index=2:2:n;

    otherwise
        
        error('Graph type not correctly chosen')

end


%search function is run for all values of n

switch atype
    
    case 'normal'
        
        data=zeros(length(index),3);
        
        for ii=1:length(index)
            [data(ii,1),data(ii,2),data(ii,3)]=search(index(ii),gtype,'normal',0);
        end
        
    case 'chiral'
        
        %Values of alpha for all graphs
        alphatot=100;
        alphavec=zeros(alphatot,1);

        for ii=1:length(alphavec)
            alphavec(ii)=pi/alphatot;
        end
        
        data=zeros(length(index),4);
        datatemp=zeros(length(alphavec),4);
        
        for ii=1:length(index)
    
            for jj=1:length(alphavec)

                [datatemp(jj,1),datatemp(jj,2),datatemp(jj,3)]=search(index(ii),gtype,'chiral',alphavec(jj));

            end
            
            %alpha corresponding to maximum success probability is chosen
            %as optimal alpha
            [max_p,max_ind]=max(datatemp(:,1));


            data(ii,1)=max_p;
            data(ii,2)=datatemp(max_ind,2);
            data(ii,3)=datatemp(max_ind,3);
            data(ii,4)=alphavec(max_ind);

        end
end

end
