%Calculates Laplacian Matrix of graph

%Output/s: 

%lapl - Laplacian matrix of graph

%Input/s: 

%n - Used to calculate the number of vertices in the graph; N
%N - Number of vertices of the graph
%gtype - Graph type input as string ('comp','hyper','bi')
%atype - Algorithm type input as string ('normal','chiral')
%alpha - Phase angle in CQW

function lapl = graphlapl(n,N,gtype,atype,alpha)


switch atype
    
    %Laplacian matrices for standard QSSA
    case 'normal'
        switch gtype

            case 'comp'
                adj=ones(N)-eye(N);
                lapl=adj-(N-1)*eye(N);
            case 'hyper'

                x=[0,1;1,0];
                adj=x;
                for ii=1:n-1 
                   adj=kron(adj,eye(2))+kron(eye(2^(ii)),x); 
                end
                lapl=adj-n*eye(N);

            case 'bi'
                adj=[zeros(n),ones(n);ones(n),zeros(n)];

                lapl=adj-n*eye(N);

            otherwise

                error('Graph type not correctly chosen')

        end
        
    %Laplacian matrices for chiral QSSA    
    case 'chiral'

        switch gtype

            case 'comp'

                adj=zeros(N);

                crow=[0,exp(-1i*alpha)*ones(1,(N-1)/2),exp(1i*alpha)*ones(1,(N-1)/2)];

                for jj=1:N
                    adj(jj,:)=circshift(crow,jj-1);
                end

                lapl=adj-(N-1)*cos(alpha)*eye(N);

            case 'hyper'

                r4=[0,exp(1i*alpha),exp(-1i*alpha),0;exp(-1i*alpha),0,0,exp(1i*alpha);exp(1i*alpha),0,0,exp(-1i*alpha);0,exp(-1i*alpha),exp(1i*alpha),0];
                adj=r4;

                for jj=1:n/2-1 
                   adj=kron(adj,eye(4))+kron(eye(4^(jj)),r4); 
                end

                lapl=adj-n*cos(alpha)*eye(N);

            case 'bi'
                birow = zeros(1,n);
                biM=zeros(n);
                for jj=1:n
                    birow(jj) = exp(-1i*((-1)^jj)*alpha);
                end

                for jj=1:n
                    biM(jj,:)=circshift(birow,jj-1);
                end

                adj = [zeros(n),biM;biM',zeros(n)];

                lapl = adj-n*cos(alpha)*eye(N);

                otherwise

                error('Graph type not correctly chosen')
        end
end


end

