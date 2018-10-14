%Computes successprobability and run time of normal QSSA or chiral QSSA

function [pmax,tmax,gamma]=search(n,gtype,atype,alpha)

%Output/s: 

%pmax - Success probability
%tmax - Run time 
%gamma - Optimal gamma 

%Input/s: 

%n - Used to calculate the number of vertices in the graph; N
%gtype - Graph type input as string ('comp','hyper','bi')
%atype - Algorithm type input as string ('normal','chiral')
%alpha - Phase angle in CQW


%Define number of vertices
if strcmp('comp',gtype) || strcmp('cycle',gtype)
    N=n; 
elseif strcmp('bi',gtype)
    N=2*n;
elseif strcmp('hyper',gtype)
    N=2^n;
end

Fs=4;       %Sampling frequency
T=1/Fs;     %Sampling period
L=N*Fs;     %Length of Signal

%Probability data
p1 = zeros(L,1);
    
%Laplacian Matrix
lapl=graphlapl(n,N,gtype,atype,alpha);

%Oracle
orac=zeros(N);
orac(1,1)=1;

%Optimal gamma
gamma=gcalc(lapl);

%Values of alpha for which gamma>10^5 are rejected
if(abs(gamma)>10^5)
    gamma=0;
    pmax=0;
    tmax=0;
else
    %Hamiltonian
    h=-gamma*lapl-orac;

    %eigenvectors and eigenvalues
    [V,D] = eig(h);
    D=diag(D);

    %Probability Function
    for ii=1:L

        t=(ii-1)*T;

        exD=exp(-1i*D*t);

        upart=bsxfun(@times,V',exD);

        u1=V(1,:)*upart;

        psi1=u1*ones(N,1);

        psi1=psi1/sqrt(N);

        %Probability matrix
        p1(ii)=(abs(psi1))^2;
    end

    %Power spectrum
    M=repmat(mean(p1),L,1);
    pmean=p1-M;

    Y=fft(pmean);
    P2=abs(Y./L);
    P1=P2(1:L/2+1);
    P1(2:end-1)=2*P1(2:end-1);
    f=Fs*(0:(L/2))/L;

    [maxiM,fmaxrow]=max(P1);
    
    %Maximum/periodic frequency is calculated
    fmax=f(fmaxrow);
    
    
    %New sampling frequency
    Fs=50; 
    T=1/Fs; 

    ttrim=1/fmax;

    pmax=0;

    count=1;
    
    %Success probability and run time in truncated data is found
    while true
        t=(count-1)*T;

        if(t>ttrim)
            break;
        end

        exD=exp(-1i*D*t);
        upart = bsxfun(@times,V',exD);
        u1=V(1,:)*upart;

        psi1=u1*ones(N,1);

        psi1=psi1/sqrt(N);

        pmaxnew=abs(psi1)^2;

        if pmaxnew>pmax
            pmax=pmaxnew;
            tmax=t;
        end

        count=count+1;

    end
    
end

end