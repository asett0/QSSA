%Calculates gamma using sum approximation formula

%Output/s: 

%gamma - Optimal gamma

%Input/s: 

%lapl - Laplacian matrix of graph

function gamma = gcalc(lapl)

%Initialisation
gamma=0;
N=size(lapl,1);
    
%Eigenvalues are calculated and sorted   
eigvals=eig(-lapl);

eigvals = sort(eigvals);

%Gamma is calculated
for ii=2:N
    gamma=gamma+1/eigvals(ii);
end

gamma=gamma/N;
    
end
