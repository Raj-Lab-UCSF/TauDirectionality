function Xt = spreadModel(Xo,t,L,b)
%This is a function to model the progression of pathology, Xt, via the
%graph diffusion equation X(t) = e^(-beta*L*t) * X(0)
%here the graph diffusion also models accumulation, so the equation
%becomes: X(t) = e^(-beta*L*t) * X(i) + X(i) for each iteration of the
%model, i.  

%Xt = the outcome of the model, a matrix of the outcome of all iterations
%Xo = the initial seeding of the model, a vector of all nodes with val = 0
%or 1
%t = a vector of times, if evenly spaced all = same #, if not, make jumps
%proportionate
%L = Normalized Laplacian (Adjacency) Matrix
%b = ease of diffusion constant, joins with t to make the bt time diffusion
%constant

%L = genLplcns(mat);

%ease of diffusion constant
beta = b;

%Xi = X outcome vector at each iteration, Xt = matrix of outcome vectors
Xi = Xo;
Xt = zeros(size(Xo,1),length(t));

%number of iterations of model = length of vector of times
for i = 1:length(t)
%     Xi = expm(-beta*L*t(i)) * Xo + Xi;
%     Xi = expm(-beta*L*t(i)) * Xo;
    Xi = expm(-beta*L*t(i)) * Xi + Xi;
    Xt(:,i) = Xi;
end

end



