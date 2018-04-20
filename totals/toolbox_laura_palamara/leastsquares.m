function y = leastsquares(tm, om, nconst, u, npts)
% Calculates the solution of a 
% least squares system of equations.

dim = 2*nconst+1;
A = zeros(dim, dim);
temp = zeros(1, dim);
c = zeros(1, dim);
for i = 1:1:npts;
   temp(1) = 1.0;
   j = 1:1:nconst;
   temp(2*j)= cos(om(j)*tm(i));
   temp(2*j+1) = sin(om(j)*tm(i));
   for k = 1:1:dim;
      for l = 1:1:dim;
         A(k,l) = temp(k)*temp(l)+A(k,l);
      end
   end
   c(1) = u(i);
   m = 1:1:nconst;
   c(2*m) = u(i)*temp(2*m) + c(2*m);
   c(2*m+1) = u(i)*temp(2*m+1) + c(2*m+1);
end
for i = 1:1:dim;
    j = 1:1:dim;
    A(i, j) = A(i, j)./npts;
    c(i) = c(i)./npts;
end   
y = A\c';
