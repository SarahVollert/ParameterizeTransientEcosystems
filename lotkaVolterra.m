function dndt = lotkaVolterra(~,n,A,r)
%This is the Lokta-Volterra equations.

dndt = n.*r +( A*n).*n;

end