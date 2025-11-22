function cons = f_consumption(aprime,a,z,alpha,delta,r,tau_a)
% Action space: (aprime,a,z)
% State variables: (a,z)

w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));

income = w*z+r*a;
taxes  = tau_a*a;

cons  = a + income - taxes - aprime; 

end %end function