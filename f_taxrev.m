function taxes = f_taxrev(aprime,a,z,alpha,delta,r,tau_a)
% f_taxrev computes total taxes paid by household with state (a,z)
% Action space: (aprime,a,z)
% State variables: (a,z)

w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));

income = w*z+r*a;
taxes  = tau_a*a;

end %end function