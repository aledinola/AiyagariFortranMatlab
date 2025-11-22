function F = f_returnfn(aprime,a,z,alpha,delta,r,tau_a,crra)
% Action space: (aprime,a,z)
% State variables: (a,z)

w=(1-alpha)*((r+delta)/alpha)^(alpha/(alpha-1));

income = w*z+r*a;
taxes = tau_a*a;

cons  = a + income - taxes - aprime; 

F=-Inf;
if cons>0
    if crra==1
        F = log(cons);
    else
        F = cons^(1-crra)/(1-crra);
    end
end

end %end function