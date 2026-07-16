function ops = stevens_ops(J)
%STEVENS_OPS Angular momentum and Stevens operator matrices for given J.
% Basis ordered mJ = J..-J. Conventions: Hutchings; O64s carries the 1/(4i) factor.
dim = round(2*J+1);
mJ  = (J:-1:-J).';
X   = J*(J+1);
Jz  = diag(mJ);
jp  = sqrt(X - mJ(2:end).*(mJ(2:end)+1));  % <m+1|J+|m>, m = mJ(2:end)
Jplus  = diag(jp, 1);
Jminus = Jplus';
E = eye(dim);
ops.Jz = Jz; ops.Jplus = Jplus; ops.Jminus = Jminus;
ops.Jx = (Jplus + Jminus)/2;
ops.Jy = (Jplus - Jminus)/(2i);
ops.O20 = 3*Jz^2 - X*E;
ops.O40 = 35*Jz^4 - (30*X - 25)*Jz^2 + (3*X^2 - 6*X)*E;
ops.O44 = 0.5*(Jplus^4 + Jminus^4);
ops.O44s = -0.5i*(Jplus^4 - Jminus^4);   % = (Jp^4-Jm^4)/(2i), sine partner of O44
ops.O60 = 231*Jz^6 - (315*X - 735)*Jz^4 + (105*X^2 - 525*X + 294)*Jz^2 ...
          + (-5*X^3 + 40*X^2 - 60*X)*E;
Cm = 11*Jz^2 - (X + 38)*E;
P4 = Jplus^4 + Jminus^4;  M4 = Jplus^4 - Jminus^4;
ops.O64c = 0.25 *(P4*Cm + Cm*P4);
ops.O64s = -0.25i*(M4*Cm + Cm*M4);
end
