% Calculates the displacements velocity and acceleration at time t+dt and
% returns them along with global slopes at t+dt

function [nu,dn1,vn1,an1] = Newmark(I1,I2,I3,I4,I5,I6,dn,vn,an,Fnet,M,K_inv,C)

%calculate effective load at t + dt
F_eff = Fnet + M*(I1*dn + I2*vn + I3*an) + C*(I4*dn - I5*vn - I6*an);       %%Fnet should be at t+dt

%%%%%%% Fnet on RHS is at next time step %%% CHECK %%%%


% solve for displacement at t+dt
dn1 = K_inv*F_eff;

% calculate accln and velocity at t+dt
vn1 = I4*(dn1-dn) + I5*vn + I6*an;
an1 = I1*(dn1 - dn) - I2*vn - I3*an;



%% extracting global slopes from the displacement dn1 @t+dt
nu = [];
n = 4;
for i = 1:length(dn1)
    nu = [nu;dn1(n:n+2)];
    n = n+6;
    if n > length(dn1)
        break;
    end
end

