% Calculates value of all the shape functions and its derivatives at a
% specified "value"
function [N_atValue] = sf_atValue(x1,x2,value)

le = x2-x1;

%% Axial u
Nu1 = 1 - (value/le); %u
Nu2 = 0;           %v
Nu3 = 0;            %w
Nu4 = 0;    %theta_x
Nu5 = 0;            %theta_y
Nu6 = 0;            %theta_z
Nu7 = (value/le);      %u_i+1
Nu8 = 0;
Nu9 = 0;
Nu10 = 0;
Nu11 = 0;
Nu12 = 0;
Nu = [Nu1,Nu2,Nu3,Nu4,Nu5,Nu6,Nu7,Nu8,Nu9,Nu10,Nu11,Nu12]';

%% transverse v
Nv1 = 0;
Nv2 = 1-3*(value/le)^2+2*(value/le)^3;
Nv3 = 0;
Nv4 = 0;
Nv5 = 0;
Nv6 = -value*(1-(value/le))^2;
Nv7 = 0;
Nv8 = 3*(value/le)^2-2*(value/le)^3;
Nv9 = 0;
Nv10 = 0;
Nv11 = 0;
Nv12 = -value*((value/le)^2-(value/le));
Nv = [Nv1,Nv2,Nv3,Nv4,Nv5,Nv6,Nv7,Nv8,Nv9,Nv10,Nv11,Nv12]';

%% lateral w
Nw1 = 0;
Nw2 = 0;
Nw3 = 1-3*(value/le)^2+2*(value/le)^3;
Nw4 = 0;
Nw5 = -value*(1-(value/le))^2;
Nw6 = 0;
Nw7 = 0;
Nw8 = 0;
Nw9 = 3*(value/le)^2-2*(value/le)^3;
Nw10 = 0;
Nw11 = -value*((value/le)^2-(value/le));
Nw12 = 0;
Nw = [Nw1,Nw2,Nw3,Nw4,Nw5,Nw6,Nw7,Nw8,Nw9,Nw10,Nw11,Nw12]';

%% torsional theta x
Ntx1 = 0;
Ntx2 = 0;
Ntx3 = 0;
Ntx4 = 1 - (value/le);
Ntx5 = 0;
Ntx6 = 0;
Ntx7 = 0;
Ntx8 = 0;
Ntx9 = 0;
Ntx10 = (value/le);
Ntx11 = 0;
Ntx12 = 0;
Ntx = [Ntx1,Ntx2,Ntx3,Ntx4,Ntx5,Ntx6,Ntx7,Ntx8,Ntx9,Ntx10,Ntx11,Ntx12]';

%% slope theta y
Nty1 = 0;
Nty2 = 0;
Nty3 = ((-6*value)/(le^2))*(1-(value/le));
Nty4 = 0;
Nty5 = -(1+3*(value/le)^2-4*(value/le));
Nty6 = 0;
Nty7 = 0;
Nty8 = 0;
Nty9 = ((6*value)/(le^2))*(1-(value/le));
Nty10 = 0;
Nty11 = -(value/le)*(3*(value/le)-2);
Nty12 = 0;
Nty = [Nty1,Nty2,Nty3,Nty4,Nty5,Nty6,Nty7,Nty8,Nty9,Nty10,Nty11,Nty12]';

%% slope theta z
Ntz1 = 0;
Ntz2 = ((-6*value)/(le^2))*(1-(value/le));
Ntz3 = 0;
Ntz4 = 0;
Ntz5 = 0;
Ntz6 = -(1+3*(value/le)^2-4*(value/le));
Ntz7 = 0;
Ntz8 = ((6*value)/(le^2))*(1-(value/le));
Ntz9 = 0;
Ntz10 = 0;
Ntz11 = 0;
Ntz12 = -(value/le)*(3*(value/le)-2);
Ntz = [Ntz1,Ntz2,Ntz3,Ntz4,Ntz5,Ntz6,Ntz7,Ntz8,Ntz9,Ntz10,Ntz11,Ntz12]';



%% Derivatives of shape function
%% Axial u
Bu1 = -(1/le);
Bu2 = 0;
Bu3 = 0;
Bu4 = 0;
Bu5 = 0;
Bu6 = 0;
Bu7 = (1/le);
Bu8 = 0;
Bu9 = 0;
Bu10 = 0;
Bu11 = 0;
Bu12 = 0;
Bu = [Bu1,Bu2,Bu3,Bu4,Bu5,Bu6,Bu7,Bu8,Bu9,Bu10,Bu11,Bu12]';

%% transverse v
Bv1 = 0;
Bv2 = ((-6*value)/(le^2))*(1-(value/le));
Bv3 = 0;
Bv4 = 0;
Bv5 = 0;
Bv6 = -(1+3*(value/le)^2-4*(value/le));
Bv7 = 0;
Bv8 = ((6*value)/(le^2))*(1-(value/le));
Bv9 = 0;
Bv10 = 0;
Bv11 = 0;
Bv12 = -(value/le)*(3*(value/le)-2);
Bv = [Bv1,Bv2,Bv3,Bv4,Bv5,Bv6,Bv7,Bv8,Bv9,Bv10,Bv11,Bv12]';

%% transverse w
Bw1 = 0;
Bw2 = 0;
Bw3 = ((-6*value)/(le^2))*(1-(value/le));
Bw4 = 0;
Bw5 = -(1+3*(value/le)^2-4*(value/le));
Bw6 = 0;
Bw7 = 0;
Bw8 = 0;
Bw9 = ((6*value)/(le^2))*(1-(value/le));
Bw10 = 0;
Bw11 = -(value/le)*(3*(value/le)-2);
Bw12 = 0;
Bw = [Bw1,Bw2,Bw3,Bw4,Bw5,Bw6,Bw7,Bw8,Bw9,Bw10,Bw11,Bw12]';

%% torsional theta x
Btx1 = 0;
Btx2 = 0;
Btx3 = 0;
Btx4 = -(1/le);
Btx5 = 0;
Btx6 = 0;
Btx7 = 0;
Btx8 = 0;
Btx9 = 0;
Btx10 = (1/le);
Btx11 = 0;
Btx12 = 0;
Btx = [Btx1,Btx2,Btx3,Btx4,Btx5,Btx6,Btx7,Btx8,Btx9,Btx10,Btx11,Btx12]';

%% slope theta y
Bty1 = 0;
Bty2 = 0;
Bty3 = -(6/le^2)*(1-2*(value/le));
Bty4 = 0;
Bty5 = -(2/le)*(3*(value/le)-2);
Bty6 = 0;
Bty7 = 0;
Bty8 = 0;
Bty9 = (6/le^2)*(1-2*(value/le));
Bty10 = 0;
Bty11 = (2/le)*(1-3*(value/le));
Bty12 = 0;
Bty = [Bty1,Bty2,Bty3,Bty4,Bty5,Bty6,Bty7,Bty8,Bty9,Bty10,Bty11,Bty12]';

%% slope theta z
Btz1 = 0;
Btz2 = -(6/le^2)*(1-2*(value/le));
Btz3 = 0;
Btz4 = 0;
Btz5 = 0;
Btz6 = -(2/le)*(3*(value/le)-2);
Btz7 = 0;
Btz8 = (6/le^2)*(1-2*(value/le));
Btz9 = 0;
Btz10 = 0;
Btz11 = 0;
Btz12 = (2/le)*(1-3*(value/le));
Btz = [Btz1,Btz2,Btz3,Btz4,Btz5,Btz6,Btz7,Btz8,Btz9,Btz10,Btz11,Btz12]';
%


N_atValue = [Nu,Nv,Nw,Ntx,Nty,Ntz];
