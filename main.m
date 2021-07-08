%% Start
clear variables;
close all;

%% Number of Gauss points for gauss quadrature

% p+1 = 2Ngp
% Ngp : Number of gauss points
% p : degree of the polynomial upto which gauss quad gives correct values
Ngp = 3;

%% Geometry specifications
% specify the domain limits: x0 - initial coordinate, L - length of domain

Nnodes = 10;    %number of nodes
x0 = 0;
L  = 51.07; %material prop 1
% L = 1;    %material prop 2
% L = 10;   %material prop 3
% L =0.5;     %material prop 4
Le = L/(Nnodes-1); %Length of element

%% Create a mesh for the specified domain
% want to input number of nodes
% want to get nodal coordinates and element connectivity

Nelems = Nnodes - 1;    %number of elements
Ndof = 6*Nnodes;        %degrees of freedom (6 dof @each node: u,v,w,theta_x,theta_y,theta_z)
%  xs: Nnodesx1 vector containing nodal coordinates in sequence
%  eConn: (Nnodes-1) x 12 matrix containing element connectivity
%  xConn: Connects eConn to cordinates correpondingly.
[xs, eConn,xConn] = createMesh(x0, L, Nnodes);

%% Specification of material properties

%%% material properties 1 (From UT austin theses)
ro = 0.1715;
ri = 0.0714;
I = (pi/4)*(ro^4 - ri^4);
A = pi*(ro^2 - ri^2);
J = (pi/2)*(ro^4 - ri^4);
rho = 7830;
E = 206.8*10^9;
mu = 0.3;
G = E/(2*(1+mu));
rho_A = rho*A;
E_I = E*I;
G_J = G*J;
rho_J = rho*J;
E_A = E*A;


%%% material prop 2
% b = 0.01;
% h = 0.01;
% I = (b*h^3)/12;
% A = b*h;
% J = 2.25*b^4;
% rho = 7850;
% E = 215*10^9;
% mu = 0.3;
% G = E/(2*(1+mu));
% rho_A = rho*A;
% E_I = E*I;
% G_J = G*J;
% rho_J = rho*J;
% E_A = E*A;

% %% material prop 3 Pai
% b = 0.1;            %breadth of beam
% h = 0.01;           %height of beam
% E = 10^7;           %Young's modulus
% I = (b*h^3)/12;     %Second moment of area
% A = b*h;            %Area of cross section
% J = (1/12)*b*h*(h^2+b^2);       %polar moment of inertia
% rho = 7850;         %density of beam
% mu = 0.3;           %poisson's ratio
% G = E/(2*(1+mu));   %shear modulus from young's modulus
%
% rho_A = rho*A;
% E_I = E*I;
% G_J = G*J;
% rho_J = rho*J;
% E_A = E*A;


% %% material prop 4 Sanjeevi
% D = 0.016;
% I = (pi*(D^4))/64;
% A = (pi*(D^2))/4;
% J = (pi*(D^4))/32;
% rho = 1150;
% E = 3e9;
% mu = 41/100;
% G = E/(2*(1+mu));
% rho_A = rho*A;
% E_I = E*I;
% G_J = G*J;
% rho_J = rho*J;
% E_A = E*A;



%% External force

%external force : fext = fc*dirac(x)
%fext : Distribution of external force
% fc : centralized external force of order 6 X 1
% fc = [Fx  Fy  Fz  T_(theta_x)  M_(theta_y)  M_(theta_z)]
% dirac(x) : dirac delta function

n = 1;  %constant for various magnitudes of force/moment

fc = zeros(6,1);  %initializing centralized external force

% fc(2) = 50000;
fc(2) =(n*pi^2*E_I)/(4*L^2);        %applying transverse force Fy
% fc(6) = (n*2*pi*E_I)/L;           %applying moment M_(theta_z)

%Force ==> (n*pi^2*E_I)/(4*L^2);
%Moment ==> (n*2*pi*E_I)/L;
%this is sufficient moment/force to cause large deflection

%% Rayleigh damping constants
% d1*M + d2*K

d1 = 5; %2      %%mass proportional constant
d2 = 0.1; %0.1  %%stiffness proportional constant

%% formulate mass property and stiffness prop matrix
prop1 = [E_A,0,0,G_J,E_I,E_I];
prop2 = [rho_A,rho_A,rho_A,rho_J,0,0];
stiff = diag(prop1);
mass = diag(prop2);


%% Computing global matrices
[M,K,C,F,R,elemMat] = createGlobalMatrices(eConn,xConn,E_I,G_J,E_A,Ngp,fc,Le,Nelems,Ndof,d1,d2,mass,stiff);



%% Eigenvalues and eigenvectors
%
% %creating copies of matrices
% %mass matrices
% M_ax = M;   %axial
% M_tr = M;   %transverse
% M_lt = M;   %lateral
% M_to = M;   %torsional
%
% %stiffness matrices
% K_ax = K;
% K_tr = K;
% K_lt = K;
% K_to = K;
%
%
% %mass and stiffness matrix for finding axial natural frequencies
% for i = 2:Nnodes
%     M_ax(:,i:i+4) = [];
%     M_ax(i:i+4,:) = [];
%     K_ax(:,i:i+4) = [];
%     K_ax(i:i+4,:) = [];
% end
%
%
%
% %mass and stiffness matrix for finding transverse natural frequencies
% ui = 1;
% for i = 1:Nnodes
%     M_tr(:,ui) = [];
%     M_tr(ui,:) = [];
%     K_tr(:,ui) = [];
%     K_tr(ui,:) = [];
%     ui = ui + 5;
%     if ui>size(M_tr,2)
%         break;
%     end
% end
% ui = 2;
% for i = 1:Nnodes
%     M_tr(:,ui:ui+2) = [];
%     M_tr(ui:ui+2,:) = [];
%     K_tr(:,ui:ui+2) = [];
%     K_tr(ui:ui+2,:) = [];
%     ui = ui + 2;
%     if ui>size(M_tr,2)
%         break;
%     end
% end
%
%
% %mass and stiffness matrix for finding lateral natural frequencies
%
% ui = 1;
% for i = 1:Nnodes
%     M_lt(:,ui) = [];
%     M_lt(ui,:) = [];
%     K_lt(:,ui) = [];
%     K_lt(ui,:) = [];
%     ui = ui + 5;
%     if ui>size(M_lt,2)
%         break;
%     end
% end
% ui = 2;
% for i = 1:Nnodes
%     M_lt(:,ui:ui+2) = [];
%     M_lt(ui:ui+2,:) = [];
%     K_lt(:,ui:ui+2) = [];
%     K_lt(ui:ui+2,:) = [];
%     ui = ui + 2;
%     if ui>size(M_lt,2)
%         break;
%     end
% end
%
%
% %%mass and stiffness matrix for finding torsional natural frequencies
%
% M_to(:,1:3)=[]; M_to(1:3,:)=[];
% K_to(:,1:3)=[]; K_to(1:3,:)=[];
% for i = 2:Nnodes
%     if i==Nnodes
%         M_to(:,i:i+1) = [];
%         M_to(i:i+1,:) = [];
%         K_to(:,i:i+1) = [];
%         K_to(i:i+1,:) = [];
%         break;
%     end
%     M_to(:,i:i+4) = [];
%     M_to(i:i+4,:) = [];
%     K_to(:,i:i+4) = [];
%     K_to(i:i+4,:) = [];
%
% end
%
%
%
% a = linspace(0,L,Nnodes); %vector of discretized beam length
% a=a/L;
%
% %Eigenvalues for axial
% eval_ax = eig(M_ax\K_ax);
% ef_ax = sort(eval_ax);
% nf_ax = sort(sqrt(real(eval_ax)));      %Natural frequency axial
% [vec_ax,val_ax] = eig(M_ax\K_ax);
% ms1_ax = vec_ax(:,find(eval_ax == ef_ax(1),1));
% ms2_ax = vec_ax(:,find(eval_ax == ef_ax(2),1));
% ms3_ax = vec_ax(:,find(eval_ax == ef_ax(3),1));
% % ms4_ax = vec_ax(:,find(eval_ax == ef_ax(4),1));
% figure(1);
% plot(a,[0;ms1_ax/max(abs(ms1_ax))],'LineWidth',2);
% hold on
% plot(a,[0;ms2_ax/max(abs(ms2_ax))],'LineWidth',2);
% plot(a,[0;ms3_ax/max(abs(ms3_ax))],'LineWidth',2);
% % plot(a,[0;ms4_ax/max(abs(ms4_ax))],'LineWidth',2);
% title('Axial vibration modes');
% xlabel('Length')
% ylabel('Displacement');
% grid minor
% legend({'Mode 1','Mode 2','Mode 3'},'Location','best')
% saveas(gcf,'AxialvibrationModes.jpg')
%
%
%
%
% %Eigenvalues for transverse
% eval_tr = eig(M_tr\K_tr);
% ef_tr = sort(eval_tr);
% [vec_tr,val_tr] = eig(M_tr\K_tr);
% nf_tr = sort(sqrt(real(eval_tr)));
% ms1_tr = vec_tr(:,find(eval_tr == ef_tr(1),1));
% ms2_tr = vec_tr(:,find(eval_tr == ef_tr(2),1));
% ms3_tr = vec_tr(:,find(eval_tr == ef_tr(3),1));
% % ms4_tr = vec_tr(:,find(eval_tr == ef_tr(4),1));
% figure(2);
% plot(a,[0;ms1_tr(1:2:end)/max(abs(ms1_tr(1:2:end)))],'LineWidth',2);
% hold on
% plot(a,[0;ms2_tr(1:2:end)/max(abs(ms2_tr(1:2:end)))],'LineWidth',2);
% plot(a,[0;ms3_tr(1:2:end)/max(abs(ms3_tr(1:2:end)))],'LineWidth',2);
% % plot(a,[0;ms4_tr(1:2:end)/max(abs(ms4_tr(1:2:end)))],'LineWidth',2);
% title('Transverse vibration modes');
% xlabel('Length')
% ylabel('Displacement');
% grid minor
% legend({'Mode 1','Mode 2','Mode 3'},'Location','best')
% saveas(gcf,'TransverseVibrationModes.jpg')
%
%
%
% %Eigenvalues for lateral
% eval_lt = eig(M_lt\K_lt);
% ef_lt = sort(eval_lt);
% [vec_lt,val_lt] = eig(M_lt\K_lt);
% nf_lt = sort(sqrt(real(eval_lt)));
% ms1_lt = vec_lt(:,find(eval_lt == ef_lt(1),1));
% ms2_lt = vec_lt(:,find(eval_lt == ef_lt(2),1));
% ms3_lt = vec_lt(:,find(eval_lt == ef_lt(3),1));
% % ms4_lt = vec_lt(:,find(eval_lt == ef_lt(4),1));
% figure(3);
% plot(a,[0;ms1_lt(1:2:end)/max(abs(ms1_lt(1:2:end)))],'LineWidth',2);
% hold on
% plot(a,[0;ms2_lt(1:2:end)/max(abs(ms2_lt(1:2:end)))],'LineWidth',2);
% plot(a,[0;ms3_lt(1:2:end)/max(abs(ms3_lt(1:2:end)))],'LineWidth',2);
% % plot(a,[0;ms4_lt(1:2:end)/max(ms4_lt(1:2:end))],'LineWidth',2);
% title('Lateral vibration modes');
% xlabel('Length')
% ylabel('Displacement');
% grid minor
% legend({'Mode 1','Mode 2','Mode 3'},'Location','best')
% saveas(gcf,'LateralVibrationModes.jpg')
%
%
%
% %Eigenvalues for torsional
% eval_to = eig(M_to\K_to);
% ef_to = sort(eval_to);
% [vec_to,val_to] = eig(M_to\K_to);
% nf_to = sort(sqrt(real(eval_to)));
% ms1_to = vec_to(:,find(eval_to == ef_to(1),1));
% ms2_to = vec_to(:,find(eval_to == ef_to(2),1));
% ms3_to = vec_to(:,find(eval_to == ef_to(3),1));
% % ms4_to = vec_to(:,find(eval_to == ef_to(4),1));
% figure(4);
% plot(a,[0;ms1_to/max(abs(ms1_to))],'LineWidth',2);
% hold on
% plot(a,[0;ms2_to/max(abs(ms2_to))],'LineWidth',2);
% plot(a,[0;ms3_to/max(abs(ms3_to))],'LineWidth',2);
% % plot(a,[0;ms4_to/max(ms4_to)],'LineWidth',2);
% title('Torsional vibration modes');
% xlabel('Length')
% ylabel('Displacement');
% grid minor
% legend({'Mode 1','Mode 2','Mode 3'},'Location','best')
% saveas(gcf,'TorsionalVibrationModes.jpg')
%

%% Newmark method

% Initial conditions
dn = zeros(Ndof-6,1);       %displacement at time t     order: u;v;w;theta_x;theta_y;theta_z
vn = zeros(Ndof-6,1);       %velocity at time t
% an = zeros(Ndof-6,1);
Fnet = F - R;
an = M\Fnet;                  %accleration at time t



dt = 1e-3;     %timestep
time = 0:dt:40;  %simulation time


% Newmark method parameters beta=b,gamma=g
b = 0.25;
g = 0.5;

%integration constants for newmark method
I1 = 1/(b*dt^2);
I2 = 1/(b*dt);
I3 = (1/(2*b) -1);
I4 = g*dt*I1;
I5 = 1 - (g*dt*I2);
I6 = dt*(1 - g*I3 -g);

%effective stiffness matrix
K_eff = K + I1*M + I4*C;
K_inv = inv(K_eff);

%Triangularization of K_eff
% [LT,Diag] = ldl(K_eff);
% K_inv = inv(LT*Diag*LT'); %triangularization and inverting


displ = dn;   %displacements of all nodes for all time initialize
for dt = time
    dt
    [nu,dn1,vn1,an1] = Newmark(I1,I2,I3,I4,I5,I6,dn,vn,an,Fnet,M,K_inv,C);
    
    displ = [displ,dn1]; %appending displacements
    
    %     fprintf('dt = %0.6f\t end point displ = %0.2f\n',dt,displ(end-4,end));
    
    [M,K,C,F,R] = updateEqn(nu,Ngp,eConn,xConn,E_I,G_J,E_A,dn1,Nelems,Ndof,elemMat);
    
    %updating initial conditions for next timestep
    dn = dn1;
    vn = vn1;
    an = an1;
    Fnet = F - R;
    %effective stiffness matrix
    K_eff = K + I1*M + I4*C;
    K_inv = inv(K_eff);
    
    %Triangularization of K_eff
    %     [LT,Diag] = ldl(K_eff);
    %     K_inv = inv(LT*Diag*LT'); %triangularization and inverting
    %
end
ax_disp = displ(1:6:end,:); %axial disp of all nodes
tr_disp = displ(2:6:end,:); %transverse disp of all nodes
lat_disp = displ(3:6:end,:);  %lateral disp
theta_x = displ(4:6:end,:);
theta_y = displ(5:6:end,:);
theta_z = displ(6:6:end,:);


%% finding coordinates by keeping length of the element constant
%
Xt = [];
Yt = [];
for p = 1:Nelems %loop over each element
    X = [];
    Y = [];
    for k = 1:length(time) %loop over all time for pth elem
        yim1 = 0;
        xim1 = 0;
        for i = 1:p  %loop to cal x and y variation for pth elem
            yi = yim1 + Le*sin(theta_z(i,k));
            xi = xim1 + Le*cos(theta_z(i,k));
            yim1 = yi;
            xim1 = xi;
            %             yi = i.*Le*sin(theta_z(i,k));      <==   This gives different plot
            %             xi = i.*Le*cos(theta_z(i,k));     <==
        end
        X = [X;xi]; %variation in x for pth elem
        Y = [Y;yi]; %variation in y for pth elem
    end
    Xt = [Xt,X]; %appending variation in x for all elem
    Yt = [Yt,Y]; %appending variation in y for all elem
end

%normalizing
Xt = Xt/L;
Yt = Yt/L;


%% Animation

figure(1);
title('Beam displacement')
xlabel('X')
ylabel('Y')
grid on

for i = 1:100:size(Xt,1)
    for j = 1:p-1
        X1 = Xt(i,j);
        X2 = Xt(i,j+1);
        Y1 = Yt(i,j);
        Y2 = Yt(i,j+1);
        if j == 1
            plot([0,X1],[0,Y1],'-o','LineWidth',2);
            axis([-1 1 -1 1]);
        end
        hold on;
        plot([X1,X2],[Y1,Y2],'-o','LineWidth',2);
    end
    hold off;
    drawnow;
%     pause
    
end

