% Forms the global mass, stiffness, damping, external force and reaction
% force matrices
% M : Global assembled mass matrix
% K : Global assembled stiffness matrix
% C : Global assembled damping matrix
% F : Global assembled external force matrix
% R : Global assembled reaction force matrix

function [M,K,C,F,R,elemMat] = createGlobalMatrices(eConn,xConn,E_I,G_J,E_A,Ngp,fc,Le,Nelems,Ndof,d1,d2,mass,stiff)

%% global assembled matrices initialization
M = zeros(Ndof, Ndof);
K = zeros(Ndof, Ndof);
C = zeros(Ndof, Ndof);
F = zeros(Ndof,1);
R = zeros(Ndof,1);


%% Gauss quadrature data
% qs : abscissas
% ws : weights
[qs, ws] = quad_data(Ngp);

%%

elemMat.Me = {};
elemMat.Ke = {};
elemMat.Ce = {};
elemMat.Fe = {};

%% Forming elemental matrices and assembly
for e = 1:Nelems  %loop over all elements
    
    %%Typical element:
    %     *------------------*
    %     1,2,3,            7,8,9
    %     4,5,6             10,11,12
    
    
    n1 = eConn(e,1);        %left end number of element
    n2 = eConn(e,12);       %right end number of element
    x1 = xConn(e,1);        %left end coordinate of element
    x2 = xConn(e,12);       %right end coordinate of element
    f0 = zeros(6,1);        %internal force caused due to large deformation
    
    
    %%element mass matrices for K,M,R
    Ke_elem = zeros(12,12);
    Me_elem = zeros(12,12);
    Re_elem = zeros(12,1);
    Fe_elem = zeros(12,1);
    
    
    %%inital position of the beam from which f0 is to be calculated
    %%% no initial displacement, therefore f0 = N*ue_elem = 0;
    ue_elem = zeros(12,1);
    
    
    %%integrating local elemental matrices using gauss quadrature
    for q = 1:Ngp
        xi = qs(q);     %xi : value of the abscissa at which shape function
        %data is to be calculated in natural coordinate system
        
        [Nu,Nv,Nw,Ntx,Nty,Ntz,Bu,Bv,Bw,Btx,Bty,Btz,j] = sf_data(x1, x2, xi);
        
        N = [Nu,Nv,Nw,Ntx,Nty,Ntz];     %N - shape function matrix of order 12 X 6
        B = [Bu,Bv,Bw,Btx,Bty,Btz];     %B - derivative of shape function matrix 12 X 6
        
        %         Nu : Axial
        %         Nv : Transverse
        %         Nw : Lateral
        %         Ntx : N_theta_x(torsional)
        %         Nty : N_theta_y(lateral bending)
        %         Ntz : N_theta_z(transverse bending)
        %
        %         Similarly B are derivatives of corresponding shape functions
        
        %%calculating components of f0 based on initial well profile/beam position
        %%% no initial displacement, therefore f0 = 0;
        f0(1) = E_A*Bu'*ue_elem;        %EA*epsilon ==> strain: epsilon = (\partial u)/(\partial x); In place of u substitute N^T*ue
        f0(2) = 0;                      %similarly for other terms
        f0(3) = 0;
        f0(4) = G_J*Btx'*ue_elem;       %GJ*tau
        f0(5) = E_I*Bty'*ue_elem;       %EI*kappa_y
        f0(6) = E_I*Btz'*ue_elem;       %EI*kappa_z
        
        %%forming matrices
        Ke_elem = Ke_elem  + (B*stiff*B')*j*ws(q);   %equation of elemental stiffness : (B*stiff*B'); Multiplied by weights of gauss quad and jacobian
        Me_elem = Me_elem + (N*mass*N')*j*ws(q);     %eq of elemental mass matrix : (N*mass*N')
        Re_elem = Re_elem + B*f0*j*ws(q);                   %eq of elemental internal force : B*f0;
        
    end
    
    %%damping matrix
    Ce_elem = d1*Me_elem + d2*Ke_elem;
    
    %%elemental force matrix fe
    
    %     point force/moment is applied on the last node. (right node)
    %     we want to find the distribution of force/moment over the beam
    %     So, we calculate the value of force at each node.
    %
    %     elemental force, fe = \int_(0^le)(N*fext)dx
    %     where, fext = fc*diracdelta(x).
    %     Using the property of dirac delta function we can say:
    %     \int_(0^le){ (N*fext)dx } = \int_(0^le) { N*fc*diracdelta(x-le)dx } = N( @le )*fc
    %     {x-le becuase we need the force at le for each element}
    
    if e==Nelems
        value = Le; %right node of element
        [N_atValue] = sf_atValue(x1,x2,value); %gives the 12 X 6 values of shape function @value
        Fe_elem = N_atValue*fc;     %N( @le )*fc
    end
    
    
    
    elemMat.Me{e} = Me_elem;
    elemMat.Ke{e} = Ke_elem;
    elemMat.Ce{e} = Ce_elem;
    elemMat.Fe{e} = Fe_elem;
    
    
    %%Global element matrices
    
    %as Tr = I in initial case,
    Ke_g = Ke_elem;
    Me_g = Me_elem;
    Fe_g = Fe_elem;
    Ce_g = Ce_elem;
    Re_g = Re_elem;
    
    %%Assembly of matrices
    K((n1:n2), (n1:n2)) = K((n1:n2), (n1:n2)) + Ke_g;
    M((n1:n2), (n1:n2)) = M((n1:n2), (n1:n2)) + Me_g;
    F(n1:n2) = F(n1:n2) + Fe_g;
    R(n1:n2) = R(n1:n2) + Re_g;
    C((n1:n2), (n1:n2)) = C((n1:n2), (n1:n2)) + Ce_g;
end

%% Boundary conditions

% deleting 1st six rows and coloums as left end is fixed
M(1:6,:)=[];
M(:,1:6)=[];

K(1:6,:)=[];
K(:,1:6)=[];

C(1:6,:)=[];
C(:,1:6)=[];

F(1:6,:)=[];

R(1:6,:)=[];


