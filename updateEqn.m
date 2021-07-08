% Forms the elemental mass and stiffness matrix for each element and
% transforms them using Rotation matrix which is calculated using the
% nu(slope) values
% Also forms F and R and transforms them.
function [M,K,C,F,R] = updateEqn(nu,Ngp,eConn,xConn,E_I,G_J,E_A,dn1,Nelems,Ndof,elemMat)
nu = [zeros(3,1);nu];   %adding the initial zeros to nu which were deleted by boundary conditions for simplicity of code further
% dn1 = [zeros(6,1);dn1]; %adding the initial zeros to dn1 which were deleted to complete ue_elem (12 X 1)


%% global assembled matrices initialize
K = zeros(Ndof, Ndof);
M = zeros(Ndof, Ndof);
C = zeros(Ndof, Ndof);
F = zeros(Ndof, 1);
R = zeros(Ndof, 1);



%% Connectivity matrix for Rot
%to transform the particular node
nuConn = zeros(Nelems,2);
index = 1;
for i = 1:Nelems
    nuConn(i,:) = [index,index+3];
    index = index+3;
end

%% Gauss quadrature data
% qs :  abscissas
% ws : weights
[qs, ws] = quad_data(Ngp);

%% Formulate element matrices
for e = 1:Nelems
    n1 = eConn(e,1);        %left end number of element
    n2 = eConn(e,12);       %right end number of element
    x1 = xConn(e,1);        %left end coordinate of element
    x2 = xConn(e,12);       %right end coordinate of element
    f0 = zeros(6,1);        %internal force caused due to large deformation
    
    %element matrices initialize
    Re_elem = zeros(12,1);
    
    %%position of the beam at t+dt from which f0 is to be calculated
    ue_elem = zeros(12,1);
    %     ue_elem = dn1(n1:n2);
    
    %integrating local elemental stiffness and mass matrices
    for q = 1:Ngp
        xi = qs(q);     %xi : value of the abscissa at which shape function
        %data is to be calculated in natural coordinate system
        [~,~,~,~,~,~,Bu,Bv,Bw,Btx,Bty,Btz,j] = sf_data(x1, x2, xi);
        
        B = [Bu,Bv,Bw,Btx,Bty,Btz];     %B - derivative of shape function matrix 12 X 6
        
        %%calculating components of f0 based on initial well profile/beam position
        
        f0(1) = E_A*Bu'*ue_elem;        %EA*epsilon ==> strain: epsilon = (\partial u)/(\partial x); In place of u substitute N^T*ue
        f0(2) = 0;                      %similarly for other terms
        f0(3) = 0;
        f0(4) = G_J*Btx'*ue_elem;       %GJ*tau
        f0(5) = E_I*Bty'*ue_elem;       %EI*kappa_y
        f0(6) = E_I*Btz'*ue_elem;       %EI*kappa_z
        
        %%forming matrices
        Re_elem = Re_elem + B*f0*j*ws(q);                   %eq of elemental internal force : B*f0;
        
    end
    
    %% calculating Tr for each node and Rot for each element
    nuC1 = nuConn(e,1); %values to help calculate at respective nodes
    nuC2 = nuConn(e,2);
    
    kx=0; ky=0; kz=1;
    theta = [0,-kz,ky;
        kz,0,-kx;
        -ky,kx,0];
    
    for j = [nuC1,nuC2]
        nuvec = [nu(j);nu(j+1);nu(j+2)]; %nux = nu(j); slope in x
        %nuy = nu(j+1); slope in y
        %nuz=nu(j+2);  slope in z
        
        
        if j == nuC2 %%Tr for 2nd node (right)
            
            %rodrigues formula
            Tr2 = eye(3)+sin(nuvec(3))*theta+(1-cos(nuvec(3)))*theta^2;
            
            if  isnan(Tr2)      %if Tr2 is giving NaN as nuvec is zero; Set Tr1 = Identity and break
                Tr2 = eye(3);
                break;
            end
            break;
        end
        
        %Tr for 1st node (left)
        %rodrigues formula
        Tr1 = eye(3)+sin(nuvec(3))*theta+(1-cos(nuvec(3)))*theta^2;
        %         Tr1 = eye(3) + (sin(norm(nuvec))/norm(nuvec))*theta...
        %             + ((1-cos(norm(nuvec)))/norm(nuvec))*theta^2;
        
        if  isnan(Tr1)  %if Tr1 is giving NaN as nuvec is zero; Set Tr1 = Identity
            Tr1 = eye(3);
        end
        
    end
    
    %rotation matrix
    Rot = [Tr1,zeros(3),zeros(3),zeros(3);
        zeros(3),Tr1,zeros(3),zeros(3);
        zeros(3),zeros(3),Tr2,zeros(3);
        zeros(3),zeros(3),zeros(3),Tr2];
    
    Ke_elem = elemMat.Ke{1,e};
    Me_elem = elemMat.Me{1,e};
    Ce_elem = elemMat.Ce{1,e};
    Fe_elem = elemMat.Fe{1,e};
    
    %% Global element matrices
    Ke_g = Rot'*Ke_elem*Rot;
    Me_g = Rot'*Me_elem*Rot;
    Ce_g = Rot'*Ce_elem*Rot;
    %     Fe_g = Rot'*Fe_elem;
    %     Re_g = Rot'*Re_elem;
    Fe_g = Fe_elem;
    Re_g = Re_elem;
    
    %% Assembly of matrices
    K((n1:n2), (n1:n2)) = K((n1:n2), (n1:n2)) + Ke_g;
    M((n1:n2), (n1:n2)) = M((n1:n2), (n1:n2)) + Me_g;
    C((n1:n2), (n1:n2)) = C((n1:n2), (n1:n2)) + Ce_g;
    F(n1:n2) = F(n1:n2) + Fe_g;
    R(n1:n2) = R(n1:n2) + Re_g;
    
end

%% Boundary conditions

% deleting 1st six rows and coloums as left end is fixed
M(1:6,:)=[];
M(:,1:6)=[];
% M(1,1)=1;M(2,2)=1;M(3,3)=1;M(4,4)=1;M(5,5)=1;M(6,6)=1;

K(1:6,:)=[];
K(:,1:6)=[];

C(1:6,:)=[];
C(:,1:6)=[];

F(1:6,:)=[];
R(1:6,:)=[];