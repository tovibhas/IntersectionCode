%%%% AN 88 LINE TOPOLOGY OPTIMIZATION CODE Nov, 2010 %%%%
%%%% Intersection Indicator, Threshold Stress, Intersection Objective and Sensitivity
%%%  by Vibhas Mishra

clear all
close all

%% Design space and ring domain definition - Fixed
nelx=300;                                              % number of elements in x-direction
nely= 100;                                             % number of element in y-direction
volfrac= 0.5;                                          % volume fraction needed in the design space
penal= 3;                                              % penality for SIMP interpolation
rmin= 2.5;                                             % rmin radius for filtering density
ft= 2;                                                 % functionality ft = 1(No filtering), ft = 2(With filtering)
ElementNumber = reshape([1:nelx*nely],nely,nelx);      % the number of the elements
pI = 1;                                                % Penalization on the intersection objective
lambdaN = 0.02;                                        % lambda for normalizing objective function
%% MATERIAL PROPERTIES
E0 = 210000;                 % youngs modulus of material
Emin = E0*1e-9;              % youngs modulus of void
nu = 0.3;                    % poisons ration

%% Variables to be sweeped
ratio = 2.5/10;                 % ratio = rmin/rmin1,
rmin1 = rmin/ratio;             % rmin1 = radius of the circle for intersection size determination
vt = 1;                         % The penalization on weight which is fucntion of size, default = 1
wt = 0.2;                       % Weight on the intersection objective, ranges from 0.2 to 1.0.
IepsF = 1;                      % Factor to reduce the threshold stress, default = 1



%% PATH definition for saving files

pathtop = "./result/" ...
    + num2str(nelx) + "X" + num2str(nely) + "/"...
    + "vol=" + num2str(volfrac) + "/"...
    + "SIMP=" + num2str(penal) + "/"...
    + "rmin=" + num2str(rmin)  + "/"...
    + "ratio=" + num2str(ratio)  + "/"...
    + "vt=" + num2str(vt)  + "/"...
    + "wt=" + num2str(wt) + "/"...
    + "IepsF=" + num2str(IepsF) + "/";
mkdir(char(pathtop))

%% PREPARE FINITE ELEMENT ANALYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);                          % Element stiffness matrix definition
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);                                   % Node number definition
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);         % Element dof matrix definition
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
% DEFINE LOADS AND SUPPORTS (HALF MBB-BEAM)
F = sparse(2*(nely+1)*(nelx+1)-nely,1,-1000,2*(nely+1)*(nelx+1),1);
U = zeros(2*(nely+1)*(nelx+1),1);
fixeddofs = [1:2*(nely+1)];
alldofs = [1:2*(nely+1)*(nelx+1)];
freedofs = setdiff(alldofs,fixeddofs);


%% PREPARE FILTER - DENSITY FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);      % size of the filter based on rmin
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);           % weights have been assigned
Hs = sum(H,2);                  % sum of the weights

%% PREPARE FILTER - FOR WEIGHT CALCULATION DEPENDING ON INTERSECTION SIZE
iDI = ones(nelx*nely*(2*(ceil(rmin1)-1)+1)^2,1);
jDI = ones(size(iDI));
sDI = zeros(size(iDI));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin1)-1),1):min(i1+(ceil(rmin1)-1),nelx)
            for j2 = max(j1-(ceil(rmin1)-1),1):min(j1+(ceil(rmin1)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iDI(k) = e1;
                jDI(k) = e2;
                sDI(k) = max(0,rmin1-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
DI = sparse(iDI,jDI,sDI);                                    % weights have been assigned
DIs = max(sum(DI,2)).*ones(nelx*nely,1);                     % sum of the weights

%% WEIGHT CALCULATION - for normalization of multi-objective formulation

pathtopbaseline = "./baseline/" ...
    + num2str(nelx) + "X" + num2str(nely) + "/"...
    + "vol=" + num2str(volfrac) + "/"...
    + "SIMP=" + num2str(penal) + "/"...
    + "rmin=" + num2str(rmin)  + "/";

load(pathtopbaseline+'baseline.mat');

% Intersection size determination
xPhys1 = reshape(((DI*xPhys(:))./DIs),nely,nelx);

sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
K = sparse(iK,jK,sK); K = (K+K')/2;
U(freedofs) = K(freedofs,freedofs)\F(freedofs);
ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));

% B-matrix : displacement and strain relation
B = [-0.5,0,0.5,0,0.5,0,-0.5,0;...
    0,-0.5,0,-0.5,0,0.5,0,0.5;...
    -0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,-0.5];

% D-matrix : stress and strain relation
D = ((E0-Emin)/(1-nu^2)).*[1,nu,0;...
    nu,1,0;...
    0,0,0.5*(1-nu)];

% stress calculation
Istress = (xPhys(:).^penal).*(D*B*U(edofMat'))';

% principal stress calculation
Istress1 = ((0.5*(Istress(:,1)+Istress(:,2)))+((0.25*(Istress(:,1)-Istress(:,2)).^2) + (Istress(:,3).^2)).^0.5);
Istress2 = ((0.5*(Istress(:,1)+Istress(:,2)))-((0.25*(Istress(:,1)-Istress(:,2)).^2) + (Istress(:,3).^2)).^0.5);
Istressp = [Istress1,Istress2];

% Threshold stress calculation
IstressVM = (0.5.*((Istress1 - Istress2).^2 + Istress1.^2 + Istress2.^2)).^0.5;
IstressVMStr = sort(IstressVM(find(xPhys(:)>0.1 & xPhys(:)<0.5)),'descend');
Percent90VM = IstressVMStr(floor(0.1*length(IstressVMStr(:))));
Iepsilon = (Percent90VM/IepsF).^2;
Iepsilon1 = Iepsilon;

% Intersection Indicator calculation
Istresspabs = (Istressp.^2 + Iepsilon1);
Istresspratio = [Istresspabs(:,1)./Istresspabs(:,2),Istresspabs(:,2)./Istresspabs(:,1)];
T = 1./log2((Istresspratio(:,1)+Istresspratio(:,2)));
O = (xPhys(:).^penal).*T;
Oavg = DI*(O(:)./DIs);

% weight calculation based on size of intersection
w = (1-xPhys1.^vt).^(1/vt);
I = (sum((w(:).*O(:)).^pI)^(1/pI))/ (mean(xPhys(:)));

% Weights for normalization
wcn = c;
win = I;

%% INITIALIZE ITERATION
if wt < 1
    x = (volfrac).*ones(nely,nelx);
else
    x = load(pathtopbaseline+'baseline.mat').xPhys;
end
xPhys = reshape((H*x(:))./Hs,nely,nelx);
xPhys1 = reshape(((DI*xPhys(:))./DIs),nely,nelx);
xHeavi = x;
loop = 0;
change = 1;
m = 1;
n = nelx*nely;
xmin = (1e-6)*ones(n,1);
xmax = ones(n,1);
a0mma = 1;
amma = zeros(m,1);
cmma = 100000*ones(m,1);
dmma = zeros(m,1);
xold1(1:nely,1:nelx) = volfrac;
xold2(1:nely,1:nelx) = volfrac;
low = (1e-6)*ones(nelx*nely,1);
upp = ones(nelx*nely,1);
% resetting the compliance and intersection objective values
c =0;
I =0;

%% START ITERATION
while change > 0.01 && loop <1000
    
    loop = loop + 1;
    move = 0.2;
    
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K = sparse(iK,jK,sK); K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs);
    
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    ce = reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    c = sum(sum((Emin+xPhys.^penal*(E0-Emin)).*ce));
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    dv_max = ones(nely,nelx);
    dv_min = ones(nely,nelx);
    
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv_max(:) = H*(dv_max(:)./Hs);
        dv_min(:) = H*(dv_min(:)./Hs);
    end
    
    %% INTERSECTION OBJECTIVE CALCULATION
    % DEFINITON OF B and D MATRIX
    B = [-0.5,0,0.5,0,0.5,0,-0.5,0;...
        0,-0.5,0,-0.5,0,0.5,0,0.5;...
        -0.5,-0.5,-0.5,0.5,0.5,0.5,0.5,-0.5];
    
    D = ((E0-Emin)/(1-nu^2)).*[1,nu,0;...
        nu,1,0;...
        0,0,0.5*(1-nu)];
    
    % CALCULATION OF INTERSECTION OBJECTIVE
    Istress = (xPhys(:).^penal).*(D*B*U(edofMat'))';
    Istress1 = ((0.5*(Istress(:,1)+Istress(:,2)))+((0.25*(Istress(:,1)-Istress(:,2)).^2) + (Istress(:,3).^2)).^0.5);
    Istress2 = ((0.5*(Istress(:,1)+Istress(:,2)))-((0.25*(Istress(:,1)-Istress(:,2)).^2) + (Istress(:,3).^2)).^0.5);
    Istressp = [Istress1,Istress2];
    IstressVM = (0.5.*((Istress1 - Istress2).^2 + Istress1.^2 + Istress2.^2)).^0.5;
    Istresspabs = (Istressp.^2 + Iepsilon1);
    Istresspratio = [Istresspabs(:,1)./Istresspabs(:,2),Istresspabs(:,2)./Istresspabs(:,1)];
    T = 1./log2((Istresspratio(:,1)+Istresspratio(:,2)));
    O = ((xPhys(:).^penal).*T);
    w = (1-xPhys1.^vt).^(1/vt);
    I = (sum((w(:).*O(:)).^pI)^(1/pI))/ (mean(xPhys(:)));
    
    Istrain = (B*U(edofMat'))';
    Istrain1 = ((0.5*(Istrain(:,1)+Istrain(:,2)))+((0.25*(Istrain(:,1)-Istrain(:,2)).^2) + (Istrain(:,3).^2)).^0.5);
    Istrain2 = ((0.5*(Istrain(:,1)+Istrain(:,2)))-((0.25*(Istrain(:,1)-Istrain(:,2)).^2) + (Istrain(:,3).^2)).^0.5);
    
    % calculation of sensitivity of T with respect to principal stresses
    dTdsp1 = (2.*(Istress1./((Istress2.^2) + Iepsilon))) - ((2.*(Istress1).*((Istress2.^2) + Iepsilon))./(((Istress1.^2) + Iepsilon).^2));
    dTdsp2 = (2.*(Istress2./((Istress1.^2) + Iepsilon))) - ((2.*(Istress2).*((Istress1.^2) + Iepsilon))./(((Istress2.^2) + Iepsilon).^2));
    dTdsp = [dTdsp1,dTdsp2];
    
    % calculation of sensitivity of principal stress with respect to stress components
    dsp1ds1 = 0.5 + 0.25*((Istress(:,1)-Istress(:,2))./((0.25*(Istress(:,1)-Istress(:,2)).^2) + Istress(:,3).^2).^0.5);
    dsp1ds2 = 0.5 + 0.25*((Istress(:,2)-Istress(:,1))./((0.25*(Istress(:,1)-Istress(:,2)).^2) + Istress(:,3).^2).^0.5);
    dsp1ds12 = ((Istress(:,3))./((0.25*(Istress(:,1)-Istress(:,2)).^2) + Istress(:,3).^2).^0.5);
    dsp2ds1 = 0.5 - 0.25*((Istress(:,1)-Istress(:,2))./((0.25*(Istress(:,1)-Istress(:,2)).^2) + Istress(:,3).^2).^0.5);
    dsp2ds2 = 0.5 - 0.25*((Istress(:,2)-Istress(:,1))./((0.25*(Istress(:,1)-Istress(:,2)).^2) + Istress(:,3).^2).^0.5);
    dsp2ds12 = -((Istress(:,3))./((0.25*(Istress(:,1)-Istress(:,2)).^2) + Istress(:,3).^2).^0.5);
    
    
    
    dIt1 = penal.*(xPhys(:).^(penal-1)).*w(:).*(1./log2((Istresspratio(:,1)+Istresspratio(:,2))));
    dIt21 = -(log(2)).*(1./(Istresspratio(:,1)+Istresspratio(:,2))).*...
        ((1./log((Istresspratio(:,1)+Istresspratio(:,2)))).^2).*...
        (penal*(xPhys(:).^(2*penal-1)).*w(:)).*...
        sum(([(dTdsp1.*dsp1ds1 + dTdsp2.*dsp2ds1),(dTdsp1.*dsp1ds2 + dTdsp2.*dsp2ds2),(dTdsp1.*dsp1ds12 + dTdsp2.*dsp2ds12)]).*(D*B*U(edofMat'))',2);
    
    
    % Adjoint equation
    adj11 = B'*D'*[(dTdsp1.*dsp1ds1 + dTdsp2.*dsp2ds1),(dTdsp1.*dsp1ds2 + dTdsp2.*dsp2ds2),(dTdsp1.*dsp1ds12 + dTdsp2.*dsp2ds12)]';
    adj12 = (log(2)).*(1./(Istresspratio(:,1)+Istresspratio(:,2))).*((1./log((Istresspratio(:,1)+Istresspratio(:,2)))).^2).*w(:).*(xPhys(:).^(2*penal)).*adj11';
    adj13 = (sparse(edofMat(:),ones(length(edofMat(:)),1),adj12(:)));
    lambda22 = zeros(2*(nelx+1)*(nely+1),1);
    lambda22(freedofs) = K(freedofs,freedofs)\adj13(freedofs);
    
    dIt22t = (penal*(E0-Emin)*(xPhys.^(penal-1))).*(reshape(sum((lambda22(edofMat)*KE).*U(edofMat),2),nely,nelx));
    dIt22 = dIt22t(:);
    
    % Total sensitivity
    dI = (nelx*nely)*(((sum(xPhys(:))*((dIt1 + dIt21 + dIt22)...
        - ((DI*(((((1-xPhys1(:).^vt).^((1/vt)-1)).*(xPhys1(:).^(vt-1))).*O)./DIs)))))...
        -(sum((w(:).*O(:)).^pI)^(1/pI)) )/ ((sum(xPhys(:)))^2));
    dI = reshape(dI,nely,nelx);
    dI(:) = H*(dI(:)./Hs);
    
    %% MMA IMPLEMENTATION
    
    if wt < 2 && loop > 0
        f0val = (wt)*(c/(wcn*lambdaN)) + (1-wt)*(I/(win*lambdaN)) ;
        df0dx = ((wt)*(dc(:)./(wcn*lambdaN)) + (1-wt)*(dI(:)./(win*lambdaN)));
        df0dx2 = 0*df0dx;
    else
        f0val = c;
        df0dx = dc(:);
        df0dx2 = 0*df0dx;
    end
    
    fval01 = sum(xPhys(:))/(n*volfrac) -1;
    df01dx = dv_max(:)./(nelx*nely*volfrac);
    df01dx2 = 0*df01dx;
    
    fval = [fval01];
    dfdx = [df01dx'];
    dfdx2 = 0*dfdx;
    
    xval = x(:);
    
    xmax = min(1,x(:)+move);
    xmin = max(0,x(:)-move);
    
    [xmma,ymma,zmma,lam,xsi,eta,mu1,zet,s1,low,upp] = ...
        mmasub(m,n,loop,xval,xmin,xmax,xold1,xold2, ...
        f0val,df0dx,df0dx2,fval,dfdx,dfdx2,low,upp,a0mma,amma,cmma,dmma);
    
    
    xold2 = xold1;
    xold1 = xval;
    xnew = reshape(xmma,nely,nelx);
    
    if ft == 1
        xPhys = xnew;
    elseif ft == 2
        xPhys(:) = (H*xnew(:))./Hs;
        xPhys1(:) = ((DI*xPhys(:))./DIs);
        xPhys_old(:)= (H*xval(:))./Hs;
    end
    
    if loop < 0
        change = 1;
    else
        change = max(abs(xPhys(:)-xPhys_old(:)));
    end
    x = xnew;
    
    %% FIGURES
    
    figure(1)
    colormap(gray); imagesc(1-xPhys); caxis([0 1]); axis equal;axis tight manual; axis off; drawnow;daspect([1 1 1])
    
    figure(2)
    imagesc(reshape(T,nely,nelx)); axis equal,axis tight manual,axis off,colormap jet;caxis([0 1]);drawnow;daspect([1 1 1]); colorbar;
    
    figure(3)
    imagesc(reshape(w(:),nely,nelx)); axis equal,axis tight manual,axis off,colormap jet;caxis([0 1]);drawnow;daspect([1 1 1]); colorbar;
    
    figure(4)
    imagesc(reshape(w(:).*O(:),nely,nelx)); axis equal,axis tight manual,axis off,colormap jet;caxis([0 1]);drawnow;daspect([1 1 1]); colorbar;
    
    %% PRINT RESULTS
    fprintf(' It.:%5i ObjC.:%11.4f ObjI.:%11.4f Vol.:%7.3f ch.:%7.3f\n',loop,c,I, ...
        mean(xPhys(:)),change);
    %
    
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by E. Andreassen, A. Clausen, M. Schevenels,%
% B. S. Lazarov and O. Sigmund,  Department of Solid  Mechanics,           %
%  Technical University of Denmark,                                        %
%  DK-2800 Lyngby, Denmark.                                                %
% Please sent your comments to: sigmund@fam.dtu.dk                         %
%                                                                          %
% The code is intended for educational purposes and theoretical details    %
% are discussed in the paper                                               %
% "Efficient topology optimization in MATLAB using 88 lines of code,       %
% E. Andreassen, A. Clausen, M. Schevenels,                                %
% B. S. Lazarov and O. Sigmund, Struct Multidisc Optim, 2010               %
% This version is based on earlier 99-line code                            %
% by Ole Sigmund (2001), Structural and Multidisciplinary Optimization,    %
% Vol 21, pp. 120--127.                                                    %
%                                                                          %
% The code as well as a postscript version of the paper can be             %
% downloaded from the web-site: http://www.topopt.dtu.dk                   %
%                                                                          %
% Disclaimer:                                                              %
% The authors reserves all rights but do not guaranty that the code is     %
% free from errors. Furthermore, we shall not be liable in any event       %
% caused by the use of the program.                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

