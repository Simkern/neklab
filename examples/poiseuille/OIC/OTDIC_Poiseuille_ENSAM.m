% Optimal Disturbance
%
% compute the Orr-Sommerfeld matrix for three-
% dimensional Poiseuille or Couette flows and 
% compute the energy weight matrix 
%
% compute and displays the optimal initial condition and
% the corresponding flow response: inpout/output
%
% INPUT 
%
% Re        = Reynolds number
% alpha     = alpha (streamwise wave number)
% beta      = beta  (spanwise wave number)
% iflow     = type of flow  (Poiseuille=1, Couette=2)  
% N         = total number of modes for normal velocity
% T         = time of optimal growth
%
    clear
    close all
    clc
    addpath('/scratch/skern/matlab/TutorialCodes')
    addpath('/scratch/skern/git/matlab/nekmatlab/')
    addpath('/scratch/skern/nek/v17/02_OTD/_OptimalIC_Poiseuille/')
    set(groot, 'defaultAxesTickLabelInterpreter','latex');
    set(groot, 'defaultLegendInterpreter','latex');
    set(groot, 'defaultTextInterpreter','latex');

    global D0 D1 D2 D4 
    global qb
    
    zi = sqrt(-1);
    %...input data
    iflow  = 1;          % Poiseuille (1) or Couette flow (2)
    N      = 96;        % number of Chebyshev polynomials
    Re     = 5000;       % Reynolds number
    alpha0 = 1;          % alpha
    beta0  = 1;          % beta
    Tmax   = 20.26605; %25.14365;          % Time of maximum energy
    T      = [0,2*Tmax]; % T range for amplification envelope 
    if3d   = 1;
    outnek = 1;
    r      = 1;          % number of modes to retain;
    xdom   = 2*pi;	 % domain extent in x direction (alpha)
    if if3d==1
    	zdom = 2*pi;	 % domain extent in z direction (beta)
    else
	zdom = 0;
    end

    % 3D
    meshfname2d='poiseuille2D.fld'
    % 2D
    meshfname3d='poiseuille3D.fld'
    
    fprintf('\nOptimal initial condition for Poiseuille flow:\n\n')
    fprintf('\t# cheb  = %i\n',N)
    fprintf('\tRe      = %i\n',Re)
    fprintf('\talpha   = %f\n',alpha0)
    fprintf('\tbeta    = %f\n',beta0)
    fprintf('\tTmax    = %f\n',Tmax)
    fprintf('\t# modes = %i\n',r)
    fprintf('\nMesh to be read: %s\n',meshfname2d)
    fprintf('\nMesh to be read: %s\n',meshfname3d)
    fprintf('\noutnek : %i',outnek)
    %fprintf('\nsavef  : %i\n',savef)
    %fprintf('\nDirectory suffix: %s\n\n',dirsuffix)

    %...generate Chebyshev differentiation matrices
    [D0,D1,D2,D4] = ChebMat(N);

    % choose all possible alpha,beta combinations
    alpha = alpha0; beta = beta0;
    %...set up Orr-Sommerfeld matrices A and B 
    if (iflow == 1)
        [A,B] = PoiseuilleMatrix(N,alpha,beta,Re);
    else
        [A,B] = CouetteMatrix(N,alpha,beta,Re);
    end

    %...generate energy weight matrix
    k2 = alpha^2 + beta^2;
    M  = EnergyMatrix(N+1,N+1,k2);

    %...compute the Orr-Sommerfeld matrix (by inverting B)
    OS = inv(B)*A;

%    %...compute the numerical range  
    EM    = GetMatrix(OS,M,k2);
    e     = eig(EM);
    outeigs2=sort(unique(round(imag(e),8)),'descend');
    
     %...compute the growth rate envelope
    fprintf('\nCompute energy growth envelope for t = [ %4.2f, %4.2f ]\n',T(1),T(end))
    [~,~,gg] = OptimalIC(OS,T,M,k2,1,1);
    TGenv = [gg(:,1)';gg(:,2)'];
   %[~,idx] = max(gg(:,2));
   %Tmax = gg(idx,1);
    GTmax = interp1(gg(:,1),gg(:,2),Tmax);
    
    fprintf('\nCompute Optimal IC and OC for T = %f\n',Tmax)
    [flowin,flowot,gg] = OptimalIC(OS,Tmax,M,k2,2,r);
    
    %...visualize the optimal perturbations
    vin    = D0*flowin(1:N+1,:);
    dvin   = D1*flowin(1:N+1,:);
    etain  = D0*flowin(N+2:2*(N+1),:);
    vout   = D0*flowot(1:N+1,:); 
    dvout  = D1*flowot(1:N+1,:);
    etaout = D0*flowot(N+2:2*(N+1),:); 
    ycoord = D0(:,2);
    
    uinvar = 1i*alpha/k2*dvin-1i*beta/k2*etain;
    duinvar = 1i*alpha*uinvar;
    vinvar = vin;
    winvar = 1i*beta/k2*dvin+1i*alpha/k2*etain;
    uoutvar = 1i*alpha/k2*dvout-1i*beta/k2*etaout;
    voutvar = vout;
    woutvar = 1i*beta/k2*dvout+1i*alpha/k2*etaout;
    
    xcoord = 0:pi/10:40*pi;
    zcoord = 0:pi/20:2*pi;
    
    %% read mesh
    disp(' ') 
    disp(['Reading mesh file: ' meshfname2d]);
    [data2d,lr12d,elmap2d,~,~,fields2d,emode,wdsz,etag,header2d,status] = readnek(meshfname2d);
    disp(['Reading mesh file: ' meshfname3d]);
    [data3d,lr13d,elmap3d,~,~,fields3d,emode,wdsz,etag,header3d,status] = readnek(meshfname3d);
    if status ~= 0
        disp('Error in readnek.m. Aborting.')
        return
    end
    disp('Done.')
  
    %% to choose the right fields
    if if3d == 1
        ndim = 3;
    else
        ndim = 2;
    end
    
    nel2d = numel(data2d(:,:,1));
    nel3d = numel(data3d(:,:,1));
    xx2d = reshape(data2d(:,:,1),nel2d,1);
    yy2d = reshape(data2d(:,:,2),nel2d,1);

    xx3d = reshape(data3d(:,:,1),nel3d,1);
    yy3d = reshape(data3d(:,:,2),nel3d,1);
    zz3d = reshape(data3d(:,:,3),nel3d,1);
    
    utmp = interp1(ycoord,vinvar(:,1),yy2d);
    utmp2 = chebyshev_interp_1d(length(ycoord),ycoord,vinvar(:,1),length(yy2d),yy2d);

    save('OIC_mat.mat','uinvar','vinvar','winvar','uoutvar','voutvar','woutvar','alpha','beta','ycoord')
%% OPTIMAL INITIAL CONDITION
    time = 0;
    istep = 0;
    disp(['Optimal initial condition for maximum growth at T= ' num2str(Tmax) 's:' ])
    for nic=1:1
        utmp2 = interp1(ycoord,uinvar(:,nic),yy2d);
        vtmp2 = interp1(ycoord,vinvar(:,nic),yy2d);
        wtmp2 = interp1(ycoord,winvar(:,nic),yy2d);
        utmp3 = interp1(ycoord,uinvar(:,nic),yy3d);
        vtmp3 = interp1(ycoord,vinvar(:,nic),yy3d);
        wtmp3 = interp1(ycoord,winvar(:,nic),yy3d);
            
        utmp2dr = real(utmp2.*exp(-1i*xx2d*alpha));
        vtmp2dr = real(vtmp2.*exp(-1i*xx2d*alpha));
        wtmp2dr = real(wtmp2.*exp(-1i*xx2d*alpha));
        utmp2di = imag(utmp2.*exp(-1i*xx2d*alpha));
        vtmp2di = imag(vtmp2.*exp(-1i*xx2d*alpha));
        wtmp2di = imag(wtmp2.*exp(-1i*xx2d*alpha));
            
        ndim = 2
        data2d(:,:,ndim+1) = reshape(utmp2dr,size(data2d(:,:,1)))/sqrt(2);
        data2d(:,:,ndim+2) = reshape(vtmp2dr,size(data2d(:,:,1)))/sqrt(2);
        
        fname2D=[ 'OptimalIC2D_' num2str(nic,'%02i') '_uv_r.fld'];
        status = writenek(fname2D,data2d,lr12d,elmap2d,time,istep,fields2d,emode,wdsz,etag);
        fprintf('  2D Optimal IC field (r=%i) written: %s\n',nic,fname2D)
        
        data2d(:,:,ndim+1) = reshape(utmp2dr,size(data2d(:,:,1)))/sqrt(2);
        data2d(:,:,ndim+2) = reshape(wtmp2dr,size(data2d(:,:,1)))/sqrt(2);
        fname2D=[ 'OptimalIC2D_' num2str(nic,'%02i') '_uw_r.fld'];
        status = writenek(fname2D,data2d,lr12d,elmap2d,time,istep,fields2d,emode,wdsz,etag);
        fprintf('  2D Optimal IC field (r=%i) written: %s\n',nic,fname2D)
            
        data2d(:,:,ndim+1) = reshape(utmp2di,size(data2d(:,:,1)))/sqrt(2);
        data2d(:,:,ndim+2) = reshape(vtmp2di,size(data2d(:,:,1)))/sqrt(2);
        
        fname2D=[ 'OptimalIC2D_' num2str(nic,'%02i') '_uv_i.fld'];
        status = writenek(fname2D,data2d,lr12d,elmap2d,time,istep,fields2d,emode,wdsz,etag);
        fprintf('  2D Optimal IC field (r=%i) written: %s\n',nic,fname2D)
        
        data2d(:,:,ndim+1) = reshape(utmp2di,size(data2d(:,:,1)))/sqrt(2);
        data2d(:,:,ndim+2) = reshape(wtmp2di,size(data2d(:,:,1)))/sqrt(2);
        fname2D=[ 'OptimalIC2D_' num2str(nic,'%02i') '_uw_i.fld'];
        status = writenek(fname2D,data2d,lr12d,elmap2d,time,istep,fields2d,emode,wdsz,etag);
        fprintf('  2D Optimal IC field (r=%i) written: %s\n',nic,fname2D)
        
        utmp3d = real(utmp3.*exp(-1i*(xx3d*alpha + zz3d*beta)));
        vtmp3d = real(vtmp3.*exp(-1i*(xx3d*alpha + zz3d*beta)));
        wtmp3d = real(wtmp3.*exp(-1i*(xx3d*alpha + zz3d*beta)));
        
        ndim = 3
        data3d(:,:,ndim+1) = reshape(utmp3d,size(data3d(:,:,1)));
        data3d(:,:,ndim+2) = reshape(vtmp3d,size(data3d(:,:,1)));
        data3d(:,:,ndim+3) = reshape(wtmp3d,size(data3d(:,:,1)));
             
        fname3D=[ 'OptimalIC3D_' num2str(nic,'%02i') '.fld'];
        status = writenek(fname3D,data3d,lr13d,elmap3d,time,istep,fields3d,emode,wdsz,etag);
        fprintf('  3D Optimal IC field (r=%i) written: %s\n',nic,fname3D)
    end
    disp('Done.')
%% OPTIMAL OUTPUT
    disp(['Optimal output at T= ' num2str(Tmax) 's:' ])
    for nic=1:r
        utmp2 = interp1(ycoord,uoutvar(:,nic),yy2d);
        vtmp2 = interp1(ycoord,voutvar(:,nic),yy2d);
        wtmp2 = interp1(ycoord,woutvar(:,nic),yy2d);
        utmp3 = interp1(ycoord,uoutvar(:,nic),yy3d);
        vtmp3 = interp1(ycoord,voutvar(:,nic),yy3d);
        wtmp3 = interp1(ycoord,woutvar(:,nic),yy3d);
            
        utmp2dr = real(utmp2.*exp(-1i*xx2d*alpha));
        vtmp2dr = real(vtmp2.*exp(-1i*xx2d*alpha));
        wtmp2dr = real(wtmp2.*exp(-1i*xx2d*alpha));
        utmp2di = imag(utmp2.*exp(-1i*xx2d*alpha));
        vtmp2di = imag(vtmp2.*exp(-1i*xx2d*alpha));
        wtmp2di = imag(wtmp2.*exp(-1i*xx2d*alpha));
            
        ndim = 2
        data2d(:,:,ndim+1) = reshape(utmp2dr,size(data2d(:,:,1)))/sqrt(2);
        data2d(:,:,ndim+2) = reshape(vtmp2dr,size(data2d(:,:,1)))/sqrt(2);
        
        fname2D=[ 'OptimalOC2D_' num2str(nic,'%02i') '_uv_r.fld'];
        status = writenek(fname2D,data2d,lr12d,elmap2d,time,istep,fields2d,emode,wdsz,etag);
        fprintf('  2D Optimal OC field (r=%i) written: %s\n',nic,fname2D)
        
        data2d(:,:,ndim+1) = reshape(utmp2dr,size(data2d(:,:,1)))/sqrt(2);
        data2d(:,:,ndim+2) = reshape(wtmp2dr,size(data2d(:,:,1)))/sqrt(2);
        fname2D=[ 'OptimalOC2D_' num2str(nic,'%02i') '_uw_r.fld'];
        status = writenek(fname2D,data2d,lr12d,elmap2d,time,istep,fields2d,emode,wdsz,etag);
        fprintf('  2D Optimal OC field (r=%i) written: %s\n',nic,fname2D)
            
        data2d(:,:,ndim+1) = reshape(utmp2di,size(data2d(:,:,1)))/sqrt(2);
        data2d(:,:,ndim+2) = reshape(vtmp2di,size(data2d(:,:,1)))/sqrt(2);
        
        fname2D=[ 'OptimalOC2D_' num2str(nic,'%02i') '_uv_i.fld'];
        status = writenek(fname2D,data2d,lr12d,elmap2d,time,istep,fields2d,emode,wdsz,etag);
        fprintf('  2D Optimal IC field (r=%i) written: %s\n',nic,fname2D)
        
        data2d(:,:,ndim+1) = reshape(utmp2di,size(data2d(:,:,1)))/sqrt(2);
        data2d(:,:,ndim+2) = reshape(wtmp2di,size(data2d(:,:,1)))/sqrt(2);
        fname2D=[ 'OptimalOC2D_' num2str(nic,'%02i') '_uw_i.fld'];
        status = writenek(fname2D,data2d,lr12d,elmap2d,time,istep,fields2d,emode,wdsz,etag);
        fprintf('  2D Optimal OC field (r=%i) written: %s\n',nic,fname2D)
        
        utmp3d = real(utmp3.*exp(-1i*(xx3d*alpha + zz3d*beta)));
        vtmp3d = real(vtmp3.*exp(-1i*(xx3d*alpha + zz3d*beta)));
        wtmp3d = real(wtmp3.*exp(-1i*(xx3d*alpha + zz3d*beta)));
        
        ndim = 3
        data3d(:,:,ndim+1) = reshape(utmp3d,size(data3d(:,:,1)));
        data3d(:,:,ndim+2) = reshape(vtmp3d,size(data3d(:,:,1)));
        data3d(:,:,ndim+3) = reshape(wtmp3d,size(data3d(:,:,1)));
             
        fname3D=[ 'OptimalOC3D_' num2str(nic,'%02i') '.fld'];
        status = writenek(fname3D,data3d,lr13d,elmap3d,time,istep,fields3d,emode,wdsz,etag);
        fprintf('  3D Optimal OC field (r=%i) written: %s\n',nic,fname3D)
    end    
    disp('Done.')
