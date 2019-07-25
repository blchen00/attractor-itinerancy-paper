%% two bistable units with constant stimulus
%  update: 08-13-18
%  1. Varying weight matrix to estimate reachable final states
%  2. Show transient state-transitions in rich regions in I-tau phase
%  diagram.

clear
sympref('HeavisideAtOrigin', 1);  %  heaviside(0)=1

%% parameters
n = 2;
taur = 0.01*ones(n,1);
taus = 0.05*ones(n,1);
taud = 0.25*ones(n,1);
ts = taus./taur;    % dimensionless synaptic time constant
td = taud./taur;    % dimensionless depression time constant

alpha = 1./ts;
epsilon = 1./td;

rmax = 50*ones(n,1);
p0 = 0.5*ones(n,1);
alp = ones(n,1);

a = p0.*taud.*rmax;         % dimensionless depression parameter
b = alp.*p0.*taus.*rmax;    % dimensionless synpatic parameter

% parameters for graphics
fsize = 32;                 % set fontsize for plot labels
lw=2.5;
aw=1.5;
axisrange2d = [-0.05 1.05 -0.05 1.05];    
axisrange3d = [-0.05 1.05 -0.05 1.05 -0.05 1.05];


%% bistable region in the w-theta plane
% location of the cusp point
wc = 4*(a+b+1)./b;
thc = 2+log(a+b+1);

% choose self-excitation and threshold in the bistable region
% already in dimensionless form
w = 40*ones(n,1);     % self-excitation weight
W = diag(w);          % diagonal matrix

U = 2*(rand(n)-0.5);  % uniformly-random connection matrix \in (-1,1)
U = U-diag(diag(U));  % only added to off-diagnal elements

W = W + 0*U;
w = diag(W);          % self-coupling strength as diagnal elements

theta = 5*ones(n,1);
theta_noise = 2*(rand(n,1)-0.5);    % uniformly-random noise \in (-1,1) 
theta = theta + 0*theta_noise;      % added to theshoulds 
% 
% W = [47,-1.2;-0.4,54];
% % % W = eye(2).*(diag(W));
% theta = [5.6;6.4]
% W = [40,-1;-1,40];
% W = [40,0.25;0.25,40];


%% solve ODEs
f =@(x) 1./(1+exp(-x));
g =@(x) log(x./(1-x));

sinf =@(r,d) b.*r.*d./(1+b.*r.*d);  % steady solution when sdot = 0
dinf =@(r) 1./(1+a.*r);             % steady solution when ddot = 0
sigma =@(r) sinf(r,dinf(r));        % steady solution when sdot = ddot =0

% dimensionless parameters
ti = 0;
tf = 1000;
dt = 0.01;
trange = ti:dt:tf;
nt = length(trange);
tspan = [ti tf];

% odes with time-varying stimuli
t0 = 1;    % onset time of stimuli
tau_fixed = 60;       % duration of stimuli
Iapp_fixed = 0*0.45;   % amplitude of stimuli

options = odeset('RelTol',1e-13,'InitialStep',1e-6,'MaxStep',1);


% % % create sub-directories for data and fig storage
current_path = pwd();
new_fig_dir = '/phase_diagram_figs/';
new_data_dir = '/phase_diagram_data/';
mkdir(current_path,new_fig_dir);
mkdir(current_path,new_data_dir);


%% find all 3^n fxpts with an initial guess r0
% random search for fxpts
r0 = 0.5*ones(n,1);
R = fxpt_solver(r0,W,theta,0*Iapp_fixed,a,b);   % R contains all 3^n fxpts

rc = zeros(3^n,n);
rc_norm = zeros(3^n,1);
for i=1:3^n
    rc(i,:) = R(:,:,i+1);           % save fxpts to rc
    rc_norm(i) = norm(rc(i,:));     % compute the norm of a fxpt for sorting
end

[rcnorm_sort,rcnorm_sort_idx] = sort(rc_norm);
rc_sort = rc(rcnorm_sort_idx,:);    % sort rc according to norms, small->large

rc_vec = rc_sort';                      % get rc vectors
sc_vec = b.*rc_vec./(1+(a+b).*rc_vec);  % get sc vectors
dc_vec = 1./(1+a.*rc_vec);              % get dc vectors
Xc_vec = zeros(3*n,(n+1)^2);            
for i=1:3^n
Xc_vec(:,i) = [rc_vec(:,i);sc_vec(:,i);dc_vec(:,i)];    % get Xc vectors
end

% compute the eigenvalues of all found fxpts
Lambda = linearization(rc_vec,W,a,b,alpha,epsilon);
re_Lambda = real(Lambda);           % real parts of eigenvalues
im_Lambda = imag(Lambda);           % imaginary parts of eigenvalues
sgn_re_Lambda = sign(re_Lambda);    % positive real parts of eigenvalues 
if any(any(sgn_re_Lambda==0))
    disp('found zero real eigenvalues => non-hyperbolic fxpt. EXIT!!!')
else
    fxpts_idx = sum(sgn_re_Lambda==1); % get Poincare indices for each fxpt
end

% color coding attractors
attractors_idx =find(fxpts_idx==0);
saddles_idx = find(fxpts_idx>0);
rc_attractors = rc_vec(:,attractors_idx);       % attractor vectors
rc_saddles = rc_vec(:,saddles_idx);             % saddle vectors

% get digital representation of attractors: if rc<median => '0'; else '1'
temp_list = num2str(transpose(rc_attractors>median(rc_attractors,2)));
% convert digitals to binary characters for legends
legend_list = dec2bin(bin2dec(temp_list));

% my_colormap = brewermap(length(attractors_idx),'YlGnBu');
number_attractors = length(attractors_idx);     % number of stable nodes
number_saddles = length(saddles_idx);           % number of saddles



%% vectors for varying amplitudes and durations
ntrial = 30;
Iapp_vec = linspace(0,5,ntrial); 
tau_vec = linspace(1,200,ntrial); 

% duration and amplitude mesh
[taum,Iappm] = meshgrid(tau_vec,Iapp_vec);
zm = zeros(size(taum));
ngrid = numel(taum);

parfor itrial=1:ngrid
 
% varying durations and amplitudes
Iapp = Iappm(itrial);             tau = taum(itrial);

Iext =@(t) Iapp*(heaviside(t-t0)-heaviside(t-t0-tau)).*ones(n,1);

rdot =@(t,r,s,d)    -r + f(W*s - theta + Iext(t));
sdot =@(r,s,d)   (-s + b.*r.*d.*(1-s))./ts;
ddot =@(r,s,d)   (1 - d - a.*r.*d)./td;

deqns =@(t,X) [
               rdot(t, X(1:n),X(n+1:2*n),X(2*n+1:3*n) );
               sdot( X(1:n),X(n+1:2*n),X(2*n+1:3*n) );
               ddot( X(1:n),X(n+1:2*n),X(2*n+1:3*n) )
              ];

% set the initial state 
X0 = Xc_vec(:,1);       % initial state as (OFF,OFF)

X0 = Xc_vec(:,attractors_idx(1));  
% integrate with a specific initial condition
    [~,X] = ode45(deqns,tspan,X0,options);

    dX_vec = X(end,:)' - Xc_vec;
    dist_vec = vecnorm(dX_vec);
    
    [~,dist_min_idx] = min(dist_vec);
    asymp_state_idx = dist_min_idx;
    
    zm(itrial) = asymp_state_idx;
    disp([itrial,asymp_state_idx])

   X=[];dx_vec=[];dist_vec=[];dist_min_idx=[];asymp_state_idx=[];
end


figure 
zm_colorcode = zm;
for i=1:length(attractors_idx)
    zm_colorcode(zm_colorcode==attractors_idx(i))=i;
end
basin_plt = pcolor(taum,Iappm,zm_colorcode);
hold on

colormap(lines(number_attractors))
set(basin_plt,'edgecolor','none')
caxis([1 number_attractors])
colorbar('Ticks',1:length(attractors_idx),'TickLabels',legend_list)

axis square
set(gca,'fontsize',fsize,'fontname','times')
xlabel('$\tau_{\rm{dur}}$','FontSize',fsize,'Interpreter','latex')
ylabel('$I_{\rm{app}}$','FontSize',fsize,'Interpreter','latex')


%% %% save data to file
data_pathname = './phase_diagram_data/';
data_filename = sprintf('phase_n%d.txt', n);
dlmwrite(fullfile(data_pathname,data_filename),zm_colorcode,'delimiter','\t');
% 

% save figs
fig_pathname = './phase_diagram_figs/';
fig_filename = sprintf('phase_diagram_n_%d', n);
saveas(gcf,fullfile(fig_pathname,fig_filename),'epsc');
saveas(gcf,fullfile(fig_pathname,fig_filename),'fig');

%% save data to file
data_filename1 = sprintf('w_theta.txt');
data_filename2 = sprintf('numAttractors.txt');
data_filename3 = sprintf('attractorIndices.txt');
% 
% data format
dlmwrite(fullfile(data_pathname,data_filename1),W, ...
    '-append','delimiter','\t');

dlmwrite(fullfile(data_pathname,data_filename1),theta, ...
    '-append','delimiter','\t');

dlmwrite(fullfile(data_pathname,data_filename2),[n, number_attractors, attractors_idx], ...
    '-append','delimiter','\t');

dlmwrite(fullfile(data_pathname,data_filename3),[n, attractors_idx], ...
    '-append','delimiter','\t');



%% Evaluating all eigenvalues of all fxpts
function eigvals_vec = linearization(rc_vec,W,a,b,alpha,epsilon)
    [n,m] = size(rc_vec);  % n=# of units; m=# of fxpts = 3^n
    krondel =@(i,j) double(i==j);
    
    eigvals_vec = zeros(3*n,m);
    
    for k = 1:m
        r = rc_vec(:,k);
        L = [];
        for i = 1:n
            l = [];
            for j = 1:n

                l11 = -krondel(i,j);
                l12 = W(i,j)*r(i)*(1-r(i));
                l13 = 0;
                l21 = krondel(i,j)*alpha(i)*b(i)/(1+(a(i)+b(i))*r(i));
                l22 = -krondel(i,j)*alpha(i)*(1+(a(i)+b(i))*r(i))/(1+a(i)*r(i));
                l23 = krondel(i,j)*alpha(i)*b(i)*r(i)*(1+a(i)*r(i))/(1+(a(i)+b(i))*r(i));
                l31 = -krondel(i,j)*epsilon(i)*a(i)/(1+a(i)*r(i));
                l32 = 0;
                l33 = -krondel(i,j)*epsilon(i)*(1+a(i)*r(i));

                l_ij =   [l11,l12,l13;
                          l21,l22,l23;
                          l31,l32,l33];

                l = [l,l_ij];
            end
            L = [L;l];
        end
%         eigvals_vec(:,k) = eig(L);
        eigvals_vec(:,k) = vpa(eig(sym(L)),30);
    end
  
end

%% Finding all 3^n fxpts
function R = fxpt_solver(R0,W,theta,Iext,a,b)
% function R = fxpt_solver(R0,W,theta,Iext,a,b,rmin,rmax)

    tolerance = 1e-15;
%     itermax = 1e3;
    [n,m,i]=size(R0);
    
    rmin=0; rmax=1;   
    
    g =@(x) log(x./(1-x));
    sinf =@(x) b.*x./(1+(a+b).*x);
    F =@(x) g(x) - (W*sinf(x) - theta + Iext);
    
    opt=optimset; opt.Display='off'; opt1=opt;
    opt.TolX=tolerance;  opt.TolFun=tolerance; 
    
    R = R0;
%     iter = 0;
    while i<3^n+1   %&& iter<itermax
    r0 = rmin + (rmax-rmin)*rand(n,m);
    [r,fr,key] = fsolve(@(r) F(r),r0,opt1);
    
    
    if key>0
        N=size(R,3);
        for j=1:N
            if norm(R(:,:,j)-r)<1e-5
                key = 0; 
                break; 
            end
        end
        if key>0
            [r1,fr,key] = fsolve(@(r) F(r),r,opt);
            if norm(r-r1)<1e-5 & key>0
                R(:,:,i+1) = r1;
                i = i+1;
            end
        end
    end
%     iter = iter+1;
    end

end