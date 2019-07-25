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
lw=5;
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
W = [47,-1.2;-0.4,54];
theta = [5.6;6.4];



%% solve ODEs
f =@(x) 1./(1+exp(-x));
g =@(x) log(x./(1-x));

sinf =@(r,d) b.*r.*d./(1+b.*r.*d);  % steady solution when sdot = 0
dinf =@(r) 1./(1+a.*r);             % steady solution when ddot = 0
sigma =@(r) sinf(r,dinf(r));        % steady solution when sdot = ddot =0

% dimensionless parameters
ti = 0;
tf = 400;
dt = 0.01;
trange = ti:dt:tf;
nt = length(trange);
tspan = [ti tf];

% odes with time-varying stimuli
t0 = 20;    % onset time of stimuli
tau_fixed = 60;       % duration of stimuli
Iapp_fixed = 0*0.45;   % amplitude of stimuli


%% find all 3^n fxpts with an initial guess r0

% random search for fxpts
r0 = 0.5*ones(n,1);
tic
R = fxpt_solver(r0,W,theta,0*Iapp_fixed,a,b);   % R contains all 3^n fxpts
toc

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
    disp('found zero real eigenvalues => non-hyperbolic fxpt. Exit')
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

% % % -----
% % a few sample points in Iapp_vec and tau_vec
% % with given W and theta:
% % W = [47,-1.2;-0.4,54];
% % theta = [5.6;6.4];

% % initial state = (OFF,OFF)
tau_vec = [5,11,23,40];  
Iapp_fixed = 2;
ntrial=length(tau_vec);

% Iapp_vec = [1.5,1.758,2,3];
% tau_fixed = 27;
% ntrial=length(Iapp_vec);
% % % -------------------









%% ========

for itrial=1:ntrial
% % % fixed amplitude and duration
% Iapp = Iapp_fixed;      tau = tau_fixed;

% % % fixed amplitude and varying duration
Iapp = Iapp_fixed;      tau = tau_vec(itrial);

% % fixed duration and varying amplitude
% Iapp = Iapp_vec(itrial);      tau = tau_fixed;



Iext =@(t) Iapp*(heaviside(t-t0)-heaviside(t-t0-tau)).*ones(n,1);

rdot =@(t,r,s,d)    -r + f(W*s - theta + Iext(t));
sdot =@(r,s,d)   (-s + b.*r.*d.*(1-s))./ts;
ddot =@(r,s,d)   (1 - d - a.*r.*d)./td;

deqns =@(t,X) [
               rdot(t, X(1:n),X(n+1:2*n),X(2*n+1:3*n) );
               sdot( X(1:n),X(n+1:2*n),X(2*n+1:3*n) );
               ddot( X(1:n),X(n+1:2*n),X(2*n+1:3*n) )
              ];

options = odeset('RelTol',1e-13,'InitialStep',1e-6,'MaxStep',1);

% set the initial state
% X0 = rand(3*n,1);
X0 = Xc_vec(:,1);       % initial state as (OFF,OFF)

% integrate from ti to tf
[t,X] = ode23s(deqns,tspan,X0,options);
r = X(:,1:n);  s = X(:,n+1:2*n);  d = X(:,2*n+1:3*n);


% % % -----------------
% % time evolutions of units
figure

hold on

if Iapp>0
    rectangle('Position',[t0,0.01,tau,Iapp],'FaceColor', [0.8 0.8 0.8],'linestyle','none')
else
    rectangle('Position',[t0,Iapp-0.01,tau,abs(Iapp)],'FaceColor', [0.8 0.8 0.8],'linestyle','none')
end

plot(t,r,'linewidth',lw)
ylim([0 1.5])
yticks([0 0.5 1])

set(gca,'fontsize',fsize,'fontname','times','linewidth',aw)
axis square
box on
xlabel('$t$','FontSize',fsize,'Interpreter','latex')
ylabel('$r$','FontSize',fsize,'Interpreter','latex')
text(200,0.3,sprintf('$\\tau_{\\rm{dur}} = %.0f$', tau),'FontSize',fsize,'Interpreter','latex')
text(200,0.5,sprintf('$I_{\\rm{app}} = %.0f$', Iapp),'FontSize',fsize,'Interpreter','latex')









% fig=figure;
% 
% left_color = [0 0 0];
% right_color = [0.8 0.8 0.8];
% set(fig,'defaultAxesColorOrder',[left_color; right_color]);
% hold on
% 
% yyaxis right
% if Iapp>0
%     rectangle('Position',[t0,0.01,tau,Iapp],'FaceColor', [0.8 0.8 0.8 0.6],'linestyle','none')
% else
%     rectangle('Position',[t0,Iapp-0.01,tau,abs(Iapp)],'FaceColor', [0.8 0.8 0.8 0.6],'linestyle','none')
% end
% ylim([0 2.5])
% yticks([0 1 2])
% ylabel('$I_{\rm{app}}$','FontSize',fsize,'Interpreter','latex')
% 
% yyaxis left
% plot(t,r(:,1),'b-','linewidth',lw)
% plot(t,r(:,2),'r-','linewidth',lw)
% set(gca,'fontsize',fsize,'fontname','times','linewidth',aw)
% axis square
% box on
% xlabel('$t$','FontSize',fsize,'Interpreter','latex')
% ylabel('$r$','FontSize',fsize,'Interpreter','latex')
% text(200,0.3,sprintf('$\\tau_{\\rm{dur}} = %.0f$', tau),'FontSize',fsize,'Interpreter','latex')
end

% % % % ---------------------


% % % %----------------------------------
% % % trajectories in r1-r2 plane
% figure
% hold on
% plot(r(:,1),r(:,2),'k','linewidth',lw)
% plot(r(1,1),r(1,2),'k>','markersize',12)
% % fixed points on r1-r2 plane
% for i=1:length(attractors_idx)
%     plot(rc_attractors(1,i),rc_attractors(2,i),'.','markersize',35)
% %     plot(rc_attractors(1,i),rc_attractors(2,i),'.','markersize',35,...
% %     'Color',[1,1-0.9/(2^n)*i,1-0.9/(2^n)*i])
% %     'Color',my_colormap(i,:))
%     hold on
% end
%     plot(rc_saddles(1,:),rc_saddles(2,:),'kx','markersize',10)
%     axis(axisrange2d)
%     axis square
%     xticks([0 0.5 1])
%     yticks([0 0.5 1])
%     set(gca,'fontsize',fsize,'fontname','times','linewidth',aw)
%     xlabel('$r_1$','FontSize',fsize,'Interpreter','latex')
%     ylabel('$r_2$','FontSize',fsize,'Interpreter','latex')
% % 	legend(legend_list,'Location','bestoutside') 
% % % %-------------------------------------------



% % % % -----------------------
% % % s-d phase plane & nullcines
% figure(22)
% hold on
% plot(s(:,1),d(:,1),'b','linewidth',1.5)
% plot(s(:,2),d(:,2),'r','linewidth',1.5)
% % % plot the initial state
% plot(s(1,1),d(1,1),'b>','markersize',12)
% plot(s(1,2),d(1,2),'r>','markersize',12)
% % % plot the final state
% plot(s(end,1),d(end,1),'b^','markersize',12)
% plot(s(end,2),d(end,2),'r^','markersize',12)
% 
% 
% % % nullclines for unit 1
% fimplicit(@(s1,d1) -s1+b(1)*f(W(1,1)*s1-theta(1)+0*Iapp)*d1*(1-s1), [0 1 0 1],'k','linewidth',1)
% fimplicit(@(s1,d1) 1-d1-a(1)*f(W(1,1)*s1-theta(1)+0*Iapp)*d1, [0 1 0 1],'k','linewidth',1)
% 
% % % % nullclines for unit 2
% % fimplicit(@(s2,d2) -s2+b(2)*f(W(2,2)*s2-theta(2)+0*Iapp)*d2*(1-s2), [0 1 0 1],'k-.','linewidth',1)
% % fimplicit(@(s2,d2) 1-d2-a(2)*f(W(2,2)*s2-theta(2)+0*Iapp)*d2, [0 1 0 1],'k-.','linewidth',1)
% 
% axis([0 0.4 0 1])
% axis square
% set(gca,'fontsize',fsize,'fontname','times','linewidth',aw)
% xlabel('$s$','FontSize',fsize,'Interpreter','latex')
% ylabel('$d$','FontSize',fsize,'Interpreter','latex')
% % % % -----------------------------


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