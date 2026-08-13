function [ic, S, V] = solve_brusselator_1d(ic, par, T, flag)
% solves brusselator equation on [0,L]
% flag 	=0 no plot
%		=1 space-time plot

%% load initial data
	%ic = load('ic.txt');

%% setup
	x  = ic(:,1);
	u  = [ic(:,2); ic(:,3)];
	N  = length(x);
	dx = x(2)-x(1);
	L  = (N+1)*dx;

%% parameters
	% De Wit
	a	= 2.5;
	d1	= 4.11;
	d2	= 9.73;

	sigma = par.sigma;
	b     = par.b;
	d1    = sigma*d2;
	
	% Tzou
	%a	= 1.4;
	%b	= 3.21;
	%d1	= 0.2666;
	%d2	= 1;
	
%% stepsize
	t  = 0;
	dt = 1;
	%disp(['Time stepsize:    ',num2str(dt)]);
	%disp(['Spatial stepsize: ',num2str(dx)]);
	
%% compute finite difference approximation D2 of second derivative
	e  = ones(N,1);
	D2 = sparse(1:N-1,[2:N-1 N],ones(N-1,1),N,N) - sparse(1:N,[1:N],e,N,N);
	D2 = D2+D2';
	D2(1,2) = 2; D2(N,N-1) = 2; % Neumann boundary conditions
								% (comment out for Dirichlet conditions)
	D2 = D2/dx^2;
	
%% define right-hand side of partial differential equation
	function f = rhs(t,u)
		u1 = u(1:N);
		u2 = u(N+1:2*N);
		f  = [d1*D2*u1 + a - (b+1)*u1 + u1.^2.*u2; ...
			  d2*D2*u2 + b*u1 - u1.^2.*u2];
	end
	
%% compute solution
	[S, UV] = ode45(@rhs, [0:0.1:T], u);
	V = UV(:, N+1:end);
	ic = [x, UV(end,1:N)', UV(end,N+1:end)'];

%% space-time plot
	if flag == 1
		figure(1); clf;
		image(V, 'CDataMapping', 'scaled');
		ax = gca; ax.YDir = 'normal';
%		colormap(ax, brewermap([],'RdYlBu'));
		colormap(ax, 'sky');
		colorbar(ax);
		drawnow;
	end

end 

