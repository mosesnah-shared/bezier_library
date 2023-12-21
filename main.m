%% Initialization 
clear all; close all; clc;

%% Draw the trajectory and save the Bezier Parameters
Runit_interactiveBezier

%% Check the parameters and plot the trajectory
% Parameters saved as Bezier structure
% The number of trajectories 
Ntraj = length( Bez );
Nd    = size( Bez( 1 ).Q, 2 );
Ns    = 200;
t     = linspace( 0, 1, Ns );

f = figure( ); a = axes( 'parent', f );
axis equal
hold( a, 'on' )

% Whole x y data and its time array 
x_arr = zeros( 1, Ns*Ntraj );
y_arr = zeros( 1, Ns*Ntraj );
t_arr = zeros( 1, Ns*Ntraj );

for i = 1 : Ntraj
    [ C, theta, K, L]  = BezierEval( Bez( i ).Q, t, Bez( i ).w ); %rational
    plot( a, C( 1, : ), C( 2, : ) ) 
    t_arr( (i-1)*Ns+1: i*Ns ) = t+(i-1);
    x_arr( (i-1)*Ns+1: i*Ns ) = C( 1, : );
    y_arr( (i-1)*Ns+1: i*Ns ) = C( 2, : );
end

% Find the best-fit curve
% For discrete movement
fx = fit( t_arr', x_arr', 'poly9' );
fy = fit( t_arr', y_arr', 'poly9' );

% Double-check with the plot
plot( a, fx( t_arr ), fy( t_arr ), 'linewidth', 1, 'color', 'k' );

% Saving also the derivatives
[ dfx, ddfx ] = differentiate( fx, t_arr );
[ dfy, ddfy ] = differentiate( fy, t_arr );

% Saving the data for DMP MATLAB
  p_data = [   fx( t_arr ),   fy( t_arr ) ]';
 dp_data = [           dfx,         dfy   ]';
ddp_data = [          ddfx,        ddfy   ]';

data = struct( );
data.t_arr    =    t_arr;
data.p_data   =   p_data;
data.dp_data  =  dp_data;
data.ddp_data = ddp_data;

save( 'O_distorted.mat', 'data' );

%%
% Find the best-fit curve
fx = fit( t_arr', x_arr', 'poly9' );
fy = fit( t_arr', y_arr', 'poly9' );

% Get the derivatives
[ dfx, ddfx ] = differentiate( fx, t_arr );
[ dfy, ddfy ] = differentiate( fy, t_arr );

subplot( 2, 1, 1 );
plot( t_arr, x_arr ) 
hold on
plot( t_arr, dfx' )
plot( fx )

subplot( 2, 1, 2 );
plot( t_arr, y_arr ) 
hold on
plot( fy )


%%
for i = 1 : Ntraj
    traj = 0;
    for j = 1 : Nd
        P = Bez( i ).Q( :, j );        
        traj = traj + nchoosek( Nd-1, j-1 ) * (1-t_arr).^(Nd-j) .* t_arr.^(j-1) .* P ;
    end
    plot( a , traj( 1, : ), traj( 2, : ), 'linewidth', 5 )
end

%%
% Also, we lift up the order 
Bez2 = RaiseBezier( Bez, 4 );
Nd2  = size( Bez2( 1 ).Q, 2 );

for i = 1 : Ntraj
    traj = 0;
    for j = 1 : Nd2
        P = Bez2( i ).Q( :, j );        
        traj = traj + nchoosek( Nd2-1, j-1 ) * (1-t_arr).^(Nd2-j) .* t_arr.^(j-1) .* P ;
    end
    plot( a , traj( 1, : ), traj( 2, : ), 'linewidth', 5 )
end

%% We also plot the symbolic form of the trajectory
syms tsym

  traj_sym  = cell( 1, Ntraj );
 dtraj_sym  = cell( 1, Ntraj );
ddtraj_sym  = cell( 1, Ntraj );

  traj_func = cell( 1, Ntraj );
 dtraj_func = cell( 1, Ntraj );
ddtraj_func = cell( 1, Ntraj );

  traj_sym2  = cell( 1, Ntraj );
 dtraj_sym2  = cell( 1, Ntraj );
ddtraj_sym2  = cell( 1, Ntraj );

  traj_func2 = cell( 1, Ntraj );
 dtraj_func2 = cell( 1, Ntraj );
ddtraj_func2 = cell( 1, Ntraj );


for i = 1 : Ntraj
    traj_sym{ i } = 0;

    for j = 1 : Nd
        P = Bez( i ).Q( :, j );        
        traj_sym{ i } = traj_sym{ i } + nchoosek( Nd-1, j-1 ) * (1-tsym).^(Nd-j) .* tsym.^(j-1) .* P ;
    end

    traj_sym2{ i } = 0;
    for j = 1 : Nd2
        P = Bez2( i ).Q( :, j );        
        traj_sym2{ i } = traj_sym2{ i } + nchoosek( Nd2-1, j-1 ) * (1-tsym).^(Nd2-j) .* tsym.^(j-1) .* P ;
    end

     dtraj_sym{ i } = diff(  traj_sym{ i }, tsym );
    ddtraj_sym{ i } = diff( dtraj_sym{ i }, tsym );

     dtraj_sym2{ i } = diff(  traj_sym2{ i }, tsym );
    ddtraj_sym2{ i } = diff( dtraj_sym2{  i }, tsym );    

    % Save these to MATLAB functions
      traj_func{ i } = matlabFunction(   traj_sym{ i } );
     dtraj_func{ i } = matlabFunction(  dtraj_sym{ i } );
    ddtraj_func{ i } = matlabFunction( ddtraj_sym{ i } );

    % Save these to MATLAB functions
      traj_func2{ i } = matlabFunction(   traj_sym2{ i } );
     dtraj_func2{ i } = matlabFunction(  dtraj_sym2{ i } );
    ddtraj_func2{ i } = matlabFunction( ddtraj_sym2{ i } );    
end

%% Plotting the velocity and acceleration
f = figure( ); a = axes( 'parent', f );
hold( a, 'on' ); 
t0 = linspace( 0, 1, 100); 
for i = 1 : Ntraj
    t = linspace( i-1, i, 100 );
    plot( a, t,  dtraj_func{ i }( t0 ), 'linewidth', 4, 'color', 'k' )
    plot( a, t, ddtraj_func{ i }( t0 ), 'linewidth', 4, 'color', 'r' )

    plot( a, t,  dtraj_func2{ i }( t0 ), 'linewidth', 4, 'color', 'k' )
    plot( a, t, ddtraj_func2{ i }( t0 ), 'linewidth', 4, 'color', 'g' )    
    
end
