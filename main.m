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


% Once we get the whole x_arr and y_arr, we set this as origin 0 
x_arr = x_arr - x_arr( 1 );
y_arr = y_arr - y_arr( 1 );

% Fitting the function, which was manually discovered
myfittype = fittype("( 6*(t/a)^5-15*(t/a)^4+10*(t/a)^3 ) *( b10*(t-a)^10 + b9*(t-a)^9 + b8*(t-a)^8 + b7*(t-a)^7 + b6 * (t-a)^6 + b5 * (t-a)^5 + b4 * (t-a)^4 + b3*(t-a)^3 + c )", ...
                'independent',"t", 'coefficients', [ "b3" "b4" "b5" "b6" "b7" "b8" "b9" "b10"  ], 'problem' ,["a", "c"] );

% Specify the value for a_fix, c_fix
a_fix = max( t_arr );
c_fix = x_arr( end );

% Perform the fit, passing the fixed 'b' as a 'problem' parameter
f = figure( ); a = axes( 'parent', f );
[fitresult_x, ~] = fit( t_arr', x_arr', myfittype, 'problem', { a_fix, c_fix } );
plot( a, t_arr, fitresult_x( t_arr ), t_arr, x_arr )

% Specify the value for 'b'
a_fix = max( t_arr );
c_fix = y_arr( end );

% Perform the fit, passing the fixed 'b' as a 'problem' parameter
f = figure( ); a = axes( 'parent', f );
[fitresult_y, ~] = fit( t_arr', y_arr', myfittype, 'problem', { a_fix, c_fix } );
plot( a, t_arr, fitresult_y( t_arr ), t_arr, y_arr )


% Double-check with the plot
f = figure( ); a = axes( 'parent', f );
plot( a, fitresult_x( t_arr ), fitresult_y( t_arr ), 'linewidth', 1, 'color', 'k' );
hold on
plot( a, x_arr, y_arr, 'linewidth', 1, 'color', 'c' );

% Saving also the derivatives
[ dfx, ddfx ] = differentiate( fitresult_x, t_arr );
[ dfy, ddfy ] = differentiate( fitresult_y, t_arr );

% Saving the data for DMP MATLAB
  p_data = [   fitresult_x( t_arr ),   fitresult_y( t_arr )]';
 dp_data = [           dfx,         dfy   ]';
ddp_data = [          ddfx,        ddfy   ]';

data = struct( );
data.t_arr    =    t_arr;
data.p_data   =   p_data;
data.dp_data  =  dp_data;
data.ddp_data = ddp_data;


save( 'C.mat', 'data' );
