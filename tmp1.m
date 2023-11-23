% Generate curve
t=linspace(0,1,100);

f = figure( ); a = axes( 'parent', f );
axis equal
hold( a, 'on' )
for i = 1 : length( Bez )
    [ C, theta, K, L]  = BezierEval( Bez(i).Q, t, Bez(i).w); %rational
    plot( a, C( 1, : ), C( 2, : ) ) 
end

% Drawing the Bezier Curve
t2 = t;

points = Bez(1).Q;

P1 = points( :, 1 );
P2 = points( :, 2 );
P3 = points( :, 3 );
P4 = points( :, 4 );
P5 = points( :, 5 );
P6 = points( :, 6 );

traj = nchoosek(5,0) * (1-t2).^5 .* t2.^0 .* P1 + ...
     + nchoosek(5,1) * (1-t2).^4 .* t2.^1 .* P2 ...
     + nchoosek(5,2) * (1-t2).^3 .* t2.^2 .* P3+ ...
     + nchoosek(5,3) * (1-t2).^2 .* t2.^3 .* P4+ ...
     + nchoosek(5,4) * (1-t2).^1 .* t2.^4 .* P5+ ...
     + nchoosek(5,5) * (1-t2).^0 .* t2.^5 .* P6;

plot( a, traj( 1, : ), traj( 2, : ), 'color', 'k', 'linewidth', 10 )

points = Bez(2).Q;

P1 = points( :, 1 );
P2 = points( :, 2 );
P3 = points( :, 3 );
P4 = points( :, 4 );
P5 = points( :, 5 );
P6 = points( :, 6 );

traj = nchoosek(5,0) * (1-t2).^5 .* t2.^0 .* P1 + ...
     + nchoosek(5,1) * (1-t2).^4 .* t2.^1 .* P2 ...
     + nchoosek(5,2) * (1-t2).^3 .* t2.^2 .* P3+ ...
     + nchoosek(5,3) * (1-t2).^2 .* t2.^3 .* P4+ ...
     + nchoosek(5,4) * (1-t2).^1 .* t2.^4 .* P5+ ...
     + nchoosek(5,5) * (1-t2).^0 .* t2.^5 .* P6;

plot( a, traj( 1, : ), traj( 2, : ), 'color', 'k', 'linewidth', 10 )


points = Bez(3).Q;

P1 = points( :, 1 );
P2 = points( :, 2 );
P3 = points( :, 3 );
P4 = points( :, 4 );
P5 = points( :, 5 );
P6 = points( :, 6 );

traj = nchoosek(5,0) * (1-t2).^5 .* t2.^0 .* P1 + ...
     + nchoosek(5,1) * (1-t2).^4 .* t2.^1 .* P2 ...
     + nchoosek(5,2) * (1-t2).^3 .* t2.^2 .* P3+ ...
     + nchoosek(5,3) * (1-t2).^2 .* t2.^3 .* P4+ ...
     + nchoosek(5,4) * (1-t2).^1 .* t2.^4 .* P5+ ...
     + nchoosek(5,5) * (1-t2).^0 .* t2.^5 .* P6;

plot( a, traj( 1, : ), traj( 2, : ), 'color', 'k', 'linewidth', 10 )


points = Bez(4).Q;

P1 = points( :, 1 );
P2 = points( :, 2 );
P3 = points( :, 3 );
P4 = points( :, 4 );
P5 = points( :, 5 );
P6 = points( :, 6 );

traj = nchoosek(5,0) * (1-t2).^5 .* t2.^0 .* P1 + ...
     + nchoosek(5,1) * (1-t2).^4 .* t2.^1 .* P2 ...
     + nchoosek(5,2) * (1-t2).^3 .* t2.^2 .* P3+ ...
     + nchoosek(5,3) * (1-t2).^2 .* t2.^3 .* P4+ ...
     + nchoosek(5,4) * (1-t2).^1 .* t2.^4 .* P5+ ...
     + nchoosek(5,5) * (1-t2).^0 .* t2.^5 .* P6;

plot( a, traj( 1, : ), traj( 2, : ), 'color', 'k', 'linewidth', 10 )


points = Bez(5).Q;

P1 = points( :, 1 );
P2 = points( :, 2 );
P3 = points( :, 3 );
P4 = points( :, 4 );
P5 = points( :, 5 );
P6 = points( :, 6 );

traj = nchoosek(5,0) * (1-t2).^5 .* t2.^0 .* P1 + ...
     + nchoosek(5,1) * (1-t2).^4 .* t2.^1 .* P2 ...
     + nchoosek(5,2) * (1-t2).^3 .* t2.^2 .* P3+ ...
     + nchoosek(5,3) * (1-t2).^2 .* t2.^3 .* P4+ ...
     + nchoosek(5,4) * (1-t2).^1 .* t2.^4 .* P5+ ...
     + nchoosek(5,5) * (1-t2).^0 .* t2.^5 .* P6;

plot( a, traj( 1, : ), traj( 2, : ), 'color', 'k', 'linewidth', 10 )


%% Symbolic Form
% Get the symbolic form too
syms tsym

% Drawing the Bezier Curve
t2 = t;

points = Bez(1).Q;

P1 = points( :, 1 );
P2 = points( :, 2 );
P3 = points( :, 3 );
P4 = points( :, 4 );
P5 = points( :, 5 );
P6 = points( :, 6 );


traj1_sym = nchoosek(5,0) * (1-tsym)^5 * tsym^0 * P1 + ...
          + nchoosek(5,1) * (1-tsym)^4 * tsym^1 * P2 +...
          + nchoosek(5,2) * (1-tsym)^3 * tsym^2 * P3 + ...
          + nchoosek(5,3) * (1-tsym)^2 * tsym^3 * P4 + ...
          + nchoosek(5,4) * (1-tsym)^1 * tsym^4 * P5 + ...
          + nchoosek(5,5) * (1-tsym)^0 * tsym^5 * P6;

 dtraj1_sym = diff(  traj1_sym, tsym );
ddtraj1_sym = diff( dtraj1_sym, tsym );

  traj1 = matlabFunction(   traj1_sym );
 dtraj1 = matlabFunction(  dtraj1_sym );
ddtraj1 = matlabFunction( ddtraj1_sym );

f = figure( ); a = axes( 'parent', f );
plot( t2, dtraj1( t2 ) ) 
hold on

points = Bez(2).Q;

P1 = points( :, 1 );
P2 = points( :, 2 );
P3 = points( :, 3 );
P4 = points( :, 4 );
P5 = points( :, 5 );
P6 = points( :, 6 );


traj2_sym = nchoosek(5,0) * (1-tsym)^5 * tsym^0 * P1 + ...
          + nchoosek(5,1) * (1-tsym)^4 * tsym^1 * P2 +...
          + nchoosek(5,2) * (1-tsym)^3 * tsym^2 * P3 + ...
          + nchoosek(5,3) * (1-tsym)^2 * tsym^3 * P4 + ...
          + nchoosek(5,4) * (1-tsym)^1 * tsym^4 * P5 + ...
          + nchoosek(5,5) * (1-tsym)^0 * tsym^5 * P6;

 dtraj2_sym = diff(  traj2_sym, tsym );
ddtraj2_sym = diff( dtraj2_sym, tsym );

  traj2 = matlabFunction(   traj2_sym );
 dtraj2 = matlabFunction(  dtraj2_sym );
ddtraj2 = matlabFunction( ddtraj2_sym );

plot( t2 + 1, dtraj2( t2 ) )  