function [fitresult, gof] = createFit(bn, hst)
%CREATEFIT(BN,HST)
%  Create a fit.
%
%  Data for 'untitled fit 2' fit:
%      X Input: bn
%      Y Output: hst
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 27-Jun-2022 15:41:22


%% Fit: 'untitled fit 2'.
[xData, yData] = prepareCurveData( bn, hst );

% Set up fittype and options.
ft = fittype( 'a*exp(-(x^2)/(2*var2))+b*exp(-((x-xplus)^2)/(2*var2))+c*exp(-((x+xminus)^2)/(2*var2))', 'independent', 'x', 'dependent', 'y' );
excludedPoints = abs(xData) > 40;
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 0 0 0 0];
opts.StartPoint = [0.94707022979172 0.994312868974457 0.698307887979348 0.48383621174715 0.932513710090837 0.750090183634553];
opts.Exclude = excludedPoints;

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 2' );
h = plot( fitresult, xData, yData, excludedPoints, 'predobs' );
legend( h, 'hst vs. bn', 'Excluded hst vs. bn', 'untitled fit 2', 'Lower bounds (untitled fit 2)', 'Upper bounds (untitled fit 2)', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'bn', 'Interpreter', 'none' );
ylabel( 'hst', 'Interpreter', 'none' );
grid on


