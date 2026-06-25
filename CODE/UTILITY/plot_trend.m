function[ k , r2 , h , r, p , m] = plot_trend(x,y,col,fig , line_width)
%% Edited 2023 - 07 -20

mdl = fitlm(x,y ) ;

r2    = mdl.Rsquared.Ordinary;
r = sqrt(r2);
m = table2array(mdl.Coefficients('(Intercept)', 'Estimate' ));
k = table2array(mdl.Coefficients('x1', 'Estimate' ) );

p = coefTest(mdl);

%disp(['Equation is y = ' num2str(m) '*x + ' num2str(k)]);

hold on
%plot(x,y_est, 'color', col, 'LineWidth',1)
h = plot(mdl,'HandleVisibility','off');
delete(h(1))

delete(findall(fig, 'type', 'legend'));

cbHandles = findobj(h,'DisplayName','Confidence bounds');
cbHandles = findobj(h,'LineStyle',cbHandles.LineStyle, 'Color', cbHandles.Color);
fitHandle = findobj(h,'DisplayName','Fit');

fitHandle.Color = col; %orange 

if nargin > 4  
set(fitHandle, 'LineWidth', line_width)
end

set(cbHandles, 'Color', col, 'LineWidth', 1)

lgd = findobj('type', 'legend');
delete(lgd)

title('')
ylabel('')
xlabel('')
end

