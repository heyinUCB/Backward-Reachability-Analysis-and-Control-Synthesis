%%
y=linspace(-pi/12,pi/12);y=y(:);
csin = y.^3\(sin(y)-y);
ccos = y.^2\(cos(y)-1);

if 1
    % Plot LS fits and compare w/ sine and cosine
    x = linspace(-pi/12,pi/12);
    figure
    subplot(211)
    plot(x,sin(x),'xr',x,csin*x.^3+x,'ob')
    xlabel('x');ylabel('sin');
    legend('sin','LS fit','location','best');
    subplot(212)
    plot(x,cos(x),'xr',x,ccos*x.^2+1,'ob')
    xlabel('x');ylabel('cos');
    legend('cos','LS fit','location','best');
end
% sin(x) is approximated by -0.166*x^3 + x for x in [-pi/12, pi/12].
% cos(x) is approximated by -0.498*x^2 + 1 for x in [-pi/12, pi/12].