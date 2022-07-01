t = linspace(0,2*pi);
R=0.1/N*2;
hold on 
for i = 0:0
    for n = 1:N
        plot(c(n,1)+R*cos(t)+i, c(n,2)+R*sin(t),'k')  
    end
end
daspect([1 1 1])
