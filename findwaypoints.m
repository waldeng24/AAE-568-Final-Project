TPBVP2
xpoints(1) = y(1,1);
ypoints(1) = y(3,1);
zpoints(1) = y(5,1);
for i = 2:length(y)
    xpoints(i) = y(1,i);
    ypoints(i) = y(3,i);
    gpoints(i) = 200*griewank([y(1,i) y(3,i)]);
    zpoints(i) = y(5,i);
    if zpoints(i)<gpoints(i)
        zpoints(i)=gpoints(i);
    end
end
figure(4)
plot3(xpoints,ypoints,zpoints,'--c',linewidth=3)