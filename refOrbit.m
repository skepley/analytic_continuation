lorenz = @(t,y)vectorizeinitdata(t,y);
ICs = [1;1;1];
T = [0,100];
[t,sol] = ode45(lorenz,T,ICs);
PX = sol(:,1:3:end)';
PY = sol(:,2:3:end)';
PZ = sol(:,3:3:end)';
orbIDX = 200:length(t);

orbit = [PX(orbIDX);PY(orbIDX);PZ(orbIDX)];
gcf;
hold on
plot3(orbit(1,:),orbit(2,:),orbit(3,:),'r','LineWidth',.25)

