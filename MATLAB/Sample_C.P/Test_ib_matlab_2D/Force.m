function F=Force(X)
global kp km dtheta K;

F=K*(X(kp,:)+X(km,:)-2*X)/(dtheta*dtheta);