# include <iostream>
# include <cmath>

using namespace std;

void orbit(double* q, double x, double y, double dx, double dy, const double mu){
  
  double r = sqrt((x+mu)*(x+mu)+y*y);
  double s = sqrt((x-1+mu)*(x-1+mu)+y*y);
  
  q[0] = dx;
  q[1] = dy;
  q[2] = x + 2*dy - (1-mu)*(x+mu)/(r*r*r) - mu*(x-1+mu)/(s*s*s);
  q[3] = y - 2*dx - (1-mu)*y/(r*r*r) - mu*y/(s*s*s);
}
void K(double* k1, double* k2, double* k3, double* k4, double* k5, double* k6, double* k7, double x, double y, double dx, double dy, double dt, double mu){
 double  tem_x, tem_y, tem_dx, tem_dy;
  tem_x=0.;
  tem_y=0.;
  tem_dx=0.;
  tem_dy=0.;
  
   orbit(k1,x,y,dx,dy,mu);
  
  tem_x = x + dt*(1./5.)*k1[0];
  tem_y = y + dt*(1./5.)*k1[1];
  tem_dx = dx + dt*(1./5.)*k1[2];
  tem_dy = dy + dt*(1./5.)*k1[3];
  
  orbit(k2,tem_x,tem_y,tem_dx,tem_dy,mu);
  
  tem_x = x + dt*(3./40.)*k1[0] + dt*(9./40.)*k2[0];
  tem_y = y + dt*(3./40.)*k1[1] + dt*(9./40.)*k2[1];
  tem_dx = dx + dt*(3./40.)*k1[2] + dt*(9./40.)*k2[2];
  tem_dy = dy + dt*(3./40.)*k1[3] + dt*(9./40.)*k2[3]; 
  
  orbit(k3,tem_x,tem_y,tem_dx,tem_dy,mu);
  
  tem_x = x + dt*(44./45.)*k1[0] - dt*(56./15.)*k2[0] + dt*(32./9.)*k3[0];
  tem_y = y + dt*(44./45.)*k1[1] - dt*(56./15.)*k2[1] + dt*(32./9.)*k3[1];
  tem_dx = dx + dt*(44./45.)*k1[2] - dt*(56./15.)*k2[2] + dt*(32./9.)*k3[2];
  tem_dy = dy + dt*(44./45.)*k1[3] - dt*(56./15.)*k2[3] + dt*(32./9.)*k3[3];
  
   orbit(k4,tem_x,tem_y,tem_dx,tem_dy,mu);
  
  tem_x = x + dt*(19372./6561.)*k1[0] - dt*(25360./2187.)*k2[0] + dt*(64448./6561.)*k3[0] - dt*(212./729.)*k4[0];
  tem_y = y + dt*(19372./6561.)*k1[1] - dt*(25360./2187.)*k2[1] + dt*(64448./6561.)*k3[1] - dt*(212./729.)*k4[1];
  tem_dx = dx + dt*(19372./6561.)*k1[2] - dt*(25360./2187.)*k2[2] + dt*(64448./6561.)*k3[2] - dt*(212./729.)*k4[2];
  tem_dy = dy + dt*(19372./6561.)*k1[3] - dt*(25360./2187.)*k2[3] + dt*(64448./6561.)*k3[3] - dt*(212./729.)*k4[3];
  
    orbit(k5,tem_x,tem_y,tem_dx,tem_dy,mu);
  
  tem_x = x + dt*(9017./3168.)*k1[0] - dt*(355./33.)*k2[0] + dt*(46732./5247.)*k3[0] + dt*(49./176.)*k4[0] - dt*(5103./18656.)*k5[0];
  tem_y = y + dt*(9017./3168.)*k1[1] - dt*(355./33.)*k2[1] + dt*(46732./5247.)*k3[1] + dt*(49./176.)*k4[1] - dt*(5103./18656.)*k5[1];
  tem_dx = dx + dt*(9017./3168.)*k1[2] - dt*(355./33.)*k2[2] + dt*(46732./5247.)*k3[2] + dt*(49./176.)*k4[2] - dt*(5103./18656.)*k5[2];
  tem_dy = dy + dt*(9017./3168.)*k1[3] - dt*(355./33.)*k2[3] + dt*(46732./5247.)*k3[3] + dt*(49./176.)*k4[3] - dt*(5103./18656.)*k5[3];
  
     orbit(k6,tem_x,tem_y,tem_dx,tem_dy,mu);
  
  tem_x = x + dt*(35./384.)*k1[0] + dt*(500./1113.)*k3[0] + dt*(125./192.)*k4[0] - dt*(2187./6784.)*k5[0] + dt*(11./84.)*k6[0];
  tem_y = y + dt*(35./384.)*k1[1] + dt*(500./1113.)*k3[1] + dt*(125./192.)*k4[1] - dt*(2187./6784.)*k5[1] + dt*(11./84.)*k6[1];
  tem_dx = dx + dt*(35./384.)*k1[2] + dt*(500./1113.)*k3[2] + dt*(125./192.)*k4[2] - dt*(2187./6784.)*k5[2] + dt*(11./84.)*k6[2];
  tem_dy = dy + dt*(35./384.)*k1[3] + dt*(500./1113.)*k3[3] + dt*(125./192.)*k4[3] - dt*(2187./6784.)*k5[3] + dt*(11./84.)*k6[3];
  
      orbit(k7,tem_x,tem_y,tem_dx,tem_dy,mu);
}


void max(double& c, double a, double b, double f, double g){
  double c1, c2;
  c1=b;
  c2=f;
 if (a>b) c1 = a;
 if (g>f) c2 = g;
 if (c1>c2) c = c1;
 else c=c2;
}

int main(){
  const double mu=0.012277471;
  double x, y, dx, dy;
  double var;
  
  double k1[4], k2[4], k3[4], k4[4], k5[4], k6[4], k7[4], tem_x, tem_y, tem_dx, tem_dy;
  
  double tol=1e-7;
  
  double dt=1e-6;
  
  cout << "0.0" << "\t" << dt << "\t"<< "0.994" << "\t" << "0.0" << "\t" << "0.0" << "\t" << "-2.00158510637908" << endl;
  
  x=0.0994;
  y=0.0;
  dx=0.0;
  dy=-2.00158510637908;
  
  tem_x=0.;
  tem_y=0.;
  tem_dx=0.;
  tem_dy=0.;
  
  for(double t=dt; t<=1.; t+=dt){
    
   K(k1,k2,k3,k4,k5,k6,k7,x,y,dx,dy,dt,mu);
      
   tem_x = x + dt*(35./384.)*k1[0] + dt*(500./1113.)*k3[0] + dt*(125./192.)*k4[0] - dt*(2187./6784.)*k5[0] + dt*(11./84.)*k6[0];
   tem_y = y + dt*(35./384.)*k1[1] + dt*(500./1113.)*k3[1] + dt*(125./192.)*k4[1] - dt*(2187./6784.)*k5[1] + dt*(11./84.)*k6[1];
   tem_dx = dx + dt*(35./384.)*k1[2] + dt*(500./1113.)*k3[2] + dt*(125./192.)*k4[2] - dt*(2187./6784.)*k5[2] + dt*(11./84.)*k6[2];
   tem_dy = dy + dt*(35./384.)*k1[3] + dt*(500./1113.)*k3[3] + dt*(125./192.)*k4[3] - dt*(2187./6784.)*k5[3] + dt*(11./84.)*k6[3]; 
      
   x += dt*(5179./57600.)*k1[0] + dt*(7571./16695.)*k3[0] + dt*(393./640.)*k4[0] - dt*(92097./339200.)*k5[0] + dt*(187./2100.)*k6[0] + dt*(1./40.)*k7[0];
   y += dt*(5179./57600.)*k1[1] + dt*(7571./16695.)*k3[1] + dt*(393./640.)*k4[1] - dt*(92097./339200.)*k5[1] + dt*(187./2100.)*k6[1] + dt*(1./40.)*k7[1];
   dx += dt*(5179./57600.)*k1[2] + dt*(7571./16695.)*k3[2] + dt*(393./640.)*k4[2] - dt*(92097./339200.)*k5[2] + dt*(187./2100.)*k6[2] + dt*(1./40.)*k7[2];
   dy += dt*(5179./57600.)*k1[3] + dt*(7571./16695.)*k3[3] + dt*(393./640.)*k4[3] - dt*(92097./339200.)*k5[3] + dt*(187./2100.)*k6[3] + dt*(1./40.)*k7[3]; 
   
   max(var,abs(x - tem_x),abs(y - tem_y),abs(dx - tem_dx),abs(dy - tem_dy));
   
   dt=0.1*dt*pow((tol/var),1./5.);
   
   cout << t << "\t" << dt << "\t" << x << "\t" << y << "\t" << dx << "\t" << dy << endl;
  }
 
  return 0;
}