#include "UserGravity.h"
#include "UserRadGravity.h"
#include "UserZGravity.h"
#include "UserPWGravity.h"
#include "TurbulentDrivingSource.h"
#include "physconst.h"
//include "parameters.h"
#include <stdio.h>
#include <math.h>

UserGravity* UserGravity::instance = 0;

UserGravity* UserGravity::gravInstances[MAX_INSTANCES] = {0};
int UserGravity::numInstances = 0;

double UserGravity::secUnit = 1;
double UserGravity::cmUnit = 1;
double  UserGravity::xPeriod = 0;
double  UserGravity::yPeriod = 0;
double  UserGravity::zPeriod = 0;

int UserGravity::type = 0;


void UserGravity::initUserGravityParameters( Parameters param) {
  UserGravity::gravInstances[ UserGravity::numInstances++] = new UserRadGravity();
  UserGravity::gravInstances[ UserGravity::numInstances++] = new TurbulentDrivingSource();
  UserGravity::gravInstances[ UserGravity::numInstances++] = new UserZGravity();
  UserGravity::gravInstances[ UserGravity::numInstances++] = new UserPWGravity();

  for( int i = 0; i < UserGravity::numInstances; i++) { 
    UserGravity::gravInstances[i] ->initializeParameters(param);
  }
}

void UserGravity::initUserGravity() {
  if( instance != 0 ) instance->initialize();
}

void UserGravity::initConstants( int gType, double sec_unit, double cm_unit, double xp, double yp, double zp) {
  UserGravity::secUnit = sec_unit;
  UserGravity::cmUnit = cm_unit;
  UserGravity::xPeriod = xp;
  UserGravity::yPeriod = yp;
  UserGravity::zPeriod = zp;
  UserGravity::type = gType;

  for( int i = 0; i < UserGravity::numInstances; i++) {
    if( gType == UserGravity::gravInstances[i]->gravType()) {
      UserGravity::instance = UserGravity::gravInstances[i];
      break;
    }
  }
}


void UserRadGravity::initializeParameters( Parameters param){
  strcpy( filename, param.rGravFile);
}


void UserRadGravity::initialize( ) {
  FILE *f;
  double radius, mass;
  f = fopen(filename,"r");
  while( fscanf(f, "%le %le", &radius, &mass) != EOF) {
    r.push_back(radius); m.push_back(mass); g.push_back(radius > 0. ? mass/(radius*radius): 0.);
  }  
}

void UserRadGravity::getGravity( double t, double x, double y, double z, double& gx, double& gy, double& gz, double &dtg) {
  double rad = sqrt(x*x + y*y + z*z);
  int rsize = r.size();
  for( int i = 0; i < rsize-1; i++) { 
    if( rad >= r[i] && rad < r[i+1]) {
      double dr = r[i+1] - r[i]; 
      double delta = rad - r[i];
      double dg = g[i+1] - g[i]; 
      double grav = g[i] + dg/dr*delta;
      gx = -grav*x/rad;
      gy = -grav*y/rad;
      gz = -grav*z/rad;
      dtg = sqrt(rad/sqrt(gx*gx+gy*gy+gz*gz));
      return;
    } 
  } 
  // treat as a point mass outside
  double grav = m[rsize-1]/(rad*rad);  
  gx = -grav*x/rad;
  gy = -grav*y/rad;
  gz = -grav*z/rad;
}

void UserZGravity::initializeParameters( Parameters param){
  gravLengthScale = param.dgravL;
  zGravity = param.dgrav;
}

void UserZGravity::getGravity( double t, double x, double y, double z, double& gx, double& gy, double& gz, double &dtg) {
  double znorm = 1.;
  gx = 0.; gy = 0.; gz = 0.;

  double w, *g;
  if( direction == 0) {
    w=x; g=&gx;
  }
  if( direction == 1) {
    w=y; g=&gy;
  }
  if( direction == 2) {
    w=z; g=&gz;
  }

  if( gravLengthScale > 0.) {
    znorm = fabs(w/gravLengthScale);
  }
  
  double sgn = w>0. ? 1. : -1.;
  double norm = znorm*znorm;
  
  if( bGravityFixedAtOne) norm = fmin(norm, 1.);

  *g = -zGravity*norm*sgn;
  dtg = sqrt(fabs(z/ *g));

  return;
}

void UserPWGravity::initializeParameters( Parameters param){
  Mbh = param.dMbh;
  rsoft = param.drbhSoft;
  xbh = param.dxbh; ybh = param.dybh; zbh = param.dzbh;
  vxbh = param.dvxbh; vybh = param.dvybh; vzbh = param.dvzbh;
  cspeed = LIGHTSPEED*param.dSecUnit/(3.0857e21*param.dKpcUnit);
}

void UserPWGravity::getGravity( double t, double x, double y, double z, double& gx, double& gy, double& gz, double &dtg) {
  double dx = x - (xbh + vxbh*t);
  double dy = y - (ybh + vybh*t);
  double dz = z - (zbh + vzbh*t);
  double r = fmax(sqrt(dx*dx + dy*dy + dz*dz),1e-10);
  double rprime = r;
  double rg = 0.; //2.*Mbh/(cspeed*cspeed);
  double rMinusRg = fmax(rprime-rg,0.);
  double g = -Mbh/(rMinusRg*rMinusRg + rsoft*rsoft);
  gx = g*dx/r; gy = g*dy/r; gz = g*dz/r;
  dtg = sqrt(r/sqrt(gx*gx+gy*gy+gz*gz));

  return;
}
