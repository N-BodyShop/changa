#include "cell_typedef.h"

void nodeBucketForce(CellMultipoleMoments *mm, int numMM, CellGravityParticle *part, int numPart) {
  CellVector3D r;
  cellSPEtype rsq;
  cellSPEtype twoh, a, b, c, d;
  cellSPEtype dir, qirx, qiry, qirz, qir, tr, qir3;
  cellSPEtype u,dih;
  cellSPEtype idt2;

  for (i=0; i<numPart; ++i) {
    for (j=0; j<numMM; ++j) {
      r.x = part[i].core.position.x - mm[j].cm.x;
      r.y = part[i].core.position.y - mm[j].cm.y;
      r.z = part[i].core.position.z - mm[j].cm.z;
      rsq = r.x * r.x + r.y * r.y + r.z * r.z;
      twoh = mm[j].soft + part[i].core.soft;
      dir = 1.0/sqrt(rsq);

      if ((rsq) < (twoh)*(twoh)) {
        dih = 2.0/(twoh);
        u = dih/dir;
        if (u < 1.0) {
          a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
              - 1.0/10.0*u*u*u*u*u);
          b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
          c = dih*dih*dih*dih*dih*(12.0/5.0 - 3.0/2.0*u);
          d = 3.0/2.0*dih*dih*dih*dih*dih*dih*dir;
        }
        else {
          a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u + u*u*u
              - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
          b = -1.0/15.0*dir*dir*dir + dih*dih*dih*(8.0/3.0 - 3.0*u + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
          c = -1.0/5.0*dir*dir*dir*dir*dir + 3.0*dih*dih*dih*dih*dir
          + dih*dih*dih*dih*dih*(-12.0/5.0 + 1.0/2.0*u);
          d = -dir*dir*dir*dir*dir*dir*dir
          + 3.0*dih*dih*dih*dih*dir*dir*dir
          - 1.0/2.0*dih*dih*dih*dih*dih*dih*dir;
        }
      }
      else {
        a = dir;
        b = a*a*a;
        c = 3.0*b*a*a;
        d = 5.0*c*a*a;
      }

      qirx = mm[j].xx*r.x + mm[j].xy*r.y + mm[j].xz*r.z;
      qiry = mm[j].xy*r.x + mm[j].yy*r.y + mm[j].yz*r.z;
      qirz = mm[j].xz*r.x + mm[j].yz*r.y + mm[j].zz*r.z;
      qir = 0.5*(qirx*r.x + qiry*r.y + qirz*r.z);
      tr = 0.5*(mm[j].xx + mm[j].yy + mm[j].zz);
      qir3 = b*mm[j].totalMass + d*qir - c*tr;
      part[i].potential -= mm[j].totalMass * a + c*qir - b*tr;
      part[i].treeAcceleration.x -= qir3*r.x - c*qirx;
      part[i].treeAcceleration.y -= qir3*r.y - c*qiry;
      part[i].treeAcceleration.z -= qir3*r.z - c*qirz;
      idt2 = (part[i].core.mass + mm[j].totalMass)*b;
      if(idt2 > part[i].dtGrav)
        part[i].dtGrav = idt2;
    }
  }
}

void particleBucketForce(CellExternalGravityParticle *ext, int numExt, CellGravityParticle *part, int numPart) {
  CellVector3D r;
  cellSPEtype rsq;
  cellSPEtype twoh, a, b;
  cellSPEtype idt2;
  int i,j;
  
  for (i=0; i<numPart; ++i) {
    for (j=0; j<numExt; ++j) {
      r.x = ext[j].position.x - part[i].core.position.x;
      r.y = ext[j].position.y - part[i].core.position.y;
      r.z = ext[j].position.z - part[i].core.position.z;
      rsq = r.x * r.x + r.y * r.y + r.z * r.z;
      twoh = ext[j].soft + part[i].core.soft;
      if(rsq != 0) {

        cellSPEtype r2, u,dih,dir;
        r2 = sqrt(rsq);
        if (r2 < (twoh)) {
          dih = 2.0/(twoh);
          u = r2*dih;
          if (u < 1.0) {
            a = dih*(7.0/5.0 - 2.0/3.0*u*u + 3.0/10.0*u*u*u*u
                - 1.0/10.0*u*u*u*u*u);
            b = dih*dih*dih*(4.0/3.0 - 6.0/5.0*u*u + 1.0/2.0*u*u*u);
          }
          else {
            dir = 1.0/r2;
            a = -1.0/15.0*dir + dih*(8.0/5.0 - 4.0/3.0*u*u + u*u*u
                - 3.0/10.0*u*u*u*u + 1.0/30.0*u*u*u*u*u);
            b = -1.0/15.0*dir*dir*dir + dih*dih*dih*(8.0/3.0 - 3.0*u + 6.0/5.0*u*u - 1.0/6.0*u*u*u);
          }
        }
        else {
          a = 1.0/r2;
          b = a*a*a;
        }
        
        idt2 = (part[i].core.mass + ext[j].mass)*b; // (timescale)^-2
        // of interaction
        part[i].treeAcceleration.x += r.x * (b * ext[j].mass);
        part[i].treeAcceleration.y += r.y * (b * ext[j].mass);
        part[i].treeAcceleration.z += r.z * (b * ext[j].mass);
        part[i].potential -= ext[j].mass * a;
        if(idt2 > part[i].dtGrav)
          part[i].dtGrav = idt2;
      }
    }
  }
}

void funcLookup(int funcIndex, void* readWritePtr, int readWriteLen, void* readOnlyPtr, int readOnlyLen,
                void* writeOnlyPtr, int writeOnlyLen, DMAListEntry* dmaList ) {
  if (funcIndex==1) {
    // node interaction
    nodeBucketForce(readOnlyPtr, readOnlyLen/sizeof(CellMultipoleMoments),
                    readWritePtr, readWriteLen/sizeof(CellGravityParticle));
  } else if (funcIndex==2) {
    // particle interaction
    particleBucketForce(readOnlyPtr, readOnlyLen/sizeof(CellExternalGravityParticle),
                        readWritePtr, readWriteLen/sizeof(CellGravityParticle));
  }
}