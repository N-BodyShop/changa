#include "cell_typedef.h"
#include "stdio.h"
#include "math.h"
#include "spu_intrinsics.h"
#include "simdmath.h"

/* Original version */
/*
void nodeBucketForce(CellMultipoleMoments *mm, int numMM, CellGravityParticle *part, int numPart) {
  CellVector4D r;
  cellSPEtype rsq;
  cellSPEtype twoh, a, b, c, d;
  cellSPEtype dir, qirx, qiry, qirz, qir, tr, qir3;
  cellSPEtype u,dih;
  cellSPEtype idt2;
  int i, j;

  for (i=0; i<numPart; ++i) {
    CellGravityParticle *particle = &part[i];
    vec_float4 *particlePosition_vec = (vec_float4*)&particle->core.position;
    for (j=0; j<numMM; ++j) {
      //printf("extnode %d: %f %f %f (%f %f %f) %f %f %f %f %f %f\n",j,mm[j].radius,mm[j].soft,mm[j].totalMass,mm[j].cm.x,mm[j].cm.y,mm[j].cm.z,mm[j].xx,mm[j].xy,mm[j].xz,mm[j].yy,mm[j].yz,mm[j].zz);
      CellMultipoleMoments *moments = &mm[j];
      CellVector3D *momentsCM = &moments->cm;
      vec_float4 *momentsCM_vec = (vec_float4*)momentsCM;
      vec_float4 r_vec = *particlePosition_vec - *momentsCM_vec;
      r.x = spu_extract(r_vec, 0);
      r.y = spu_extract(r_vec, 1);
      r.z = spu_extract(r_vec, 2);
      //r.x = particle->core.position.x - momentsCM->x;
      //r.y = particle->core.position.y - momentsCM->y;
      //r.z = particle->core.position.z - momentsCM->z;
      rsq = r.x * r.x + r.y * r.y + r.z * r.z;
      twoh = moments->soft + particle->core.soft;
      //dir = 1.0f/sqrt(rsq);
      //__asm__ ("FRSQEST dir, rsq");
      //__asm__ ("frsqest %0, %1" : "=r"(dir) : "r"(rsq));
      qword vec_rsq = si_from_float((cellSPEtype)rsq);
      vec_float4 vec_dir = spu_rsqrte((vec_float4)vec_rsq);
      dir = (cellSPEtype)si_to_float((qword)vec_dir);

      if ((rsq) < (twoh)*(twoh)) {
	printf("never occuring case\n");
        dih = 2.0f/(twoh);
        u = dih/dir;
        if (u < 1.0f) {
	  printf("case 1-1\n");
          a = dih*(7.0f/5.0f - 2.0f/3.0f*u*u + 3.0f/10.0f*u*u*u*u
              - 1.0f/10.0f*u*u*u*u*u);
          b = dih*dih*dih*(4.0f/3.0f - 6.0f/5.0f*u*u + 1.0f/2.0f*u*u*u);
          c = dih*dih*dih*dih*dih*(12.0f/5.0f - 3.0f/2.0f*u);
          d = 3.0f/2.0f*dih*dih*dih*dih*dih*dih*dir;
        }
        else {
	  printf("case 1-2\n");
          a = -1.0f/15.0f*dir + dih*(8.0f/5.0f - 4.0f/3.0f*u*u + u*u*u
              - 3.0f/10.0f*u*u*u*u + 1.0f/30.0f*u*u*u*u*u);
          b = -1.0f/15.0f*dir*dir*dir + dih*dih*dih*(8.0f/3.0f - 3.0f*u + 6.0f/5.0f*u*u - 1.0f/6.0f*u*u*u);
          c = -1.0f/5.0f*dir*dir*dir*dir*dir + 3.0f*dih*dih*dih*dih*dir
          + dih*dih*dih*dih*dih*(-12.0f/5.0f + 1.0f/2.0f*u);
          d = -dir*dir*dir*dir*dir*dir*dir
          + 3.0f*dih*dih*dih*dih*dir*dir*dir
          - 1.0f/2.0f*dih*dih*dih*dih*dih*dih*dir;
        }
      }
      else {
	//printf("case 2\n");
        a = dir;
        b = a*a*a;
        c = 3.0f*b*a*a;
        d = 5.0f*c*a*a;
      }

      qirx = moments->xx*r.x + moments->xy*r.y + moments->xz*r.z;
      qiry = moments->xy*r.x + moments->yy*r.y + moments->yz*r.z;
      qirz = moments->xz*r.x + moments->yz*r.y + moments->zz*r.z;
      qir = 0.5f*(qirx*r.x + qiry*r.y + qirz*r.z);
      tr = 0.5f*(moments->xx + moments->yy + moments->zz);
      qir3 = b*moments->totalMass + d*qir - c*tr;
      particle->potential -= moments->totalMass * a + c*qir - b*tr;
      particle->treeAcceleration.x -= qir3*r.x - c*qirx;
      particle->treeAcceleration.y -= qir3*r.y - c*qiry;
      particle->treeAcceleration.z -= qir3*r.z - c*qirz;
      idt2 = (particle->core.mass + moments->totalMass)*b;
      //if(idt2 > particle->dtGrav)
      //  particle->dtGrav = idt2;
      particle->dtGrav = (cellSPEtype)si_to_float((qword)fmaxf4((vec_float4)si_from_float((float)idt2), (vec_float4)si_from_float((float)(particle->dtGrav))));
      //printf("intpart %d: %f (%f %f %f) %f / (%f %f %f) %f %f\n",i,particle->core.mass,particle->core.position.x,particle->core.position.y,particle->core.position.z,particle->core.soft,particle->treeAcceleration.x,particle->treeAcceleration.y,particle->treeAcceleration.z,particle->potential,particle->dtGrav);
    }
  }
}
*/
/* Version vectorized */
void nodeBucketForce(CellMultipoleMoments *mm, int numMM, CellGravityParticle *part, int numPart) {
  //CellVector4D r;
  vec_float4 rx, ry, rz;
  vec_float4 rsq;
  vec_float4 twoh, a, b, c, d;
  vec_float4 a1, a2, b1, b2, c1, c2, d1, d2;
  vec_float4 dir, qirx, qiry, qirz, qir, qir3;
  float tr;
  vec_float4 u,dih;
  vec_float4 idt2;
  int i, j, k;
  vec_float4 partPosX, partPosY, partPosZ, partSoft, partMass, partAccX, partAccY, partAccZ, partDtGrav, partPot;

  for (i=0; i<numPart; i+=4) {
    for (k=0; k<4; k++) {
      CellExternalGravityParticle *particle = &part[i+k].core;
      partPosX = spu_insert(particle->position.x, partPosX, k);
      partPosY = spu_insert(particle->position.y, partPosY, k);
      partPosZ = spu_insert(particle->position.z, partPosZ, k);
      partSoft = spu_insert(particle->soft, partSoft, k);
      partMass = spu_insert(particle->mass, partMass, k);
    }
    partAccX = partAccY = partAccZ = partDtGrav = partPot = spu_splats(0.0f);
    //CellGravityParticle *particle = &part[i];
    //vec_float4 *particlePosition_vec = (vec_float4*)&particle->core.position;
    for (j=0; j<numMM; ++j) {
      //printf("extnode %d: %f %f %f (%f %f %f) %f %f %f %f %f %f\n",j,mm[j].radius,mm[j].soft,mm[j].totalMass,mm[j].cm.x,mm[j].cm.y,mm[j].cm.z,mm[j].xx,mm[j].xy,mm[j].xz,mm[j].yy,mm[j].yz,mm[j].zz);
      CellMultipoleMoments *moments = &mm[j];
      CellVector3D *momentsCM = &moments->cm;
      //vec_float4 *momentsCM_vec = (vec_float4*)momentsCM;
      //vec_float4 r_vec = *particlePosition_vec - *momentsCM_vec;
      //r.x = spu_extract(r_vec, 0);
      //r.y = spu_extract(r_vec, 1);
      //r.z = spu_extract(r_vec, 2);
      rx = partPosX - spu_splats(momentsCM->x);
      ry = partPosY - spu_splats(momentsCM->y);
      rz = partPosZ - spu_splats(momentsCM->z);
      rsq = rx * rx + ry * ry + rz * rz;
      twoh = spu_splats(moments->soft) + partSoft;
      //dir = 1.0f/sqrt(rsq);
      //__asm__ ("FRSQEST dir, rsq");
      //__asm__ ("frsqest %0, %1" : "=r"(dir) : "r"(rsq));
      //qword vec_rsq = si_from_float((cellSPEtype)rsq);
      //vec_float4 vec_dir = spu_rsqrte((vec_float4)vec_rsq);
      //dir = (cellSPEtype)si_to_float((qword)vec_dir);
      dir = spu_rsqrte(rsq);

      // comparison has 0xFFFFFFFF if rsq<twoh*twoh, 0x00000000 otherwise
      vec_uint4 comparison = spu_cmpgt(twoh*twoh, rsq);
      // the following vector has all its element equal, and can be used in the following if
      unsigned int comparator = spu_extract(spu_gather(comparison), 0);
      if (comparator) {
	//printf("never occuring case\n");
        dih = spu_splats(2.0f)/(twoh);
        u = dih/dir;
	vec_uint4 comparison2 = spu_cmpgt(spu_splats(1.0f), u);
	unsigned int comparator2 = spu_extract(spu_gather(comparison2), 0);
        if (comparator2) {
	  //printf("case 1-1\n");
          a = dih*(spu_splats(7.0f/5.0f) - spu_splats(2.0f/3.0f)*u*u
		   + spu_splats(3.0f/10.0f)*u*u*u*u
		   - spu_splats(1.0f/10.0f)*u*u*u*u*u);
          b = dih*dih*dih*(spu_splats(4.0f/3.0f) - spu_splats(6.0f/5.0f)*u*u
			   + spu_splats(1.0f/2.0f)*u*u*u);
          c = dih*dih*dih*dih*dih*(spu_splats(12.0f/5.0f) - spu_splats(3.0f/2.0f)*u);
          d = spu_splats(3.0f/2.0f)*dih*dih*dih*dih*dih*dih*dir;
        }
        if ((~comparator2) & 0xF) {
	  //printf("case 1-2\n");
          a2 = spu_splats(-1.0f/15.0f)*dir + dih*(spu_splats(8.0f/5.0f)
	      - spu_splats(4.0f/3.0f)*u*u + u*u*u
	      - spu_splats(3.0f/10.0f)*u*u*u*u + spu_splats(1.0f/30.0f)*u*u*u*u*u);
          b2 = spu_splats(-1.0f/15.0f)*dir*dir*dir + dih*dih*dih*(spu_splats(8.0f/3.0f) - spu_splats(3.0f)*u + spu_splats(6.0f/5.0f)*u*u - spu_splats(1.0f/6.0f)*u*u*u);
          c2 = spu_splats(-1.0f/5.0f)*dir*dir*dir*dir*dir + spu_splats(3.0f)*dih*dih*dih*dih*dir
	    + dih*dih*dih*dih*dih*(spu_splats(-12.0f/5.0f) + spu_splats(1.0f/2.0f)*u);
          d2 = -dir*dir*dir*dir*dir*dir*dir
	    + spu_splats(3.0f)*dih*dih*dih*dih*dir*dir*dir
	    - spu_splats(1.0f/2.0f)*dih*dih*dih*dih*dih*dih*dir;
        }
	a1 = spu_sel(a2, a, comparison2);
	b1 = spu_sel(b2, b, comparison2);
	c1 = spu_sel(c2, c, comparison2);
	d1 = spu_sel(d2, d, comparison2);
      }
      if ((~comparator) & 0xF) {
	//printf("case 2\n");
        a2 = dir;
        b2 = a2*a2*a2;
        c2 = spu_splats(3.0f)*b2*a2*a2;
        d2 = spu_splats(5.0f)*c2*a2*a2;
      }
      a = spu_sel(a2, a1, comparison);
      b = spu_sel(b2, b1, comparison);
      c = spu_sel(c2, c1, comparison);
      d = spu_sel(d2, d1, comparison);

      qirx = spu_splats(moments->xx)*rx + spu_splats(moments->xy)*ry + spu_splats(moments->xz)*rz;
      qiry = spu_splats(moments->xy)*rx + spu_splats(moments->yy)*ry + spu_splats(moments->yz)*rz;
      qirz = spu_splats(moments->xz)*rx + spu_splats(moments->yz)*ry + spu_splats(moments->zz)*rz;
      qir = spu_splats(0.5f)*(qirx*rx + qiry*ry + qirz*rz);
      tr = 0.5f*(moments->xx + moments->yy + moments->zz);
      qir3 = b*spu_splats(moments->totalMass) + d*qir - c*spu_splats(tr);
      partPot -= spu_splats(moments->totalMass) * a + c*qir - b*spu_splats(tr);
      partAccX -= qir3*rx - c*qirx;
      partAccY -= qir3*ry - c*qiry;
      partAccZ -= qir3*rz - c*qirz;
      idt2 = (partMass + spu_splats(moments->totalMass))*b;
      //if(idt2 > particle->dtGrav)
      //  particle->dtGrav = idt2;
      //particle->dtGrav = (cellSPEtype)si_to_float((qword)fmaxf4((vec_float4)si_from_float((float)idt2), (vec_float4)si_from_float((float)(particle->dtGrav))));
      partDtGrav = fmaxf4(idt2, partDtGrav);
      //printf("intpart %d: %f (%f %f %f) %f / (%f %f %f) %f %f\n",i,particle->core.mass,particle->core.position.x,particle->core.position.y,particle->core.position.z,particle->core.soft,particle->treeAcceleration.x,particle->treeAcceleration.y,particle->treeAcceleration.z,particle->potential,particle->dtGrav);
    }
    for (k=0; k<4; ++k) {
      CellGravityParticle *particle = &part[i+k];
      particle->potential = spu_extract(partPot, k);
      particle->treeAcceleration.x = spu_extract(partAccX, k);
      particle->treeAcceleration.y = spu_extract(partAccY, k);
      particle->treeAcceleration.z = spu_extract(partAccZ, k);
      particle->dtGrav = spu_extract(partDtGrav, k);
    }
  }
}
/* Version manually unrolled */
/*
void nodeBucketForce(CellMultipoleMoments *mm, int numMM, CellGravityParticle *part, int numPart) {
  CellVector4D r1, r2;
  cellSPEtype rsq1, rsq2;
  cellSPEtype twoh1, twoh2, a1, a2, b1, b2, c1, c2, d1, d2;
  cellSPEtype dir1, dir2, qirx1, qirx2, qiry1, qiry2, qirz1, qirz2, qir1, qir2, tr, qir31, qir32;
  cellSPEtype u1,u2,dih1,dih2;
  cellSPEtype idt21,idt22;
  int i, j;

  for (i=0; i<numPart; i+=2) {
    CellGravityParticle *particle1 = &part[i];
    vec_float4 *particlePosition_vec1 = (vec_float4*)&particle1->core.position;
    CellGravityParticle *particle2 = &part[i+1];
    vec_float4 *particlePosition_vec2 = (vec_float4*)&particle2->core.position;
    for (j=0; j<numMM; ++j) {
      //printf("extnode %d: %f %f %f (%f %f %f) %f %f %f %f %f %f\n",j,mm[j].radius,mm[j].soft,mm[j].totalMass,mm[j].cm.x,mm[j].cm.y,mm[j].cm.z,mm[j].xx,mm[j].xy,mm[j].xz,mm[j].yy,mm[j].yz,mm[j].zz);
      CellMultipoleMoments *moments = &mm[j];
      CellVector3D *momentsCM = &moments->cm;
      vec_float4 *momentsCM_vec = (vec_float4*)momentsCM;
      vec_float4 r_vec1 = *particlePosition_vec1 - *momentsCM_vec;
      r1.x = spu_extract(r_vec1, 0);
      r1.y = spu_extract(r_vec1, 1);
      r1.z = spu_extract(r_vec1, 2);
      vec_float4 r_vec2 = *particlePosition_vec2 - *momentsCM_vec;
      r2.x = spu_extract(r_vec2, 0);
      r2.y = spu_extract(r_vec2, 1);
      r2.z = spu_extract(r_vec2, 2);
      //r.x = particle->core.position.x - momentsCM->x;
      //r.y = particle->core.position.y - momentsCM->y;
      //r.z = particle->core.position.z - momentsCM->z;
      rsq1 = r1.x * r1.x + r1.y * r1.y + r1.z * r1.z;
      rsq2 = r2.x * r2.x + r2.y * r2.y + r2.z * r2.z;
      twoh1 = moments->soft + particle1->core.soft;
      twoh2 = moments->soft + particle2->core.soft;
      //dir = 1.0f/sqrt(rsq);
      //__asm__ ("FRSQEST dir, rsq");
      __asm__ ("frsqest %0, %1" : "=r"(dir1) : "r"(rsq1));
      __asm__ ("frsqest %0, %1" : "=r"(dir2) : "r"(rsq2));
      //qword vec_rsq = si_from_float((cellSPEtype)rsq);
      //vec_float4 vec_dir = spu_rsqrte((vec_float4)vec_rsq);
      //dir = (cellSPEtype)si_to_float((qword)vec_dir);

      if ((rsq) < (twoh)*(twoh)) {
	printf("never occuring case\n");
        dih = 2.0f/(twoh);
        u = dih/dir;
        if (u < 1.0f) {
	  printf("case 1-1\n");
          a = dih*(7.0f/5.0f - 2.0f/3.0f*u*u + 3.0f/10.0f*u*u*u*u
              - 1.0f/10.0f*u*u*u*u*u);
          b = dih*dih*dih*(4.0f/3.0f - 6.0f/5.0f*u*u + 1.0f/2.0f*u*u*u);
          c = dih*dih*dih*dih*dih*(12.0f/5.0f - 3.0f/2.0f*u);
          d = 3.0f/2.0f*dih*dih*dih*dih*dih*dih*dir;
        }
        else {
	  printf("case 1-2\n");
          a = -1.0f/15.0f*dir + dih*(8.0f/5.0f - 4.0f/3.0f*u*u + u*u*u
              - 3.0f/10.0f*u*u*u*u + 1.0f/30.0f*u*u*u*u*u);
          b = -1.0f/15.0f*dir*dir*dir + dih*dih*dih*(8.0f/3.0f - 3.0f*u + 6.0f/5.0f*u*u - 1.0f/6.0f*u*u*u);
          c = -1.0f/5.0f*dir*dir*dir*dir*dir + 3.0f*dih*dih*dih*dih*dir
          + dih*dih*dih*dih*dih*(-12.0f/5.0f + 1.0f/2.0f*u);
          d = -dir*dir*dir*dir*dir*dir*dir
          + 3.0f*dih*dih*dih*dih*dir*dir*dir
          - 1.0f/2.0f*dih*dih*dih*dih*dih*dih*dir;
        }
      }
      else {
	//printf("case 2\n");
        a1 = dir1;
        b1 = a1*a1*a1;
        c1 = 3.0f*b1*a1*a1;
        d1 = 5.0f*c1*a1*a1;
        a2 = dir2;
        b2 = a2*a2*a2;
        c2 = 3.0f*b2*a2*a2;
        d2 = 5.0f*c2*a2*a2;
      }

      qirx1 = moments->xx*r1.x + moments->xy*r1.y + moments->xz*r1.z;
      qiry1 = moments->xy*r1.x + moments->yy*r1.y + moments->yz*r1.z;
      qirz1 = moments->xz*r1.x + moments->yz*r1.y + moments->zz*r1.z;
      qir1 = 0.5f*(qirx1*r1.x + qiry1*r1.y + qirz1*r1.z);
      tr = 0.5f*(moments->xx + moments->yy + moments->zz);
      qir31 = b1*moments->totalMass + d1*qir1 - c1*tr;
      qirx2 = moments->xx*r2.x + moments->xy*r2.y + moments->xz*r2.z;
      qiry2 = moments->xy*r2.x + moments->yy*r2.y + moments->yz*r2.z;
      qirz2 = moments->xz*r2.x + moments->yz*r2.y + moments->zz*r2.z;
      qir2 = 0.5f*(qirx2*r2.x + qiry2*r2.y + qirz2*r2.z);
      qir32 = b2*moments->totalMass + d2*qir2 - c2*tr;
      particle1->potential -= moments->totalMass * a1 + c1*qir1 - b1*tr;
      particle1->treeAcceleration.x -= qir31*r1.x - c1*qirx1;
      particle1->treeAcceleration.y -= qir31*r1.y - c1*qiry1;
      particle1->treeAcceleration.z -= qir31*r1.z - c1*qirz1;
      idt21 = (particle1->core.mass + moments->totalMass)*b1;
      particle2->potential -= moments->totalMass * a2 + c2*qir2 - b2*tr;
      particle2->treeAcceleration.x -= qir32*r2.x - c2*qirx2;
      particle2->treeAcceleration.y -= qir32*r2.y - c2*qiry2;
      particle2->treeAcceleration.z -= qir32*r2.z - c2*qirz2;
      idt22 = (particle2->core.mass + moments->totalMass)*b2;
      //if(idt2 > particle->dtGrav)
      //  particle->dtGrav = idt2;
      particle1->dtGrav = (cellSPEtype)si_to_float((qword)fmaxf4((vec_float4)si_from_float((float)idt21), (vec_float4)si_from_float((float)(particle1->dtGrav))));
      particle1->dtGrav = (cellSPEtype)si_to_float((qword)fmaxf4((vec_float4)si_from_float((float)idt22), (vec_float4)si_from_float((float)(particle2->dtGrav))));
      //printf("intpart %d: %f (%f %f %f) %f / (%f %f %f) %f %f\n",i,particle->core.mass,particle->core.position.x,particle->core.position.y,particle->core.position.z,particle->core.soft,particle->treeAcceleration.x,particle->treeAcceleration.y,particle->treeAcceleration.z,particle->potential,particle->dtGrav);
    }
  }
}
*/

/* Vectorized version */
void particleBucketForce(CellExternalGravityParticle *ext, int numExt, CellGravityParticle *part, int numPart) {
  //CellVector3D r;
  vec_float4 rx, ry, rz;
  vec_float4 rsq;
  vec_float4 twoh, a, b, a1, a2, b1, b2;
  vec_float4 idt2;
  int i,j,k;
  vec_float4 partPosX, partPosY, partPosZ, partSoft, partMass, partAccX, partAccY, partAccZ, partDtGrav, partPot;
  
  for (i=0; i<numPart; ++i) {
    for (k=0; k<4; k++) {
      CellExternalGravityParticle *particle = &part[i+k].core;
      partPosX = spu_insert(particle->position.x, partPosX, k);
      partPosY = spu_insert(particle->position.y, partPosY, k);
      partPosZ = spu_insert(particle->position.z, partPosZ, k);
      partSoft = spu_insert(particle->soft, partSoft, k);
      partMass = spu_insert(particle->mass, partMass, k);
    }
    partAccX = partAccY = partAccZ = partDtGrav = partPot = spu_splats(0.0f);
    for (j=0; j<numExt; ++j) {
      CellExternalGravityParticle *extPart = &ext[j];
      vec_float4 extPartMass = spu_splats(extPart->mass);
      //printf("extpart %d: %f (%f %f %f) %f\n",j,ext[j].mass,ext[j].position.x,ext[j].position.y,ext[j].position.z,ext[j].soft);
      rx = spu_splats(extPart->position.x) - partPosX;
      ry = spu_splats(extPart->position.y) - partPosY;
      rz = spu_splats(extPart->position.z) - partPosZ;
      rsq = rx * rx + ry * ry + rz * rz;
      twoh = spu_splats(extPart->soft) + partSoft;
      //if(rsq != 0) {

        vec_float4 r2_inv, u,dih,dir,r2twoh;
	//qword vec_rsq = si_from_float((cellSPEtype)rsq);
	//vec_float4 vec_r2 = spu_rsqrte((vec_float4)vec_rsq);
	//r2_inv = (cellSPEtype)si_to_float((qword)vec_r2);
	r2_inv = spu_rsqrte(rsq);
        //r2 = sqrt(rsq);
	r2twoh = r2_inv * twoh;
	vec_uint4 comparison = spu_cmpgt(r2twoh, spu_splats(1.0f));
	unsigned int comparator = spu_extract(spu_gather(comparison), 0);
        if (comparator) {
          dih = spu_splats(2.0f)/twoh;
          u = spu_splats(2.0f)/r2twoh;
	  vec_uint4 comparison2 = spu_cmpgt(spu_splats(1.0f), u);
	  unsigned int comparator2 = spu_extract(spu_gather(comparison2), 0);
          if (comparator2) {
            a = dih*(spu_splats(7.0f/5.0f) - spu_splats(2.0f/3.0f)*u*u + spu_splats(3.0f/10.0f)*u*u*u*u
		     - spu_splats(1.0f/10.0f)*u*u*u*u*u);
            b = dih*dih*dih*(spu_splats(4.0f/3.0f) - spu_splats(6.0f/5.0f)*u*u + spu_splats(1.0f/2.0f)*u*u*u);
	  }
	  if ((~comparator2) & 0xF) {
            dir = r2_inv;
            a2 = spu_splats(-1.0f/15.0f)*dir + dih*(spu_splats(8.0f/5.0f) - spu_splats(4.0f/3.0f)*u*u + u*u*u
                - spu_splats(3.0f/10.0f)*u*u*u*u + spu_splats(1.0f/30.0f)*u*u*u*u*u);
            b2 = spu_splats(-1.0f/15.0f)*dir*dir*dir + dih*dih*dih*(spu_splats(8.0f/3.0f) - spu_splats(3.0f)*u + spu_splats(6.0f/5.0f)*u*u - spu_splats(1.0f/6.0f)*u*u*u);
          }
	  a1 = spu_sel(a2, a, comparison2);
	  b1 = spu_sel(b2, b, comparison2);
	}
	if ((~comparator) & 0xF) {
          a2 = r2_inv;
          b2 = a2*a2*a2;
        }
	a = spu_sel(a2, a1, comparison);
	b = spu_sel(b2, b1, comparison);
        
	// if rsq is zero for a particle we must delete the force just calculated, as it refers to self force
	vec_uint4 zerorsq = spu_cmpeq(rsq, spu_splats(0.0f));
	extPartMass = spu_sel(extPartMass, spu_splats(0.0f), zerorsq);
	idt2 = spu_sel((partMass + extPartMass)*b, spu_splats(0.0f), zerorsq); // (timescale)^-2
        // of interaction
	
        partAccX += rx * (b * extPartMass);
        partAccY += ry * (b * extPartMass);
        partAccZ += rz * (b * extPartMass);
        partPot -= extPartMass * a;
        //if(idt2 > part[i].dtGrav)
        //  part[i].dtGrav = idt2;
	partDtGrav = fmaxf4(idt2, partDtGrav);
      //}
    }
    for (k=0; k<4; ++k) {
      CellGravityParticle *particle = &part[i+k];
      particle->potential = spu_extract(partPot, k);
      particle->treeAcceleration.x = spu_extract(partAccX, k);
      particle->treeAcceleration.y = spu_extract(partAccY, k);
      particle->treeAcceleration.z = spu_extract(partAccZ, k);
      particle->dtGrav = spu_extract(partDtGrav, k);
    }
    //printf("intpart %d: %f (%f %f %f) %f / (%f %f %f) %f %f\n",i,part[i].core.mass,part[i].core.position.x,part[i].core.position.y,part[i].core.position.z,part[i].core.soft,part[i].treeAcceleration.x,part[i].treeAcceleration.y,part[i].treeAcceleration.z,part[i].potential,part[i].dtGrav);
  }
}
/* Original version */
/*
void particleBucketForce(CellExternalGravityParticle *ext, int numExt, CellGravityParticle *part, int numPart) {
  CellVector3D r;
  cellSPEtype rsq;
  cellSPEtype twoh, a, b;
  cellSPEtype idt2;
  int i,j;
  
  for (i=0; i<numPart; ++i) {
    for (j=0; j<numExt; ++j) {
      //printf("extpart %d: %f (%f %f %f) %f\n",j,ext[j].mass,ext[j].position.x,ext[j].position.y,ext[j].position.z,ext[j].soft);
      r.x = ext[j].position.x - part[i].core.position.x;
      r.y = ext[j].position.y - part[i].core.position.y;
      r.z = ext[j].position.z - part[i].core.position.z;
      rsq = r.x * r.x + r.y * r.y + r.z * r.z;
      twoh = ext[j].soft + part[i].core.soft;
      if(rsq != 0) {

        cellSPEtype r2_inv, u,dih,dir,r2twoh;
	qword vec_rsq = si_from_float((cellSPEtype)rsq);
	vec_float4 vec_r2 = spu_rsqrte((vec_float4)vec_rsq);
	r2_inv = (cellSPEtype)si_to_float((qword)vec_r2);
        //r2 = sqrt(rsq);
	r2twoh = r2_inv * twoh;
        if (1.0f < r2twoh) {
          dih = 2.0f/(twoh);
          u = 2.0f/r2twoh;
          if (u < 1.0f) {
            a = dih*(7.0f/5.0f - 2.0f/3.0f*u*u + 3.0f/10.0f*u*u*u*u
                - 1.0f/10.0f*u*u*u*u*u);
            b = dih*dih*dih*(4.0f/3.0f - 6.0f/5.0f*u*u + 1.0f/2.0f*u*u*u);
          }
          else {
            dir = r2_inv;
            a = -1.0f/15.0f*dir + dih*(8.0f/5.0f - 4.0f/3.0f*u*u + u*u*u
                - 3.0f/10.0f*u*u*u*u + 1.0f/30.0f*u*u*u*u*u);
            b = -1.0f/15.0f*dir*dir*dir + dih*dih*dih*(8.0f/3.0f - 3.0f*u + 6.0f/5.0f*u*u - 1.0f/6.0f*u*u*u);
          }
        }
        else {
          a = r2_inv;
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
    //printf("intpart %d: %f (%f %f %f) %f / (%f %f %f) %f %f\n",i,part[i].core.mass,part[i].core.position.x,part[i].core.position.y,part[i].core.position.z,part[i].core.soft,part[i].treeAcceleration.x,part[i].treeAcceleration.y,part[i].treeAcceleration.z,part[i].potential,part[i].dtGrav);
  }
}
*/

void BucketEwald(cellSPEtype *output, CellEwaldContainer *data) {
  //CellGravityParticle *p;
  int numPart = data->numPart;
  CellMultipoleMoments *mom = &data->rootMoments;
  vec_float4 Q2;
  vec_float4 L,fEwCut2,fInner2,alpha,alpha2,alphan_1,alphan_2,k1,ka,twoalpha2;
  vec_float4 fPot,ax,ay,az;
  vec_float4 dx,dy,dz,x,y,z,r2,r2_1,dir,dir2,a;
  vec_float4 Q2mirx,Q2miry,Q2mirz,Q2mir,Qta;
  vec_float4 g0,g1,g2,g3;
  vec_float4 g01,g11,g21,g31;
  vec_float4 g02,g12,g22,g32;
  vec_float4 hdotx,s,c;
  int i,j,ix,iy,iz,nEwReps,bInHole,bInHolex,bInHolexy;
  int offset = (numPart+3) & ~0x3;
  vec_float4 *positionX = (vec_float4*)(((char*)data)+sizeof(CellEwaldContainer));
  vec_float4 *positionY = (vec_float4*)(((char*)positionX)+offset*sizeof(cellSPEtype));
  vec_float4 *positionZ = (vec_float4*)(((char*)positionY)+offset*sizeof(cellSPEtype));
  CellEWT *ewt = (CellEWT*)(((char*)positionZ)+offset*sizeof(cellSPEtype));
  vec_float4 *pot = (vec_float4*)output;
  vec_float4 *accX = (vec_float4*)&output[offset];
  vec_float4 *accY = (vec_float4*)&output[2*offset];
  vec_float4 *accZ = (vec_float4*)&output[3*offset];
  vec_float4 zero = spu_splats(0.0f);

  Q2 = spu_splats(0.5f*(mom->xx + mom->yy + mom->zz));

  //n = req->lastParticle - req->firstParticle + 1;
  //p = &myParticles[req->firstParticle];
  nEwReps = (int) ceil(data->fEwCut);
  L = spu_splats(data->fPeriod);
  fEwCut2 = spu_splats(data->fEwCut)*spu_splats(data->fEwCut)*L*L;
  fInner2 = spu_splats(3.0e-3f)*L*L;
  nEwReps = nEwReps > data->nReps ? nEwReps : data->nReps;
  alpha = spu_splats(2.0f)*spu_re(L);
  alpha2 = alpha*alpha;
  k1 = spu_splats((float)M_PI)*spu_re(alpha2*L*L*L);
  ka = spu_splats(2.0f)*alpha*spu_rsqrte(spu_splats((float)M_PI));
  for(j=0; j<(offset>>2); ++j) {
    //if (p[j].rung < activeRung) continue;
    fPot = spu_splats(mom->totalMass)*k1;
    ax = spu_splats(0.0f);
    ay = spu_splats(0.0f);
    az = spu_splats(0.0f);
    dx = positionX[j] - spu_splats(mom->cm.x);
    dy = positionY[j] - spu_splats(mom->cm.y);
    dz = positionZ[j] - spu_splats(mom->cm.z);
    for (ix=-nEwReps;ix<=nEwReps;++ix) {
      bInHolex = (ix >= -data->nReps && ix <= data->nReps);
      x = dx + spu_splats((float)ix)*L;
      for(iy=-nEwReps;iy<=nEwReps;++iy) {
	bInHolexy = (bInHolex && iy >= -data->nReps && iy <= data->nReps);
	y = dy + spu_splats((float)iy)*L;
	for(iz=-nEwReps;iz<=nEwReps;++iz) {
	  bInHole = (bInHolexy && iz >= -data->nReps && iz <= data->nReps);
	  z = dz + spu_splats((float)iz)*L;
	  r2 = x*x + y*y + z*z;
	  vec_uint4 comparison = spu_cmpgt(r2, fEwCut2);
	  if (bInHole) comparison = spu_splats((unsigned int)0);
	  unsigned int comparator = spu_extract(spu_gather(comparison), 0);
	  if (comparator == 0xF) continue;
	  vec_uint4 comparison2 = spu_cmpgt(r2, fInner2);
	  unsigned int comparator2 = spu_extract(spu_gather(comparison2), 0);
	  twoalpha2 = spu_splats((float)2)*alpha2;
	  if (comparator2) {
	    dir = spu_rsqrte(r2);
	    dir2 = dir*dir;
	    a = expf4(-r2*alpha2);
	    a *= ka*dir2;
	    if (bInHole) g01 = -erff4(alpha/dir);
	    else g01 = erfcf4(alpha/dir);
	    g01 *= dir;
	    g11 = g01*dir2 + a;
	    alphan_1 = twoalpha2;
	    g21 = spu_splats((float)3)*g11*dir2 + alphan_1*a;
	    alphan_1 *= twoalpha2;
	    g31 = spu_splats((float)5)*g21*dir2 + alphan_1*a;
	    //alphan *= 2*alpha2;
	    //g4 = 7*g3*dir2 + alphan*a;
	    //alphan *= 2*alpha2;
	    //g5 = 9*g4*dir2 + alphan*a;
	  }
	  else if ((~comparator2) & 0xF) {
	    alphan_2 = ka;
	    r2 *= alpha2;
	    g02 = alphan_2*(spu_splats(1.0f/3.0f)*r2 - spu_splats(1.0f));
	    alphan_2 *= alpha2;
	    g12 = alphan_2*(spu_splats(1.0f/5.0f)*r2 - spu_splats(1.0f/3.0f));
	    alphan_2 *= twoalpha2;
	    g22 = alphan_2*(spu_splats(1.0f/7.0f)*r2 - spu_splats(1.0f/5.0f));
	    alphan_2 *= twoalpha2;
	    g32 = alphan_2*(spu_splats(1.0f/9.0f)*r2 - spu_splats(1.0f/7.0f));
	    //alphan *= 2*alpha2;
	    //g4 = alphan*(spu_splats(1.0f/11.0f)*r2 - spu_splats(1.0f/9.0f));
	    //alphan *= 2*alpha2;
	    //g5 = alphan*(spu_splats(1.0f/13.0f)*r2 - spu_splats(1.0f/11.0f));
	  }
	  g0 = spu_sel(g02, g01, comparison2);
	  g1 = spu_sel(g12, g11, comparison2);
	  g2 = spu_sel(g22, g21, comparison2);
	  g3 = spu_sel(g32, g31, comparison2);

	  Q2mirx = spu_splats(mom->xx)*x + spu_splats(mom->xy)*y + spu_splats(mom->xz)*z;
	  Q2miry = spu_splats(mom->xy)*x + spu_splats(mom->yy)*y + spu_splats(mom->yz)*z;
	  Q2mirz = spu_splats(mom->xz)*x + spu_splats(mom->yz)*y + spu_splats(mom->zz)*z;
	  Q2mir = spu_splats(0.5f)*(Q2mirx*x + Q2miry*y + Q2mirz*z);
	  Qta = g1*spu_splats(mom->totalMass) - g2*Q2 + g3*Q2mir;
	  fPot -= spu_sel(g0*spu_splats(mom->totalMass) - g1*Q2 + g2*Q2mir, zero, comparison);
	  ax += spu_sel(g2*Q2mirx - x*Qta, zero, comparison);
	  ay += spu_sel(g2*Q2miry - y*Qta, zero, comparison);
	  az += spu_sel(g2*Q2mirz - z*Qta, zero, comparison);
	}
      }
    }
    vec_float4 vec_hx, vec_hy, vec_hz, vec_hCfac, vec_hSfac, tmp;
    for (i=0; i<data->nEwhLoop; ++i) {
      vec_hx = spu_splats(ewt[i].hx);
      vec_hy = spu_splats(ewt[i].hy);
      vec_hz = spu_splats(ewt[i].hz);
      vec_hCfac = spu_splats(ewt[i].hCfac);
      vec_hSfac = spu_splats(ewt[i].hSfac);
      hdotx = vec_hx*dx + vec_hy*dy + vec_hz*dz;
      sincosf4(hdotx, &s, &c);
      fPot += vec_hCfac*c + vec_hSfac*s;
      tmp = vec_hCfac*s - vec_hSfac*c;
      ax += vec_hx*tmp;
      ay += vec_hy*tmp;
      az += vec_hz*tmp;
    }
    pot[j] = fPot;
    accX[j] = ax;
    accY[j] = ay;
    accZ[j] = az;
  }
  return;
}

void funcLookup(int funcIndex, void* readWritePtr, int readWriteLen, void* readOnlyPtr, int readOnlyLen,
                void* writeOnlyPtr, int writeOnlyLen, DMAListEntry* dmaList ) {
  CellContainer *container = (CellContainer*)readOnlyPtr;
  //printf("funcLookup %d: %d int %d ext, %d %d\n",funcIndex,container->numInt,container->numExt,readWriteLen,readOnlyLen);
  if (funcIndex==1) {
    // node interaction
    nodeBucketForce(&container->data, container->numExt,
                    readWritePtr, container->numInt);
  } else if (funcIndex==2) {
    // particle interaction
    particleBucketForce(&container->data, container->numExt,
                        readWritePtr, container->numInt);
  } else if (funcIndex==3) {
    // ewald computation
    BucketEwald((cellSPEtype*)writeOnlyPtr, (CellEwaldContainer*)readOnlyPtr);
  }
}
