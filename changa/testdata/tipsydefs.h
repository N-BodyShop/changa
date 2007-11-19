#define MAXDIM 3
#define forever for(;;)
#include <rpc/types.h>
#include <rpc/xdr.h>

struct gas_particle {
    float mass;
    float pos[MAXDIM];
    float vel[MAXDIM];
    float rho;
    float temp;
    float hsmooth;
    float metals ;
    float phi ;
} ;

struct dark_particle {
    float mass;
    float pos[MAXDIM];
    float vel[MAXDIM];
    float eps;
    float phi ;
} ;

struct star_particle {
    float mass;
    float pos[MAXDIM];
    float vel[MAXDIM];
    float metals ;
    float tform ;
    float eps;
    float phi ;
} ;

struct dump {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
} ;

struct dump header ;

struct dump_pad {
    double time ;
    int nbodies ;
    int ndim ;
    int nsph ;
    int ndark ;
    int nstar ;
    int pad;
} ;
struct dump_pad header_pad ;

int xdr_header(XDR *pxdrs,struct dump *ph)
{
        int pad = 0;

        if (!xdr_double(pxdrs,&ph->time)) return 0;
        if (!xdr_int(pxdrs,&ph->nbodies)) return 0;
        if (!xdr_int(pxdrs,&ph->ndim)) return 0;
        if (!xdr_int(pxdrs,&ph->nsph)) return 0;
        if (!xdr_int(pxdrs,&ph->ndark)) return 0;
        if (!xdr_int(pxdrs,&ph->nstar)) return 0;
        if (!xdr_int(pxdrs,&pad)) return 0;
        return 1;
        }
  
int xdr_gas(XDR *pxdrs,struct gas_particle *ph)
{ 
        int i;
        if (!xdr_float(pxdrs,&ph->mass)) return 0;
        for(i=0; i<MAXDIM; i++) if (!xdr_float(pxdrs,&ph->pos[i])) return 0;
        for(i=0; i<MAXDIM; i++) if (!xdr_float(pxdrs,&ph->vel[i])) return 0;
        if (!xdr_float(pxdrs,&ph->rho)) return 0;
        if (!xdr_float(pxdrs,&ph->temp)) return 0;
        if (!xdr_float(pxdrs,&ph->hsmooth)) return 0;
        if (!xdr_float(pxdrs,&ph->metals)) return 0;
        if (!xdr_float(pxdrs,&ph->phi)) return 0;
        return 1;
        }

int xdr_dark(XDR *pxdrs,struct dark_particle *ph)
{
        int i;
        if (!xdr_float(pxdrs,&ph->mass)) return 0;
        for(i=0; i<MAXDIM; i++) if (!xdr_float(pxdrs,&ph->pos[i])) return 0;
        for(i=0; i<MAXDIM; i++) if (!xdr_float(pxdrs,&ph->vel[i])) return 0;
        if (!xdr_float(pxdrs,&ph->eps)) return 0;
        if (!xdr_float(pxdrs,&ph->phi)) return 0;
        return 1;
        }

int xdr_star(XDR *pxdrs,struct star_particle *ph)
{
        int i;
        if (!xdr_float(pxdrs,&ph->mass)) return 0;
        for(i=0; i<MAXDIM; i++) if (!xdr_float(pxdrs,&ph->pos[i])) return 0;
        for(i=0; i<MAXDIM; i++) if (!xdr_float(pxdrs,&ph->vel[i])) return 0;
        if (!xdr_float(pxdrs,&ph->metals)) return 0;
        if (!xdr_float(pxdrs,&ph->tform)) return 0;
        if (!xdr_float(pxdrs,&ph->eps)) return 0;
        if (!xdr_float(pxdrs,&ph->phi)) return 0;
        return 1;
        }

