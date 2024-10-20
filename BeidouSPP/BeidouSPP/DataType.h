#pragma once
#include <math.h>
#include <time.h>
#include <stdint.h>

#define PI          3.1415926535897932  /* pi */
#define CLIGHT      299792458.0         /* speed of light (m/s) */
#define GM          3.986004418E14      /*  standard gravitational parameter of earth   */
#define  OMGe         7.2921151467E-5    /* rad of earth's rotation*/
#define RE_WGS84    6378137.0           /* earth semimajor axis (WGS84) (m) */
#define FE_WGS84    (1.0/298.257223563) /* earth flattening (WGS84) */
#define FREQ1_CMP   1.561098E9          /* BDS B1I     frequency (Hz) */
#define FREQ1       1.57542E9           /* L1/E1/B1C  frequency (Hz) */

#define NFREQ        1                   /* number of carrier frequencies */

/* time struct */
struct gtime_t
{
    time_t time;        /* time (s) expressed by standard time_t */
    double sec;         /* fraction of second under 1 s */
};

struct gpst_t
{
    int week;               /* gps week*/
    double sec;            /* time of week*/
};

/* observation data record */
struct obsd_t
{
    bool svh;                   /* sv health flag*/
    int EpochIndex;        /* index of epoch */
    CString sat;              /* satellite/receiver number */
    gtime_t time;       /* receiver sampling time (GPST) */
    double P[NFREQ]; /* observation data pseudorange (m) */
};

/* observation data */
struct obs_t
{
    int*validIndex;     /* the index of valid data in all data */
    int n,nmax;         /* number of obervation data/allocated*/
    obsd_t* data;       /* observation data records */
};

/* BeiDou broadcast ephemeris type */
struct eph_t
{
    CString sat;            /* satellite number */
    //int flag;           /* GPS/QZS: L2 P data flag */
    /* BDS: nav type (0:unknown,1:IGSO/MEO,2:GEO) */
    gtime_t toc, ttr; /* Toe,Toc,T_trans */
    double toe;         /* Toe (s) in week */
    int sva;               /* SV accuracy (URA index) */
    /* SV orbit parameters */
    double sqrtA, e, i0, OMG0, omg, M0, deln, OMGd, idot;
    double crc, crs, cuc, cus, cic, cis;        
    double fit;         /* fit interval (h) */
    double f0, f1, f2;    /* SV clock parameters (af0,af1,af2) */
    double tgd[6];      /* group delay parameters */
    /* CMP:tgd[0]=TGD_B1I ,tgd[1]=TGD_B2I/B2b,tgd[2]=TGD_B1Cp */
    /*     tgd[3]=TGD_B2ap,tgd[4]=ISC_B1Cd   ,tgd[5]=ISC_B2ad */
};

/* navigation data type */
struct  nav_t
{
    int n;         /* number of broadcast ephemeris */
    eph_t* eph;         /* BDS ephemeris */
    double utc_cmp[8];  /* BeiDou UTC parameters {A0,A1,Tot,WNt,dt_LS,WN_LSF,DN,dt_LSF} */
    double ion_cmp[8];  /* BeiDou iono model parameters {a0,a1,a2,a3,b0,b1,b2,b3} */
};

/* solution type */
struct sol_t
{
    int chisqrTest;   /* Chi-square Test */
    gtime_t time;       /* time (GPST) */
    double rr[3];       /* position (m) */
    /* {x,y,z} or {e,n,u} */
    double  qr[6];       /* position variance/covariance (m^2) */
    /* {c_xx,c_yy,c_zz,c_xy,c_yz,c_zx} or */
    /* {c_ee,c_nn,c_uu,c_en,c_nu,c_ue} */
    double dtr;      /* receiver clock bias to time systems (s) */
    double gdop, pdop;  /* gdop/pdop */
    int ns;         /* number of valid satellites */
};

/* solution buffer type */
struct solbuf_t
{
    int n, nmax;         /* number of solution/max number of buffer */
    int start, end;      /* start/end index */
    gtime_t time;       /* current solution time */
    sol_t* data;        /* solution data */
    double rb[3];       /* reference position {x,y,z} (ecef) (m) */
};
