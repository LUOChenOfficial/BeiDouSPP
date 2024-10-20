#pragma once
#include "..//src/Matrix.h"
#include "..//src/Angle.h"
#include <math.h>
#include <time.h>
#include "DataType.h"

class CBdSPP
{
private:
    obs_t*obss;          /* observation data of all epochs */
    nav_t navs;          /* navigation data */
    solbuf_t sol;          /* postioning solution*/
    gtime_t ts, te;     /* time of start and end */
    double ti;              /* time interval */
    int EpochCount;             /*total count of epochs */
    double maxGDOP;         /* threshold of GDOP */
    CAngle minElevAngle;    /* minium of satellite elevation angle */
    double rrInit[4];               /* approximate receiver pos from obs data */
    double rrTmp[4];               /* receiver pos of current time */
public:
    CBdSPP();
    ~CBdSPP();
public:
    /* main functions */
    /* read obs data */
    void ReadObsData(CString ObsFile);
    /* read nav data */
    void ReadNavData(CString NavFile);
    /* calculate satellite position of every epoch */
    CMatrix CalSatPosClkEachEpoch(int epoIndex);
    /* form design matrix and constant matrix (B&L) */
    void FormDesConsMatEachEpoch(int epoIndex, CMatrix SatPosClkEpoch,CMatrix& B, CMatrix& L, double* varIon, double* varTrp);
    /* form weight matrix (P) */
    void FormWeightMatEachEpoch(int epoIndex, CMatrix SatPosClkEpoch, CMatrix& P, double*varIon, double*varTrp);
    /* output solutions */
    void OutputSol(CString SolFile, CString strObsFile, CString strNavFile);

    void mainfunc();
private:
    /* tool functions*/
    /* split string*/
    CString* SplitString(CString str, char split, int& iSubStrs);
    /* convert string ("yyyy mm dd hh mm ss") to gtime_t struct*/
    void str2time(CString* strTime, gtime_t&time);
    /* convert gtime_t struct to calendar day */
    void time2epoch(gtime_t t, double* ep);
    /* convert gtime_t struct to gpst*/
    gpst_t time2gpst(gtime_t t);
    /* compute the diff of gtime_t struct */
    double timediff(gtime_t t1, gtime_t t2);
    /* compute the sum of gtime_t struct and double */
    gtime_t timeadd(gtime_t t, double sec);
    /* select ephemeris */
    eph_t seleph(obsd_t data);
    /* compute geometric distance and receiver-to-satellite vector *///including Sagnac effect correction
    double geodist(double* rr, double* rs, double* vec, double& SagEffectCor);
    /* compute satellite azimuth/elevation */
    void satazel(double* rr, double* vec,double* azel);
    /* convert ecef to blh*/
    void ecef2blh(double* ecef, double* blh);
    /* compute tropospheric delay by standard atmosphere and saastamoinen model */
    double tropmodel(double* pos,double el, double humi);
    /* compute ionospheric delay by broadcast ionosphere model (klobuchar model) */
    double ionmodel(gtime_t t, double* ion, double* pos,double* azel);
};

