#include "pch.h"
#include "BdSPP.h"

CBdSPP::CBdSPP()
{
	obss = { 0 };
	navs = { 0 };
	sol = { 0 };
	maxGDOP = 30.0;
	minElevAngle(DEG) = 10;
}

CBdSPP::~CBdSPP()
{
}

CString* CBdSPP::SplitString(CString str, char split, int& iSubStrs)
{
	int iPos = 0; //分割符位置
	int iNums = 0; //分割符的总数
	CString strTemp = str;
	CString strRight;
	//先计算子字符串的数量
	while (iPos != -1)
	{
		iPos = strTemp.Find(split);
		if (iPos == -1)
		{
			break;
		}
		strRight = strTemp.Mid(iPos + 1, str.GetLength());
		strTemp = strRight;
		iNums++;
	}
	if (iNums == 0) //没有找到分割符
	{
		//子字符串数就是字符串本身
		iSubStrs = 1;
		return NULL;
	}
	//子字符串数组
	iSubStrs = iNums + 1; //子串的数量 = 分割符数量 + 1
	CString* pStrSplit;
	pStrSplit = new CString[iSubStrs];
	strTemp = str;
	CString strLeft;
	for (int i = 0; i < iNums; i++)
	{
		iPos = strTemp.Find(split);
		//左子串
		strLeft = strTemp.Left(iPos);
		//右子串
		strRight = strTemp.Mid(iPos + 1, strTemp.GetLength());
		strTemp = strRight;
		pStrSplit[i] = strLeft;
	}
	pStrSplit[iNums] = strTemp;
	return pStrSplit;
}

void CBdSPP::str2time(CString* strTime, gtime_t& time)
{
	int year, mon, day, hour, min;
	double sec;
	int day_all;
	int doy[12] = { 1,32,60,91,121,152,182,213,244,274,305,335 };
	year = _ttoi(strTime[0]);
	mon = _ttoi(strTime[1]);
	day = _ttoi(strTime[2]);
	hour = _ttoi(strTime[3]);
	min = _ttoi(strTime[4]);
	sec = _tstof(strTime[5]);
	day_all = (year - 1970) * 365 + (year - 1969) / 4 + doy[int(mon - 1)] + day - 2;
	if (((year % 4 == 0 && year % 100 != 0) || year % 400 == 0) && mon >= 3)
		day_all = day_all + 1;
	time.time = day_all * 86400 + hour * 3600 + min * 60 + int(sec);
	time.sec = sec - int(sec);
}

void CBdSPP::time2epoch(gtime_t t, double* ep)
{
	const int mday[] = { /* # of days in a month */
		31,28,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31,
		31,29,31,30,31,30,31,31,30,31,30,31,31,28,31,30,31,30,31,31,30,31,30,31
	};
	int days, sec, mon, day;

	/* leap year if year%4==0 in 1901-2099 */
	days = (int)(t.time / 86400);
	sec = (int)(t.time - (time_t)days * 86400);
	for (day = days % 1461, mon = 0; mon < 48; mon++) {
		if (day >= mday[mon]) day -= mday[mon]; else break;
	}
	ep[0] = 1970 + days / 1461 * 4 + mon / 12; ep[1] = mon % 12 + 1; ep[2] = day + 1;
	ep[3] = sec / 3600; ep[4] = sec % 3600 / 60; ep[5] = sec % 60 + t.sec;
}

gpst_t CBdSPP::time2gpst(gtime_t t)
{
	gpst_t gpst;
	CString* gpst0 = new CString[6];
	gpst0[0] = _T("1980");
	gpst0[1] = _T("1");
	gpst0[2] = _T("6");
	gpst0[3] = _T("0");
	gpst0[4] = _T("0");
	gpst0[5] = _T("0");
	gtime_t t0;
	str2time(gpst0, t0);
	time_t sec = t.time - t0.time;
	int week = (int)(sec / (86400 * 7));
	gpst.week = week;
	gpst.sec = (double)(sec - (double)week * 86400 * 7) + t.sec;
	return gpst;
}

double CBdSPP::timediff(gtime_t t1, gtime_t t2)
{
	return difftime(t1.time, t2.time) + t1.sec - t2.sec;
}

gtime_t CBdSPP::timeadd(gtime_t t, double sec)
{
	double tt;

	t.sec += sec; tt = floor(t.sec); t.time += (int)tt; t.sec -= tt;
	return t;
}

eph_t CBdSPP::seleph(obsd_t data)
{
	int isFound = 0;
	eph_t ephTmp = navs.eph[0];
	gpst_t time = time2gpst(data.time);
	double dt= fabs(time.sec-navs.eph[0].toe);
	for (int i = 0; i < navs.n; i++)
	{
		if (navs.eph[i].sat != data.sat)continue;
		isFound = 1;
		double dtTmp;
		dtTmp = fabs(time.sec - navs.eph[i].toe);
		if (dtTmp > dt)continue;
		ephTmp = navs.eph[i];
		dt = dtTmp;
	}
	if (!isFound)return { 0 };
	return ephTmp;
}

double CBdSPP::geodist(double* rr, double* rs, double* vec,double& SagEffectCor)
{
	double dist,distSquared=0;
	/* check satellite postion */
	if (sqrt(pow(rs[0], 2) + pow(rs[1], 2) + pow(rs[2], 2)) < RE_WGS84)return-1;
	for (int i = 0; i < 3; i++)
	{
		vec[i] = rs[i] - rr[i];
		distSquared += pow(rr[i] - rs[i], 2);
	}
	dist = sqrt(distSquared);
	/* progress sagnac effect correction*/
	SagEffectCor = OMGe * (rs[0] * rr[1] - rs[1] * rr[0]) / CLIGHT;
	return dist;
}

void CBdSPP::satazel(double* rr, double* vec,double* azel)
{
	/* convert rr(ecef) to rr(blh) */
	double rr_blh[3];
	ecef2blh(rr, rr_blh);
	/* compute vec(ecef) to vec(enu) rotate matrix */
	CMatrix M(3,3);
	double B = rr_blh[0], L = rr_blh[1], H = rr_blh[2];
	M(0, 0) = -sin(L);				   M(0, 1) = cos(L);					M(0, 2) = 0.0;
	M(1, 0) = -sin(B) * cos(L);  M(1, 1) = -sin(B) * sin(L);  M(1, 2) = cos(B);
	M(2, 0) = cos(B) * cos(L);   M(2, 1) = cos(B) * sin(L);  M(2, 2) = sin(B);
	/* convert vec(ecef) to vec(enu) */
	CMatrix VecMat(3, 1),VecTransMat(3,1);
	VecMat(0, 0) = vec[0];	VecMat(1, 0) = vec[1]; VecMat(2, 0) = vec[2];
	VecTransMat = M * VecMat;
	double e = VecTransMat(0, 0);
	double n = VecTransMat(1, 0);
	double u = VecTransMat(2, 0);
	double norm = sqrt(pow(VecTransMat(0, 0), 2) + pow(VecTransMat(1, 0), 2) + pow(VecTransMat(2, 0), 2));
	double az = (pow(VecTransMat(0, 0), 2) + pow(VecTransMat(1, 0), 2)) < 1E-12 ? 0.0 : atan2(e/norm, n/norm);
	if (az < 0.0) az += 2 * PI;
	azel[0] = az;
	azel[1]=asin(u / norm);
}

void CBdSPP::ecef2blh(double* ecef, double* blh)
{
	double e2 = FE_WGS84 * (2.0 - FE_WGS84);
	double r2 = pow(ecef[0], 2) + pow(ecef[1], 2) , z, zk, v = RE_WGS84, sinp;
	for (z = ecef[2], zk = 0.0; fabs(z - zk) >= 1E-4;) {
		zk = z;
		sinp = z / sqrt(r2 + z * z);
		v = RE_WGS84 / sqrt(1.0 - e2 * sinp * sinp);
		z = ecef[2] + v * e2 * sinp;
	}
	blh[0] = r2 > 1E-12 ? atan(z / sqrt(r2)) : (ecef[2] > 0.0 ? PI / 2.0 : -PI / 2.0);
	blh[1] = r2 > 1E-12 ? atan2(ecef[1], ecef[0]) : 0.0;
	blh[2] = sqrt(r2 + z * z) - v;
}

double CBdSPP::tropmodel(double* pos,double el, double humi)
{
	const double temp0 = 15.0; /* temparature at sea level */
	double hgt, pres, temp, e, z, trph, trpw;

	if (pos[2] < -100.0 || 1E4 < pos[2] || el <= 0) return 0.0;

	/* standard atmosphere */
	hgt = pos[2] < 0.0 ? 0.0 : pos[2];

	pres = 1013.25 * pow(1.0 - 2.2557E-5 * hgt, 5.2568);
	temp = temp0 - 6.5E-3 * hgt + 273.16;
	e = 6.108 * humi * exp((17.15 * temp - 4684.0) / (temp - 38.45));

	/* saastamoninen model */
	z = PI / 2.0 - el;
	trph = 0.0022768 * pres / (1.0 - 0.00266 * cos(2.0 * pos[0]) - 0.00028 * hgt / 1E3) / cos(z);
	trpw = 0.002277 * (1255.0 / temp + 0.05) * e / cos(z);
	return trph + trpw;
}

double CBdSPP::ionmodel(gtime_t t, double* ion, double* pos,double* azel)
{
	double ion_default[] = { /* 2004/1/1 */
		0.1118E-07,-0.7451E-08,-0.5961E-07, 0.1192E-06,
		0.1167E+06,-0.2294E+06,-0.1311E+06, 0.1049E+07
	};
	double tt, f, psi, phi, lam, amp, per, x;
	int week;
	double norm = 0;
	for (int i = 0; i < 8; i++)norm += pow(ion[i], 2);

	if (pos[2] < -1E3 || azel[1] <= 0) return 0.0;
	if (norm <= 0.0) ion = ion_default;

	/* earth centered angle (semi-circle) */
	psi = 0.0137 / (azel[1] / PI + 0.11) - 0.022;

	/* subionospheric latitude/longitude (semi-circle) */
	phi = pos[0] / PI + psi * cos(azel[0]);
	if (phi > 0.416) phi = 0.416;
	else if (phi < -0.416) phi = -0.416;
	lam = pos[1] / PI + psi * sin(azel[0]) / cos(phi * PI);

	/* geomagnetic latitude (semi-circle) */
	phi += 0.064 * cos((lam - 1.617) * PI);

	//double sec = time2gpst(t).sec;
	/* local time (s) */
	tt = 43200.0 * lam + time2gpst(t).sec;
	tt -= floor(tt / 86400.0) * 86400.0; /* 0<=tt<86400 */

	/* slant factor */
	f = 1.0 + 16.0 * pow(0.53 - azel[1] / PI, 3.0);
	/* ionospheric delay */
	amp = ion[0] + phi * (ion[1] + phi * (ion[2] + phi * ion[3]));
    per = ion[4] + phi * (ion[5] + phi * (ion[6] + phi * ion[7]));
	amp = amp < 0.0 ? 0.0 : amp;
	per = per < 72000.0 ? 72000.0 : per;
	x = 2.0 * PI * (tt - 50400.0) / per;
	//double result=CLIGHT * f * (fabs(x) < 1.57 ? 5E-9 + amp * (1.0 + x * x * (-0.5 + x * x / 24.0)) : 5E-9);
	return CLIGHT * f * (fabs(x) < 1.57 ? 5E-9 + amp * (1.0 + x * x * (-0.5 + x * x / 24.0)) : 5E-9);
}

void CBdSPP::ReadObsData(CString ObsFile)
{
	CStdioFile sf;
	if (!sf.Open(ObsFile, CFile::modeRead)) return;
	CString strTmp;
	int iSubStrs;

	/* read 1st part of header until  'APPROX POSITION XYZ' */
	int HeadCount1 = 9;
	while (HeadCount1--)
	{
		sf.ReadString(strTmp);
	}

	/* read 'APPROX POSITION XYZ' */
	sf.ReadString(strTmp);
	rrInit[0] = _tstof(strTmp.Mid(1, 13));
	rrInit[1] = _tstof(strTmp.Mid(16, 12));
	rrInit[2] = _tstof(strTmp.Mid(30, 12));
	rrInit[3] = 0;
	memcpy(rrTmp, rrInit, sizeof(rrInit));
	iSubStrs = 0;

	/* read 2st part of header until  'TIME INTERVAL' */
	int HeadCount2 = 13;
	while (HeadCount2--)
	{
		sf.ReadString(strTmp);
	}

	/* read 'TIME INTERVAL' */
	sf.ReadString(strTmp);
	ti = _tstof(strTmp.Mid(4,6));
	
	iSubStrs = 0;

	/* read 'TIME OF FIRST OBS'*/
	sf.ReadString(strTmp);
	CString* strTimeStart = new CString[6];
	strTimeStart[0] = strTmp.Mid(2, 4);
	strTimeStart[1] = strTmp[11];
	strTimeStart[2] = strTmp[17];
	strTimeStart[3] = strTmp[23];
	strTimeStart[4] = strTmp[29];
	strTimeStart[5] = strTmp.Mid(34, 9);
	str2time(strTimeStart, ts);
	delete[]strTimeStart;
	strTimeStart = NULL;
	iSubStrs = 0;

	/* read 'TIME OF LAST OBS'*/
	sf.ReadString(strTmp);
	CString* strTimeEnd = new CString[6];
	strTimeEnd[0] = strTmp.Mid(2, 4);
	strTimeEnd[1] = strTmp[11];
	strTimeEnd[2] = strTmp[17];
	strTimeEnd[3] = strTmp.Mid(22, 2);
	strTimeEnd[4] = strTmp.Mid(28, 2);
	strTimeEnd[5] = strTmp.Mid(33, 10);
	str2time(strTimeEnd, te);
	delete[]strTimeEnd;
	strTimeEnd = NULL;
	iSubStrs = 0;

	/* read rest header*/
	strTmp.Trim();
	CString strEndFlag("END OF HEADER");
	int num = 0;
	
	while (strTmp != strEndFlag)
	{
		sf.ReadString(strTmp);
		strTmp.Trim();
	}
	
	/* compute epochs and allocate space */
	EpochCount = int(difftime(te.time, ts.time) / ti) + 1;
	sol.data = new sol_t[EpochCount];
	obss = new obs_t[EpochCount];
	
	/* read body*/
	for (int index = 0; index < EpochCount; index++)
	{
		sf.ReadString(strTmp);
		CString*strData=new CString[6];
		strData[0] = strTmp.Mid(2, 4);
		strData[1] = strTmp.Mid(7, 2);
		strData[2] = strTmp.Mid(10, 2);
		strData[3] = strTmp.Mid(13, 2);
		strData[4] = strTmp.Mid(16, 2);
		strData[5] = strTmp.Mid(19, 2);
		for (int i = 0; i < iSubStrs; i++)
		{
			strData[i].Trim();
		}
		CString* strTime = new CString[6];
		strTime = strData;

		/* read gtime_t of each Epoch*/
		gtime_t TimeEpoch;
		str2time(strTime, TimeEpoch);
		int DataEpochCount = _ttoi(strTmp.Mid(33,3).Trim());
		iSubStrs = 0;
		obs_t obssTmp;
		int BdObsDataCount = 0;

		obssTmp.data = new obsd_t[DataEpochCount];
		/* read data records of each Epoch*/
		for (int i = 0; i < DataEpochCount; i++)
		{
			sf.ReadString(strTmp);
			//strTmp.Trim();
			obssTmp.data[i].sat = strTmp.Mid(0, 3);
			obssTmp.data[i].P[0] = _tstof(strTmp.Mid(21, 12));
			obssTmp.data[i].time = TimeEpoch;
			if (strTmp.Left(1) == 'C')BdObsDataCount++;
		}

		/* allocate space extract BeiDou records*/
		obss[index].nmax = BdObsDataCount;
		obss[index].data = new obsd_t[BdObsDataCount];
		for (int m = 0, n = 0; n < BdObsDataCount; m++, n++)
		{
			obss[index].data[n].P[0] = obssTmp.data[m].P[0];
			obss[index].data[n].time = obssTmp.data[m].time;
			obss[index].data[n].sat = obssTmp.data[m].sat;
			obss[index].data[n].EpochIndex = index;
			obss[index].data[n].svh = 1;
		}
	}
	
		sf.Close();

}

void CBdSPP::ReadNavData(CString NavFile)
{
	CStdioFile sf;
	if (!sf.Open(NavFile, CFile::modeRead)) return;
	CString strTmp;
	/* calculate BeiDou eph data count*/
	int BdEphDataCount = 0;
	while (sf.ReadString(strTmp))
	{
		if (strTmp.Left(1) == 'C')BdEphDataCount++;
	}
	sf.Close();

	navs.n = BdEphDataCount;
	navs.eph = new eph_t[BdEphDataCount];
	CStdioFile sfTmp;
	if (!sfTmp.Open(NavFile, CFile::modeRead)) return;

	while (sfTmp.ReadString(strTmp))
	{
		if (strTmp.Left(4) == "GPSA")
		{
			navs.ion_cmp[0] = _tstof(strTmp.Mid(6, 11).Trim());
			navs.ion_cmp[1] = _tstof(strTmp.Mid(18, 11).Trim());
			navs.ion_cmp[2] = _tstof(strTmp.Mid(30, 11).Trim());
			navs.ion_cmp[3] = _tstof(strTmp.Mid(42, 11).Trim());
			sfTmp.ReadString(strTmp);
			navs.ion_cmp[4] = _tstof(strTmp.Mid(6, 11).Trim());
			navs.ion_cmp[5] = _tstof(strTmp.Mid(18, 11).Trim());
			navs.ion_cmp[6] = _tstof(strTmp.Mid(30, 11).Trim());
			navs.ion_cmp[7] = _tstof(strTmp.Mid(42, 11).Trim());
			break;
		}
	}
	int i = 0;
	while (sfTmp.ReadString(strTmp))
	{
		if (strTmp.Left(1) == 'C')
		{
			navs.eph[i].sat = strTmp.Mid(0, 3);
			CString* strTime = new CString[6];
			gtime_t toc;
			strTime[0] = strTmp.Mid(4, 4);
			strTime[1] = strTmp.Mid(9, 2);
			strTime[2] = strTmp.Mid(12, 2);
			strTime[3] = strTmp.Mid(15, 2);
			strTime[4] = strTmp.Mid(18, 2);
			strTime[5] = strTmp.Mid(21, 2);
			str2time(strTime, toc);
			navs.eph[i].toc = toc;
			/* BDST to GPST*/
			navs.eph[i].toc.time += 14;
			navs.eph[i].f0 = _tstof(strTmp.Mid(23, 19).Trim());
			navs.eph[i].f1 = _tstof(strTmp.Mid(42, 19).Trim());
			navs.eph[i].f2 = _tstof(strTmp.Mid(61, 19).Trim());
			sfTmp.ReadString(strTmp);
			navs.eph[i].crs = _tstof(strTmp.Mid(23, 19).Trim());
			navs.eph[i].deln = _tstof(strTmp.Mid(42, 19).Trim());
			navs.eph[i].M0 = _tstof(strTmp.Mid(61, 19).Trim());
			sfTmp.ReadString(strTmp);
			navs.eph[i].cuc = _tstof(strTmp.Mid(4, 19).Trim());
			navs.eph[i].e= _tstof(strTmp.Mid(23, 19).Trim());
			navs.eph[i].cus = _tstof(strTmp.Mid(42, 19).Trim());
			navs.eph[i].sqrtA = _tstof(strTmp.Mid(61, 19).Trim());
			sfTmp.ReadString(strTmp);
			navs.eph[i].toe= _tstof(strTmp.Mid(4, 19).Trim());
			navs.eph[i].cic = _tstof(strTmp.Mid(23, 19).Trim());
			navs.eph[i].OMG0 = _tstof(strTmp.Mid(42, 19).Trim());
			navs.eph[i].cis = _tstof(strTmp.Mid(61, 19).Trim());
			sfTmp.ReadString(strTmp);
			navs.eph[i].i0 = _tstof(strTmp.Mid(4, 19).Trim());
			navs.eph[i].crc = _tstof(strTmp.Mid(23, 19).Trim());
			navs.eph[i].omg = _tstof(strTmp.Mid(42, 19).Trim());
			navs.eph[i].OMGd = _tstof(strTmp.Mid(61, 19).Trim());
			sfTmp.ReadString(strTmp);
			navs.eph[i].idot = _tstof(strTmp.Mid(4, 19).Trim());
			sfTmp.ReadString(strTmp);
			navs.eph[i].sva = _ttoi(strTmp.Mid(4, 19).Trim());
			navs.eph[i].tgd[0] = _tstof(strTmp.Mid(42, 19).Trim());
			sfTmp.ReadString(strTmp);
			i++;
		}
	}
	sfTmp.Close();
}

CMatrix CBdSPP::CalSatPosClkEachEpoch(int epoIndex)
{
		CMatrix SatPosClkMax(obss[epoIndex].nmax, 4);
		//double*dtsMax = new double[obss[epoIndex].nmax];
		//double*TGDMax= new double[obss[epoIndex].nmax];
		/* calculate detlta clock and position of every satellite */
		for (int satIndex = 0; satIndex < obss[epoIndex].nmax; satIndex++)
		{
			/* detlta clock */
			gtime_t ts, tr;
			double P, f0, f1, f2;
			double dt0,dt,dts;
			eph_t eph = seleph(obss[epoIndex].data[satIndex]);
			if (eph.e == 0)
			{
				SatPosClkMax(satIndex, 0) = 0;
				SatPosClkMax(satIndex, 1) = 0;
				SatPosClkMax(satIndex, 2) = 0;
				SatPosClkMax(satIndex, 3) = 0;
				continue;
			}
			double tk;
			gpst_t ts_gpst;
			double A = pow(eph.sqrtA, 2);
			double M;
			double e = eph.e, E = 0, ETmp = 10;

			tr = obss[epoIndex].data[satIndex].time;
			P = obss[epoIndex].data[satIndex].P[0];
			f0 = eph.f0;
			f1 = eph.f1;
			f2 = eph.f2;
			ts = timeadd(tr, -P / CLIGHT);
			dt0 = timediff(ts, eph.toc);
			double dtTmp = dt0;
			for (int i = 0; i < 2; i++)dt0 = dtTmp - (f0 + f1 * dt0 + f2 * dt0 * dt0);
			dt = f0 + f1 * dt0 + f2 * dt0 * dt0;

			/* sat postion */
			ts = timeadd(ts, -dt);
			ts_gpst = time2gpst(ts);
			tk = ts_gpst.sec-(eph.toe+14);

			M = eph.M0 + (sqrt(GM / pow(A, 3)) + eph.deln) * tk;
			while (fabs(ETmp - E) > 1e-13)
			{
				ETmp = E;
				E -= (E - e * sin(E) - M) / (1 - e * cos(E));
			};
			CString sat = obss[epoIndex].data[satIndex].sat;

			double u, r, i,u0;
			u0 = atan2(sqrt(1 - pow(eph.e, 2)) * sin(E) , (cos(E) - eph.e)) + eph.omg;
			r = A * (1 - eph.e * cos(E));
			i = eph.i0 + eph.idot * tk;
			u = u0 + eph.cus * sin(2 * u0) + eph.cuc * cos(2 * u0);
			r = r + eph.crs * sin(2 * u0) + eph.crc * cos(2 * u0);
			i = i + eph.cis * sin(2 * u0) + eph.cic * cos(2 * u0);
			double x = r * cos(u);
			double y = r * sin(u);
			/* MEO/IGSO satellite*/
			double L = eph.OMG0 + tk * (eph.OMGd - OMGe) - OMGe * eph.toe;
			SatPosClkMax(satIndex, 0) = x * cos(L) - y * cos(i) * sin(L);
			SatPosClkMax(satIndex, 1) = x * sin(L) + y * cos(i) * cos(L);
			SatPosClkMax(satIndex, 2) = y * sin(i);
			/* GEO satellite*/
			if (obss[epoIndex].data[satIndex].sat == _T("C01") ||
				obss[epoIndex].data[satIndex].sat == _T("C02") ||
				obss[epoIndex].data[satIndex].sat == _T("C03") ||
				obss[epoIndex].data[satIndex].sat == _T("C04") ||
				obss[epoIndex].data[satIndex].sat == _T("C05") ||
				obss[epoIndex].data[satIndex].sat == _T("C59") ||
				obss[epoIndex].data[satIndex].sat == _T("C60"))
			{
				double f = -5.0 / 180 * PI;
				double p = OMGe * tk;
				double L = eph.OMG0 + tk * eph.OMGd - OMGe * eph.toe;
				double X = x * cos(L) - y * cos(i) * sin(L);
				double Y = x * sin(L) + y * cos(i) * cos(L);
				double Z = y * sin(i);
				CMatrix Rx(3, 3), Rz(3, 3);
				Rx(0, 0) = 1;				Rx(0, 1) = 0;				Rx(0, 2) = 0;
				Rx(1, 0) = 0;				Rx(1, 1) = cos(f);	    Rx(1, 2) = sin(f);
				Rx(2, 0) = 0;				Rx(2, 1) = -sin(f);		Rx(2, 2) = cos(f);

				Rz(0, 0) = cos(p);		Rz(0, 1) = sin(p);		Rz(0, 2) = 0;
				Rz(1, 0) = -sin(p);		Rz(1, 1) = cos(p);		Rz(1, 2) = 0;
				Rz(2, 0) = 0;				Rz(2, 1) = 0;				Rz(2, 2) = 1;

				CMatrix GeoRealPos(3, 1), GeoPos(3, 1);
				GeoPos(0, 0) = X;
				GeoPos(1, 0) = Y;
				GeoPos(2, 0) = Z;
				GeoRealPos = Rz * Rx * GeoPos;
				SatPosClkMax(satIndex, 0) = GeoRealPos(0, 0);
				SatPosClkMax(satIndex, 1) = GeoRealPos(1, 0);
				SatPosClkMax(satIndex, 2) = GeoRealPos(2, 0);
				/*
				CString sat = obss[epoIndex].data[satIndex].sat;
				double X1 = GeoRealPos(0, 0);
				double Y1 = GeoRealPos(1, 0);
				double Z1 = GeoRealPos(2, 0);
				Z1;
				*/
			}
			
			/* relative effec correction */
			double dtre = (-2 * e * sqrt(GM * A) / (CLIGHT * CLIGHT) * sin(E)) * CLIGHT;
			tk = timediff(ts, eph.toc);
			dts = f0 + f1 * tk + f2 * tk * tk;

			SatPosClkMax(satIndex, 3) = dts * CLIGHT;
			SatPosClkMax(satIndex, 3) -=dtre;
		}

		/* exclude satellites by min elevation*/
		obss[epoIndex].n = obss[epoIndex].nmax;
		for (int satIndex = 0; satIndex < obss[epoIndex].nmax; satIndex++)
		{
			double rr[3], rs[3], vec[3], el, dist,azel[2];
			memcpy(rr, rrTmp, sizeof(rrTmp));
			rs[0] = SatPosClkMax(satIndex, 0);
			rs[1] = SatPosClkMax(satIndex, 1);
			rs[2] = SatPosClkMax(satIndex, 2);
			//CString sat = obss[epoIndex].data[satIndex].sat;
			double SagEffectCor;
			dist = geodist(rr, rs, vec, SagEffectCor);
			satazel(rr, vec, azel);
			el = azel[1];
			if (fabs(el) < minElevAngle(RAD) || (rs[0] == 0 && rs[1] == 0 && rs[2] == 0))
			{
				obss[epoIndex].data[satIndex].svh = 0;
				obss[epoIndex].n = obss[epoIndex].n - 1;
			}
		}
		/* extract valid satellites */
		CMatrix SatPosClkValid(obss[epoIndex].n, 4);
		obss[epoIndex].validIndex = new int[obss[epoIndex].n];
		for (int i = 0, j = 0; i < obss[epoIndex].nmax; i++)
		{
			if (obss[epoIndex].data[i].svh == 0)continue;
			//CString sat = obss[epoIndex].data[i].sat;
			SatPosClkValid(j, 0) = SatPosClkMax(i, 0);
			SatPosClkValid(j, 1) = SatPosClkMax(i, 1);
			SatPosClkValid(j, 2) = SatPosClkMax(i, 2);
			SatPosClkValid(j, 3) = SatPosClkMax(i, 3) ;
			obss[epoIndex].validIndex[j] = i;
			//double t = SatPosClkValid(j, 3);
			j++;
		}
		return SatPosClkValid;
}

void CBdSPP::FormDesConsMatEachEpoch(int epoIndex, CMatrix SatPosClkEpoch,CMatrix& B, CMatrix& L,double*varIon,double*varTrp)
{
		obs_t obs = obss[epoIndex];
		int nValid = SatPosClkEpoch.Row();
		B.SetSize(nValid, 4);
		L.SetSize(nValid, 1);
		double rr[4];
		memcpy(rr, rrTmp, sizeof(rrTmp));
		double rr_blh[3];
		ecef2blh(rr, rr_blh);
		/* form design matrix(B) */
		for(int i=0;i<obs.n;i++)
		{
			int ind = obs.validIndex[i];
			eph_t eph = seleph(obs.data[ind]);
			double rs[3],dts;
			rs[0] = SatPosClkEpoch(i, 0);
			rs[1] = SatPosClkEpoch(i, 1);
			rs[2] = SatPosClkEpoch(i, 2);
			dts= SatPosClkEpoch(i, 3);
			double vec[3];
			CString sat = obs.data[ind].sat;
			/* ρ(geodetic distance) including Sagnac effect correction*/
			double SagEffectCor;
			double rou = geodist(rr, rs, vec,SagEffectCor);
			/* form design matrix(B) */
			B(i, 0) = -vec[0]/ rou;			B(i, 1) = -vec[1] / rou;			B(i, 2) = -vec[2] / rou;
			B(i, 3) = 1;
			rou += SagEffectCor;
			/* form constant matrix(L) */
			double P, TGD, dtrp, dion;
			/* Pseudorange including TGD correction */
			TGD = eph.tgd[0] * CLIGHT;
			P = obs.data[ind].P[0]- TGD;
			/* tropospheric delay */
			double azel[2];
			satazel(rr, vec,azel);
			double el = azel[1];
			dtrp = tropmodel(rr_blh, el, 0.70);
			varTrp[i] = dtrp * 0.3;
			/* Ionosphere delay */
			gtime_t time = obs.data[ind].time;
			double ion[8];
			memcpy(ion, navs.ion_cmp, sizeof(navs.ion_cmp));
			dion = ionmodel(time, ion, rr_blh, azel);
			dion *= pow((FREQ1 / FREQ1_CMP), 2);
			varIon[i] = dion * 0.5;
			double dtr = rr[3];

			double l = P - (rou + dtr - dts + dion + dtrp);
			L(i, 0) = l;
		}
}

void CBdSPP::FormWeightMatEachEpoch(int epoIndex, CMatrix SatPosClkEpoch, CMatrix& P, double* varIon, double* varTrp)
{
	double ura[16] = { 2.4,3.4,4.85,6.85,9.65,13.65,24.0,48.0,96.0,192.0,384.0,768.0,1536.0,3072.0,6144.0,0.0 };
	obs_t obs = obss[epoIndex];
	int nValid = SatPosClkEpoch.Row();
	P.SetSize(nValid, nValid);
	P.Unit();
	double rr[3];
	memcpy(rr, rrTmp, sizeof(rrTmp));
	double rr_blh[3];
	ecef2blh(rr, rr_blh);
	/* form weight matrix(P) */
	for (int i = 0; i < obs.n; i++)
	{
		int ind = obs.validIndex[i];
		eph_t eph = seleph(obs.data[ind]);
		double rs[3], dts;
		rs[0] = SatPosClkEpoch(i, 0);
		rs[1] = SatPosClkEpoch(i, 1);
		rs[2] = SatPosClkEpoch(i, 2);
		dts = SatPosClkEpoch(i, 3);
		double vec[3];
		//CString sat = obs.data[ind].sat;
		double SagEffectCor;
		double rou = geodist(rr, rs, vec, SagEffectCor);
		double azel[2],el;
		satazel(rr, vec, azel);
		el = azel[1];
		double p = 0.3, q = 0.3;
		double sigmaURA, sigmaION, sigmaTRP, sigmaELE;
		sigmaELE = 1.0/(p * p + q * q / sin(el));
		sigmaION = varIon[i];
		sigmaTRP = varTrp[i];
		sigmaURA = ura[eph.sva-1];
		P(i, i) = 1 / (sigmaURA + sigmaION + sigmaTRP + sigmaELE);
	}
}

void CBdSPP::OutputSol(CString SolFile, CString strObsFile, CString strNavFile)
{
	CStdioFile sf;
	if (!sf.Open(SolFile, CFile::modeCreate|| CFile::modeNoTruncate)) return;
	//if (!sf.Open(SolFile, CFile::modeWrite)) return;
	CString strLine;
	// 查找从右向左第一个遇到的反斜杠 '\'
	int pos1 = strObsFile.ReverseFind(_T('\\'));
	CString strObsFileNotPath = strObsFile.Mid(pos1 + 1);
	int pos2 = strNavFile.ReverseFind(_T('\\'));
	CString strNavFileNotPath = strNavFile.Mid(pos2 + 1);

	/* write header*/
	
	strLine.Format(_T("%% obs file  : "));
	strLine += strObsFileNotPath;
	strLine += _T("\r");
	sf.WriteString(strLine);
	strLine.Empty();
	
	strLine.Format(_T("%% nav file  : "));
	strLine += strNavFileNotPath;
	strLine += _T("\r");
	sf.WriteString(strLine);
	strLine.Empty();

	strLine.Format(_T("%% elev mask : %.3f deg\r"), minElevAngle(DEG));
	sf.WriteString(strLine);
	strLine.Empty();
	sf.WriteString(_T("% ionos opt : broadcast\r"));
	sf.WriteString(_T("% tropo opt : saastamoinen\r"));
	sf.WriteString(_T("% ephemeris : broadcast\r"));
	sf.WriteString(_T("% navi sys  : beidou\r"));
	sf.WriteString(_T("%                       GPST      x-ecef(m)      y-ecef(m)      z-ecef(m)     clock bias  ns   sdx(m)   sdy(m)   sdz(m)  sdxy(m)  sdyz(m)  sdzx(m)   pdop   gdop chisqrTest\r"));

	/* write body*/
	for (int epoIndex = 0; epoIndex < EpochCount; epoIndex++)
	{
		CString strTmp;
		double ep[6];
		time2epoch(sol.data[epoIndex].time, ep);
		strLine.Format(_T("%4.0f\t%02.0f\t%02.0f\t%02.0f\t%02.0f\t%02.3f\t"), ep[0], ep[1], ep[2], ep[3], ep[4], ep[5]);
		strTmp.Format(_T("%.4f\t%.4f\t%.4f\t%.4f\t  %d   "), sol.data[epoIndex].rr[0], sol.data[epoIndex].rr[1], sol.data[epoIndex].rr[2], sol.data[epoIndex].dtr, sol.data[epoIndex].ns);
		strLine += strTmp;
		strTmp.Format(_T("%.4f   %.4f\t%.4f\t%.4f   %.4f  %.4f    "), sol.data[epoIndex].qr[0], sol.data[epoIndex].qr[1], sol.data[epoIndex].qr[2], sol.data[epoIndex].qr[3], sol.data[epoIndex].qr[4], sol.data[epoIndex].qr[5]);
		strLine += strTmp;
		strTmp.Format(_T("%.2f   %.2f\t%d\r"), sol.data[epoIndex].pdop, sol.data[epoIndex].gdop, sol.data[epoIndex].chisqrTest);
		strLine += strTmp;
		sf.WriteString(strLine);
		strTmp.Empty();
		strLine.Empty();
	}
	sf.Close();
}

void CBdSPP::mainfunc()
{
	sol.nmax = EpochCount;
	sol.data = new sol_t[EpochCount];
	double chisqr[100] = 
	{ 
		10.8,13.8,16.3,18.5,20.5,22.5,24.3,26.1,27.9,29.6,
		31.3,32.9,34.5,36.1,37.7,39.3,40.8,42.3,43.8,45.3,
		46.8,48.3,49.7,51.2,52.6,54.1,55.5,56.9,58.3,59.7,
		61.1,62.5,63.9,65.2,66.6,68.0,69.3,70.7,72.1,73.4,
		74.7,76.0,77.3,78.6,80.0,81.3,82.6,84.0,85.4,86.7,
		88.0,89.3,90.6,91.9,93.3,94.7,96.0,97.4,98.7,100 ,
		101 ,102 ,103 ,104 ,105 ,107 ,108 ,109 ,110 ,112 ,
		113 ,114 ,115 ,116 ,118 ,119 ,120 ,122 ,123 ,125 ,
		126 ,127 ,128 ,129 ,131 ,132 ,133 ,134 ,135 ,137 ,
		138 ,139 ,140 ,142 ,143 ,144 ,145 ,147 ,148 ,149 
	};
	int epoIndex ;
	for ( epoIndex = 0; epoIndex < EpochCount; epoIndex++)
	{
		sol.data[epoIndex].chisqrTest = 0;
		memcpy(rrTmp, rrInit, sizeof(rrInit));
		CMatrix SatPosClk;
		CMatrix B, BT, L, P, detX, N;
		double norm = 1;
		SatPosClk = CalSatPosClkEachEpoch(epoIndex);
		while (norm > 1e-6)
		{
			double* varTrp = new double[SatPosClk.Row()];
			double* varIon = new double[SatPosClk.Row()];
			FormDesConsMatEachEpoch(epoIndex, SatPosClk, B, L, varIon, varTrp);
			FormWeightMatEachEpoch(epoIndex, SatPosClk, P, varIon, varTrp);
			BT = ~B;
			N = BT * P * B;
			detX = N.Inv() * BT * P * L;
			rrTmp[0] += detX(0, 0);
			rrTmp[1] += detX(1, 0);
			rrTmp[2] += detX(2, 0);
			rrTmp[3] += detX(3, 0);
			norm = pow(detX(0, 0), 2) + pow(detX(1, 0), 2) + pow(detX(2, 0), 2) + pow(detX(3, 0), 2);
		}
		/* evaluate result */
		CMatrix Q = N.Inv();
		double gdop = sqrt(Q(0, 0) + Q(1, 1) + Q(2, 2) + Q(3, 3));
		double pdop = sqrt(Q(0, 0) + Q(1, 1) + Q(2, 2));
		CMatrix V = B * detX - L;
		double sigma0 = sqrt((~V * P * V)(0, 0) / (V.Row() - 4));
		double chisqrTest = (~V * P * V)(0, 0) / (V.Row() - 5);
		/* Chi-square Test */
		if (chisqrTest < chisqr[V.Row() - 6])
			sol.data[epoIndex].chisqrTest = 1;
		double sdx, sdy, sdz, sdxy, sdyz, sdzx;
		//PrintMatrix(Q, strMatFile);
		sdx = sigma0 * sigma0 * Q(0, 0);
		sdy = sigma0 * sigma0 * Q(1, 1);
		sdz = sigma0 * sigma0 * Q(2, 2);
		sdxy = sigma0 * sigma0 * Q(0, 1);
		sdyz = sigma0 * sigma0 * Q(1, 2);
		sdzx = sigma0 * sigma0 * Q(2, 0);
		/* restore result */
		for (int i = 0; i < 3; i++)
			sol.data[epoIndex].rr[i] = rrTmp[i];
		sol.data[epoIndex].dtr = rrTmp[3] ;
		sol.data[epoIndex].qr[0] = sdx;
		sol.data[epoIndex].qr[1] = sdy;
		sol.data[epoIndex].qr[2] = sdz;
		sol.data[epoIndex].qr[3] = sdxy;
		sol.data[epoIndex].qr[4] = sdyz;
		sol.data[epoIndex].qr[5] = sdzx;
		sol.data[epoIndex].time = obss[epoIndex].data[0].time;
		sol.data[epoIndex].gdop = gdop;
		sol.data[epoIndex].pdop = pdop;
		sol.data[epoIndex].ns = obss[epoIndex].n;
	}
	delete[]obss;
	obss = NULL;
	delete[]navs.eph;
	navs.eph = NULL;
}