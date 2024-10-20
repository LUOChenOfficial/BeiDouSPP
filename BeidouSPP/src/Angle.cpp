#include "pch.h"
#include "Angle.h"
#include "math.h"
const double EPSILON = 1.0E-12;
const double PI = 4.0 * atan(1.0);

//���ع��캯������ȱʡֵ
CAngle::CAngle(double value, AngleStyle style)
{
	dValue = value;
	nCurStyle = style;
}

CAngle::~CAngle(void)
{
}
//���أ�������
double& CAngle::operator() (AngleStyle style) //ָ�������ͻ�ȡ�Ƕ�ֵ
{
	//double dAngleValue;
	if (style == DMS)
	{
		if (nCurStyle == DEG)
		{
			dValue = Dms(dValue);
		}
		else if (nCurStyle == RAD)
		{
			dValue = Dms(dValue * 180.0 / PI);
		}
		nCurStyle = DMS;

	}
	else if (style == DEG)
	{
		if (nCurStyle == DMS)
		{
			dValue = Deg(dValue);
		}
		else if (nCurStyle == RAD)
		{
			dValue = dValue * 180.0 / PI;
		}
		nCurStyle = DEG;
	}
	else
	{
		if (nCurStyle == DMS)
		{
			dValue = Deg(dValue) * PI / 180;
		}
		else if (nCurStyle == DEG)
		{
			dValue = dValue * PI / 180;
		}
		nCurStyle = RAD;
	}
	return dValue;
}
//���أ����������ú����ǳ�������ֻ�ܱ���CAngle����ʹ��
double CAngle::operator() (AngleStyle style) const //ָ�������ͻ�ȡ�Ƕ�ֵ
{
	double dAngleValue;
	if (style == DMS)
	{
		if (nCurStyle == DEG)
		{
			dAngleValue = Dms(dValue);
		}
		else if (nCurStyle == RAD)
		{
			dAngleValue = Dms(dValue * 180.0 / PI);
		}
		else
		{
			dAngleValue = dValue;
		}

	}
	else if (style == DEG)
	{
		if (nCurStyle == DMS)
		{
			dAngleValue = Deg(dValue);
		}
		else if (nCurStyle == RAD)
		{
			dAngleValue = dValue * 180.0 / PI;
		}
		else
		{
			dAngleValue = dValue;
		}
	}
	else
	{
		if (nCurStyle == DMS)
		{
			dAngleValue = Deg(dValue) * PI / 180;
		}
		else if (nCurStyle == DEG)
		{
			dAngleValue = dValue * PI / 180;
		}
		else
		{
			dAngleValue = dValue;
		}
	}
	return dAngleValue;
}


//˽�г�Ա���ȷ�����ʮ���ƶ�ת��
double CAngle::Deg(double dDms) const
{
	int iDeg, iMin;
	double dSec;

	iDeg = int(dDms + EPSILON);//��//��һ����С�������Է�ֹȡ��ʱ�ĳ���
	iMin = int((dDms - iDeg) * 100 + EPSILON);//��
	dSec = ((dDms - iDeg) * 100 - iMin) * 100;//��
	return iDeg + (double)iMin / 60 + dSec / 3600;
}

//˽�г�Ա��ʮ���ƶ���ȷ���ת��
double CAngle::Dms(double dDeg) const
{
	int iDeg, iMin;
	double dSec;
	double dTmp;

	iDeg = int(dDeg + EPSILON);//�������ֶ�
	dTmp = (dDeg - iDeg) * 60;//С������ת���ɷ�
	iMin = int(dTmp + EPSILON);//ȡ�ֵ���������
	dSec = (dTmp - iMin) * 60;//��ȡ��

	return iDeg + (double)iMin / 100 + dSec / 10000;
}

//��Ԫ����+����
CAngle operator + (const CAngle& m1, const CAngle& m2)
{
	CAngle addAngle(0, RAD);
	addAngle(RAD) = m1(RAD) + m2(RAD);
	return addAngle;
}
//��Ԫ����-����
CAngle operator - (const CAngle& m1, const CAngle& m2)
{
	CAngle subAngle(0, RAD);
	subAngle(RAD) = m1(RAD) - m2(RAD);
	return subAngle;
}