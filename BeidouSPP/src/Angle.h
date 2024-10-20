#pragma once

//ö���������ͣ����ڴ���Ƕ���ʽ
enum AngleStyle
{
	DEG,
	DMS,
	RAD
};
class CAngle
{
public:
	CAngle(double value = 0, AngleStyle style = DMS);
	~CAngle(void);
private:
	double dValue;//�Ƕ�ֵ
	AngleStyle  nCurStyle;//��ǰ�Ƕ�ֵ����
private:
	//���ó���Ա���������ã�1.���Ա���ᱻ�ı�
	//2.���Ա������������
	double Deg(double dDms) const;
	double Dms(double dDeg) const;

public:
	//��ȡָ�������ͻ�ȡ�Ƕ�ֵ��
	//���ڷ��ص���dValue�����ã����Ը�ֵ��С���Ըı䣬�����Խ��и�ֵ
	double& operator() (AngleStyle style);   //CAanle ang1,   ang1.operator(DEG); angl1(DEG)

	//���أ���ȡָ�������ͻ�ȡ�Ƕ�ֵ����ֵ���ɸı䣬const CAngle���ͱ�������
	double operator() (AngleStyle style) const;
	//���������+/-
	friend CAngle operator + (const CAngle& m1, const CAngle& m2);
	friend CAngle operator - (const CAngle& m1, const CAngle& m2);
};


