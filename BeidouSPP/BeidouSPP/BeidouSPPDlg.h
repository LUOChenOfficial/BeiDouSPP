
// BeidouSPPDlg.h: 头文件
//

#pragma once
#include "BdSPP.h"

// CBeidouSPPDlg 对话框
class CBeidouSPPDlg : public CDialogEx
{
// 构造
public:
	CBeidouSPPDlg(CWnd* pParent = nullptr);	// 标准构造函数

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_BEIDOUSPP_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;
	CBdSPP m_CBdSPP;
	CString ObsFile, NavFile;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnBnClickedReadobsdata();
	afx_msg void OnBnClickedReadnavdata();
//	afx_msg void OnBnClickedCalsatpos();
	afx_msg void OnBnClickedPosfunc();
	afx_msg void OnBnClickedOutputsol();
};
