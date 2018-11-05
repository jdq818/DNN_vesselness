
#include "vtkPoints.h"
#include "vtkLine.h"
#include "vtkPolyLine.h"


#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkSmartPointer.h"


#include "vtkUnstructuredGrid.h"

// for read
#include "vtkUnstructuredGridReader.h"

// for write
#include "vtkUnstructuredGridWriter.h"

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "zxhImageData.h"
#include "zxhImageGipl.h"
#include "zxhImageNifti.h"
#include "miiMinHeap.h"

#define TAB_CHAR	9
#define SPA_CHAR	1
using namespace std;

typedef struct
{
	float x;
	float y;
	float z;
	void cvPoint(float px,float py,float pz);
	int n;
}PointCordTypeDef;
void PointCordTypeDef::cvPoint(float px,float py,float pz)
{
	x=px;
	y=py;
	z=pz;
};
typedef struct
{
	float fMean;
	float fStdVar;
	float fMaxErr;
}ErrTypeDef;

void ReadVtk(char *chFileName, vector<PointCordTypeDef> &PointCord)
{
	if (chFileName == NULL)
	{
		cout << "Cannot find VTK-file!" << endl;
		return;
	}
	if (!PointCord.empty())
	{
		PointCord.clear();
	}

	vtkSmartPointer<vtkUnstructuredGridReader> iVtkReader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	iVtkReader->SetFileName( chFileName );
	iVtkReader->Update();

	vtkSmartPointer<vtkUnstructuredGrid> iGridRead = iVtkReader->GetOutput();

	int nPointNum = iGridRead->GetMaxCellSize();

	double dCord[3];
	PointCordTypeDef strctTempPoint;

	for (int i = 0; i < nPointNum; i++)
	{
		iGridRead->GetPoint(i, dCord);
		strctTempPoint.x = dCord[0];
		strctTempPoint.y = dCord[1];
		strctTempPoint.z = dCord[2];
		PointCord.push_back(strctTempPoint);
	}
}
bool ReadTxt(char *chFileName,vector<short>&sMolPontInts)
{
	ifstream inMolIntdata(chFileName);
	double data;string strLine;
	bool bsizeeque=true;
	if(!inMolIntdata)
	{
		cout << "Unable to open txt-file: " << endl;
		cout << chFileName << endl;
		return false; //exit(1); // terminate with error
	}
	else
	{

		while(getline(inMolIntdata,strLine,'\n'))
		{ 
			data = atof(strLine.c_str());
			sMolPontInts.push_back(data);
		}
		return true;

	}

}
bool ReadTxt_Pdef(char *chFileName,vector<PointCordTypeDef>&pPont)
{
	string strNum;
	int nStart = 0, nEnd = 0;
	PointCordTypeDef strctTempPoint;
	ifstream iFileIn(chFileName);
	if(!iFileIn)
	{
		cout << "Unable to open txt-file: " << endl;
		cout << chFileName << endl;
		return false; //exit(1); // terminate with error
	}
	string strLine;
	while(getline(iFileIn, strLine))
	{
		nStart=0;
		// read 1st number
		nEnd = strLine.find_first_of(' ',nStart);
		strNum.assign(strLine, nStart, nEnd - nStart);
		strctTempPoint.x =atof(strNum.c_str());
		//fMean=1024*fMean/1896;//Add by JDQ
		// read 2nd number
		nStart = nEnd + 1;
		nEnd = strLine.find_first_of(' ', nStart);
		strNum.assign(strLine, nStart, nEnd - nStart);
		strctTempPoint.y =atof(strNum.c_str());
		// read 3rd number
		nStart = nEnd + 1;
		nEnd = strLine.find_first_of(' ', nStart);
		strNum.assign(strLine, nStart, nEnd - nStart);
		strctTempPoint.z =atof(strNum.c_str());
		strctTempPoint.x=-1*strctTempPoint.x;
		strctTempPoint.y=-1*strctTempPoint.y;
		pPont.push_back(strctTempPoint);
	}
	return true;
}

void WriteVtk(vector< PointCordTypeDef > PointCord, char* chFileName)
{
	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();

	/*	int nPointNum = PointCord.size();*/

	int nPointNum = PointCord.size();

	for (int i = 0; i < nPointNum; i++)
	{
		iPoints->InsertNextPoint(PointCord[i].x, PointCord[i].y, PointCord[i].z);
	}	

	vtkSmartPointer<vtkPolyLine> iLine = vtkSmartPointer<vtkPolyLine>::New();
	iLine->GetPointIds()->SetNumberOfIds(nPointNum);
	for (int i = 0; i < nPointNum; i++)
	{
		iLine->GetPointIds()->SetId(i, i);
	}

	vtkSmartPointer<vtkUnstructuredGrid> iGrid = vtkUnstructuredGrid::New();
	iGrid->Allocate(1, 1);	
	iGrid->InsertNextCell(iLine->GetCellType(), iLine->GetPointIds());
	iGrid->SetPoints(iPoints);

	vtkSmartPointer<vtkUnstructuredGridWriter> iVtkWriter = vtkUnstructuredGridWriter::New();
	iVtkWriter->SetInput(iGrid);
	iVtkWriter->SetFileName(chFileName);
	iVtkWriter->Write();
}
void WriteTxt(vector< PointCordTypeDef > PointCord, char* chFileName)
{
	ofstream WriteFileTxt(chFileName);
	int nPointNum = PointCord.size();
	float fImgPixel[3];	
	for (int i = 0; i < nPointNum; i++)
	{
		float fx=-1*PointCord[i].x;
		float fy=-1*PointCord[i].y;
		float fz=PointCord[i].z;
		WriteFileTxt <<right<<fixed<<setfill('4')<<setprecision(4)<<fx<<' '<<fy<<' '<<fz<<'\n';
	}	
}
void WriteCAIntTxt(vector<short> sMolPontInts,char* chFileName)
{
	ofstream WriteFileTxt(chFileName,ios::out);
	int nPointNum = sMolPontInts.size();
	for (int i = 0; i < nPointNum; i++)
	{
		WriteFileTxt<<right<<fixed<<setfill('0')<<setprecision(4) <<sMolPontInts[i] <<"\n";

	}	
}
bool ReadPointTxt(char *filename,vector< miiCNode<double, float> > &cl)
{
	string strNum;
	int nStart = 0, nEnd = 0;
	miiCNode<double, float> strctTempPoint;
	ifstream iFileIn(filename);
	if(!iFileIn)
	{
		cout << "Unable to open txt-file: " << endl;
		cout << filename << endl;
		return false; //exit(1); // terminate with error
	}
	string strLine;
	if (!getline(iFileIn, strLine))
	{
		cout << "The format of Mean and StdDev file is wrong!" << endl;
		return false;
	}
	//// read 1st number
	//nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	//strNum.assign(strLine, nStart, nEnd - nStart);
	//strctTempPoint.x =atof(strNum.c_str());
	////fMean=1024*fMean/1896;//Add by JDQ
	//// read 2nd number
	//nStart = nEnd + 1;
	//nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	//strNum.assign(strLine, nStart, nEnd - nStart);
	//strctTempPoint.y =atof(strNum.c_str());
	//// read 3rd number
	//nStart = nEnd + 1;
	//nEnd = strLine.find_first_of(TAB_CHAR, nStart);
	//strNum.assign(strLine, nStart, nEnd - nStart);
	//strctTempPoint.z =atof(strNum.c_str());
	//cl.push_back(strctTempPoint);

	ifstream ifs(filename);

	char *buffer=new char[100];int ibuffersize=100;
	float s;
	std::string sComment,sContent,sLine;
	while(ifs.eof()==false&&ifs.fail()==false)
	{
		ifs.getline(buffer,ibuffersize);
		sLine=buffer;
		zxh::ParseStringLine(sContent,sComment,sLine);
		zxh::trim_both(sContent);

		if(sContent.length()!=0)
		{
			std::istringstream istr;
			istr.str(sContent);
			float a,b,c;
			istr>>a>>b>>c;
			strctTempPoint.x=-a;
			strctTempPoint.y=-b;
			strctTempPoint.z=c;
		};
	}
	cl.push_back(strctTempPoint);
	return true;
}
float bezier2funcX(float uu,PointCordTypeDef *controlP){  
	float part0 = controlP[0].x * uu * uu;  
	float part1 = 2 * controlP[1].x * uu *(1 - uu);  
	float part2 = controlP[2].x * (1 - uu) * (1 - uu);     
	return part0 + part1 + part2 ;   
} 
float bezier2funcY(float uu,PointCordTypeDef *controlP){  
	float part0 = controlP[0].y * uu * uu;  
	float part1 = 2 * controlP[1].y * uu *(1 - uu);  
	float part2 = controlP[2].y * (1 - uu) * (1 - uu);     
	return part0 + part1 + part2 ;   
}  
float bezier2funcZ(float uu,PointCordTypeDef *controlP){  
	float part0 = controlP[0].z * uu * uu;  
	float part1 = 2 * controlP[1].z * uu *(1 - uu);  
	float part2 = controlP[2].z * (1 - uu) * (1 - uu);     
	return part0 + part1 + part2 ;   
}  
float bezier3funcX(float uu,PointCordTypeDef *controlP){  
	float part0 = controlP[0].x * uu * uu * uu;  
	float part1 = 3 * controlP[1].x * uu * uu * (1 - uu);  
	float part2 = 3 * controlP[2].x * uu * (1 - uu) * (1 - uu);  
	float part3 = controlP[3].x * (1 - uu) * (1 - uu) * (1 - uu);     
	return part0 + part1 + part2 + part3;   
}      
float bezier3funcY(float uu,PointCordTypeDef *controlP){  
	float part0 = controlP[0].y * uu * uu * uu;  
	float part1 = 3 * controlP[1].y * uu * uu * (1 - uu);  
	float part2 = 3 * controlP[2].y * uu * (1 - uu) * (1 - uu);  
	float part3 = controlP[3].y * (1 - uu) * (1 - uu) * (1 - uu);     
	return part0 + part1 + part2 + part3;   
}
float bezier3funcZ(float uu,PointCordTypeDef *controlP){  
	float part0 = controlP[0].z * uu * uu * uu;  
	float part1 = 3 * controlP[1].z * uu * uu * (1 - uu);  
	float part2 = 3 * controlP[2].z * uu * (1 - uu) * (1 - uu);  
	float part3 = controlP[3].z * (1 - uu) * (1 - uu) * (1 - uu);     
	return part0 + part1 + part2 + part3;   
}  
void createCurve0(PointCordTypeDef originPoint[4],int originCount,vector<PointCordTypeDef> &curvePoint)
{  
	//生成4控制点，产生贝塞尔曲线  
	PointCordTypeDef controlPoint2[4];
	int i=0;
	controlPoint2[0] = originPoint[i];  
	controlPoint2[1] = originPoint[i+1];
	controlPoint2[2] = originPoint[i+2];
	controlPoint2[3] =originPoint[i+3];
	float u = 1;  
	while(u >= 0){  
		float px = bezier3funcX(u,controlPoint2);  
		float py = bezier3funcY(u,controlPoint2);
		float pz= bezier3funcZ(u,controlPoint2);
		//u的步长决定曲线的疏密  
		u -= 0.005;  
		PointCordTypeDef tempP;
		tempP.cvPoint(px,py,pz);
		//存入曲线点   
		curvePoint.push_back(tempP);  
	}     
}  
void createCurve1(PointCordTypeDef originPoint[4],int originCount,vector<PointCordTypeDef> &curvePoint)
{  
	//控制点收缩系数 ，经调试0.6较好
	float scale = 0.6;  
	PointCordTypeDef *midpoints=new PointCordTypeDef[originCount];  
	//生成中点       
	for(int i = 0 ;i < originCount-1 ; i++){      
		int nexti = i + 1;  
		midpoints[i].x = (originPoint[i].x + originPoint[nexti].x)/2.0;  
		midpoints[i].y = (originPoint[i].y + originPoint[nexti].y)/2.0;  
		midpoints[i].z = (originPoint[i].z + originPoint[nexti].z)/2.0;
	}      
	//平移中点  
	PointCordTypeDef *extrapoints=new PointCordTypeDef[2 * originCount];   
	for(int i = 1 ;i < originCount-1 ; i++)
	{  
		int nexti = i + 1;  
		int backi = i - 1;  
		PointCordTypeDef midinmid;  
		midinmid.x = (midpoints[i].x + midpoints[backi].x)/2.0;  
		midinmid.y = (midpoints[i].y + midpoints[backi].y)/2.0;  
		midinmid.z = (midpoints[i].z + midpoints[backi].z)/2.0;
		float offsetx = originPoint[i].x - midinmid.x;  
		float offsety = originPoint[i].y - midinmid.y;
		float offsetz = originPoint[i].z - midinmid.z;
		int extraindex = 2 * i;  
		extrapoints[extraindex].x = midpoints[backi].x + offsetx;  
		extrapoints[extraindex].y = midpoints[backi].y + offsety;
		extrapoints[extraindex].z = midpoints[backi].z + offsetz;  
		//朝 originPoint[i]方向收缩   
		float addx = (extrapoints[extraindex].x - originPoint[i].x) * scale;  
		float addy = (extrapoints[extraindex].y - originPoint[i].y) * scale;
		float addz = (extrapoints[extraindex].z - originPoint[i].z) * scale;  
		extrapoints[extraindex].x = originPoint[i].x + addx;  
		extrapoints[extraindex].y = originPoint[i].y + addy;  
		extrapoints[extraindex].z = originPoint[i].z + addz;  
		int extranexti = extraindex + 1;  
		extrapoints[extranexti].x = midpoints[i].x + offsetx;  
		extrapoints[extranexti].y = midpoints[i].y + offsety;
		extrapoints[extranexti].z = midpoints[i].z + offsetz;
		//朝 originPoint[i]方向收缩   
		addx = (extrapoints[extranexti].x - originPoint[i].x) * scale;  
		addy = (extrapoints[extranexti].y - originPoint[i].y) * scale;
		addz = (extrapoints[extranexti].z - originPoint[i].z) * scale;  
		extrapoints[extranexti].x = originPoint[i].x + addx;  
		extrapoints[extranexti].y = originPoint[i].y + addy;  
		extrapoints[extranexti].z = originPoint[i].z + addz;  
	}      

	PointCordTypeDef controlPoint1[3];  
	//生成3控制点，产生第一段贝塞尔曲线  
	int i1=0;
	controlPoint1[0] = originPoint[i1];  
	int extraindex1 = i1+2;  
	controlPoint1[1] = extrapoints[extraindex1];  
	int nexti1 = i1 + 1;  
	controlPoint1[2] = originPoint[nexti1];      
	float u1 = 1;  
	while(u1 >= 0)
	{  
		float px = bezier2funcX(u1,controlPoint1);  
		float py = bezier2funcY(u1,controlPoint1);
		float pz= bezier2funcZ(u1,controlPoint1);
		//u的步长决定曲线的疏密  
		u1 -= 0.005;  
		PointCordTypeDef tempP;
		tempP.cvPoint(px,py,pz);
		//存入曲线点   
		curvePoint.push_back(tempP);  
	}      

	PointCordTypeDef controlPoint2[4];  
	//生成4控制点，产生贝塞尔曲线  
	for(int i = 1 ;i < originCount-2 ; i++){  
		controlPoint2[0] = originPoint[i];  
		int extraindex = 2 * i;  
		controlPoint2[1] = extrapoints[extraindex + 1];  
		int extranexti = extraindex + 2;  
		controlPoint2[2] = extrapoints[extranexti];  
		int nexti = i + 1;  
		controlPoint2[3] = originPoint[nexti];      
		float u = 1;  
		while(u >= 0){  
			float px = bezier3funcX(u,controlPoint2);  
			float py = bezier3funcY(u,controlPoint2);
			float pz= bezier3funcZ(u,controlPoint2);
			//u的步长决定曲线的疏密  
			u -= 0.005;  
			PointCordTypeDef tempP;
			tempP.cvPoint(px,py,pz);
			//存入曲线点   
			curvePoint.push_back(tempP);  
		}      
	}  
	PointCordTypeDef controlPoint3[3];  
	//生成3控制点，产生最后一段贝塞尔曲线  
	int i3=originCount-2;
	controlPoint3[0] = originPoint[i3];  
	int extraindex3 = i3+3;  
	controlPoint3[1] = extrapoints[extraindex3];  
	int nexti3 = i3 + 1;  
	controlPoint3[2] = originPoint[nexti3];      
	float u3 = 1;  
	while(u3 >= 0){  
		float px = bezier2funcX(u3,controlPoint3);  
		float py = bezier2funcY(u3,controlPoint3);
		float pz= bezier2funcZ(u3,controlPoint3);
		//u的步长决定曲线的疏密  
		u3 -= 0.005;  
		PointCordTypeDef tempP;
		tempP.cvPoint(px,py,pz);
		//存入曲线点   
		curvePoint.push_back(tempP);  
	}    
	delete[]midpoints;
	delete[]extrapoints;
}  
void createCurve2(PointCordTypeDef originPoint[4],int originCount,vector<PointCordTypeDef> &vcurvePoint,vector<vector<PointCordTypeDef>> &vvcurvePoint)//这个函数在1的基础上分段保存曲线
{  
	//控制点收缩系数 ，经调试0.6较好
	float scale = 0.6;  
	PointCordTypeDef *midpoints=new PointCordTypeDef[originCount];  
	//生成中点       
	for(int i = 0 ;i < originCount-1 ; i++){      
		int nexti = i + 1;  
		midpoints[i].x = (originPoint[i].x + originPoint[nexti].x)/2.0;  
		midpoints[i].y = (originPoint[i].y + originPoint[nexti].y)/2.0;  
		midpoints[i].z = (originPoint[i].z + originPoint[nexti].z)/2.0;
	}      
	//平移中点  
	PointCordTypeDef *extrapoints=new PointCordTypeDef[2 * originCount];   
	for(int i = 1 ;i < originCount-1 ; i++)
	{  
		int nexti = i + 1;  
		int backi = i - 1;  
		PointCordTypeDef midinmid;  
		midinmid.x = (midpoints[i].x + midpoints[backi].x)/2.0;  
		midinmid.y = (midpoints[i].y + midpoints[backi].y)/2.0;  
		midinmid.z = (midpoints[i].z + midpoints[backi].z)/2.0;
		float offsetx = originPoint[i].x - midinmid.x;  
		float offsety = originPoint[i].y - midinmid.y;
		float offsetz = originPoint[i].z - midinmid.z;
		int extraindex = 2 * i;  
		extrapoints[extraindex].x = midpoints[backi].x + offsetx;  
		extrapoints[extraindex].y = midpoints[backi].y + offsety;
		extrapoints[extraindex].z = midpoints[backi].z + offsetz;  
		//朝 originPoint[i]方向收缩   
		float addx = (extrapoints[extraindex].x - originPoint[i].x) * scale;  
		float addy = (extrapoints[extraindex].y - originPoint[i].y) * scale;
		float addz = (extrapoints[extraindex].z - originPoint[i].z) * scale;  
		extrapoints[extraindex].x = originPoint[i].x + addx;  
		extrapoints[extraindex].y = originPoint[i].y + addy;  
		extrapoints[extraindex].z = originPoint[i].z + addz;  
		int extranexti = extraindex + 1;  
		extrapoints[extranexti].x = midpoints[i].x + offsetx;  
		extrapoints[extranexti].y = midpoints[i].y + offsety;
		extrapoints[extranexti].z = midpoints[i].z + offsetz;
		//朝 originPoint[i]方向收缩   
		addx = (extrapoints[extranexti].x - originPoint[i].x) * scale;  
		addy = (extrapoints[extranexti].y - originPoint[i].y) * scale;
		addz = (extrapoints[extranexti].z - originPoint[i].z) * scale;  
		extrapoints[extranexti].x = originPoint[i].x + addx;  
		extrapoints[extranexti].y = originPoint[i].y + addy;  
		extrapoints[extranexti].z = originPoint[i].z + addz;  
	}      

	PointCordTypeDef controlPoint1[3];  
	vvcurvePoint.clear();
	//生成3控制点，产生第一段贝塞尔曲线  
	int i1=0;
	controlPoint1[0] = originPoint[i1];  
	int extraindex1 = i1+2;  
	controlPoint1[1] = extrapoints[extraindex1];  
	int nexti1 = i1 + 1;  
	controlPoint1[2] = originPoint[nexti1];      
	float u1 = 1;  
	while(u1 >= 0)
	{  
		float px = bezier2funcX(u1,controlPoint1);  
		float py = bezier2funcY(u1,controlPoint1);
		float pz= bezier2funcZ(u1,controlPoint1);
		//u的步长决定曲线的疏密  
		u1 -= 0.005;  
		PointCordTypeDef tempP;
		tempP.cvPoint(px,py,pz);
		//存入曲线点   
		vcurvePoint.push_back(tempP);  
	}   
	vvcurvePoint.push_back(vcurvePoint);
	vcurvePoint.clear();
	PointCordTypeDef controlPoint2[4];  
	//生成4控制点，产生贝塞尔曲线  
	for(int i = 1 ;i < originCount-2 ; i++){  
		controlPoint2[0] = originPoint[i];  
		int extraindex = 2 * i;  
		controlPoint2[1] = extrapoints[extraindex + 1];  
		int extranexti = extraindex + 2;  
		controlPoint2[2] = extrapoints[extranexti];  
		int nexti = i + 1;  
		controlPoint2[3] = originPoint[nexti];      
		float u = 1;  
		while(u >= 0){  
			float px = bezier3funcX(u,controlPoint2);  
			float py = bezier3funcY(u,controlPoint2);
			float pz= bezier3funcZ(u,controlPoint2);
			//u的步长决定曲线的疏密  
			u -= 0.005;  
			PointCordTypeDef tempP;
			tempP.cvPoint(px,py,pz);
			//存入曲线点   
			vcurvePoint.push_back(tempP);  
		}      
	}  
	vvcurvePoint.push_back(vcurvePoint);
		vcurvePoint.clear();


	PointCordTypeDef controlPoint3[3];  
	//生成3控制点，产生最后一段贝塞尔曲线  
	int i3=originCount-2;
	controlPoint3[0] = originPoint[i3];  
	int extraindex3 = i3+3;  
	controlPoint3[1] = extrapoints[extraindex3];  
	int nexti3 = i3 + 1;  
	controlPoint3[2] = originPoint[nexti3];      
	float u3 = 1;  
	while(u3 >= 0){  
		float px = bezier2funcX(u3,controlPoint3);  
		float py = bezier2funcY(u3,controlPoint3);
		float pz= bezier2funcZ(u3,controlPoint3);
		//u的步长决定曲线的疏密  
		u3 -= 0.005;  
		PointCordTypeDef tempP;
		tempP.cvPoint(px,py,pz);
		//存入曲线点   
		vcurvePoint.push_back(tempP);  
	}    
	vvcurvePoint.push_back(vcurvePoint);
		vcurvePoint.clear();
	delete[]midpoints;
	delete[]extrapoints;
}  
//void createCurve(PointCordTypeDef originPoint[4],int originCount,vector<PointCordTypeDef> &curvePoint)
//{  
//    //控制点收缩系数 ，经调试0.6较好，CvPoint是opencv的，可自行定义结构体(x,y)  
//    float scale = 0.6;  
//    PointCordTypeDef *midpoints=new PointCordTypeDef[originCount];  
//    //生成中点       
//    for(int i = 0 ;i < originCount ; i++){      
//        int nexti = (i + 1) % originCount;  
//        midpoints[i].x = (originPoint[i].x + originPoint[nexti].x)/2.0;  
//        midpoints[i].y = (originPoint[i].y + originPoint[nexti].y)/2.0;  
//    }      
//      
//    //平移中点  
//    PointCordTypeDef *extrapoints=new PointCordTypeDef[2 * originCount];   
//    for(int i = 0 ;i < originCount ; i++)
//	{  
//         int nexti = (i + 1) % originCount;  
//         int backi = (i + originCount - 1) % originCount;  
//         PointCordTypeDef midinmid;  
//         midinmid.x = (midpoints[i].x + midpoints[backi].x)/2.0;  
//         midinmid.y = (midpoints[i].y + midpoints[backi].y)/2.0;  
//         int offsetx = originPoint[i].x - midinmid.x;  
//         int offsety = originPoint[i].y - midinmid.y;  
//         int extraindex = 2 * i;  
//         extrapoints[extraindex].x = midpoints[backi].x + offsetx;  
//         extrapoints[extraindex].y = midpoints[backi].y + offsety;  
//         //朝 originPoint[i]方向收缩   
//         int addx = (extrapoints[extraindex].x - originPoint[i].x) * scale;  
//         int addy = (extrapoints[extraindex].y - originPoint[i].y) * scale;  
//         extrapoints[extraindex].x = originPoint[i].x + addx;  
//         extrapoints[extraindex].y = originPoint[i].y + addy;  
//           
//         int extranexti = (extraindex + 1)%(2 * originCount);  
//         extrapoints[extranexti].x = midpoints[i].x + offsetx;  
//         extrapoints[extranexti].y = midpoints[i].y + offsety;  
//         //朝 originPoint[i]方向收缩   
//         addx = (extrapoints[extranexti].x - originPoint[i].x) * scale;  
//         addy = (extrapoints[extranexti].y - originPoint[i].y) * scale;  
//         extrapoints[extranexti].x = originPoint[i].x + addx;  
//         extrapoints[extranexti].y = originPoint[i].y + addy;  
//           
//    }      
//      
//    PointCordTypeDef controlPoint[4];  
//    //生成4控制点，产生贝塞尔曲线  
//    for(int i = 0 ;i < originCount ; i++){  
//           controlPoint[0] = originPoint[i];  
//           int extraindex = 2 * i;  
//           controlPoint[1] = extrapoints[extraindex + 1];  
//           int extranexti = (extraindex + 2) % (2 * originCount);  
//           controlPoint[2] = extrapoints[extranexti];  
//           int nexti = (i + 1) % originCount;  
//           controlPoint[3] = originPoint[nexti];      
//           float u = 1;  
//           while(u >= 0){  
//               int px = bezier3funcX(u,controlPoint);  
//               int py = bezier3funcY(u,controlPoint);  
//               //u的步长决定曲线的疏密  
//               u -= 0.005;  
//               PointCordTypeDef tempP = cvPoint(px,py,);  
//               //存入曲线点   
//               curvePoint.push_back(tempP);  
//           }      
//    }  
//}  
//三次贝塞尔曲线  
void fitbesscurvUSE4points(PointCordTypeDef originPoint[4],vector<PointCordTypeDef>&vcurvePoint)
{
	//order the A B
	float fPointS[3]={originPoint[0].x,originPoint[0].y,originPoint[0].z};
	float fPointA[3]={originPoint[1].x,originPoint[1].y,originPoint[1].z};
	float fPointB[3]={originPoint[2].x,originPoint[2].y,originPoint[2].z};
	float fSADist=zxh::VectorOP_Distance(fPointS,fPointA,3);
	float fSBDist=zxh::VectorOP_Distance(fPointS,fPointB,3);
	if(fSADist>fSBDist)
	{
		PointCordTypeDef PTemp=originPoint[1];
		originPoint[1]=originPoint[2];
		originPoint[2]=PTemp;
		cout<<"Current Order is S B A E;"<<endl;
	}
	int originCount=4;
	//createCurve1(originPoint,originCount,curvePoint);
	
	//createCurve2(originPoint,originCount,vcurvePoint);


	//createCurve0(originPoint,originCount,curvePoint);
	//vector<miiCNode<double,float>>vCurvePointWorld;
	//miiCNode<double,float> ptemp;
	//for (int i=0;i<curvePoint.size();i++)
	//{
	//	ptemp.x=curvePoint[i].x;
	//	ptemp.y=curvePoint[i].y;
	//	ptemp.z=curvePoint[i].z;
	//	vCurvePointWorld.push_back(ptemp);
	//}

}
void fitbesscurvUSE4points1(PointCordTypeDef originPoint[4],vector<PointCordTypeDef>&vcurvePoint,vector<vector<PointCordTypeDef>> &vvcurvePoint)
{
	//order the A B
	float fPointS[3]={originPoint[0].x,originPoint[0].y,originPoint[0].z};
	float fPointA[3]={originPoint[1].x,originPoint[1].y,originPoint[1].z};
	float fPointB[3]={originPoint[2].x,originPoint[2].y,originPoint[2].z};
	float fSADist=zxh::VectorOP_Distance(fPointS,fPointA,3);
	float fSBDist=zxh::VectorOP_Distance(fPointS,fPointB,3);
	if(fSADist>fSBDist)
	{
		PointCordTypeDef PTemp=originPoint[1];
		originPoint[1]=originPoint[2];
		originPoint[2]=PTemp;
		cout<<"Current Order is S B A E;"<<endl;
	}
	int originCount=4;
	//createCurve1(originPoint,originCount,curvePoint);
	
	createCurve2(originPoint,originCount,vcurvePoint,vvcurvePoint);


	//createCurve0(originPoint,originCount,curvePoint);
	//vector<miiCNode<double,float>>vCurvePointWorld;
	//miiCNode<double,float> ptemp;
	//for (int i=0;i<curvePoint.size();i++)
	//{
	//	ptemp.x=curvePoint[i].x;
	//	ptemp.y=curvePoint[i].y;
	//	ptemp.z=curvePoint[i].z;
	//	vCurvePointWorld.push_back(ptemp);
	//}

}
void GropInd1(int nsize,vector<vector<int>>&vv_nid)
{
	vector<int>v_nid;
	v_nid.clear();
	int nlastp=0;
	int nst=0;
	while (nlastp<=nsize-4)//分组下标
	{
		for(int k=nst;k<=nst+3;k++)
		{
			v_nid.push_back(k);

		}
			vv_nid.push_back(v_nid);
             int vnidsize=v_nid.size();
		     nlastp=v_nid[vnidsize-1];
			v_nid.clear();
			nst=nlastp;
		}
		

	}
void GropInd2(int nsize,vector<vector<int>>&vv_nid)
{
	vector<int>v_nid;
	v_nid.clear();
	int nlastp=0;
	int nst=0;
	while (nlastp<=nsize-2)//分组下标
	{
		for(int k=nst;k<=nst+3;k++)
		{
			v_nid.push_back(k);

		}
			vv_nid.push_back(v_nid);
             int vnidsize=v_nid.size();
		     nlastp=v_nid[vnidsize-1];
			v_nid.clear();
			nst=nlastp-2;
		}
		

	}
void Create_Bessel_inGro(vector< PointCordTypeDef> vREFPointsWold,vector<vector<int>> vv_nid,vector< PointCordTypeDef>&whol_curPonts)
{
	vector<vector<PointCordTypeDef>>vvcurvePoint;
	vvcurvePoint.clear();
	for(int g=0;g<vv_nid.size();g++)
	{
		int r=vv_nid[g][0];//index
		int s=vv_nid[g][1];
		int c=vv_nid[g][2];
		int q=vv_nid[g][3];
		PointCordTypeDef originPoint[4];
		originPoint[0]=vREFPointsWold[r];//pints
		originPoint[1]=vREFPointsWold[s];
		originPoint[2]=vREFPointsWold[c];
		originPoint[3]=vREFPointsWold[q];
		vector<PointCordTypeDef>curvePoint;
		curvePoint.clear();
		fitbesscurvUSE4points(originPoint,curvePoint);

		vvcurvePoint.push_back(curvePoint);	
		if (g>=2)
		{
			int nsize=vvcurvePoint[g-1].size();
			float fPointS[3]={vvcurvePoint[g-1][nsize-1].x,vvcurvePoint[g-1][nsize-1].y,vvcurvePoint[g-1][nsize-1].z};
			float fPointA[3]={vvcurvePoint[g][0].x,vvcurvePoint[g][0].y,vvcurvePoint[g][0].z};
			float fSADist=zxh::VectorOP_Distance(fPointS,fPointA,3);
			float fxxx=0;
		}

	}

	//变成一个点集

	for (int j=0;j<vvcurvePoint.size();j++)
	{
		vector<PointCordTypeDef>tmpvep;
		tmpvep=vvcurvePoint[j];
		for (int k=0;k<tmpvep.size();k++)
		{
			PointCordTypeDef tmpPon=tmpvep[k];
			whol_curPonts.push_back(tmpPon);
		}
	}

}
void Putinto_wholvec(vector<vector<PointCordTypeDef>>vvcurvePoint,vector<vector<vector<PointCordTypeDef>>>&vvWhcurvePoint,PointCordTypeDef originPoint[4])
{
	int nXdnum=originPoint[0].n;
	vvWhcurvePoint[nXdnum].push_back(vvcurvePoint[0]);
	vvWhcurvePoint[nXdnum+1].push_back(vvcurvePoint[1]);
	vvWhcurvePoint[nXdnum+2].push_back(vvcurvePoint[2]);
}
void Aver_loc_Ponts(vector<vector<vector<PointCordTypeDef>>> vvWhcurvePoint,vector<PointCordTypeDef>&vcurvePoint)
{
	vcurvePoint.clear();
	for(int i=0;i<vvWhcurvePoint.size();i++)
	{
		vector<vector<PointCordTypeDef>> vvptmp=vvWhcurvePoint[i];
		int npointsNUM=vvptmp[0].size();
		for(int k=0;k<npointsNUM;k++)
		{
			float fsumx=0;
			float fsumy=0;
			float fsumz=0;
			int nsumNUM=0;
			for(int j=0;j<vvptmp.size();j++)
			{
				vector<PointCordTypeDef>vptmp=vvptmp[j];
				fsumx=fsumx+vptmp[k].x;
				fsumy=fsumy+vptmp[k].y;
				fsumz=fsumz+vptmp[k].z;
				nsumNUM++;
			}
			float favex=fsumx/nsumNUM;
			float favey=fsumy/nsumNUM;
			float favez=fsumz/nsumNUM;
			//
			PointCordTypeDef Pnttmp;
			Pnttmp.x=favex;
			Pnttmp.y=favey;
			Pnttmp.z=favez;
			vcurvePoint.push_back(Pnttmp);
		}
	}
}
void Create_Bessel_inGro1(vector< PointCordTypeDef> vREFPointsWold,vector<vector<int>> vv_nid,vector< PointCordTypeDef>&whol_curPonts,vector<vector<vector<PointCordTypeDef>>> &vvWhcurvePoint)
{
	vector<vector<PointCordTypeDef>>vvcurvePoint;//保存四个点生成的一组（三条）线段
	vvcurvePoint.clear();
	for(int g=0;g<vv_nid.size();g++)
	{
		int r=vv_nid[g][0];//index
		int s=vv_nid[g][1];
		int c=vv_nid[g][2];
		int q=vv_nid[g][3];
		PointCordTypeDef originPoint[4];
		originPoint[0]=vREFPointsWold[r];//pints
		originPoint[1]=vREFPointsWold[s];
		originPoint[2]=vREFPointsWold[c];
		originPoint[3]=vREFPointsWold[q];
		vector<PointCordTypeDef>vcurvePoint;
		vcurvePoint.clear();
		fitbesscurvUSE4points1(originPoint,vcurvePoint,vvcurvePoint);
		//vvcurvePoint.push_back(vcurvePoint);	
	/*	if (g>=2)
		{
			int nsize=vvcurvePoint[g-1].size();
			float fPointS[3]={vvcurvePoint[g-1][nsize-1].x,vvcurvePoint[g-1][nsize-1].y,vvcurvePoint[g-1][nsize-1].z};
			float fPointA[3]={vvcurvePoint[g][0].x,vvcurvePoint[g][0].y,vvcurvePoint[g][0].z};
			float fSADist=zxh::VectorOP_Distance(fPointS,fPointA,3);
			float fxxx=0;
		}*/
		Putinto_wholvec(vvcurvePoint,vvWhcurvePoint,originPoint);
	}
	//对每个局部的点，分别做平均
	//vector<PointCordTypeDef> vave_whocur;

	Aver_loc_Ponts( vvWhcurvePoint,whol_curPonts);
	//变成一个点集

	//for (int j=0;j<vvcurvePoint.size();j++)
	//{
	//	vector<PointCordTypeDef>tmpvep;
	//	tmpvep=vvcurvePoint[j];
	//	for (int k=0;k<tmpvep.size();k++)
	//	{
	//		PointCordTypeDef tmpPon=tmpvep[k];
	//		whol_curPonts.push_back(tmpPon);
	//	}
	//}

}
void Create_Bessel_inGro_OVlap(vector< PointCordTypeDef> vREFPointsWold,vector<vector<int>> vv_nid,vector< PointCordTypeDef>&whol_curPonts)
{
	vector<vector<PointCordTypeDef>>vvcurvePoint;
	vvcurvePoint.clear();
	for(int g=0;g<vv_nid.size();g++)
	{
		int r=vv_nid[g][0];//index
		int s=vv_nid[g][1];
		int c=vv_nid[g][2];
		int q=vv_nid[g][3];
		PointCordTypeDef originPoint[4];
		originPoint[0]=vREFPointsWold[r];//pints
		originPoint[1]=vREFPointsWold[s];
		originPoint[2]=vREFPointsWold[c];
		originPoint[3]=vREFPointsWold[q];
		vector<PointCordTypeDef>curvePoint;
		curvePoint.clear();
		//bessel fitting
		fitbesscurvUSE4points(originPoint,curvePoint);
		vvcurvePoint.push_back(curvePoint);	
		if (g>=2)
		{
			int nsize=vvcurvePoint[g-1].size();
			float fPointS[3]={vvcurvePoint[g-1][nsize-1].x,vvcurvePoint[g-1][nsize-1].y,vvcurvePoint[g-1][nsize-1].z};
			float fPointA[3]={vvcurvePoint[g][0].x,vvcurvePoint[g][0].y,vvcurvePoint[g][0].z};
			float fSADist=zxh::VectorOP_Distance(fPointS,fPointA,3);
			float fxxx=0;
		}

	}

	//变成一个点集

	for (int j=0;j<vvcurvePoint.size();j++)
	{
		vector<PointCordTypeDef>tmpvep;
		tmpvep=vvcurvePoint[j];
		for (int k=0;k<tmpvep.size();k++)
		{
			PointCordTypeDef tmpPon=tmpvep[k];
			whol_curPonts.push_back(tmpPon);
		}
	}

}
void GropXd(vector<PointCordTypeDef>&vREFPoints)
{
	for (int k=0;k<vREFPoints.size()-1;k++)
	{
		vREFPoints[k].n=k;
	}
	vREFPoints[vREFPoints.size()-1].n=-1;
}
bool StorUnique(vector<PointCordTypeDef> &vUnorgaPointsWorld)
{
	bool bunique=true;
	for(int i=0;i<vUnorgaPointsWorld.size();i++)
	{
		float fi[3]={vUnorgaPointsWorld[i].x,vUnorgaPointsWorld[i].y,vUnorgaPointsWorld[i].z};
		for(int j=i+1;j<vUnorgaPointsWorld.size();j++)
		{
			float fj[3]={vUnorgaPointsWorld[j].x,vUnorgaPointsWorld[j].y,vUnorgaPointsWorld[j].z};
			float dist=zxh::VectorOP_Distance(fi,fj,3);
			if(dist<0.00001)
			{
				vUnorgaPointsWorld.erase(vUnorgaPointsWorld.begin()+j);
				j--;
			}
		}
	}
	return true;
}
int main(int argc, char *argv[])

{
	//debug
	//char *chpontsFile ="E:/work_jdq/CTA_data/RCAAEF/training_jdq/01/OrdP_image01_R0_1.txt";
	//char *chResultPathName ="E:/work_jdq/CTA_data/RCAAEF/training_jdq/01/Bess_image01_R0_1";
	//debug
	//release
	if( argc < 2 )
	{
		cerr << "Usage: " << endl;
		cerr << "jdqBesselCurve	line(.txt)	Results" << endl;
		return -1;
	}
	char *chpontsFile =argv[1];
	char *chResultPathName =argv[2];
	//release

	//读取手工选取的坐标点
	vector< PointCordTypeDef> vREFPoints;
	ReadTxt_Pdef(chpontsFile,vREFPoints);
	//将手工点的下标分组，每四个一组
	int nGe=4;
	int nsize=vREFPoints.size();
	vector<vector<int>>vv_nid;
	vv_nid.clear();
	//GropInd1(nsize,vv_nid);
	GropInd2(nsize,vv_nid);
	//将手工点按照将要生成的线段分组，比如第一个起点属于第一段
	GropXd(vREFPoints);
	//以组为单位生成多组Bessel曲线
	PointCordTypeDef originPoint[4];
	vector<PointCordTypeDef> whol_curPonts;
	whol_curPonts.clear();
	vector<vector<vector<PointCordTypeDef>>> vvWhcurvePoint(nsize-1);
	Create_Bessel_inGro1(vREFPoints,vv_nid,whol_curPonts,vvWhcurvePoint);
	cout << "Bessel curve fitting is down!" << endl;
	//删除重复点
	StorUnique(whol_curPonts);
	// save
	int iCurveNameLen = strlen(chResultPathName) + strlen(".txt");
	char *CurveNamevtk = (char*)malloc(iCurveNameLen);
	strcpy(CurveNamevtk, chResultPathName);
	strcat(CurveNamevtk, ".vtk");
	WriteVtk(whol_curPonts, CurveNamevtk);
	//
	char *CurveNametxt = (char*)malloc(iCurveNameLen);
	strcpy(CurveNametxt, chResultPathName);
	strcat(CurveNametxt, ".txt");
	WriteTxt(whol_curPonts, CurveNametxt);
	cout << "Save successfully: Bessel curve !" << endl;
}
