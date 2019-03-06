

/*=========================================================================

Program:   ZXH CoronaryArteryExtraction Software
Author:	 Dengqiang Jia
Module:    $RCSfle: zxhcaeDMP.cpp    $
Language:  C++
Date:      $Date: From  2011-03 $
Version:   $Revision: 2.2.1 $

=========================================================================*/
/// \brief
/// Spatially encoded mutual information + free-form deformation registration
/// Linear
/// gradient descent
/// save all registration info into zxhRegistrationStruct StructRegistration which set to optimiser
/// For concatenated transforms all (for FFD and LocallyAffines in optimization regularisation step)
///                                 First update zxhRegistrationStruct::m_pConcatenatedTransformsByRegridding(spacing 1mm)
///                                 then transform all images in ref_space using zxhRegistrationStruct.m_CopyRefXXXOrg
///                                 finally set current transform to identity
/// For preset transformation, unless same spacing FFD for the first Reg, otherwise would be treated as concatenation
///

//
//void Help()
//{
//	std::cout<<" An simple example for registration of images, target.nii.gz and source.nii.gz, result save as res: \n" ;
//	std::cout<<" zxhsemi0 -test target.nii.gz -ref source.nii.gz -o result0 -ffd 20 20 20 -bending 0.001\n";
//	std::cout<<" zxhsemi -test target.nii.gz -ref source.nii.gz -o result -pre 0 result0.FFD -ffd 20 20 20 -ffd 10 10 10 -Reg 2 -sub 2 2 2 -sub 1 1 1 -bending 0.001\n";
//	std::cout<<"OR \n";
//	std::cout<<" zxhsemi -test target.nii.gz -ref source.nii.gz -o result -hierarchy 3 -bending 0.0031\n\n";
//
//	std::cout<<"  <-test/target img.>     (test or target image)\n" ;
//	std::cout<<"  <-ref/source img.>      (reference or source image)\n" ;
//	std::cout<<"  <-o savename>           (string for saving transformed -ref/-source image, file names prefix)\n" ;
//	std::cout<<"  <-maskt img.>           (mask image on test image, use -maskr on ref image) \n" ;
//	std::cout<<"  USE -ffd: zxhsemi0 to fast init and get the .FFD for setting -pre, and then set -ffd\n" ;
//	std::cout<<"  <-ffd fx fy fz>         (spacing of free form deformation, FFD, multi-input related to -Reg)\n" ;
//	std::cout<<"  <-pre 0 s>              (pre set transformation field)\n";
//	std::cout<<"  <-sub fx fy fz [ft]>    ([3 3 3], sampling spacing; positive: millimeters interval; negative: pixels interval)\n";
//	std::cout<<"  <-Reg i>                ([1] number of multi-level registrations)\n" ;
//	std::cout<<"  OR USE -hierarchy, simple and not need to set -ffd,-sub,-Reg:\n" ; 
//	std::cout<<"  <-hierarchy n>          ([3] number of multi-level FFD registration, where\n";
//	std::cout<<"                           the first level of -ffd spacing is one forth of test image extent, and halve in subsequence level\n" ; 
//	std::cout<<"                           the final level of -sub sampling is pixel size of the test image\n" ; 	 
//	std::cout<<"\n";
//	std::cout<<"  <-bending f..f>         ([0.001]weighting for bending energy, recommend f=0.001~0.01)\n" ;
//	std::cout<<"  OPTIONS of spatially encoding scheme\n"  ;// similarity computation, default normalized mutual information\n" ; 
//	std::cout<<"  <-semiradius f...f>     (radius of local region, default set to twice ffd spacing)\n";
//	//std::cout<<"  <-semiwidth f...f>      (or -semisize, width/size of local region, default 80mm)\n";
//	//std::cout<<"  <-semikernel s>         ([0], -1='Gaussian', 0='ZeroBSpline', 3='BSpline')\n";
//
//	std::cout<<"\n" ;
//
//} 
//void HELP()
//{ 
//	Help() ;
//	std::cout<<"------------------------------------------------------\n" ;
//	std::cout<<"  OPTIONS of gradient optimization computation; use setting in previous -Reg when un-set \n" ; 
//}  
//int main(int argc, char* argv[])
//{
//	zxh::echo_zxh(argc,argv);
//	if( argc == 1 )
//	{
//		std::cout<<"  zxhsemi [options]        \n";
//		Help();
//		return 0 ;
//	}
//	if( argc == 2 && strcmp( argv[1],"-H" )==0 )
//	{
//		std::cout<<"  zxhsemi [options]        \n";
//		HELP();
//		return -1 ;
//	}
//	if( glbVerboseOutput>0 )
//	{
//		std::cout<<"\n * \n * zxhsemi, version of 2011-03  \n * \n";
//		zxh::echo_arguments( argc, argv ) ;
//	} 
//	return zxhsemi_main(argc,argv);
//} 
//int zxhsemi_main(int argc, char* argv[])
//{}
//


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
#include <math.h>
#include<stdlib.h>
#include<stdio.h>

#include "zxhImageData.h"
#include "zxhImageGipl.h"
#include "zxhImageNifti.h"

#include "jdq2017util.h"
#include "jdqdijkstra.h"
#include "jdqPoint.h"

#define TAB_CHAR	9
#define M_PI 3.14159265358979323846
#define SPHERE_RADIUS 4
using namespace std;

typedef struct
{
	float x;
	float y;
	float z;
}PointCordTypeDef;


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
		fImgPixel[0] = PointCord[i].x;
		fImgPixel[1] = PointCord[i].y;
		fImgPixel[2] = PointCord[i].z;
		WriteFileTxt<<right<<fixed<<setfill('0')<<setprecision(4)<< -fImgPixel[0] << " " << -fImgPixel[1] << " " << fImgPixel[2] << "\n ";

	}	
}



char *GetPathEXT(char *chFileName)
{  char path_buffer[_MAX_PATH];  
char drive[_MAX_DRIVE];  
char dir[_MAX_DIR];  
char fname[_MAX_FNAME];  
char ext[_MAX_EXT];  

_splitpath( chFileName, drive, dir, fname, ext );  

return ext;
}
bool ReadTxtAsPC(char *chFileName,float fresamle, vector<PointCordTypeDef> &PointCord)
{
	if (chFileName == NULL)
	{
		cout << "Cannot find txt-file!" << endl;
		return 1;
	}
	//按照jda point格式读取line文件
	std::vector<jdq2017::point3D> ref1,ref;
	std::vector<jdq2017::point3D> cl;
	if ( ! jdq2017::readCenterline(chFileName, ref1))
	{
		std::vector<double> rad;
		std::vector<double> io;
		if ( ! jdq2017::readReference(chFileName,ref1,rad,io))
		{
			std::cerr << "Error in reading line data" << std::endl;
			return 1;
		}
	}
	cout<<" Original number of points: " <<ref1.size()-1<<endl;
	//resample the points
	jdq2017::ResamplePaths<jdq2017::point3D>Respler;
	Respler(ref1,fresamle);
	ref=Respler.resultPath();
	//将ref points 的数据转存成PointCord
	if (!PointCord.empty())
	{
		PointCord.clear();
	}
	PointCordTypeDef strctTempPoint;
	for (int i = 0; i < ref.size(); i++)
	{
		strctTempPoint.x =-1*ref[i]._x;
		strctTempPoint.y =-1*ref[i]._y;
		strctTempPoint.z =ref[i]._z;
		PointCord.push_back(strctTempPoint);
	}
	return true;
}

bool ReadPointTxt_ostia(char *filename,vector< PointCordTypeDef> &cl)
{
	string strNum;
	int nStart = 0, nEnd = 0;
	PointCordTypeDef strctTempPoint;
	if(cl.size()>0)
		cl.clear();
	ifstream iFileIn(filename,ios::in);

	if(!iFileIn)
	{
		cout << "Unable to open txt-file: " << endl;
		cout << filename << endl;
		return false; //exit(1); // terminate with error
	}
	else
	{
		while(!iFileIn.eof())
		{
			float fx,fy,fz;
			iFileIn>>fx>>fy>>fz;
			strctTempPoint.x=fx;
			strctTempPoint.y=fy;
			strctTempPoint.z=fz;
			cl.push_back(strctTempPoint);
		}
	}

	return true;
}
bool ResampleCurvebyImagespacing(string strRaw,char *chCurvename,vector<PointCordTypeDef>&vPathPointsWorld)
{
	zxhImageDataT<short> imgReadRaws;//Change by JDQ
	float fNewImgSpacing[]={1,1,1,1};//Add by JDQ
	if( zxh::OpenImage( &imgReadRaws, strRaw ) == false )
	{
		std::cerr << "Raw image(nifti-file) is not found!"<<endl; 
		return -1;
	}
	imgReadRaws.GetImageSpacing(fNewImgSpacing[0],fNewImgSpacing[1],fNewImgSpacing[2],fNewImgSpacing[3] );//Add by JDQ
	float fresamle=zxh::minf(fNewImgSpacing[1],fNewImgSpacing[2])*0.5;
	////////

	char* chext=GetPathEXT(chCurvename);
	if(strcmp(chext,".vtk")==0)// read vtk file
	{
		cout<<"open successfully, format: vtk"<<endl;
		fstream DetectFile;
		DetectFile.open(chCurvename,ios::in);
		if(DetectFile)
		{
			ReadVtk(chCurvename, vPathPointsWorld);

		}
		DetectFile.close();

	}
	else if(strcmp(chext,".txt")==0)// read vtk file
	{
		cout<<"open successfully, format: txt"<<endl;
		fstream DetectFile;
		DetectFile.open(chCurvename,ios::in);
		if(DetectFile)
		{
			ReadTxtAsPC(chCurvename,fresamle, vPathPointsWorld);

		}
		DetectFile.close();

	}
}
bool ResampleCurvebyLength(float freslen,char *chCurvename,vector<PointCordTypeDef>&vPathPointsWorld)
{

	////////

	char* chext=GetPathEXT(chCurvename);
	if(strcmp(chext,".vtk")==0)// read vtk file
	{
		cout<<"open successfully, format: vtk"<<endl;
		fstream DetectFile;
		DetectFile.open(chCurvename,ios::in);
		if(DetectFile)
		{
			ReadVtk(chCurvename, vPathPointsWorld);

		}
		DetectFile.close();

	}
	else if(strcmp(chext,".txt")==0)// read vtk file
	{
		cout<<"open successfully, format: txt"<<endl;
		fstream DetectFile;
		DetectFile.open(chCurvename,ios::in);
		if(DetectFile)
		{
			ReadTxtAsPC(chCurvename,freslen, vPathPointsWorld);

		}
		DetectFile.close();

	}
	return true;
}
int main(int argc, char *argv[])
{
	//if( argc < 3 )
	//{
	//	cerr << "Usage: " << endl;
	//	cerr << "jdqresampleCurve.cpp	curve	resamplelength results  -IorL -txt(-vtk)" << endl;
	//	return -1;
	//}
	argc=4;
	string strRaw="";
	char *chCurvename ="";//  "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
	float resampllength =0.1; //"F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/reference.vtk";
	char *chResultPathName ="";
	char *chIorL =argv[4];
	char *chvtkortxt ="";
	vector<PointCordTypeDef> vPathPointsWorld;
	chvtkortxt=argv[5];
	if( argc ==4&&strcmp(chIorL,"-I")==0)// resample by image spacing
	{
		cout<<"Resample the curve by image spacing!"<<endl;
		string strRaw=string(argv[2]);
		chCurvename =argv[1];// 
		chResultPathName =argv[3];
		//string strRaw="J:/work_jdq/data_DNN/RCAAEF_32/training/dataset00/image/image.nii.gz";
		//chCurvename ="J:/work_jdq/data_DNN/RCAAEF_32/training/dataset00/centerlines/vessel0.txt";// 
		//chResultPathName ="J:/work_jdq/data_DNN/RCAAEF_32/training/dataset00/centerlines/res_vessel0";
		ResampleCurvebyImagespacing(strRaw,chCurvename,vPathPointsWorld);
	}
	if(argc ==4&&strcmp(chIorL,"-L")==0)// resample by lenth
	{
		cout<<"Resample the curve by the length!"<<endl;
		chCurvename =argv[1];//  "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
		resampllength =atof(argv[2]);//"F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/reference.vtk";
		chResultPathName =argv[3];// "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/MCLine.nii.gz";
		//chCurvename ="J:/work_jdq/data_DNN/RCAAEF_32/training/dataset00/centerlines/vessel0.txt";// 
		//chResultPathName ="J:/work_jdq/data_DNN/RCAAEF_32/training/dataset00/centerlines/res_vessel0";
		ResampleCurvebyLength(resampllength,chCurvename,vPathPointsWorld);
	}
	cout<<"New number of points: " <<vPathPointsWorld.size()-1<<endl;
	int iCurveNameLen = strlen(chResultPathName) + strlen(".txt");
	if (strcmp(chvtkortxt,"-vtk")==0)
	{
		char *CurveNamevtk = (char*)malloc(iCurveNameLen);
		strcpy(CurveNamevtk, chResultPathName);
		strcat(CurveNamevtk, ".vtk");
		WriteVtk(vPathPointsWorld, CurveNamevtk);
	}
	//
	if (strcmp(chvtkortxt,"-txt")==0)
	{
		char *CurveNametxt = (char*)malloc(iCurveNameLen);
		strcpy(CurveNametxt, chResultPathName);
		strcat(CurveNametxt, ".txt");
		WriteTxt(vPathPointsWorld, CurveNametxt);
	}

	cout << "resmaple Curve successfully!" << endl;
}

