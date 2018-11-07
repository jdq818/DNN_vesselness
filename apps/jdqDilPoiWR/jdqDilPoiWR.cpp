

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

#include <algorithm>
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

typedef struct
{
	float x;
	float y;
	float z;
	float r;
}PoiRadCordTypeDef;

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

//Read reference standard and method centerline


bool ReadTxtAsPC(char *chFileName, vector<PoiRadCordTypeDef> &PointCord)
{
	if (chFileName == NULL)
	{
		cout << "Cannot find VTK-file!" << endl;
		return 1;
	}
	//按照jda point格式读取line文件
	std::vector<jdq2017::point3D> ref;
	std::vector<double> rad;
	std::vector<double> io;
	std::vector<jdq2017::point3D> cl;
	if ( ! jdq2017::readReference(chFileName, ref, rad, io))
	{
		std::cerr << "Error in reading line data" << std::endl;
		return 1;
	}

	//将ref points 的数据转存成PointCord
	if (!PointCord.empty())
	{
		PointCord.clear();
	}
	PoiRadCordTypeDef strctTempPoint;
	for (int i = 0; i < ref.size(); i++)
	{
		strctTempPoint.x =-1*ref[i]._x;
		strctTempPoint.y =-1*ref[i]._y;
		strctTempPoint.z =ref[i]._z;
		strctTempPoint.r =rad[i];
		//cout<<"Points:" <<strctTempPoint.x<<" "<<strctTempPoint.y<<" "<<strctTempPoint.z<<" "<<endl;
		PointCord.push_back(strctTempPoint);
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
char *GetPathEXT(char *chFileName)
{  
char path_buffer[_MAX_PATH];  
char drive[_MAX_DRIVE];  
char dir[_MAX_DIR];  
char fname[_MAX_FNAME];  
char ext[_MAX_EXT];  

_splitpath( chFileName, drive, dir, fname, ext );  

return ext;
}
bool ReadPointTxt(char *filename,vector< PointCordTypeDef> &cl)
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
			iFileIn>>strctTempPoint.x>>strctTempPoint.y>>strctTempPoint.z;
			cl.push_back(strctTempPoint);
		}
	}

	return true;
}
bool BoundaryCorrect(int *PointPos,int ImgNewVslsSize[4])
{
	for (int i=0;i<3;i++)
	{
		PointPos[i]=zxh::maxf(0,PointPos[i]);
		PointPos[i]=zxh::minf(ImgNewVslsSize[i]-1,PointPos[i]);
	}
	return true;
}
void DiaPoiWBal(zxhImageDataT<short>&imgNewBlaRaw,vector<PoiRadCordTypeDef> &vPathPointsWorld)
{
	int ImgSize[4]={1};
	imgNewBlaRaw.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	for(int it=0;it<ImgSize[3];++it)//初始化空白图像
			for(int iz=0;iz<ImgSize[2];++iz)
				for(int iy=0;iy<ImgSize[1];++iy)
					for(int ix=0;ix<ImgSize[0];++ix)
					{
						imgNewBlaRaw.SetPixelByGreyscale(ix,iy,iz,it,0);
					}
	//遍历点
					float fspacing[4]={0,0,0,0};
    imgNewBlaRaw.GetImageSpacing(fspacing[0],fspacing[1],fspacing[2],fspacing[3]);
	float fmin=zxh::minf(fspacing[0],fspacing[1]);
	float fminspacing=zxh::minf(fmin,fspacing[2]);
	for (int i=0;i<vPathPointsWorld.size();i++)
	{
		float PointPosWorld[ZXH_ImageDimensionMax]={0};
		int PointPos[4]={0};//获取当前点的图像坐标
		float fR=0;
		PointPosWorld[0]=vPathPointsWorld[i].x;
		PointPosWorld[1]=vPathPointsWorld[i].y;
		PointPosWorld[2]=vPathPointsWorld[i].z;
		fR=vPathPointsWorld[i].r;
		imgNewBlaRaw.GetImageInfo()->WorldToImage(PointPosWorld);
		PointPos[0]=zxh::round(PointPosWorld[0]);
		PointPos[1]=zxh::round(PointPosWorld[1]);
		PointPos[2]=zxh::round(PointPosWorld[2]);
		BoundaryCorrect(PointPos,ImgSize);

		if(PointPos[0]==156&&PointPos[1]==58&&PointPos[2]==34)
							int xxx=0;
		//局部‘卷积',注意边界溢出
		int nkSize[3]={fR/fminspacing*1.5,fR/fminspacing*1.5,fR/fminspacing*1.5};
		int nzs=zxh::maxf(PointPos[2]-nkSize[2],0);
		int nys=zxh::maxf(PointPos[1]-nkSize[1],0);
		int nxs=zxh::maxf(PointPos[0]-nkSize[0],0);

		int nze=zxh::minf(PointPos[2]+nkSize[2],ImgSize[2]-1);
		int nye=zxh::minf(PointPos[1]+nkSize[1],ImgSize[1]-1);
		int nxe=zxh::minf(PointPos[0]+nkSize[0],ImgSize[0]-1);



		for (int nz=nzs;nz<=nze;nz++)
			for (int ny=nys;ny<=nye;ny++)
				for (int nx=nxs;nx<=nxe;nx++)
				{
					float PointWorld[ZXH_ImageDimensionMax]={0};
					int PointLocPos[4]={0};
					float PointLocWorld[4]={0};
					PointLocPos[0]=nx;
					PointLocPos[1]=ny;
					PointLocPos[2]=nz;
					BoundaryCorrect(PointLocPos,ImgSize);
					PointLocWorld[0]=PointLocPos[0];
					PointLocWorld[1]=PointLocPos[1];
					PointLocWorld[2]=PointLocPos[2];
					//获取世界坐标
					float PointCenWorld[3]={PointPos[0],PointPos[1],PointPos[2]};//中心点
					imgNewBlaRaw.GetImageInfo()->ImageToWorld(PointCenWorld);//中心点世界坐标
					imgNewBlaRaw.GetImageInfo()->ImageToWorld(PointLocWorld);//邻点
					//计算与中心点的距离
					float fdist=zxh::VectorOP_Distance(PointCenWorld,PointLocWorld,3);
					//判断是否在球内
					if(fdist<=fR)
					{
						
						float fInt=imgNewBlaRaw.GetPixelGreyscale(PointLocPos[0],PointLocPos[1],PointLocPos[2],0);
						if (fInt>0) continue;
						else
						imgNewBlaRaw.SetPixelByGreyscale(PointLocPos[0],PointLocPos[1],PointLocPos[2],0,1);
					}
				}
	}

}

int main(int argc, char *argv[])
{
	if( argc < 3 )
	{
		cerr << "Usage: " << endl;
		cerr << "jdqDilPoiWR image line result.nii.gz" << endl;
		return -1;
	}
	string strFileNameRaw =string(argv[1]);//  "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
	char *chFileName =argv[2];//"F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/reference.vtk";
	string chResultFile =string(argv[3]);// "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/MCLine.nii.gz";


	/*string strFileNameRaw ="E:/work_jdq/RCAAEF/train/rcaaef_32_img00/image.nii.gz";
	char *chFileName = "E:/work_jdq/CTA_data_preprocess/RCAAEF/training_jdq/00/linewithrad/Bess_OrdP_image00_R_0.txt";
	string chResultFile = "E:/work_jdq/RCAAEF/train/rcaaef_32_img00/lab_R_0.nii.gz";*/

	//读取图像文件
	zxhImageDataT<short> imgReadRaws,imgNewBlaRaw;//Change by JDQ
	float fNewImgSpacing[]={1,1,1,1};//Add by JDQ
	float fOldImgSpacing[]={1,1,1,1};
	if( zxh::OpenImage( &imgReadRaws, strFileNameRaw ) == false )
	{
		std::cerr << "Raw image(nifti-file) is not found!"; 
		return -1;
	}

	//判断曲线文件类型是txt or vtk，并读取
	vector<PoiRadCordTypeDef> vPathPointsWorld;
	char* chext=GetPathEXT(chFileName);
	if(strcmp(chext,".vtk")==0)// read vtk file
	{

		/*	fstream DetectFile;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
		ReadVtkAsPC(chFileName, vPathPointsWorld);

		}
		DetectFile.close();*/
	}
	else if(strcmp(chext,".txt")==0)// read vtk file
	{

		fstream DetectFile;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadTxtAsPC(chFileName, vPathPointsWorld);

		}
		DetectFile.close();
	}
	//生成空白图像->除intensity信息外，此图像其他信息来自于原始图像
	imgNewBlaRaw.NewImage(imgReadRaws.GetImageInfo() );
	//
	DiaPoiWBal(imgNewBlaRaw,vPathPointsWorld);//dilation the points with a ball
	string chFileName2(chResultFile);
	zxh::SaveImage(&imgNewBlaRaw,chResultFile);
	cout << "Line has already dialate with radius in an image!" << endl;


}

