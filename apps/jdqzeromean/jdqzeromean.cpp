

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
#include<numeric>

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
bool TransDatanum2Greynum(int np[3],int &npd,int ImgNewSize[3])
{
	npd=np[2] * ImgNewSize[1] * ImgNewSize[0] + np[1] * ImgNewSize[0] + np[0];
	return true;
}
bool TransGreynum2Datanum(int npd,int np[3],int ImgNewSize[3])
{
	np[2]=npd/(ImgNewSize[1] * ImgNewSize[0]);
	int nxy=npd%(ImgNewSize[1] * ImgNewSize[0]);
	np[1]=nxy/ImgNewSize[0];
	np[0]=nxy%ImgNewSize[0];
	return true;
}
void DiaPoiWBal(const short *sintData,zxhImageDataT<short>&imgNewBlaRaw,short *snewintData,float flen)
{
	int voxnum=imgNewBlaRaw.GetNumberOfPixels();
	int ImgSize[4]={1};
	imgNewBlaRaw.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	for(int i=0;i<voxnum;i++)//初始化空白图像
	{
		short sinten=sintData[i];	
		int PointPos[4]={0};//获取当前点的图像坐标
		TransGreynum2Datanum(i,PointPos,ImgSize);
		if(sinten!=0)
		{
			float PointPosWorld[ZXH_ImageDimensionMax]={0};

			//局部‘卷积',注意边界溢出
			int nkSize[3]={40,40,40};
			int nzs=zxh::maxf(PointPos[2]-nkSize[2],0);
			int nys=zxh::maxf(PointPos[1]-nkSize[1],0);
			int nxs=zxh::maxf(PointPos[0]-nkSize[0],0);

			int nze=zxh::minf(PointPos[2]+nkSize[2],ImgSize[2]);
			int nye=zxh::minf(PointPos[1]+nkSize[1],ImgSize[1]);
			int nxe=zxh::minf(PointPos[0]+nkSize[0],ImgSize[0]);
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
						int nnp=0;
						TransDatanum2Greynum(PointLocPos,nnp,ImgSize);
						float fInt=snewintData[nnp];	
						if (fInt>0) continue;
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
						if(fdist<=flen)
						{
							snewintData[nnp]=1;
						}
					}
		}
		else
		{
			continue;
		}
	}

}
bool Cal_meanvarise(vector<float>resultSet,float &mean,float &stdev)
{
	double sum = std::accumulate(std::begin(resultSet), std::end(resultSet), 0.0);
	 mean =  sum / resultSet.size(); //均值

	double accum  = 0.0;
	std::for_each (std::begin(resultSet), std::end(resultSet), [&](const double d) {
		accum  += (d-mean)*(d-mean);
	});
	 stdev = sqrt(accum/(resultSet.size()-1)); //方差
	return true;

}
bool zeromeanvar(zxhImageDataT<short>&imgReadRaws,zxhImageDataT<float>&imgNewBlaRaw,float fmean,float fvar)
{
	int ImgNewSize[4]={1};
	imgNewBlaRaw.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);

	for(int it=0;it<ImgNewSize[3];++it)
			for(int iz=0;iz<ImgNewSize[2];++iz)
				for(int iy=0;iy<ImgNewSize[1];++iy)
					for(int ix=0;ix<ImgNewSize[0];++ix)
					{
						short sintens=imgReadRaws.GetPixelGreyscale(ix,iy,iz,it);
						float fnewint=(sintens-fmean)/fvar;
						imgNewBlaRaw.SetPixelByGreyscale(ix,iy,iz,it,fnewint);
					}
	return true;
}

int main(int argc, char *argv[])
{
	if( argc < 3 )
	{
		cerr << "Usage: " << endl;
		cerr << "jdqzeromean image length result.nii.gz" << endl;
		return -1;
	}
	string strFileNameRaw=string(argv[1]);//  "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
	string strFileNameMaskRaw=string(argv[2]);//"F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/reference.vtk";
	string chResultFile =string(argv[3]);// "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/MCLine.nii.gz";


	/*string strFileNameRaw ="E:/work_jdq/RCAAEF/train/rcaaef_32_img00/image.nii.gz";
	string strFileNameMaskRaw ="E:/work_jdq/RCAAEF/train/rcaaef_32_img00/whs_lab_image_ROI.nii.gz";
	string chResultFile ="E:/work_jdq/RCAAEF/train/rcaaef_32_img00/zmv_whs_lab_image.nii.gz";*/

	//读取图像文件
	zxhImageDataT<short> imgReadRaws,imgMaskRaws;
	zxhImageDataT<float> imgNewBlaRaw;//Change by JDQ
	if( zxh::OpenImage( &imgReadRaws, strFileNameRaw ) == false )
	{
		std::cerr << "Raw image(nifti-file) is not found!"; 
		return -1;
	}
	if( zxh::OpenImage( &imgMaskRaws, strFileNameMaskRaw ) == false )
	{
		std::cerr << "Mask Raw image(nifti-file) is not found!"; 
		return -1;
	}
	const short *sintData = imgReadRaws.GetImageData();
	const short *sintMaskData = imgMaskRaws.GetImageData();
	//计算均值和方差
	float fmean=0;
	float fvar=0;
	vector<float>vint;
	vint.clear();
	int ImgSizeMaskraw[4]={1};
	imgMaskRaws.GetImageSize(ImgSizeMaskraw[0],ImgSizeMaskraw[1],ImgSizeMaskraw[2],ImgSizeMaskraw[3]);
	for (int i = 0; i < ImgSizeMaskraw[0]*ImgSizeMaskraw[1]*ImgSizeMaskraw[2]-1; i ++) 
	{
		short sintmask=sintMaskData[i];
		if(sintmask!=0)
		vint.push_back(sintData[i]);

	}
	Cal_meanvarise(vint,fmean,fvar);
	cout<<"mean: "<<fmean<<" , "<<"var: "<<fvar<<endl;

	//生成空白图像->除intensity信息外，此图像其他信息来自于原始图像
	imgNewBlaRaw.NewImage(imgReadRaws.GetImageInfo() );	
	int ImgSize[4]={1};
	imgNewBlaRaw.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	//
	zeromeanvar(imgReadRaws,imgNewBlaRaw,fmean,fvar);//遍历点 归一化
	string chFileName2(chResultFile);
	zxh::SaveImage(&imgNewBlaRaw,chResultFile);
	cout << "zeromeanvar: successfully!" << endl;


}

