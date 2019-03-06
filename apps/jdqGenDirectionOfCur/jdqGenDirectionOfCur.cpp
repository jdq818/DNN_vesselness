

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
#include<direct.h>
#include <io.h>
#include<algorithm>

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
	float rad;
}PointCordTypeDef;
typedef struct
{
	float dx;
	float dy;
	float dz;
}DireCordTypeDef;

typedef struct
{
	PointCordTypeDef pcor;
	vector<DireCordTypeDef>vdire;
}PointCordDireTypeDef;

bool ReadDireTxt(const char *chFileName, vector<DireCordTypeDef> &PointCord)
{
	if (chFileName == NULL)
	{
		cout << "Cannot find txt-file!" << endl;
		return 1;
	}
	//按照jda point格式读取line文件
	std::vector<jdq2017::point3D>ref;
	if ( ! jdq2017::readCenterline(chFileName, ref))
	{
		std::cerr << "Error in direction data" << std::endl;
		return 1;
	}
	//将ref points 的数据转存成PointCord
	if (!PointCord.empty())
	{
		PointCord.clear();
	}
	DireCordTypeDef strctTempPoint;
	for (int i = 0; i < ref.size(); i++)
	{
		strctTempPoint.dx =ref[i]._x;
		strctTempPoint.dy =ref[i]._y;
		strctTempPoint.dz =ref[i]._z;
		PointCord.push_back(strctTempPoint);
	}
	return true;
}
void cf_findFileFromDir2(string mainDir, vector<string> &files)
{
	files.clear();
	const char *dir = mainDir.c_str();
	_chdir(dir);
	long hFile;
	_finddata_t fileinfo;

	if ((hFile = _findfirst("*.*", &fileinfo)) != -1)
	{
		do
		{
			if (!(fileinfo.attrib & _A_SUBDIR))//找到文件
			{
				char filename[_MAX_PATH];
				strcpy_s(filename, dir);
				//strcat_s(filename, "\\");
				strcat_s(filename, fileinfo.name);
				string temfilename = filename;
				files.push_back(temfilename);
			}
		} while (_findnext(hFile, &fileinfo) == 0);
		_findclose(hFile);
	}
}
bool ReadTxtAsPC(const char *chFileName,float fresamle, vector<PointCordTypeDef> &PointCord)
{
	if (chFileName == NULL)
	{
		cout << "Cannot find txt-file!" << endl;
		return 1;
	}
	//按照jda point格式读取line文件
	std::vector<jdq2017::point3D> refm,ref;
	std::vector<double> rad;
	std::vector<double> radm;
	std::vector<double> io;
	std::vector<jdq2017::point3D> cl;
	if ( !jdq2017::readReference(chFileName, refm, radm, io))
	{
		std::cerr << "Error in reading line data" << std::endl;
		return 1;
	}
	cout << "read success: !" << chFileName<<endl;
	//resample the points
	jdq2017::ResamplePaths<jdq2017::point3D>Respler;
	Respler(refm,radm,fresamle);
	ref=Respler.resultPath();
	rad=Respler.resultRadius();
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
		strctTempPoint.rad=rad[i];
		PointCord.push_back(strctTempPoint);
	}
	return true;
}
bool Calc_divec(int i,float fdivec[3],vector<PointCordTypeDef> &ref,float clSampling)
{

	int nFBsize=1/clSampling;
	int nFPosi=zxh::minf(i+nFBsize,ref.size()-1);
	int nBPosi=zxh::maxf(i-nFBsize,0);
	PointCordTypeDef curPontW=ref[i];
	PointCordTypeDef FPontW=ref[nFPosi];
	PointCordTypeDef BPontW=ref[nBPosi];
	float fvec1[3]={curPontW.x-BPontW.x,curPontW.y-BPontW.y,curPontW.z-BPontW.z};
	float fvec2[3]={FPontW.x-curPontW.x,FPontW.y-curPontW.y,FPontW.z-curPontW.z};
	float fcor=0.5;
	if(i==nFPosi||i==nBPosi) fcor=1;
	fdivec[0]=fcor*(fvec1[0]+fvec2[0]);
	fdivec[1]=fcor*(fvec1[1]+fvec2[1]);
	fdivec[2]=fcor*(fvec1[2]+fvec2[2]);
	zxh::VectorOP_Normalise(fdivec,3);
	return true;
}
bool World2Image(float fnPontWor[3],int nPontCor[3],zxhImageDataT<short>&imgReadRaws)
{
	float fWorldCoord[3]={fnPontWor[0],fnPontWor[1],fnPontWor[2]};
	imgReadRaws.GetImageInfo()->WorldToImage(fWorldCoord);
	nPontCor[0] = zxh::round(fWorldCoord[0]);
	nPontCor[1]= zxh::round(fWorldCoord[1]);
	nPontCor[2]= zxh::round(fWorldCoord[2]);
	return true;
}
bool Image2World(int nPontCor[3],float fnPontWor[3],zxhImageDataT<short>&imgReadRaws)
{
	float fImageCoord[3]={nPontCor[0],nPontCor[1],nPontCor[2]};
	imgReadRaws.GetImageInfo()->ImageToWorld(fImageCoord);
	fnPontWor[0] =fImageCoord[0];
	fnPontWor[1]= fImageCoord[1];
	fnPontWor[2]= fImageCoord[2];
	return true;
}
bool fsavePontsDirn(int nPontCor[3],float fdivec[3],PointCordTypeDef &pcor, DireCordTypeDef &pdir)
{
	pcor.x=nPontCor[0];
	pcor.y=nPontCor[1];
	pcor.z=nPontCor[2];
	pdir.dx=fdivec[0];
	pdir.dy=fdivec[1];
	pdir.dz=fdivec[2];
	return true;
}
bool bexitPont(PointCordTypeDef tpcor,vector<PointCordDireTypeDef>&vtmpPontDir,int &ind)
{
	for (int i=0;i<vtmpPontDir.size();i++)
	{
		PointCordDireTypeDef tmpCorandDir=vtmpPontDir[i];
		PointCordTypeDef tmpcor=tmpCorandDir.pcor;
		if (tpcor.x==tmpCorandDir.pcor.x&&tpcor.y==tmpCorandDir.pcor.y&&tpcor.z==tmpCorandDir.pcor.z)
		{
			ind=i;
			return true;
		}
	}

	return false;
}
bool ReadAndSaveThePointsAndDirection(vector<string> vcurvefiles,float fresamle,zxhImageDataT<short>&imgReadRaws,vector<vector<PointCordDireTypeDef>>&vvCurvePontswithDire,vector<PointCordDireTypeDef>&vALLPontDir)
{
	vector<PointCordTypeDef> vPathPointsWorld;
	for(int i=0; i<vcurvefiles.size();i++)//遍历每个curve文件
	{	
		vPathPointsWorld.clear();
		string curvename=vcurvefiles[i];
		ReadTxtAsPC(curvename.c_str(),fresamle, vPathPointsWorld);
		vector<PointCordDireTypeDef> vtmpPontDir;//存储单个curve文件的点和方向
		vtmpPontDir.clear();
		for(int j=0; j<vPathPointsWorld.size();j++)//遍历每个curve文件
		{	
			//计算第j个点的切向量
			float fdivec[3]={0,0,0};
			Calc_divec(j,fdivec,vPathPointsWorld,fresamle);
			int nPontCor[3]={0,0,0};
			float fnPontWor[3]={vPathPointsWorld[j].x,vPathPointsWorld[j].y,vPathPointsWorld[j].z};
			World2Image(fnPontWor,nPontCor,imgReadRaws);

			DireCordTypeDef tpdir;
			PointCordTypeDef tpcor;
			fsavePontsDirn(nPontCor,fdivec,tpcor,tpdir);
			//-------------------------------------按照curve 文件存
			//如果某个点已经存在；
			int ind=0;
			if(bexitPont(tpcor,vtmpPontDir,ind))
			{
				vtmpPontDir[ind].vdire.push_back(tpdir);
			}
			else//如果不存在
			{
				PointCordDireTypeDef tempPonDire;
				vector<DireCordTypeDef>vtpcor;
				vtpcor.push_back(tpdir);
				tempPonDire.pcor=tpcor;
				tempPonDire.vdire=vtpcor;
				vtmpPontDir.push_back(tempPonDire);
			}
			//-------------------------------------按照所有点 存
			int indall=0;
			if(bexitPont(tpcor,vALLPontDir,indall))
			{
				vALLPontDir[indall].vdire.push_back(tpdir);
			}
			else//如果不存在
			{
				PointCordDireTypeDef tempPonDire1;
				vector<DireCordTypeDef>vtpcor1;
				vtpcor1.push_back(tpdir);
				tempPonDire1.pcor=tpcor;
				tempPonDire1.vdire=vtpcor1;
				vALLPontDir.push_back(tempPonDire1);
			}
		}
		vvCurvePontswithDire.push_back(vtmpPontDir);
	}

	return true;
}
bool StorUnique(vector<DireCordTypeDef> &vUnorgaPointsWorld)
{
	bool bunique=true;
	for(int i=0;i<vUnorgaPointsWorld.size();i++)
	{
		float fi[3]={vUnorgaPointsWorld[i].dx,vUnorgaPointsWorld[i].dy,vUnorgaPointsWorld[i].dz};
		for(int j=i+1;j<vUnorgaPointsWorld.size();j++)
		{
			float fj[3]={vUnorgaPointsWorld[j].dx,vUnorgaPointsWorld[j].dy,vUnorgaPointsWorld[j].dz};
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
bool MergeAllTheDirecOnSphere(vector<PointCordDireTypeDef> vALLPontDir,vector<DireCordTypeDef>vDisDir,vector<PointCordDireTypeDef>&vALLPontMergDir)
{
	vector<int>vsizedirc;
	for (int i=0;i<vALLPontDir.size();i++)
	{
		PointCordDireTypeDef tempPonDir;
		tempPonDir=vALLPontDir[i];
		vector<DireCordTypeDef>tmpvdire=tempPonDir.vdire;
		vector<DireCordTypeDef>vtmpmaxdire;
		vtmpmaxdire.clear();
		for (int k=0;k<tmpvdire.size();k++)//逐个访问一个点的所有方向
		{
			DireCordTypeDef tmpdire=tmpvdire[k];
			float fcurdier[3]={tmpdire.dx,tmpdire.dy,tmpdire.dz};
			float fmaxcos=0.00001;
			int nmaxnd=-1;
			for(int nd=0;nd<vDisDir.size();nd++)//跟球面向量比较，找出最接近的球面向量
			{
				float ftmpdire[3]={vDisDir[nd].dx,vDisDir[nd].dy,vDisDir[nd].dz};
				float fcos=zxh::VectorOP_Cosine(ftmpdire,fcurdier,3);
				if(fcos>fmaxcos)
				{
					fmaxcos=fcos;
					nmaxnd=nd;
				}
			}
			DireCordTypeDef bestdire=vDisDir[nmaxnd];
			vtmpmaxdire.push_back(bestdire);
		}
		//合并相同的向量
		StorUnique(vtmpmaxdire);
		//保存成新的坐标点与其合并后的方向
		PointCordTypeDef tmppcor=tempPonDir.pcor;

		PointCordDireTypeDef tmpmerge;
		tmpmerge.pcor=tmppcor;//坐标
		tmpmerge.vdire=vtmpmaxdire;
		int vsize=vtmpmaxdire.size();
		if(vsize==5)
		{
			int x5=0;
		}
			if(vsize==4)
		{
			int x4=0;
		}
				if(vsize==3)
		{
			int x3=0;
		}

		vsizedirc.push_back(vsize);
		vALLPontMergDir.push_back(tmpmerge);

	}
	sort(vsizedirc.begin(),vsizedirc.end(),greater<int>());
	int num5= count(vsizedirc.begin(),vsizedirc.end(),5);
	int num4= count(vsizedirc.begin(),vsizedirc.end(),4);
	int num3= count(vsizedirc.begin(),vsizedirc.end(),3);
	int num2= count(vsizedirc.begin(),vsizedirc.end(),2);
	int num1= count(vsizedirc.begin(),vsizedirc.end(),1);
	return true;
}
int main(int argc, char *argv[])
{
	//这个工具是是用来将曲线的每个点的切向量合并成一个，带有方向的点的文件。
	/*if( argc < 5 )
	{
	cerr << "Usage: " << endl;
	cerr << "jdqMapPont2Img.cpp	imageRaw(.nii)	Line or points(.vtk) results -NorR" << endl;
	return -1;
	}*/

	string strFileNameRaw ="J:/work_jdq/for_DNN_vsls_v5/data/whole_ori/dataset00/image/image.nii.gz";
	string strCurveFilefolder ="J:/work_jdq/for_DNN_vsls_v5/data/whole_ori_jdq/dataset00/centerlinesRadius/";
	string strDisDir="J:/work_jdq/for_DNN_vsls_v5/data/Otherdata/dispersedDirections/dispersedDirectionsSphere500.txt";
	string strResultName ="J:/work_jdq/for_DNN_vsls_v5/infiles";

	//读取图像
	zxhImageDataT<short> imgReadRaws;
	if( zxh::OpenImageSafe( &imgReadRaws, strFileNameRaw ) == false )
	{
		std::cerr << "Raw image(nifti-file) is not found!"<<endl; 
		return -1;
	}
	//读取单位球离散化向量
	vector<DireCordTypeDef>vDisDir;
	const char *chDisDir=strDisDir.c_str();
	ReadDireTxt(chDisDir,vDisDir);
	//读取所有curve文件的名称
	vector<string> vcurvefiles;
	cf_findFileFromDir2(strCurveFilefolder,vcurvefiles);
	//读取最小的spacing
	float fNewImgSpacing[]={1,1,1,1};//Add by JDQ
	imgReadRaws.GetImageSpacing(fNewImgSpacing[0],fNewImgSpacing[1],fNewImgSpacing[2],fNewImgSpacing[3] );//Add by JDQ
	float fresamle=zxh::minf(fNewImgSpacing[1],fNewImgSpacing[2])*0.5;
	//-----------------------1，读取所有点的坐标（按照图像坐标保存），并把属于这个点的方向全部计算后保存-----------。

	vector<vector<PointCordDireTypeDef>> vvCurvePontswithDire;//存储所有curve文件的点和其方向:按照curve文件的数目存取
	vvCurvePontswithDire.clear();
	vector<PointCordDireTypeDef> vALLPontDir;//所有文件的点放在一起
	vALLPontDir.clear();
	ReadAndSaveThePointsAndDirection(vcurvefiles,fresamle,imgReadRaws,vvCurvePontswithDire,vALLPontDir);

	//-----------------------2，将每一个点的坐标和已经确定的球面向量（500个)比对，合并后再保存-----------。
	vector<PointCordDireTypeDef>vALLPontMergeDir;
	MergeAllTheDirecOnSphere(vALLPontDir,vDisDir,vALLPontMergeDir);

	//保存成图像
	/*zxhImageDataT<short> DTIImgx,DTIImgy,DTIImgz;
	imgReadNewRaws.NewImage( imgReadRaws.GetImageInfo() );

	zxh::SaveImage(&imgReadNewRaws,chFileName2.c_str());*/

}

