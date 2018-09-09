

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

#define TAB_CHAR	9
#define M_PI 3.14159265358979323846
#define SPHERE_RADIUS 4
using namespace std;

typedef struct
{
	float x;
	float y;
	float z;
	float b;
}PointCordTypeDef;
typedef struct
{
	PointCordTypeDef p;
	PointCordTypeDef parent;
	PointCordTypeDef lchild;
	PointCordTypeDef rchild;
}BraTypeDef;

bool ConnectPonts(PointCordTypeDef &PpbE,vector<PointCordTypeDef> &vPBPont,vector<PointCordTypeDef> &vPBOrderPont,zxhImageDataT<short>&imgReadRaws)
{

	//find the position
	int nPBPPosi=0;
	for(int i=0;i<vPBPont.size();i++)
	{
		if(vPBPont[i].x==PpbE.x&&vPBPont[i].y==PpbE.y&&vPBPont[i].z==PpbE.z)
		{
			nPBPPosi=i;
		}
	}
	vPBOrderPont.push_back(PpbE);
	//find the nearest points
	vPBPont.erase(vPBPont.begin()+nPBPPosi);
	if(vPBPont.empty())
	{
		return true;	

	}
	else
	{
		float fMinDist=100000;
		float nP[3]={PpbE.x,PpbE.y,PpbE.z};
		PointCordTypeDef PpbNP;
		imgReadRaws.GetImageInfo()->ImageToWorld(nP);
		for(int i=0;i<vPBPont.size();i++)
		{
			float nPN[3]={vPBPont[i].x,vPBPont[i].y,vPBPont[i].z};
			imgReadRaws.GetImageInfo()->ImageToWorld(nPN);
			float fDist=zxh::VectorOP_Distance(nP,nPN,3);
			if(fDist<fMinDist)
			{
				fMinDist=fDist;
				PpbNP=vPBPont[i];
			}
		}
		ConnectPonts(PpbNP,vPBPont,vPBOrderPont,imgReadRaws);

	}



}
bool bnbrnotInnBr(int nbr,vector<int>&vnBr)
{
	for(int nbrannum=0;nbrannum<vnBr.size();nbrannum++)
	{
		int nBr=vnBr[nbrannum];
		if(nBr==nbr)
		{
			return false;
		}

	}

	return true;
}
void WriteModCA2Txt(char *chFileName,vector<PointCordTypeDef>vPBr)
{
	ofstream WriteFileTxt(chFileName,ios::out);
	int nPointNum = vPBr.size();
	for (int i = 0; i < nPointNum; i++)
	{
		
		WriteFileTxt <<right<<fixed<<setfill('4')<<setprecision(4)<<vPBr[i].x<<' '<<vPBr[i].y<<' '<<vPBr[i].z<<'\n';

	}	

}
bool GenBranch(vector<PointCordTypeDef> vPSEPont,vector<BraTypeDef>& vBraPonts)
{
	for(int nPSENUM=0;nPSENUM<vPSEPont.size();nPSENUM=nPSENUM+2)
	{
		for(int nPSENUM=0;nPSENUM<vPSEPont.size();nPSENUM=nPSENUM+2)
		{
			PointCordTypeDef PCur;
			PointCordTypeDef PNexCur;
			PCur=vPSEPont[nPSENUM];
			PNexCur=vPSEPont[nPSENUM+1];
			BraTypeDef BrPont;
			BrPont.p=PCur;
			BrPont.parent=PCur;
			if(PCur.b<PNexCur.b)
			{
				BrPont.lchild=PNexCur;

			}
			else
			{
				cout<<"Manual points are wrong"<<endl;
			}

		}
	}
	return true;
}
void WriteVtk(vector< PointCordTypeDef > PointCord, char* chFileName)
{
	vtkSmartPointer<vtkPoints> iPoints = vtkSmartPointer<vtkPoints>::New();

	/*	int nPointNum = PointCord.size();*/

	float fImgPixel[3];
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


void GetPathname(char *chFileName,string &strname)
{  char path_buffer[_MAX_PATH];  
char drive[_MAX_DRIVE];  
char dir[_MAX_DIR];  
char fname[_MAX_FNAME];  
char ext[_MAX_EXT];  

_splitpath( chFileName, drive, dir, fname, ext );  
strname=fname;
}
void GetPathdir(char *chFileName,char *chdir)
{  char path_buffer[_MAX_PATH];  
char drive[_MAX_DRIVE];  
char dir[_MAX_DIR];  
char fname[_MAX_FNAME];  
char ext[_MAX_EXT];  

_splitpath( chFileName, drive, dir, fname, ext );  
chdir=dir;
}

int main(int argc, char *argv[])//此程序用来排序从图像中读取的点，输出的是图像坐标
{
	//if( argc < 4 )
	//{
	//	cerr << "Usage: " << endl;
	//	cerr << "zxhcaeDMPMToNewImg.cpp	imageRaw(.nii)	imageResoRaw(.nii) MCLine(.vtk) ResultHigh-ResolutionName(.nii.gz) ResultLow-ResolutionName(.nii.gz)" << endl;
	//	return -1;
	//}
	//string strFileNameRaw =string(argv[1]);//  "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
	//char *chFileName =argv[2];//"F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/reference.vtk";
	//char *chResultName =argv[3];// "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/MCLine.nii.gz";
	//string RorN=string(argv[4]);

	char *strFileNameRaw ="E:/work_jdq/CTA_data/RCAAEF/training_jdq/01/image01_R0_1.nii.gz";
	char *chResFilefold ="E:/work_jdq/CTA_data/RCAAEF/training_jdq/01/";


	zxhImageDataT<short> imgReadRaws;

	if( zxh::OpenImage( &imgReadRaws, strFileNameRaw ) == false )
	{
		std::cerr << "Raw image(nifti-file) is not found!"; 
		return -1;
	}
	//const float consfMindist=6;//right
	const float consfMindist=4;//left
	int nsx,nsy,nsz,nst;
	vector<PointCordTypeDef> vPKPont;
	vPKPont.clear();

	imgReadRaws.GetImageSize(nsx,nsy,nsz,nst);
	PointCordTypeDef PpbE;
	//store the manul segment points
	for(int nx=0;nx<nsx;nx++)
		for(int ny=0;ny<nsy;ny++)
			for(int nz=0;nz<nsz;nz++)
			{
				short shInt=imgReadRaws.GetPixelGreyscale(nx,ny,nz,0);
				if(shInt!=0)
				{
					PointCordTypeDef PKeyPont;
					PKeyPont.x=nx;
					PKeyPont.y=ny;
					PKeyPont.z=nz;
					PKeyPont.b=shInt;
					vPKPont.push_back(PKeyPont);
					if(nx==169&&ny==192&&nz==216)
					{
						PpbE=PKeyPont;
					}
				}

			}
			//定义叶子点
			vector<PointCordTypeDef> vLePo;
			vLePo.clear();
			vLePo.push_back(PpbE);
			//定义排序后点
			vector<PointCordTypeDef>vOrderPonts;
			vOrderPonts.clear();
			//vOrderPonts.push_back(PpbE);
			//找到离叶子节点距离最近的两个点
			vector<float> vDist;
			vDist.clear();


			bool bem=true;
			while (bem)
			{
				PointCordTypeDef LePo;
				float fMinDist=100000;
				int nMinIndex=-1;
				if(vPKPont.empty())
					break;
				float nP[3]={0,0,0};
				PointCordTypeDef pmindstP;
				for(int i=0;i<vPKPont.size();i++)//找出距离当前点最近的点
				{
					nP[0]=vPKPont[i].x;
					nP[1]=vPKPont[i].y;
					nP[2]=vPKPont[i].z;
					float nPN[3]={PpbE.x,PpbE.y,PpbE.z};
					imgReadRaws.GetImageInfo()->ImageToWorld(nP);
					imgReadRaws.GetImageInfo()->ImageToWorld(nPN);
					float fDist=zxh::VectorOP_Distance(nP,nPN,3);
					if(fDist<fMinDist)
					{
						fMinDist=fDist;
						nMinIndex=i;
					
					}
				}

				//保存找到的点，并从原始容器中移除
				//cout<<"current point"<<PpbE.x<<","<<PpbE.y<<","<<PpbE.z<<" To-"<<endl;
				//cout<<pmindstP.x<<","<<pmindstP.y<<","<<pmindstP.z<<":"<<pmindstP.b<<endl;
					float nminP[3]={vPKPont[nMinIndex].x,vPKPont[nMinIndex].y,vPKPont[nMinIndex].z};
					imgReadRaws.GetImageInfo()->ImageToWorld(nminP);
					pmindstP.x=nminP[0];
						pmindstP.y=nminP[1];
						pmindstP.z=nminP[2];
						pmindstP.b=fMinDist;
				vOrderPonts.push_back(pmindstP);
				PpbE=vPKPont[nMinIndex];
				vPKPont.erase(vPKPont.begin()+nMinIndex);      
				
			}
			//write to vtk
			string strname="";
			GetPathname(strFileNameRaw,strname);
			string stRname = strname.substr(0, strname.length() - 4);
			const char *chRname=stRname.c_str();
			//if(strcmp(chext,".vtk")==0)// read vtk file

			//	char chTemp[25];
			int nFileLen = strlen(chResFilefold) + strlen("OrdP_")+strlen(chRname) + strlen(".txt") + 1;
			char *chFileName = (char*)malloc(nFileLen);
			strcpy(chFileName, chResFilefold);
			strcat(chFileName, "OrdP_");
			strcat(chFileName, chRname);
			strcat(chFileName, ".txt");
			WriteModCA2Txt(chFileName,vOrderPonts);

}




