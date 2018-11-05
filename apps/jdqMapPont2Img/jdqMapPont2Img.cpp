

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
float ArcDist(int StartPosi,int EndPosi,vector<PointCordTypeDef> vPathPoints)
{
	float fFPoint[3]={0},fBPoint[3]={0};
	float fLength=0;
	if (StartPosi=EndPosi)fLength=0;
	else
	{
	for (int j=StartPosi;j<EndPosi;j++)
	{
		fFPoint[0]=vPathPoints[j].x;
		fFPoint[1]=vPathPoints[j].y;
		fFPoint[2]=vPathPoints[j].z;
		fBPoint[0]=vPathPoints[j+1].x;
		fBPoint[1]=vPathPoints[j+1].y;
		fBPoint[2]=vPathPoints[j+1].z;
		fLength=fLength+zxh::VectorOP_Distance(fFPoint,fFPoint,3);
	}
	}
	return fLength;
}
bool InSphere(float PosWorld[3],vector<PointCordTypeDef> vPathPoints)
{
	float MinDist=0;
	for (int i=0;i<vPathPoints.size();i++)
	{
		float PosModWorld[3]={vPathPoints[i].x,vPathPoints[i].y,vPathPoints[i].z};
		float Distmm2=zxh::VectorOP_Distance(PosWorld,PosModWorld,3);
       if(Distmm2<(float)SPHERE_RADIUS) return true;
	}
		return false;
	
}

bool VesselnessMaskBaseOnDist(vector<PointCordTypeDef> vPathPointsWorld,zxhImageDataT<short>&imgReadVsls,zxhImageDataT<short>&imgReadNewVsls)
{
	int ImgNewVslsSize[4]={1};
	int ImgVslsSize[4]={1};
	imgReadNewVsls.GetImageSize(ImgNewVslsSize[0],ImgNewVslsSize[1],ImgNewVslsSize[2],ImgNewVslsSize[3]);
	imgReadVsls.GetImageSize(ImgVslsSize[0],ImgVslsSize[1],ImgVslsSize[2],ImgVslsSize[3]);
	for(int it=0;it<ImgNewVslsSize[3];++it)
	for(int iz=0;iz<ImgNewVslsSize[2];++iz)
	for(int iy=0;iy<ImgNewVslsSize[1];++iy)
	for(int ix=0;ix<ImgNewVslsSize[0];++ix)
	{
		float Pos[ImageDimensionMax]={ix,iy,iz,it};
		imgReadNewVsls.GetImageInfo()->ImageToWorld(Pos);
		float PosWorld[ImageDimensionMax]={Pos[0],Pos[1],Pos[2],Pos[3]};
		if(!(InSphere(PosWorld,vPathPointsWorld)))
		{
			imgReadNewVsls.SetPixelByGreyscale(ix,iy,iz,it,-1);
		}
		else
		{//iMinNode.z * m_nImgWY * m_nImgWX + iMinNode.y * m_nImgWX + iMinNode.x
			short m=imgReadVsls.GetImageData()[iz*ImgVslsSize[1]*ImgVslsSize[0]+iy*ImgVslsSize[0]+ix];
			imgReadNewVsls.SetPixelByGreyscale(ix,iy,iz,it,imgReadVsls.GetImageData()[iz*ImgVslsSize[1]*ImgVslsSize[0]+iy*ImgVslsSize[0]+ix]);
		}

	}
	return true;
}
bool InterpoforNewImage(vector<PointCordTypeDef> vPathPointsWorld,zxhImageDataT<short>&imgReadNewVsls)
{
	
	for(int i=1;i<vPathPointsWorld.size()-1;i++)
	{
	float PointFPosWorld[3]={0};
	float PointBPosWorld[3]={0};
	int PointFPos[3]={0};
	int PointBPos[3]={0};
	PointFPosWorld[0]=vPathPointsWorld[i-1].x;
	PointFPosWorld[1]=vPathPointsWorld[i-1].y;
	PointFPosWorld[2]=vPathPointsWorld[i-1].z;
	PointBPosWorld[0]=vPathPointsWorld[i].x;
	PointBPosWorld[1]=vPathPointsWorld[i].y;
	PointBPosWorld[2]=vPathPointsWorld[i].z;
	imgReadNewVsls.GetImageInfo()->WorldToImage(PointFPosWorld);
	PointFPos[0]=(int)(PointFPosWorld[0]+0.5);
	PointFPos[1]=(int)(PointFPosWorld[1]+0.5);
	PointFPos[2]=(int)(PointFPosWorld[2]+0.5);
	imgReadNewVsls.GetImageInfo()->WorldToImage(PointBPosWorld);
	PointBPos[0]=(int)(PointFPosWorld[0]+0.5);
	PointBPos[1]=(int)(PointFPosWorld[1]+0.5);
	PointBPos[2]=(int)(PointFPosWorld[2]+0.5);
	if (PointFPos[0]==PointBPos[0]&&PointFPos[1]==PointBPos[1]&&PointFPos[2]==PointBPos[2]) continue;
	int gNbr[6][3] = { {-1, 0, 0}, \
					   { 1, 0, 0}, \
					   { 0,-1, 0}, \
					   { 0, 1, 0}, \
					   { 0, 0,-1}, \
					   { 0, 0, 1} };
	for (int i = 0; i < 6; i++)
		{
			int nx = PointFPos[0] + gNbr[i][0];
			int ny = PointFPos[1] + gNbr[i][1];
			int nz = PointFPos[2] + gNbr[i][2];



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
bool MapModelPointsToNewImageItself(vector<PointCordTypeDef> vPathPointsWorld,zxhImageDataT<short>&imgReadNewVsls)
{
	int ImgNewVslsSize[4]={1};
	imgReadNewVsls.GetImageSize(ImgNewVslsSize[0],ImgNewVslsSize[1],ImgNewVslsSize[2],ImgNewVslsSize[3]);
	for(int it=0;it<ImgNewVslsSize[3];++it)
	for(int iz=0;iz<ImgNewVslsSize[2];++iz)
	for(int iy=0;iy<ImgNewVslsSize[1];++iy)
	for(int ix=0;ix<ImgNewVslsSize[0];++ix)
	{
		imgReadNewVsls.SetPixelByGreyscale(ix,iy,iz,it,0);
	}
	
	for (int i=0;i<vPathPointsWorld.size();i++)
	{
		float PointPosWorld[ImageDimensionMax]={0};
		int PointPos[4]={0};
		PointPosWorld[0]=vPathPointsWorld[i].x;
		PointPosWorld[1]=vPathPointsWorld[i].y;
		PointPosWorld[2]=vPathPointsWorld[i].z;
		imgReadNewVsls.GetImageInfo()->WorldToImage(PointPosWorld);
		PointPos[0]=zxh::round(PointPosWorld[0]);
		PointPos[1]=(int)(PointPosWorld[1]+0.5);
		PointPos[2]=(int)(PointPosWorld[2]+0.5);
		BoundaryCorrect(PointPos,ImgNewVslsSize);
		imgReadNewVsls.SetPixelByGreyscale(PointPos[0],PointPos[1],PointPos[2],PointPos[3],ZXH_Foreground);
	}
	return true;
}
bool MapModelPointsToNewImage(vector<PointCordTypeDef> vPathPointsWorld,zxhImageDataT<short>&imgReadNewVsls)
{
	int ImgNewVslsSize[4]={1};
	imgReadNewVsls.GetImageSize(ImgNewVslsSize[0],ImgNewVslsSize[1],ImgNewVslsSize[2],ImgNewVslsSize[3]);
	float ImgNewVslsSpacing[4]={0};
	imgReadNewVsls.GetImageSpacing(ImgNewVslsSpacing[0],ImgNewVslsSpacing[1],ImgNewVslsSpacing[2],ImgNewVslsSpacing[3]);
	float fminpixdist=0.5*sqrt(ImgNewVslsSpacing[0]*ImgNewVslsSpacing[0]+ImgNewVslsSpacing[1]*ImgNewVslsSpacing[1]+ImgNewVslsSpacing[2]*ImgNewVslsSpacing[2]);
	for(int it=0;it<ImgNewVslsSize[3];++it)
	for(int iz=0;iz<ImgNewVslsSize[2];++iz)
	for(int iy=0;iy<ImgNewVslsSize[1];++iy)
	for(int ix=0;ix<ImgNewVslsSize[0];++ix)
	{
		imgReadNewVsls.SetPixelByGreyscale(ix,iy,iz,it,0);
	}
	vector<PointCordTypeDef> vPathPointsWorldMAPT;
	vPathPointsWorldMAPT.clear();
	for (int i=0;i<vPathPointsWorld.size();i++)
	{
		float PointPosWorld[ImageDimensionMax]={0};
		int PointPos[4]={0};
		PointCordTypeDef PointMAPT;
		PointPosWorld[0]=vPathPointsWorld[i].x;
		PointPosWorld[1]=vPathPointsWorld[i].y;
		PointPosWorld[2]=vPathPointsWorld[i].z;
		imgReadNewVsls.GetImageInfo()->WorldToImage(PointPosWorld);
		PointPos[0]=zxh::round(PointPosWorld[0]);
		PointPos[1]=(int)(PointPosWorld[1]+0.5);
		PointPos[2]=(int)(PointPosWorld[2]+0.5);
		BoundaryCorrect(PointPos,ImgNewVslsSize);
		PointMAPT.x=PointPos[0];
		PointMAPT.y=PointPos[1];
		PointMAPT.z=PointPos[2];
		vPathPointsWorldMAPT.push_back(PointMAPT);
		imgReadNewVsls.SetPixelByGreyscale(PointPos[0],PointPos[1],PointPos[2],PointPos[3],ZXH_Foreground);
	}

	for (int mapNUM=1;mapNUM<vPathPointsWorldMAPT.size()-1;mapNUM++)//��ֵ
	{
		PointCordTypeDef PointMAPTS;//start point in every step
		PointCordTypeDef PointMAPTE;//End point in every step
		PointMAPTS=vPathPointsWorldMAPT[mapNUM-1];
		PointMAPTE=vPathPointsWorldMAPT[mapNUM];
		for(int jz=PointMAPTS.z-8;jz<PointMAPTS.z+8;++jz)
			for(int jy=PointMAPTS.y-8;jy<PointMAPTS.y+8;++jy)
				for(int jx=PointMAPTS.x-8;jx<PointMAPTS.x+8;++jx)
				{
					float fL[4]={0};
					float fM[4]={0};
					float fN[4]={0};
					fL[0]=PointMAPTE.x-PointMAPTS.x;
					fL[1]=PointMAPTE.y-PointMAPTS.y;
					fL[2]=PointMAPTE.z-PointMAPTS.z;
					fM[0]=jx-PointMAPTS.x;
					fM[1]=jy-PointMAPTS.y;
					fM[2]=jz-PointMAPTS.z;
					fN[0]=jx-PointMAPTE.x;
					fN[1]=jy-PointMAPTE.y;
					fN[2]=jz-PointMAPTE.z;
					float fcosineLM=zxh::VectorOP_Cosine(fL,fM,3);
					float fcosineLN=-zxh::VectorOP_Cosine(fL,fN,3);
					float fLengthM=zxh::MagnitudeOfVector(fM,3);
					float fDist=fLengthM*sqrt(1-fcosineLM*fcosineLM);
					if(fcosineLM>0.9&&fcosineLN>0.9&&fDist<=fminpixdist)
					{
						imgReadNewVsls.SetPixelByGreyscale(jx,jy,jz,0,ZXH_Foreground);
					}
				}
	}
	return true;
}
bool MapModelPointsToNewImageWithinR(vector<PointCordTypeDef> vPathPointsWorld,zxhImageDataT<short>&imgReadNewVsls,int SearchRange[4])//map line into original image in a range
{
	int ImgNewVslsSize[4]={1};
	imgReadNewVsls.GetImageSize(ImgNewVslsSize[0],ImgNewVslsSize[1],ImgNewVslsSize[2],ImgNewVslsSize[3]);
	float ImgNewVslsSpacing[4]={0};
	imgReadNewVsls.GetImageSpacing(ImgNewVslsSpacing[0],ImgNewVslsSpacing[1],ImgNewVslsSpacing[2],ImgNewVslsSpacing[3]);
	float fminpixdist=sqrt(ImgNewVslsSpacing[0]*ImgNewVslsSpacing[0]+ImgNewVslsSpacing[1]*ImgNewVslsSpacing[1]+ImgNewVslsSpacing[2]*ImgNewVslsSpacing[2]);
	for(int it=0;it<ImgNewVslsSize[3];++it)
	for(int iz=0;iz<ImgNewVslsSize[2];++iz)
	for(int iy=0;iy<ImgNewVslsSize[1];++iy)
	for(int ix=0;ix<ImgNewVslsSize[0];++ix)
	{
		imgReadNewVsls.SetPixelByGreyscale(ix,iy,iz,it,0);
	}
	vector<PointCordTypeDef> vPathPointsWorldMAPT;
	vPathPointsWorldMAPT.clear();
	for (int i=0;i<vPathPointsWorld.size();i++)
	{
		float PointPosWorld[ImageDimensionMax]={0};
		int PointPos[4]={0};
		PointCordTypeDef PointMAPT;
		PointPosWorld[0]=vPathPointsWorld[i].x;
		PointPosWorld[1]=vPathPointsWorld[i].y;
		PointPosWorld[2]=vPathPointsWorld[i].z;
		imgReadNewVsls.GetImageInfo()->WorldToImage(PointPosWorld);
		PointPos[0]=zxh::round(PointPosWorld[0]);
		PointPos[1]=zxh::round(PointPosWorld[1]);
		PointPos[2]=zxh::round(PointPosWorld[2]);
		BoundaryCorrect(PointPos,ImgNewVslsSize);
		PointMAPT.x=PointPos[0];
		PointMAPT.y=PointPos[1];
		PointMAPT.z=PointPos[2];
		vPathPointsWorldMAPT.push_back(PointMAPT);
		for(int kx=PointPos[0]-SearchRange[0];kx<PointPos[0]+SearchRange[0];kx++)
			for(int ky=PointPos[1]-SearchRange[1];ky<PointPos[1]+SearchRange[1];ky++)
				for(int kz=PointPos[2]-SearchRange[2];kz<PointPos[2]+SearchRange[2];kz++)
				{
					int kPos[4]={0};
					kPos[0]=zxh::round(kx);
					kPos[1]=zxh::round(ky);
					kPos[2]=zxh::round(kz);
					BoundaryCorrect(kPos,ImgNewVslsSize);
					imgReadNewVsls.SetPixelByGreyscale(kPos[0],kPos[1],kPos[2],kPos[3],ZXH_Foreground);
				}
	}

	for (int mapNUM=1;mapNUM<vPathPointsWorldMAPT.size()-1;mapNUM++)//��ֵ
	{
		PointCordTypeDef PointMAPTS;//start point in every step
		PointCordTypeDef PointMAPTE;//End point in every step
		PointMAPTS=vPathPointsWorldMAPT[mapNUM-1];
		PointMAPTE=vPathPointsWorldMAPT[mapNUM];
		for(int jz=PointMAPTS.z-8;jz<PointMAPTS.z+8;++jz)
			for(int jy=PointMAPTS.y-8;jy<PointMAPTS.y+8;++jy)
				for(int jx=PointMAPTS.x-8;jx<PointMAPTS.x+8;++jx)
				{
					float fL[4]={0};
					float fM[4]={0};
					float fN[4]={0};
					fL[0]=PointMAPTE.x-PointMAPTS.x;
					fL[1]=PointMAPTE.y-PointMAPTS.y;
					fL[2]=PointMAPTE.z-PointMAPTS.z;
					fM[0]=jx-PointMAPTS.x;
					fM[1]=jy-PointMAPTS.y;
					fM[2]=jz-PointMAPTS.z;
					fN[0]=jx-PointMAPTE.x;
					fN[1]=jy-PointMAPTE.y;
					fN[2]=jz-PointMAPTE.z;
					float fcosineLM=zxh::VectorOP_Cosine(fL,fM,3);
					float fcosineLN=-zxh::VectorOP_Cosine(fL,fN,3);
					float fLengthM=zxh::MagnitudeOfVector(fM,3);
					float fDist=fLengthM*sqrt(1-fcosineLM*fcosineLM);
					if(fcosineLM>0.9&&fcosineLN>0.9&&fDist<=fminpixdist)
					{
						for(int kx=jx-SearchRange[0];kx<jx+SearchRange[0];kx++)
							for(int ky=jy-SearchRange[1];ky<jy+SearchRange[1];ky++)
								for(int kz=jz-SearchRange[2];kz<jz+SearchRange[2];kz++)
								{
									int kPos[4]={0};
									kPos[0]=zxh::round(kx);
									kPos[1]=zxh::round(ky);
									kPos[2]=zxh::round(kz);
									BoundaryCorrect(kPos,ImgNewVslsSize);
									imgReadNewVsls.SetPixelByGreyscale(kPos[0],kPos[1],kPos[2],kPos[3],ZXH_Foreground);
								}
					}
				}
	}
	return true;
}
bool MapModelPointsToNewImageNew(vector<PointCordTypeDef> vPathPointsWorld,zxhImageDataT<short>&imgReadNewVsls,int SearchRange[4])//map line into original image in a range
{
	int ImgNewVslsSize[4]={1};
	imgReadNewVsls.GetImageSize(ImgNewVslsSize[0],ImgNewVslsSize[1],ImgNewVslsSize[2],ImgNewVslsSize[3]);
	float ImgNewVslsSpacing[4]={0};
	imgReadNewVsls.GetImageSpacing(ImgNewVslsSpacing[0],ImgNewVslsSpacing[1],ImgNewVslsSpacing[2],ImgNewVslsSpacing[3]);
	float fminpixdist=1*sqrt(ImgNewVslsSpacing[0]*ImgNewVslsSpacing[0]+ImgNewVslsSpacing[1]*ImgNewVslsSpacing[1]+ImgNewVslsSpacing[2]*ImgNewVslsSpacing[2]);
	vector<PointCordTypeDef> vPathPointsWorldMAPT;
	vPathPointsWorldMAPT.clear();
	for (int i=0;i<vPathPointsWorld.size();i++)
	{
		float PointPosWorld[ImageDimensionMax]={0};
		int PointPos[4]={0};
		PointCordTypeDef PointMAPT;
		PointPosWorld[0]=vPathPointsWorld[i].x;
		PointPosWorld[1]=vPathPointsWorld[i].y;
		PointPosWorld[2]=vPathPointsWorld[i].z;
		imgReadNewVsls.GetImageInfo()->WorldToImage(PointPosWorld);
		PointPos[0]=zxh::round(PointPosWorld[0]);
		PointPos[1]=zxh::round(PointPosWorld[1]);
		PointPos[2]=zxh::round(PointPosWorld[2]);
		BoundaryCorrect(PointPos,ImgNewVslsSize);
		PointMAPT.x=PointPos[0];
		PointMAPT.y=PointPos[1];
		PointMAPT.z=PointPos[2];
		vPathPointsWorldMAPT.push_back(PointMAPT);
		imgReadNewVsls.SetPixelByGreyscale(PointPos[0],PointPos[1],PointPos[2],PointPos[3],ZXH_Foreground);

	}

	for (int mapNUM=1;mapNUM<vPathPointsWorldMAPT.size()-1;mapNUM++)//��ֵ
	{
		PointCordTypeDef PointMAPTS;//start point in every step
		PointCordTypeDef PointMAPTE;//End point in every step
		PointMAPTS=vPathPointsWorldMAPT[mapNUM-1];
		PointMAPTE=vPathPointsWorldMAPT[mapNUM];
		for(int jz=PointMAPTS.z-8;jz<PointMAPTS.z+8;++jz)
			for(int jy=PointMAPTS.y-8;jy<PointMAPTS.y+8;++jy)
				for(int jx=PointMAPTS.x-8;jx<PointMAPTS.x+8;++jx)
				{
					float fL[4]={0};
					float fM[4]={0};
					float fN[4]={0};
					fL[0]=PointMAPTE.x-PointMAPTS.x;
					fL[1]=PointMAPTE.y-PointMAPTS.y;
					fL[2]=PointMAPTE.z-PointMAPTS.z;
					fM[0]=jx-PointMAPTS.x;
					fM[1]=jy-PointMAPTS.y;
					fM[2]=jz-PointMAPTS.z;
					fN[0]=jx-PointMAPTE.x;
					fN[1]=jy-PointMAPTE.y;
					fN[2]=jz-PointMAPTE.z;
					float fcosineLM=zxh::VectorOP_Cosine(fL,fM,3);
					float fcosineLN=-zxh::VectorOP_Cosine(fL,fN,3);
					float fLengthM=zxh::MagnitudeOfVector(fM,3);
					float fDist=fLengthM*sqrt(1-fcosineLM*fcosineLM);
					if(fcosineLM>0.9&&fcosineLN>0.9&&fDist<=fminpixdist)
					{
						imgReadNewVsls.SetPixelByGreyscale(jx,jy,jz,0,ZXH_Foreground);		
					}
				}
	}
	return true;
}
bool MapModelPointsToImage(vector<PointCordTypeDef> vPathPointsWorld,zxhImageDataT<short>&imgReadNewVsls,int SearchRange[4])//map line into original image in a range
{
	int ImgNewVslsSize[4]={1};
	imgReadNewVsls.GetImageSize(ImgNewVslsSize[0],ImgNewVslsSize[1],ImgNewVslsSize[2],ImgNewVslsSize[3]);
	float ImgNewVslsSpacing[4]={0};
	float fminpixdist=1*sqrt(ImgNewVslsSpacing[0]*ImgNewVslsSpacing[0]+ImgNewVslsSpacing[1]*ImgNewVslsSpacing[1]+ImgNewVslsSpacing[2]*ImgNewVslsSpacing[2]);
	vector<PointCordTypeDef> vPathPointsWorldMAPT;
	vPathPointsWorldMAPT.clear();
	for (int i=0;i<vPathPointsWorld.size();i++)//map the vetor points to the image
	{
		float PointPosWorld[ImageDimensionMax]={0};
		int PointPos[4]={0};
		PointCordTypeDef PointMAPT;
		PointPosWorld[0]=vPathPointsWorld[i].x;
		PointPosWorld[1]=vPathPointsWorld[i].y;
		PointPosWorld[2]=vPathPointsWorld[i].z;
		//cout<<"WPoints:"<<i<<":" <<PointPosWorld[0]<<" "<<PointPosWorld[1]<<" "<<PointPosWorld[2]<<" "<<endl;
		imgReadNewVsls.GetImageInfo()->WorldToImage(PointPosWorld);
		
	    PointPos[0]=zxh::round(PointPosWorld[0]);
		PointPos[1]=zxh::round(PointPosWorld[1]);
		PointPos[2]=zxh::round(PointPosWorld[2]);

		BoundaryCorrect(PointPos,ImgNewVslsSize);
		PointMAPT.x=PointPos[0];
		PointMAPT.y=PointPos[1];
		PointMAPT.z=PointPos[2];
		vPathPointsWorldMAPT.push_back(PointMAPT);
		//cout<<"Points:" <<PointPos[0]<<" "<<PointPos[1]<<" "<<PointPos[2]<<" "<<endl;
		if (PointPos[0]==0&&PointPos[1]==82&&PointPos[2]==76)
			int x678=0;
		imgReadNewVsls.SetPixelByGreyscale(PointPos[0],PointPos[1],PointPos[2],PointPos[3],ZXH_Foreground);

	}
	return true;
}

bool MapModelPointsToImage_OCC(vector<PointCordTypeDef> vPathPointsWorld,zxhImageDataT<short>&imgReadNewVsls,int SearchRange[4])//map line into original image in a range
{
	int ImgNewVslsSize[4]={1};
	imgReadNewVsls.GetImageSize(ImgNewVslsSize[0],ImgNewVslsSize[1],ImgNewVslsSize[2],ImgNewVslsSize[3]);
	float ImgNewVslsSpacing[4]={0};
	float fminpixdist=1*sqrt(ImgNewVslsSpacing[0]*ImgNewVslsSpacing[0]+ImgNewVslsSpacing[1]*ImgNewVslsSpacing[1]+ImgNewVslsSpacing[2]*ImgNewVslsSpacing[2]);
	vector<PointCordTypeDef> vPathPointsWorldMAPT;
	vPathPointsWorldMAPT.clear();
	for (int i=0;i<vPathPointsWorld.size();i++)//map the vetor points to the image
	{

		if(i>2266&&i<2424) continue;
		float PointPosWorld[ImageDimensionMax]={0};
		int PointPos[4]={0};
		PointCordTypeDef PointMAPT;
		PointPosWorld[0]=vPathPointsWorld[i].x;
		PointPosWorld[1]=vPathPointsWorld[i].y;
		PointPosWorld[2]=vPathPointsWorld[i].z;
		imgReadNewVsls.GetImageInfo()->WorldToImage(PointPosWorld);
		PointPos[0]=zxh::round(PointPosWorld[0]);
		PointPos[1]=zxh::round(PointPosWorld[1]);
		PointPos[2]=zxh::round(PointPosWorld[2]);
		BoundaryCorrect(PointPos,ImgNewVslsSize);
		PointMAPT.x=PointPos[0];
		PointMAPT.y=PointPos[1];
		PointMAPT.z=PointPos[2];
		vPathPointsWorldMAPT.push_back(PointMAPT);
		imgReadNewVsls.SetPixelByGreyscale(PointPos[0],PointPos[1],PointPos[2],PointPos[3],ZXH_Foreground);

	}

	for (int mapNUM=1;mapNUM<vPathPointsWorldMAPT.size()-1;mapNUM++)//��ֵ
	{
		PointCordTypeDef PointMAPTS;//start point in every step
		PointCordTypeDef PointMAPTE;//End point in every step
		PointMAPTS=vPathPointsWorldMAPT[mapNUM-1];
		PointMAPTE=vPathPointsWorldMAPT[mapNUM];
		for(int jz=PointMAPTS.z-8;jz<PointMAPTS.z+8;++jz)
			for(int jy=PointMAPTS.y-8;jy<PointMAPTS.y+8;++jy)
				for(int jx=PointMAPTS.x-8;jx<PointMAPTS.x+8;++jx)
				{
					float fL[4]={0};
					float fM[4]={0};
					float fN[4]={0};
					fL[0]=PointMAPTE.x-PointMAPTS.x;
					fL[1]=PointMAPTE.y-PointMAPTS.y;
					fL[2]=PointMAPTE.z-PointMAPTS.z;
					fM[0]=jx-PointMAPTS.x;
					fM[1]=jy-PointMAPTS.y;
					fM[2]=jz-PointMAPTS.z;
					fN[0]=jx-PointMAPTE.x;
					fN[1]=jy-PointMAPTE.y;
					fN[2]=jz-PointMAPTE.z;
					float fcosineLM=zxh::VectorOP_Cosine(fL,fM,3);
					float fcosineLN=-zxh::VectorOP_Cosine(fL,fN,3);
					float fLengthM=zxh::MagnitudeOfVector(fM,3);
					float fDist=fLengthM*sqrt(1-fcosineLM*fcosineLM);
					if(fcosineLM>0.9&&fcosineLN>0.9&&fDist<=fminpixdist)
					{
						imgReadNewVsls.SetPixelByGreyscale(jx,jy,jz,0,ZXH_Foreground);		
					}
				}
	}
	return true;
}
bool VesselnessMaskBaseOnStructElement(vector<PointCordTypeDef> vPathPoints,zxhImageDataT<short>&imgReadVsls,zxhImageDataT<short>&imgReadNewVsls)
{
	InterpoforNewImage(vPathPoints,imgReadNewVsls);
	//DilateforNewImage();
	return true;
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
		cout << "Cannot find VTK-file!" << endl;
		return 1;
	}
	//����jda point��ʽ��ȡline�ļ�
	std::vector<jdq2017::point3D> ref1,ref;
	std::vector<double> rad;
	std::vector<double> io;
	std::vector<jdq2017::point3D> cl;
	if ( ! jdq2017::readReference(chFileName, ref1, rad, io))
	{
		std::cerr << "Error in reading line data" << std::endl;
		return 1;
	}
	//resample the points
	jdq2017::ResamplePaths<jdq2017::point3D>Respler;
	Respler(ref1,fresamle);
	ref=Respler.resultPath();
	//��ref points ������ת���PointCord
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
bool ReadTxtAsPCWRes(char *chFileName, vector<PointCordTypeDef> &PointCord)
{
	if (chFileName == NULL)
	{
		cout << "Cannot find VTK-file!" << endl;
		return 1;
	}
	//����jda point��ʽ��ȡline�ļ�
	std::vector<jdq2017::point3D>ref;
	if ( ! jdq2017::readCenterline(chFileName, ref))
	{
		std::cerr << "Error in reading line data" << std::endl;
		return 1;
	}
	//��ref points ������ת���PointCord
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
			float fx,fy,fz;
			iFileIn>>fx>>fy>>fz;
			cout<<"WPoints:" <<fx<<" "<<fy<<" "<<fz<<" "<<endl;
			strctTempPoint.x=-1*fx;
			strctTempPoint.y=-1*fy;
			strctTempPoint.z=fz;
			cl.push_back(strctTempPoint);
		}
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
int main(int argc, char *argv[])
{
	if( argc < 5 )
	{
		cerr << "Usage: " << endl;
		cerr << "jdqMapPont2Img.cpp	imageRaw(.nii)	Line or points(.vtk) results -NorR" << endl;
		return -1;
	}
	string strFileNameRaw =string(argv[1]);//  "F:/Coronary_0/Coronary_Niessen/CoronaryMasks/DataMaskOutLung/mol_image00.nii";
	char *chFileName =argv[2];//"F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/reference.vtk";
	char *chResultName =argv[3];// "F:/Coronary_0/Coronary_Niessen/ProcessByLL/training/dataset00/vessel2/MCLine.nii.gz";
	string RorN=string(argv[4]);
	string RorS=string(argv[5]);//resample the line or not
	
	
	//***Training Data as model****//
	
	zxhImageDataT<short> imgReadRaws,imgReadResoRaws;//Change by JDQ
	float fNewImgSpacing[]={1,1,1,1};//Add by JDQ
	float fOldImgSpacing[]={1,1,1,1};
	if( zxh::OpenImage( &imgReadRaws, strFileNameRaw ) == false )
	{
		std::cerr << "Raw image(nifti-file) is not found!"<<endl; 
		return -1;
	}
	/*if( zxh::OpenImage( &imgReadResoRaws, strFileNameResoRaw ) == false )
	{
		std::cerr << "Resoluted Raw image(nifti-file) is not found!"; 
		return -1;
	}*/
	//imgReadResoRaws.GetImageSpacing(fOldImgSpacing[0],fOldImgSpacing[1],fOldImgSpacing[2],fOldImgSpacing[3] );//Add by JDQ
	//imgReadRaw.GetImageExtent( fImgExtend ) ;
	imgReadRaws.GetImageSpacing(fNewImgSpacing[0],fNewImgSpacing[1],fNewImgSpacing[2],fNewImgSpacing[3] );//Add by JDQ
	float fresamle=zxh::minf(fNewImgSpacing[1],fNewImgSpacing[2])*0.5;
	int SearchRange[4]={2};
	for(int i=0;i<4;i++)
	{
		SearchRange[i]=int(fOldImgSpacing[i]/fNewImgSpacing[i]+0.5);
	}
	zxhImageDataT<short>imgReadNewRaw,imgReadResoNewRaws;
	imgReadNewRaw.NewImage( imgReadRaws.GetImageInfo() );
	int ImgNewSize[4]={1};
	imgReadNewRaw.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	////////
	vector<PointCordTypeDef> vPathPointsWorld;
	char* chext=GetPathEXT(chFileName);
	if(strcmp(chext,".vtk")==0)// read vtk file
	{
		cout<<"open successfully, format: vtk"<<endl;
		fstream DetectFile;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			ReadVtk(chFileName, vPathPointsWorld);

		}
		DetectFile.close();
	
	}
	else if(strcmp(chext,".txt")==0)// read vtk file
	{
		cout<<"open successfully, format: txt"<<endl;
	     fstream DetectFile;
		DetectFile.open(chFileName,ios::in);
		if(DetectFile)
		{
			if (RorS=="-N")
			{
				ReadTxtAsPCWRes(chFileName, vPathPointsWorld);
			}
			else
			{
			ReadTxtAsPC(chFileName,fresamle, vPathPointsWorld);
			}

		}
		DetectFile.close();
		
	}
	
	cout<<"Number of points:" <<vPathPointsWorld.size()-1<<endl;
	if (RorN=="-N")//map to a totally new image
	{
		cout<<"The line or points wil be mapped to the totally new image."<<endl;
		for(int it=0;it<ImgNewSize[3];++it)
			for(int iz=0;iz<ImgNewSize[2];++iz)
				for(int iy=0;iy<ImgNewSize[1];++iy)
					for(int ix=0;ix<ImgNewSize[0];++ix)
					{
						imgReadNewRaw.SetPixelByGreyscale(ix,iy,iz,it,0);
					}
					MapModelPointsToImage(vPathPointsWorld,imgReadNewRaw,SearchRange);
	}



	

	
	string chFileName2(chResultName);
	zxh::SaveImage(&imgReadNewRaw,chFileName2.c_str());

	

	cout << "Map successfully!" << endl;


}

 