

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

#include <string.h>
#include <iostream> 
#include <time.h> 
#include <math.h>
#include <fstream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <algorithm>
#include <sstream>
#include <limits>
#include<cmath>

#include "zxhImageGipl.h" 
#include "zxhImageData.h"
#include "zxhImageNifti.h"

#include "PatchExtractByWorldCoordinate.h"
#include <io.h>  
#include <direct.h> //创建文件夹
using namespace std;
#define M_PI 3.14159265358979323846

#include "jdq2017util.h"
#include "jdqdijkstra.h"
#include "jdqPoint.h"

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
	float pr;
}PointCordProbTypeDef;
string num2str(double i)
{
	stringstream ss;
	ss<<i;
	return ss.str();
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
	fdivec[0]=0.5*(fvec1[0]+fvec2[0]);
	fdivec[1]=0.5*(fvec1[1]+fvec2[1]);
	fdivec[2]=0.5*(fvec1[2]+fvec2[2]);
	zxh::VectorOP_Normalise(fdivec,3);
	return true;
}

bool Prj_to_plane_new(float fcurP[3],float fdivec[3],float fcenP[3],float fnewP[3])
{
	//ax+by+cz+D=0
	float fxyz[3]={fcurP[0],fcurP[1],fcurP[2]};
	float fD=-fdivec[0]*fcenP[0]-fdivec[1]*fcenP[1]-fdivec[2]*fcenP[2];
	float fabc[4]={fdivec[0],fdivec[1],fdivec[2],fD};
	float ft=(-fabc[0]*fxyz[0]-fabc[1]*fxyz[1]-fabc[2]*fxyz[2]-fabc[3])/(fabc[0]*fabc[0]+fabc[1]*fabc[1]+fabc[2]*fabc[2]);
	fnewP[0]=ft*fabc[0]+fxyz[0];
	fnewP[1]=ft*fabc[1]+fxyz[1];
	fnewP[2]=ft*fabc[2]+fxyz[2];
	return true;
}
bool Get_dirvec_byPrj(float fdivec[3],float fcenP[3],float fdivec1[3])
{

	float fcurP[3]={1,0,0};
	float fori[3]={0,0,0};
	float fnewP[3]={0,0,0};
	Prj_to_plane_new( fcurP, fdivec, fcenP, fnewP);
	float fnewori[3]={0,0,0};
	Prj_to_plane_new( fori, fdivec, fcenP, fnewori);
	fdivec1[0]=fnewP[0]-fnewori[0];
	fdivec1[1]=fnewP[1]-fnewori[1];
	fdivec1[2]=fnewP[2]-fnewori[2];
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
bool TransDatanum2Greynum(int np[3],int &npd,zxhImageDataT<short>&LabelImage)
{	int ImgSizeraw[]={0,0,0,0};
LabelImage.GetImageSize(ImgSizeraw[0],ImgSizeraw[1],ImgSizeraw[2],ImgSizeraw[3]);

npd=np[2] * ImgSizeraw[1] * ImgSizeraw[0] + np[1] * ImgSizeraw[0] + np[0];
return true;
}
bool TransGreynum2Datanum(int npd,int np[3],zxhImageDataT<short>&LabelImage)
{
	int ImgNewSize[]={0,0,0,0};
	LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	np[2]=npd/(ImgNewSize[1] * ImgNewSize[0]);
	int nxy=npd%(ImgNewSize[1] * ImgNewSize[0]);
	np[1]=nxy/ImgNewSize[0];
	np[0]=nxy%ImgNewSize[0];
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
bool ReadDireTxt(const char *chFileName, vector<PointCordTypeDef> &PointCord)
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
	PointCordTypeDef strctTempPoint;
	for (int i = 0; i < ref.size(); i++)
	{
		strctTempPoint.x =ref[i]._x;
		strctTempPoint.y =ref[i]._y;
		strctTempPoint.z =ref[i]._z;
		PointCord.push_back(strctTempPoint);
	}
	return true;
}


bool ReadCurve(string strCurvename,std::vector<jdq2017::point3D>&vcurveponts)
{
	cout<<"open successfully, curve file format: txt"<<endl;
	const char* chCurvename=strCurvename.c_str();
	if ( ! jdq2017::readCenterline(chCurvename, vcurveponts, 3))
	{
		std::cerr << "Error in reading line data" << std::endl;
		return 1;
	}
	return true;
}
bool ResampleCurve(std::vector<jdq2017::point3D>&vcurveponts,float fresamle)
{
	//resample the points
	jdq2017::ResamplePaths<jdq2017::point3D>Respler;
	Respler(vcurveponts,fresamle);
	vcurveponts=Respler.resultPath();
	return true;
}
float Calc_probOfPont(float PointWorldCor[3],vector<PointCordProbTypeDef>&vPointCentCord)
{
	//计算到最近金标中心点的距离
	float flambda=5;
	float fmindist=100000;
	int npotidx=-1;
	for (int i=0;i<vPointCentCord.size();i++)
	{
		PointCordProbTypeDef cupont=vPointCentCord[i];
		float fCenPontWord[3]={cupont.x,cupont.y,cupont.z};
		float fdist=zxh::VectorOP_Distance(PointWorldCor,fCenPontWord,3);
		if (fdist<fmindist)
		{
			fmindist=fdist;
			npotidx=i;

		}
	}
	float fmaxrad=vPointCentCord[npotidx].pr;
	float fexp=exp((-1*fmindist+fmaxrad)*flambda);
	float fproba=1-1/(1+fexp);
	if(fproba==0)
	{
		int x=0;
	}

	return fproba;
}
float Calc_probOfPont_MedIA18(float PointWorldCor[3],vector<PointCordProbTypeDef>&vPointCentCord)
{
	//计算到最近金标中心点的距离
	float flambda=6;
	float fmindist=100000;
	int npotidx=-1;
	for (int i=0;i<vPointCentCord.size();i++)
	{
		PointCordProbTypeDef cupont=vPointCentCord[i];
		float fCenPontWord[3]={cupont.x,cupont.y,cupont.z};
		float fdist=zxh::VectorOP_Distance(PointWorldCor,fCenPontWord,3);
		if (fdist<fmindist)
		{
			fmindist=fdist;
			npotidx=i;

		}
	}
	float fmaxrad=vPointCentCord[npotidx].pr;
	float fexp=exp((1-fmindist/fmaxrad)*flambda);
	float fproba=(fexp-1)/(exp(flambda)-1);
	if(fmindist>=fmaxrad) fproba=0;
	if(fproba==0)
	{
		int x=0;
	}

	return fproba;
}
//这个函数是用来生成沿着中心线切向取patch
bool LGeneratePatchExtractByWorldWithRandOffset(int PatchNumIdex,float fdivec[3][3],float InputOutputWorldCoord[], int PatchSize[],zxhImageData &IntensityImage,zxhImageData* PatchImageArray[3])
{
	zxhImageData *PatchImage1 = PatchImageArray[0];
	zxhImageData *PatchImage2 = PatchImageArray[1];
	zxhImageData *PatchImage3 = PatchImageArray[2];

	int N = PatchSize[0];
	int HalfPatchLength = PatchSize[1]; 



	int PatchPointWorldCoord[3] = { 0 };
	//Computing the intensity of points in the tangent plane using interpolation	
	//Interpolation
	zxhImageModelingLinear InterpolationMod;
	InterpolationMod.SetImage(&IntensityImage);
	for(int s=1;s<4;s++)
	{
		int sN=s*N;
		int sHalfPatchLength=s*HalfPatchLength;
		for (int i = -sN; i < sN + 1; i=i+s)
		{
			for (int j = -sN; j < sN + 1; j=j+s)
			{
				for (int k = -sHalfPatchLength; k < sHalfPatchLength + 1;k=k+s)//extend the patch along the gradient orientation
				{ 
					PatchPointWorldCoord[0] =InputOutputWorldCoord[0] + i*fdivec[0][0] + j*fdivec[1][0] + k*fdivec[2][0];
					PatchPointWorldCoord[1] = InputOutputWorldCoord[1] + i*fdivec[0][1] + j*fdivec[1][1] + k*fdivec[2][1];
					PatchPointWorldCoord[2] =InputOutputWorldCoord[2]+  i*fdivec[0][2] + j*fdivec[1][2] + k*fdivec[2][2];
					zxhImageData *PatchImage=PatchImageArray[s-1];
					//InterpolationMod.SetImage(pLALabel);//check the orientation of gradient				
					float IntensityValue = InterpolationMod.GetPixelFloatValueWithCheckByWorld(PatchPointWorldCoord[0], PatchPointWorldCoord[1], PatchPointWorldCoord[2], 0);

					int W2I_Coor_SetPatchPixel = (i/s + N) + (j/s + N) * (2 * N + 1) + (k/s + HalfPatchLength)* (2 * N + 1)* (2 * N + 1);

					PatchImage->SetPixelByGreyscale(W2I_Coor_SetPatchPixel, PatchNumIdex, 0, 0, IntensityValue);

				}
			}
		}
	}

	return true ;

}
//以下这个函数是在x, y, z三个方向取patch
bool LGeneratePatchExtractByWorldWithRandOffset_inxyzTrain(int PatchNumIdex,float fdivec[3][3],float InputOutputWorldCoord[], int PatchSize[],zxhImageData &IntensityImage,zxhImageData* PatchImageArray[3])
{


	zxhImageData *PatchImage1 = PatchImageArray[0];
	zxhImageData *PatchImage2 = PatchImageArray[1];
	zxhImageData *PatchImage3 = PatchImageArray[2];

	int N = PatchSize[0];
	int HalfPatchLength = PatchSize[1]; 



	int PatchPointWorldCoord[3] = { 0 };
	//Computing the intensity of points in the tangent plane using interpolation	
	//Interpolation
	zxhImageModelingLinear InterpolationMod;
	InterpolationMod.SetImage(&IntensityImage);
	for(int s=1;s<4;s++)
	{
		int sN=s*N;
		int sHalfPatchLength=s*HalfPatchLength;
		for (int i = -sN; i < sN + 1; i=i+s)
		{
			for (int j = -sN; j < sN + 1; j=j+s)
			{
				for (int k = -sHalfPatchLength; k < sHalfPatchLength + 1;k=k+s)//extend the patch along the gradient orientation
				{ 
					PatchPointWorldCoord[0] =InputOutputWorldCoord[0] + i*fdivec[0][0] + j*fdivec[1][0] + k*fdivec[2][0];
					PatchPointWorldCoord[1] = InputOutputWorldCoord[1] + i*fdivec[0][1] + j*fdivec[1][1] + k*fdivec[2][1];
					PatchPointWorldCoord[2] =InputOutputWorldCoord[2]+  i*fdivec[0][2] + j*fdivec[1][2] + k*fdivec[2][2];
					zxhImageData *PatchImage=PatchImageArray[s-1];
					//InterpolationMod.SetImage(pLALabel);//check the orientation of gradient				
					float IntensityValue = InterpolationMod.GetPixelFloatValueWithCheckByWorld(PatchPointWorldCoord[0], PatchPointWorldCoord[1], PatchPointWorldCoord[2], 0);

					int W2I_Coor_SetPatchPixel = (i/s + N) + (j/s + N) * (2 * N + 1) + (k/s + HalfPatchLength)* (2 * N + 1)* (2 * N + 1);

					PatchImage->SetPixelByGreyscale(W2I_Coor_SetPatchPixel, PatchNumIdex, 0, 0, IntensityValue);

				}
			}
		}
	}

	return true ;

}
//以下这个函数是在x, y, z三个方向取patch
bool LGeneratePatchExtractByWorldWithRandOffset_inxyzTest(int PatchNumIdex,float fdivec[3][3],float InputOutputWorldCoord[], int PatchSize[],zxhImageData &IntensityImage,zxhImageData* PatchImageArray[3])
{
	int N = PatchSize[0];
	int HalfPatchLength = PatchSize[1]; 



	int PatchPointWorldCoord[3] = { 0 };
	//Computing the intensity of points in the tangent plane using interpolation	
	//Interpolation
	zxhImageModelingLinear InterpolationMod;
	InterpolationMod.SetImage(&IntensityImage);
	for(int s=1;s<4;s++)
	{
		int sN=s*N;
		int sHalfPatchLength=s*HalfPatchLength;
		for (int i = -sN; i < sN + 1; i=i+s)
		{
			for (int j = -sN; j < sN + 1; j=j+s)
			{
				for (int k = -sHalfPatchLength; k < sHalfPatchLength + 1;k=k+s)//extend the patch along the gradient orientation
				{ 
					PatchPointWorldCoord[0] =InputOutputWorldCoord[0] + i*fdivec[0][0] + j*fdivec[1][0] + k*fdivec[2][0];
					PatchPointWorldCoord[1] = InputOutputWorldCoord[1] + i*fdivec[0][1] + j*fdivec[1][1] + k*fdivec[2][1];
					PatchPointWorldCoord[2] =InputOutputWorldCoord[2]+  i*fdivec[0][2] + j*fdivec[1][2] + k*fdivec[2][2];
					zxhImageData *PatchImage=PatchImageArray[s-1];
					//InterpolationMod.SetImage(pLALabel);//check the orientation of gradient				
					float IntensityValue = InterpolationMod.GetPixelFloatValueWithCheckByWorld(PatchPointWorldCoord[0], PatchPointWorldCoord[1], PatchPointWorldCoord[2], 0);

					int W2I_Coor_SetPatchPixel = (i/s + N) + (j/s + N) * (2 * N + 1) + (k/s + HalfPatchLength)* (2 * N + 1)* (2 * N + 1);
					PatchImage->SetPixelByGreyscale(W2I_Coor_SetPatchPixel, PatchNumIdex, 0, 0, IntensityValue);

				}
			}
		}
	}

	return true ;

}
int GeneratePontsset(zxhImageData &LabelImage,vector<PointCordTypeDef>&vPointCord)
{
	int ImgNewSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short intlabinte=LabelImage.GetPixelGreyscale(ix,iy,iz,it);

					if(ix==167&&iy==382&&iz==19)
					{
						int x=0;
					}
					if(intlabinte==1)
					{
						float PointWorldCor[3]={ix,iy,iz};
						LabelImage.GetImageInfo()->ImageToWorld(PointWorldCor);
						PointCordTypeDef Ponttemp;
						Ponttemp.x=PointWorldCor[0];
						Ponttemp.y=PointWorldCor[1];
						Ponttemp.z=PointWorldCor[2];
						vPointCord.push_back(Ponttemp);
					}
				}
				return vPointCord.size();
}

int GenerateVesCentPontsPrbset_Mask(zxhImageData &LabelImage,vector<PointCordProbTypeDef>&vPointCentCord,vector<PointCordProbTypeDef>&vPointCord,short *shMask)
{
	int ImgNewSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	int ImgSize[3]={ImgNewSize[0],ImgNewSize[1],ImgNewSize[2]};
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short intlabinte=LabelImage.GetPixelGreyscale(ix,iy,iz,it);
					float PointWorldCor[3]={ix,iy,iz};	
					LabelImage.GetImageInfo()->ImageToWorld(PointWorldCor);
					PointCordProbTypeDef Ponttemp;
					Ponttemp.x=PointWorldCor[0];
					Ponttemp.y=PointWorldCor[1];
					Ponttemp.z=PointWorldCor[2];
					if (intlabinte!=1) continue;
					if(intlabinte==1)
					{
						float fbprob=Calc_probOfPont(PointWorldCor,vPointCentCord);
						Ponttemp.pr=fbprob;
						vPointCord.push_back(Ponttemp);
						//标记已经取过；
						int npd[3]={Ponttemp.x,Ponttemp.y,Ponttemp.z};
						int np=-1;
						TransDatanum2Greynum(npd,np,LabelImage);
						shMask[np]=1;
					}

				}
				return vPointCord.size();
}
int GenerateVesCentPontsPrbset(zxhImageData &LabelImage,vector<PointCordProbTypeDef>&vPointCentCord,vector<PointCordProbTypeDef>&vPointCord)
{
	int ImgNewSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	int ImgSize[3]={ImgNewSize[0],ImgNewSize[1],ImgNewSize[2]};
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short intlabinte=LabelImage.GetPixelGreyscale(ix,iy,iz,it);
					float PointWorldCor[3]={ix,iy,iz};	
					LabelImage.GetImageInfo()->ImageToWorld(PointWorldCor);
					PointCordProbTypeDef Ponttemp;
					Ponttemp.x=PointWorldCor[0];
					Ponttemp.y=PointWorldCor[1];
					Ponttemp.z=PointWorldCor[2];
					if (intlabinte!=1) continue;
					if(intlabinte==1)
					{
						float fbprob=Calc_probOfPont_MedIA18(PointWorldCor,vPointCentCord);
						Ponttemp.pr=fbprob;
						vPointCord.push_back(Ponttemp);
					}

				}
				return vPointCord.size();
}
int SelectPointsFromlabimg(zxhImageData &LabelImage,vector<PointCordProbTypeDef>&vPointCord,int nlab)
{
	int ImgNewSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short intlabinte=LabelImage.GetPixelGreyscale(ix,iy,iz,it);
					float PointWorldCor[3]={ix,iy,iz};	
					LabelImage.GetImageInfo()->ImageToWorld(PointWorldCor);
					PointCordProbTypeDef Ponttemp;
					Ponttemp.x=PointWorldCor[0];
					Ponttemp.y=PointWorldCor[1];
					Ponttemp.z=PointWorldCor[2];
					if (intlabinte!=nlab) continue;
					if(intlabinte==nlab)
					{
						Ponttemp.pr=nlab;
						vPointCord.push_back(Ponttemp);
					}

				}
				return vPointCord.size();
}
int GenerateVesCentPontsRadset(zxhImageData &LabelImage,zxhImageData &LabelradImage,vector<PointCordProbTypeDef>&vPointCord)
{
	int ImgNewSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short intlabinte=LabelImage.GetPixelGreyscale(ix,iy,iz,it);
					float PointWorldCor[3]={ix,iy,iz};	
					LabelImage.GetImageInfo()->ImageToWorld(PointWorldCor);
					PointCordProbTypeDef Ponttemp;
					Ponttemp.x=PointWorldCor[0];
					Ponttemp.y=PointWorldCor[1];
					Ponttemp.z=PointWorldCor[2];
					if (intlabinte!=1) continue;
					if(intlabinte==1)
					{
						float fp=float(LabelradImage.GetPixelGreyscale(ix,iy,iz,it));
						float fpro=fp /1000;
						if(fpro==0)
						{
							int x=0;
						}
						Ponttemp.pr=fpro;

					}
					vPointCord.push_back(Ponttemp);
				}
				return vPointCord.size();
}

int GenerateVesBoudPontsPrbset(zxhImageData &LabelImage,vector<PointCordProbTypeDef>&vPointCentCord,vector<PointCordProbTypeDef>&vPointBoundCord)
{
	float fboud=5;
	int ImgNewSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short intlabinte=LabelImage.GetPixelGreyscale(ix,iy,iz,it);
					float PointWorldCor[3]={ix,iy,iz};	
					LabelImage.GetImageInfo()->ImageToWorld(PointWorldCor);
					PointCordProbTypeDef Ponttemp;
					Ponttemp.x=PointWorldCor[0];
					Ponttemp.y=PointWorldCor[1];
					Ponttemp.z=PointWorldCor[2];
					if (intlabinte==3) 
					{
						float fbprob=Calc_probOfPont(PointWorldCor,vPointCentCord);
						Ponttemp.pr=fbprob;
						vPointBoundCord.push_back(Ponttemp);
					}
				}
				return vPointBoundCord.size();
}
int GenerateVesNearBoudPontsPrbset(zxhImageData &LabelImage,vector<PointCordProbTypeDef>&vPointCentCord,vector<PointCordProbTypeDef>&vPointBoundCord)
{
	float fboud=5;
	int ImgNewSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short intlabinte=LabelImage.GetPixelGreyscale(ix,iy,iz,it);
					float PointWorldCor[3]={ix,iy,iz};	
					LabelImage.GetImageInfo()->ImageToWorld(PointWorldCor);
					PointCordProbTypeDef Ponttemp;
					Ponttemp.x=PointWorldCor[0];
					Ponttemp.y=PointWorldCor[1];
					Ponttemp.z=PointWorldCor[2];
					if (intlabinte==2)
					{
						float fbprob=Calc_probOfPont(PointWorldCor,vPointCentCord);
						Ponttemp.pr=fbprob;
						vPointBoundCord.push_back(Ponttemp);
					}
				}
				return vPointBoundCord.size();
}
int GenerateVesPontsPrbset(zxhImageData &LabelImage,vector<PointCordProbTypeDef>&vPointCord)
{
	int ImgNewSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short intlabinte=LabelImage.GetPixelGreyscale(ix,iy,iz,it);
					float PointWorldCor[3]={ix,iy,iz};	
					LabelImage.GetImageInfo()->ImageToWorld(PointWorldCor);
					PointCordProbTypeDef Ponttemp;
					Ponttemp.x=PointWorldCor[0];
					Ponttemp.y=PointWorldCor[1];
					Ponttemp.z=PointWorldCor[2];
					if (intlabinte==0||intlabinte==3) continue;
					if(intlabinte==1)
					{
						Ponttemp.pr=1;

					}
					if(intlabinte==2)
					{

						Ponttemp.pr=0.5;

					}
					if(intlabinte==3)
					{
						Ponttemp.pr=0.2;

					}
					vPointCord.push_back(Ponttemp);
					/*				if(intlabinte==1)
					{
					Ponttemp.pr=1;
					vPointCord1.push_back(Ponttemp);
					}
					if(intlabinte==2)
					{
					Ponttemp.pr=0.5;
					vPointCord2.push_back(Ponttemp);
					}
					if(intlabinte==3)
					{
					Ponttemp.pr=0.2;
					vPointCord3.push_back(Ponttemp);
					}*/
				}
				return vPointCord.size();
}
int GenerateRandomNonVesPontsPrbset(zxhImageData &LabelImage,zxhImageData &ROIImage,vector<PointCordProbTypeDef>&vVesCentPointCord,vector<PointCordProbTypeDef>&vNonVesPointCord)
{
	int ImgNewSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short intlabinte=LabelImage.GetPixelGreyscale(ix,iy,iz,it);
					short intROIinte=ROIImage.GetPixelGreyscale(ix,iy,iz,it);
					float PointWorldCor[3]={ix,iy,iz};	
					LabelImage.GetImageInfo()->ImageToWorld(PointWorldCor);
					PointCordProbTypeDef Ponttemp;
					Ponttemp.x=PointWorldCor[0];
					Ponttemp.y=PointWorldCor[1];
					Ponttemp.z=PointWorldCor[2];
					Ponttemp.pr=0;
					if (intROIinte==0) continue;//只选ROI里的点
					if (intlabinte!=0) continue;//只选label 0的点
					float fbprob=Calc_probOfPont(PointWorldCor,vVesCentPointCord);
					Ponttemp.pr=fbprob;
					vNonVesPointCord.push_back(Ponttemp);
				}
				return vNonVesPointCord.size();
}
int GetNonVesPontsset(zxhImageData &LabelImage,zxhImageData &ROIImage,vector<PointCordProbTypeDef>&vNonVesPointCord)
{
	int ImgNewSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					short intlabinte=LabelImage.GetPixelGreyscale(ix,iy,iz,it);
					short intROIinte=ROIImage.GetPixelGreyscale(ix,iy,iz,it);
					float PointWorldCor[3]={ix,iy,iz};	
					LabelImage.GetImageInfo()->ImageToWorld(PointWorldCor);
					PointCordProbTypeDef Ponttemp;
					Ponttemp.x=PointWorldCor[0];
					Ponttemp.y=PointWorldCor[1];
					Ponttemp.z=PointWorldCor[2];
					Ponttemp.pr=0;
					if (intROIinte==0) continue;//只选ROI里的点
					if (intlabinte!=0) continue;//只选label 0的点
					vNonVesPointCord.push_back(Ponttemp);
				}
				return vNonVesPointCord.size();
}

double generateGaussianNoise(double mu, double sigma)
{
	const double epsilon = 1.17549e-038;
	const double two_pi = 2.0*3.14159265358979323846;

	static double z0, z1;
	static bool generate;
	generate = !generate;

	if (!generate)
		return z1 * sigma + mu;

	double u1, u2;
	do
	{
		u1 = rand() * (1.0 / RAND_MAX);
		u2 = rand() * (1.0 / RAND_MAX);
	}
	while ( u1 <= epsilon );

	z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
	z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
	return z0 * sigma + mu;
}
bool Generate_gaussianrad(int nptofonepoint,float frad,vector<float>&vrandrad)
{	

	float maxrad=1.5*frad;
	vrandrad.clear();
	float epsilon =0.01;
	srand((unsigned)time(NULL));
	float fprob=100000;
	if (maxrad==0)
	{
		int x=0;
	}
	while (1)
	{
		float randrad =((rand() % (400)) / 100.0)*frad;	
		fprob=zxh::Gaussian(float (randrad),float(0),float(maxrad));
		if (fprob>= epsilon)
		{
			vrandrad.push_back(0.5+frad+randrad);//超过一个spacing
			if (vrandrad.size()>=nptofonepoint)
			{
				break;
			}
		}
	}

	return true;
}
bool Generate_gaussianrad_web(int nptofonepoint,float frad,vector<float>&vrandrad)
{	
	int tmpnp=nptofonepoint;
	vector<float>vtmprandrad;
	float maxrad=30*frad;
	vtmprandrad.clear();
	float epsilon =0.01;
	srand((unsigned)time(NULL));
	float fprob=100000;
	if (maxrad==0)
	{
		int x=0;
	}
	if (nptofonepoint<10)
	{
		tmpnp=nptofonepoint*10;
	}
	while (1)
	{
		float fgrad =generateGaussianNoise(0,maxrad/3);
		vtmprandrad.push_back(0.5+abs(fgrad));//超过一个spacing
		if (vtmprandrad.size()>=tmpnp)
		{
			break;
		}

	}
	random_shuffle(vtmprandrad.begin(), vtmprandrad.end());
	vrandrad.clear();
	vrandrad.insert(vrandrad.end(),vtmprandrad.begin(),vtmprandrad.begin()+nptofonepoint);


	return true;
}
bool GetOffCentVesPontsset_GaussDistri(int nnonvptsize,vector<PointCordProbTypeDef>&vVesCentPointCord,zxhImageData &LabelImage,zxhImageData &ROIImage,vector<PointCordProbTypeDef>&vNonVesPointCord)
{
	//定义x,y,z坐标轴的方向

	float fdivec[3][3]={0};
	//x方向
	fdivec[0][0]=1;
	fdivec[0][1]=0;
	fdivec[0][2]=0;

	//y方向
	fdivec[1][0]=0;
	fdivec[1][1]=1;
	fdivec[1][2]=0;

	//z方向
	fdivec[2][0]=0;
	fdivec[2][1]=0;
	fdivec[2][2]=1;
	//将所要挑选的offcentline的点平均分布到每个中心点上
	int nptofonepoint=nnonvptsize/vVesCentPointCord.size();
	for(int i=0;i<vVesCentPointCord.size();i++)
	{
		//获取当前点
		PointCordProbTypeDef Pcurpon;
		Pcurpon=vVesCentPointCord[i];
		//
		//--------------------随机产生一个正态分布的半径长度
		float frad=Pcurpon.pr;
		vector<float>vrandrad;
		if (frad==0)
		{
			int x=0;
		}
		Generate_gaussianrad_web(nptofonepoint,frad,vrandrad);
		///------------------将每一个随机产生的长度给一个随机向量，
		for(int nrad=0;nrad<vrandrad.size();nrad++)
		{
			float randradius=vrandrad[nrad];
			srand((unsigned)time(NULL));  
			while (1)
			{

				//--------------------随机产生一个向量--------------------
				float randoffset1 =float((rand() % (200))-100) / 100.0;	
				float randoffset2 =float((rand() % (200))-100) / 100.0;	
				float randoffset3=float((rand() % (200))-100) / 100.0;	

				float fvec1[3]={randoffset1*fdivec[0][0],randoffset1*fdivec[0][1],randoffset1*fdivec[0][2]};
				float fvec2[3]={randoffset2*fdivec[1][0],randoffset2*fdivec[1][1],randoffset2*fdivec[1][2]};
				float fvec3[3]={randoffset3*fdivec[2][0],randoffset3*fdivec[2][1],randoffset3*fdivec[2][2]};
				float fsumvec[3]={fvec1[0]+fvec2[0]+fvec3[0],fvec1[1]+fvec2[1]+fvec3[1],fvec1[2]+fvec2[2]+fvec3[2]};
				//Rand shift
				zxh::VectorOP_Normalise(fsumvec,3);

				float InputWorldCoord_x = Pcurpon.x + randradius*fsumvec[0];
				float InputWorldCoord_y = Pcurpon.y + randradius*fsumvec[1];
				float InputWorldCoord_z = Pcurpon.z + randradius*fsumvec[2];
				//
				float InputWorldCoord[3]={InputWorldCoord_x,InputWorldCoord_y,InputWorldCoord_z};
				float InputNeiNodeWorldCoord[4]={InputWorldCoord_x,InputWorldCoord_y,InputWorldCoord_z,0 };
				LabelImage.GetImageInfo()->WorldToImage(InputNeiNodeWorldCoord);
				int nscx = zxh::round(InputNeiNodeWorldCoord[0]);
				int nscy = zxh::round(InputNeiNodeWorldCoord[1]);
				int nscz = zxh::round(InputNeiNodeWorldCoord[2]);
				bool bIsInsideImage = LabelImage.InsideImage(nscx, nscy, nscz, 0); // 超过图像边界的，不给予考虑，也就是说，默认为normal myo
				if (!bIsInsideImage)
				{
					std::cout << "warning: niebour node of point"<< i<< "is not inside image " << "\n"; //-----------------
					continue;
				}
				short shlabnum=LabelImage.GetPixelGreyscale(nscx, nscy, nscz, 0);
				short shroinum=ROIImage.GetPixelGreyscale(nscx, nscy, nscz, 0);
				if (shroinum==0) continue;//只选ROI里的点
				if (shlabnum!=0) continue;//不选label 1,label2里的点
				PointCordProbTypeDef Ponttemp;
				Ponttemp.x=InputWorldCoord_x;
				Ponttemp.y=InputWorldCoord_y;
				Ponttemp.z=InputWorldCoord_z;
				float fbprob=Calc_probOfPont(InputWorldCoord,vVesCentPointCord);
				Ponttemp.pr=fbprob;
				if(fbprob>1)
				{
					int x=0;
				}
				vNonVesPointCord.push_back(Ponttemp);
				break;

			}
		}

	}
	return true;
}

bool GetOffCentVesPontsset_GaussDistri_Mask(int nnonvptsize,vector<PointCordProbTypeDef>&vVesCentPointCord,zxhImageData &LabelImage,zxhImageData &ROIImage,vector<PointCordProbTypeDef>&vNonVesPointCord,vector<PointCordTypeDef>vDisDir,short *shMask)
{
	//获取图像大小
	int ImgSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	//定义x,y,z坐标轴的方向

	float fdivec[3][3]={0};
	//x方向
	fdivec[0][0]=1;
	fdivec[0][1]=0;
	fdivec[0][2]=0;

	//y方向
	fdivec[1][0]=0;
	fdivec[1][1]=1;
	fdivec[1][2]=0;

	//z方向
	fdivec[2][0]=0;
	fdivec[2][1]=0;
	fdivec[2][2]=1;

	//将所要挑选的offcentline的点平均分布到每个中心点上(已经挑选的中心点的数量)
	int nptofonepoint=nnonvptsize/vVesCentPointCord.size();
	for(int i=0;i<vVesCentPointCord.size();i++)
	{
		//获取当前点
		PointCordProbTypeDef Pcurpon;
		Pcurpon=vVesCentPointCord[i];
		//
		//--------------------随机产生一个正态分布的半径长度
		float frad=Pcurpon.pr;
		vector<float>vrandrad;
		if (frad==0)
		{
			int x=0;
		}
		Generate_gaussianrad_web(nptofonepoint,frad,vrandrad);
		///------------------将每一个随机产生的长度给一个随机向量，
		for(int nrad=0;nrad<vrandrad.size();nrad++)
		{
			float randradius=vrandrad[nrad];
			while(1)
			{
				//打乱球面向量
				random_shuffle(vDisDir.begin(), vDisDir.end());
				int numofPoints=0;
				for(int nd=0;nd<vDisDir.size();nd++)
				{
					//--------------------随机产生一个向量(单位球面上)--------------------

					float fsumvec[3]={vDisDir[nd].x,vDisDir[nd].y,vDisDir[nd].z};
					//Rand shift
					zxh::VectorOP_Normalise(fsumvec,3);

					float InputWorldCoord_x = Pcurpon.x + randradius*fsumvec[0];
					float InputWorldCoord_y = Pcurpon.y + randradius*fsumvec[1];
					float InputWorldCoord_z = Pcurpon.z + randradius*fsumvec[2];
					//
					float InputWorldCoord[3]={InputWorldCoord_x,InputWorldCoord_y,InputWorldCoord_z};
					float InputNeiNodeWorldCoord[4]={InputWorldCoord_x,InputWorldCoord_y,InputWorldCoord_z,0 };
					LabelImage.GetImageInfo()->WorldToImage(InputNeiNodeWorldCoord);
					int nscx = zxh::round(InputNeiNodeWorldCoord[0]);
					int nscy = zxh::round(InputNeiNodeWorldCoord[1]);
					int nscz = zxh::round(InputNeiNodeWorldCoord[2]);
					bool bIsInsideImage = LabelImage.InsideImage(nscx, nscy, nscz, 0); // 超过图像边界的，不给予考虑，也就是说，默认为normal myo
					if (!bIsInsideImage)
					{
						std::cout << "warning: niebour node of point"<< i<< "is not inside image " << "\n"; //-----------------
						continue;
					}
					//判断是否已经取过
					int np[3]={nscx,nscy,nscz};
					int npt=-1;
					TransDatanum2Greynum(np,npt,LabelImage);
					short shortmas=shMask[npt];
					short shortROI=ROIImage.GetPixelGreyscale(nscx,nscy,nscz,0);
					if (shortmas==1) continue;//这个点取过了就不取
					if (shortROI==0) continue;//只选ROI里的点

					PointCordProbTypeDef Ponttemp;
					Ponttemp.x=InputWorldCoord_x;
					Ponttemp.y=InputWorldCoord_y;
					Ponttemp.z=InputWorldCoord_z;
					float fbprob=Calc_probOfPont(InputWorldCoord,vVesCentPointCord);
					Ponttemp.pr=fbprob;
					if(fbprob>1)
					{
						int x=0;
					}
					vNonVesPointCord.push_back(Ponttemp);
					numofPoints++;
					break;

				}//for
				if(numofPoints==0)//整个球内的点都没有取到
				{
					cout<<"the radius is too small"<<endl;
					randradius=randradius+0.5;
				}
				else
				{
					break;//while
				}
			}//while
		}//for

	}
	return true;
}
bool GetOffCentVesPontsset_GaussDistri_Mask_Rand(int nvcptsize,float flabbal,vector<PointCordProbTypeDef>&vVesCentPointCordrad,zxhImageData &LabelImage,zxhImageData &ROIImage,vector<PointCordProbTypeDef>&vNonVesPointCord,vector<PointCordTypeDef>vDisDir,short *shMask)
{
	//获取图像大小
	int ImgSize[4]={0,0,0,0};
	LabelImage.GetImageSize(ImgSize[0],ImgSize[1],ImgSize[2],ImgSize[3]);
	//定义x,y,z坐标轴的方向

	float fdivec[3][3]={0};
	//x方向
	fdivec[0][0]=1;
	fdivec[0][1]=0;
	fdivec[0][2]=0;

	//y方向
	fdivec[1][0]=0;
	fdivec[1][1]=1;
	fdivec[1][2]=0;

	//z方向
	fdivec[2][0]=0;
	fdivec[2][1]=0;
	fdivec[2][2]=1;
	int fvoffsetbase=zxh::minf(nvcptsize*flabbal,nvcptsize);
	int nlabbal=nvcptsize*flabbal/(zxh::minf(nvcptsize*flabbal,nvcptsize));

	//打乱中心点的顺序
	vector<PointCordProbTypeDef>vVesCentPointCord;
	random_shuffle(vVesCentPointCordrad.begin(), vVesCentPointCordrad.end());//打乱所有的正样本
	vVesCentPointCord.insert(vVesCentPointCord.end(),vVesCentPointCordrad.begin(),vVesCentPointCordrad.begin()+fvoffsetbase);
	//如果负样本数量小于中心点的总数

	//如果负样本的数量大于中心点的总数
	for(int i=0;i<vVesCentPointCord.size();i++)
	{
		//获取当前点
		PointCordProbTypeDef Pcurpon;
		Pcurpon=vVesCentPointCord[i];
		//
		//--------------------随机产生一个正态分布的半径长度
		float frad=Pcurpon.pr;
		vector<float>vrandrad;
		if (frad==0)
		{
			int x=0;
		}
		Generate_gaussianrad_web(nlabbal,frad,vrandrad);
		///------------------将每一个随机产生的长度给一个随机向量，
		for(int nrad=0;nrad<vrandrad.size();nrad++)
		{
			float randradius=vrandrad[nrad];
			while(1)
			{
				//打乱球面向量
				random_shuffle(vDisDir.begin(), vDisDir.end());
				int numofPoints=0;
				for(int nd=0;nd<vDisDir.size();nd++)
				{
					//--------------------随机产生一个向量(单位球面上)--------------------

					float fsumvec[3]={vDisDir[nd].x,vDisDir[nd].y,vDisDir[nd].z};
					//Rand shift
					zxh::VectorOP_Normalise(fsumvec,3);

					float InputWorldCoord_x = Pcurpon.x + randradius*fsumvec[0];
					float InputWorldCoord_y = Pcurpon.y + randradius*fsumvec[1];
					float InputWorldCoord_z = Pcurpon.z + randradius*fsumvec[2];
					//
					float InputWorldCoord[3]={InputWorldCoord_x,InputWorldCoord_y,InputWorldCoord_z};
					float InputNeiNodeWorldCoord[4]={InputWorldCoord_x,InputWorldCoord_y,InputWorldCoord_z,0 };
					LabelImage.GetImageInfo()->WorldToImage(InputNeiNodeWorldCoord);
					int nscx = zxh::round(InputNeiNodeWorldCoord[0]);
					int nscy = zxh::round(InputNeiNodeWorldCoord[1]);
					int nscz = zxh::round(InputNeiNodeWorldCoord[2]);
					bool bIsInsideImage = LabelImage.InsideImage(nscx, nscy, nscz, 0); // 超过图像边界的，不给予考虑，也就是说，默认为normal myo
					if (!bIsInsideImage)
					{
						std::cout << "warning: niebour node of point "<< i<< " is not inside the image " << "\n"; //-----------------
						continue;
					}
					//判断是否已经取过
					int np[3]={nscx,nscy,nscz};
					int npt=-1;
					TransDatanum2Greynum(np,npt,LabelImage);
					short shortmas=shMask[npt];
					short shortROI=ROIImage.GetPixelGreyscale(nscx,nscy,nscz,0);
					if (shortmas==1) continue;//这个点取过了就不取
					if (shortROI==0) continue;//只选ROI里的点

					PointCordProbTypeDef Ponttemp;
					Ponttemp.x=InputWorldCoord_x;
					Ponttemp.y=InputWorldCoord_y;
					Ponttemp.z=InputWorldCoord_z;
					float fbprob=Calc_probOfPont_MedIA18(InputWorldCoord,vVesCentPointCord);
					Ponttemp.pr=fbprob;
					if(fbprob>1)
					{
						int x=0;
					}
					vNonVesPointCord.push_back(Ponttemp);
					numofPoints++;
					break;

				}//for
				if(numofPoints==0)//整个球内的点都没有取到
				{
					cout<<"the radius is too small"<<endl;
					randradius=randradius+0.5;
				}
				else
				{
					break;//while
				}
			}//while
		}//for

	}
	return true;
}
int GenerateSelectVesCentPontsPrbset_Mask(int nvcpselec,vector<PointCordProbTypeDef>vALLVesCentPointCord,vector<PointCordProbTypeDef>&vSeletcVesCentPointCord,zxhImageData &LabelImage,short *shMask)
{
	random_shuffle(vALLVesCentPointCord.begin(), vALLVesCentPointCord.end());//打乱所有的正样本
	vSeletcVesCentPointCord.insert(vSeletcVesCentPointCord.end(),vALLVesCentPointCord.begin(),vALLVesCentPointCord.begin()+nvcpselec);
	//标记mask
	for(int i=0;i<vSeletcVesCentPointCord.size();i++)
	{
		float PointWorldCor[3]={vSeletcVesCentPointCord[i].x,vSeletcVesCentPointCord[i].y,vSeletcVesCentPointCord[i].z};	
		int np[3]={0,0,0};
		int npd=0;
		World2Image(PointWorldCor, np,LabelImage);
		TransDatanum2Greynum(np,npd,LabelImage);
		shMask[npd]=1;
	}
	return vSeletcVesCentPointCord.size();
}

bool GetNonVesPontsset_GaussDistri(int nnonvptsize,vector<PointCordProbTypeDef>&vVesCentPointCord,zxhImageData &LabelImage,zxhImageData &ROIImage,vector<PointCordProbTypeDef>&vNonVesPointCord)
{
	//

	float fdivec[3][3]={0};

	//x方向
	fdivec[0][0]=1;
	fdivec[0][1]=0;
	fdivec[0][2]=0;

	//y方向
	fdivec[1][0]=0;
	fdivec[1][1]=1;
	fdivec[1][2]=0;

	//z方向
	fdivec[2][0]=0;
	fdivec[2][1]=0;
	fdivec[2][2]=1;

	int nptofonepoint=nnonvptsize/vVesCentPointCord.size();
	for(int i=0;i<vVesCentPointCord.size();i++)
	{
		//获取当前点
		PointCordProbTypeDef Pcurpon;
		Pcurpon=vVesCentPointCord[i];
		//
		//--------------------随机产生一个正态分布的半径长度
		float frad=Pcurpon.pr;
		vector<float>vrandrad;
		if (frad==0)
		{
			int x=0;
		}
		Generate_gaussianrad_web(nptofonepoint,frad,vrandrad);
		///------------------将每一个随机产生的长度给一个随机向量，
		for(int nrad=0;nrad<vrandrad.size();nrad++)
		{
			float randradius=vrandrad[nrad];
			srand((unsigned)time(NULL));  
			while (1)
			{

				//--------------------随机产生一个向量--------------------
				float randoffset1 =float((rand() % (200))-100) / 100.0;	
				float randoffset2 =float((rand() % (200))-100) / 100.0;	
				float randoffset3=float((rand() % (200))-100) / 100.0;	

				float fvec1[3]={randoffset1*fdivec[0][0],randoffset1*fdivec[0][1],randoffset1*fdivec[0][2]};
				float fvec2[3]={randoffset2*fdivec[1][0],randoffset2*fdivec[1][1],randoffset2*fdivec[1][2]};
				float fvec3[3]={randoffset3*fdivec[2][0],randoffset3*fdivec[2][1],randoffset3*fdivec[2][2]};
				float fsumvec[3]={fvec1[0]+fvec2[0]+fvec3[0],fvec1[1]+fvec2[1]+fvec3[1],fvec1[2]+fvec2[2]+fvec3[2]};
				//Rand shift
				zxh::VectorOP_Normalise(fsumvec,3);

				float InputWorldCoord_x = Pcurpon.x + randradius*fsumvec[0];
				float InputWorldCoord_y = Pcurpon.y + randradius*fsumvec[1];
				float InputWorldCoord_z = Pcurpon.z + randradius*fsumvec[2];
				//
				float InputWorldCoord[3]={InputWorldCoord_x,InputWorldCoord_y,InputWorldCoord_z};
				float InputNeiNodeWorldCoord[4]={InputWorldCoord_x,InputWorldCoord_y,InputWorldCoord_z,0 };
				LabelImage.GetImageInfo()->WorldToImage(InputNeiNodeWorldCoord);
				int nscx = zxh::round(InputNeiNodeWorldCoord[0]);
				int nscy = zxh::round(InputNeiNodeWorldCoord[1]);
				int nscz = zxh::round(InputNeiNodeWorldCoord[2]);
				bool bIsInsideImage = LabelImage.InsideImage(nscx, nscy, nscz, 0); // 超过图像边界的，不给予考虑，也就是说，默认为normal myo
				if (!bIsInsideImage)
				{
					std::cout << "warning: niebour node of point"<< i<< "is not inside image " << "\n"; //-----------------
					continue;
				}
				short shlabnum=LabelImage.GetPixelGreyscale(nscx, nscy, nscz, 0);
				short shroinum=ROIImage.GetPixelGreyscale(nscx, nscy, nscz, 0);
				if (shroinum==0) continue;//只选ROI里的点
				if (shlabnum!=0) continue;//不选label 1,label2里的点
				PointCordProbTypeDef Ponttemp;
				Ponttemp.x=InputWorldCoord_x;
				Ponttemp.y=InputWorldCoord_y;
				Ponttemp.z=InputWorldCoord_z;
				float fbprob=Calc_probOfPont(InputWorldCoord,vVesCentPointCord);
				Ponttemp.pr=fbprob;
				if(fbprob>1)
				{
					int x=0;
				}
				vNonVesPointCord.push_back(Ponttemp);
				break;

			}
		}

	}
	//int ImgNewSize[4]={0,0,0,0};
	//LabelImage.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	//for(int it=0;it<ImgNewSize[3];++it)
	//	for(int iz=0;iz<ImgNewSize[2];++iz)
	//		for(int iy=0;iy<ImgNewSize[1];++iy)
	//			for(int ix=0;ix<ImgNewSize[0];++ix)
	//			{
	//				short intlabinte=LabelImage.GetPixelGreyscale(ix,iy,iz,it);
	//				short intROIinte=ROIImage.GetPixelGreyscale(ix,iy,iz,it);
	//				float PointWorldCor[3]={ix,iy,iz};	
	//				LabelImage.GetImageInfo()->ImageToWorld(PointWorldCor);
	//				PointCordProbTypeDef Ponttemp;
	//				Ponttemp.x=PointWorldCor[0];
	//				Ponttemp.y=PointWorldCor[1];
	//				Ponttemp.z=PointWorldCor[2];
	//				Ponttemp.pr=0;
	//				if (intROIinte==0) continue;//只选ROI里的点
	//				if (intlabinte!=0) continue;//只选label 0的点
	//				vNonVesPointCord.push_back(Ponttemp);
	//			}
	return true;
}
bool GenOnAndOffPonts_Rand(zxhImageDataT<short> &LabelImage,zxhImageDataT<short> &LabelradImage,zxhImageDataT<short>&ROIImage,vector<PointCordTypeDef> vDisDir,vector<PointCordProbTypeDef>&vSeletcVesCentPointCord,vector<PointCordProbTypeDef>&vSelecOffSetPointCord,float flabbal,float fvraPositive)
{
	//这个脚本是随机选取中心线上的点用来提取patch
	//-----------------读取label image,并产生正负样本的patch的中心点---------------------------
	//首先产生一个mask，用来记录已经取过的点（没取过的点为0， 取过1）
	int n=0;
	int ImgSizeraw[]={0,0,0,0};
	LabelImage.GetImageSize(ImgSizeraw[0],ImgSizeraw[1],ImgSizeraw[2],ImgSizeraw[3]);
	short *shMask=new short[ImgSizeraw[0]*ImgSizeraw[1]*ImgSizeraw[2]-1];
	for (int i = 0; i < ImgSizeraw[0]*ImgSizeraw[1]*ImgSizeraw[2]-1; i ++)shMask[i]=0;
	//读取中心点的坐标
	vector<PointCordProbTypeDef> vVesCentPointCordrad,vALLVesCentPointCord;//金标准点，概率为1,online
	int nvcp=GenerateVesCentPontsRadset(LabelImage,LabelradImage,vVesCentPointCordrad);
	int nvcptsize=GenerateVesCentPontsPrbset(LabelImage,vVesCentPointCordrad,vALLVesCentPointCord);
	int nvcpselec=fvraPositive*nvcptsize;//需要挑选的正样本的数量
	int nvselecvp=GenerateSelectVesCentPontsPrbset_Mask(nvcpselec,vALLVesCentPointCord,vSeletcVesCentPointCord,LabelImage,shMask);


	//挑选offcent的点，以中心线为基准，数量呈正太分布，并在ROI内部；每点的概率符合一个设计好的函数分布

	GetOffCentVesPontsset_GaussDistri_Mask_Rand(nvcpselec/fvraPositive,flabbal,vVesCentPointCordrad,LabelImage,ROIImage,vSelecOffSetPointCord,vDisDir,shMask);

	return true;
}
bool GenOnAndOffPonts(zxhImageDataT<short> LabelImage,zxhImageDataT<short> LabelradImage,zxhImageDataT<short>ROIImage,vector<PointCordTypeDef> vDisDir,vector<PointCordProbTypeDef>&vVesCentPointCord,vector<PointCordProbTypeDef>&vSelecOffSetPointCord,float labbal)
{
	//这个脚本是将所有的中心线上的点都用来提取patch
	//-----------------读取label image,并产生正负样本的patch的中心点---------------------------
	//首先产生一个mask，用来记录已经取过的点（没取过的点为0， 取过1）
	int n=0;
	int ImgSizeraw[]={0,0,0,0};
	LabelImage.GetImageSize(ImgSizeraw[0],ImgSizeraw[1],ImgSizeraw[2],ImgSizeraw[3]);
	short *shMask=new short[ImgSizeraw[0]*ImgSizeraw[1]*ImgSizeraw[2]-1];
	for (int i = 0; i < ImgSizeraw[0]*ImgSizeraw[1]*ImgSizeraw[2]-1; i ++)shMask[i]=0;
	//读取中心点的坐标
	vector<PointCordProbTypeDef> vVesCentPointCordrad;//金标准点，概率为1,online
	int nvcp=GenerateVesCentPontsRadset(LabelImage,LabelradImage,vVesCentPointCordrad);
	int nvcptsize=GenerateVesCentPontsPrbset_Mask(LabelImage,vVesCentPointCordrad,vVesCentPointCord,shMask);



	//挑选offline
	int noffvptsize=nvcptsize*labbal;
	//挑选offcent的点，以中心线为基准，数量呈正太分布，并在ROI内部；每点的概率符合一个设计好的函数分布

	GetOffCentVesPontsset_GaussDistri_Mask(noffvptsize,vVesCentPointCordrad,LabelImage,ROIImage,vSelecOffSetPointCord,vDisDir,shMask);

	return true;
}
PointCordTypeDef RotateByVector(double old_x, double old_y, double old_z, double vx, double vy, double vz, double theta)
{
    double r = theta * M_PI / 180;
    double c = cos(r);
    double s = sin(r);
    double new_x = (vx*vx*(1 - c) + c) * old_x + (vx*vy*(1 - c) - vz*s) * old_y + (vx*vz*(1 - c) + vy*s) * old_z;
    double new_y = (vy*vx*(1 - c) + vz*s) * old_x + (vy*vy*(1 - c) + c) * old_y + (vy*vz*(1 - c) - vx*s) * old_z;
    double new_z = (vx*vz*(1 - c) - vy*s) * old_x + (vy*vz*(1 - c) + vx*s) * old_y + (vz*vz*(1 - c) + c) * old_z;
    PointCordTypeDef Point3f;
	Point3f.x=new_x;
	Point3f.y=new_y;
	Point3f.z=new_z;
	return Point3f;
}
void Gen_fdivecByRotxyz(vector<PointCordTypeDef> &vxyzori,vector<vector<PointCordTypeDef>>&vvrandDir,int nRottime)
{
	srand((unsigned)time(NULL));
	for (int nt=0;nt<nRottime;nt++)
	{
		
		for(int nxyz=0;nxyz<vxyzori.size();nxyz++)//将一个坐标系认为是三个点的组合
		{
			vector<PointCordTypeDef> vrandDire;
			PointCordTypeDef tmpxyz;
			vector<PointCordTypeDef> otherxyz;
			for(int i=0;i<3;i++)//将这三个点首先分成两组，一组存放旋转方向，另一组存放要旋转的点
			{
				PointCordTypeDef tmpixyz=vxyzori[i];
				if(i==nxyz)
					tmpxyz=tmpixyz;//旋转前的坐标系中，依次选取x,y,z作为变得坐标系

				else
				{
					otherxyz.push_back(tmpixyz);
				}


			}	
			double old_x;
			double old_y;
			double old_z;
			double vx;
			double vy;
			double vz;
			
			vrandDire.push_back(tmpxyz);
			vx=tmpxyz.x;
			vy=tmpxyz.y;
			vz=tmpxyz.z;
			for (int j=0;j<otherxyz.size();j++)//接下来，对otherxyz，进行旋转
			{
				PointCordTypeDef tmpjxyz;
				tmpjxyz=otherxyz[j];
				old_x=tmpjxyz.x;
				old_y=tmpjxyz.y;
				old_z=tmpjxyz.z;
				//随机生成一个旋转角度
				int ntehta =-180+(int)(360.0 * rand()/(RAND_MAX+1.0));
				PointCordTypeDef newtmjxyz;
				newtmjxyz=RotateByVector(old_x,  old_y,  old_z,  vx,  vy,  vz, ntehta);

				vrandDire.push_back(newtmjxyz);
			}


vvrandDir.push_back(vrandDire);


		}
		
	}
}
bool Gen_OnCentandOffCentPoints(int numAxis,vector<PointCordProbTypeDef> &vVesCentPointCord,vector<PointCordProbTypeDef>&vSelecOffSetPointCord,int SiglePatchSize,float spacing111[4],int PatchInfo[4],zxhImageDataT<short>&IntensityImage,zxhImageDataT<short>&MaskImage,zxhImageDataT<float>&PatchProMask,float fdivec[3][3],string SavePathname)
{
	
		vector<PointCordProbTypeDef>vPointCord;
		//将正负样本合并，加入到容器中
		vPointCord.insert(vPointCord.end(),vVesCentPointCord.begin(),vVesCentPointCord.end());
		
		vPointCord.insert(vPointCord.end(),vSelecOffSetPointCord.begin(),vSelecOffSetPointCord.end());
		cout<<"Number of online Points and offline Points: "<<vVesCentPointCord.size()<<", "<<vSelecOffSetPointCord.size()<<endl;
		//将所有预备取patch的点打乱
		random_shuffle(vPointCord.begin(), vPointCord.end());
		int nsize=vPointCord.size();
		////将所有的patch中心点分成多个组，每组不超过30000个点
		int ngroup=nsize/30000;

		vector<vector<PointCordProbTypeDef>>vpatchPont;
		vector<PointCordProbTypeDef> tempcurgroupponts;
		for(int ptid=0;ptid<nsize;ptid++)
		{

			PointCordProbTypeDef tempcurpont;
			tempcurpont=vPointCord[ptid];
			tempcurgroupponts.push_back(tempcurpont);
			if (ptid==0) continue;
			if(tempcurgroupponts.size()%30000==0)//取到第30000个点
			{
				vpatchPont.push_back(tempcurgroupponts);
				tempcurgroupponts.clear();
				continue;
			}

		}
		vpatchPont.push_back(tempcurgroupponts);
		int x=0;

		for(int ig=0;ig<vpatchPont.size();ig++)
		{
			//按照组遍历，组编号为ig
			vector<PointCordProbTypeDef> tempcurgroupponts=vpatchPont[ig];
			int npt=tempcurgroupponts.size();
			int P_newsize[] = { SiglePatchSize,npt , 1, 1 };
			zxhImageDataT<short> PatchImage1,PatchImage2,PatchImage3;
			PatchImage1.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
			PatchImage2.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
			PatchImage3.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
			zxhImageData* PatchImageArray[3]={&PatchImage1,&PatchImage2,&PatchImage3};
			int PatchNumIdex = 0;
			string SavePatchinfo=SavePathname+"Patchinfo_"+num2str(ig)+"_"+num2str(numAxis)+".txt";
			ofstream outfile_norm(SavePatchinfo, ios::beg);//output
			outfile_norm << npt<<" "<<0<< "\n";
			for(int ptid=0;ptid<npt;ptid++)
			{

				//-------------------Center--------------------
				//获取当前点的patch
				float InputWorldCoord[4] ={ tempcurgroupponts[ptid].x,  tempcurgroupponts[ptid].y,  tempcurgroupponts[ptid].z, 0 };
				float inputimagecoord[4]={ tempcurgroupponts[ptid].x,  tempcurgroupponts[ptid].y,  tempcurgroupponts[ptid].z, 0 };
				if (LGeneratePatchExtractByWorldWithRandOffset_inxyzTrain(PatchNumIdex,fdivec,InputWorldCoord, PatchInfo,IntensityImage, PatchImageArray) == false)
				{
					continue;//

				}
				PatchNumIdex++;
				if (tempcurgroupponts[ptid].pr>1)
				{
					int x=0;
				}
				outfile_norm << PatchNumIdex << " " << tempcurgroupponts[ptid].pr<< "\n";
					//在mask中标记取过的点
					IntensityImage.GetImageInfo()->WorldToImage(inputimagecoord);
				int nscx = zxh::round(inputimagecoord[0]);
				int nscy = zxh::round(inputimagecoord[1]);
				int nscz = zxh::round(inputimagecoord[2]);
				MaskImage.SetPixelByGreyscale(nscx,nscy,nscz,0,1);
				PatchProMask.SetPixelByGreyscale(nscx,nscy,nscz,0,tempcurgroupponts[ptid].pr);
			}
			outfile_norm.close();
			string Str_T1 = SavePathname + "Patch_s1_"+num2str(ig)+"_"+num2str(numAxis)+".nii.gz";
			string Str_T2 = SavePathname + "Patch_s2_"+num2str(ig)+"_"+num2str(numAxis)+".nii.gz";
			string Str_T3 = SavePathname + "Patch_s3_"+num2str(ig)+"_"+num2str(numAxis)+".nii.gz";
			zxh::SaveImage(PatchImageArray[0], Str_T1);
			zxh::SaveImage(PatchImageArray[1], Str_T2);
			zxh::SaveImage(PatchImageArray[2], Str_T3);
		
		}
		    string Str_Mask = SavePathname + "P_Mask.nii.gz";
		    zxh::SaveImage(&MaskImage, Str_Mask);
	 string Str_ProbMask = SavePathname + "P_ProbMask.nii.gz";
		    zxh::SaveImage(&PatchProMask, Str_ProbMask);
			return true;
}
bool Gen_OnCentPointsUsingRot(int numAxis,vector<PointCordProbTypeDef> &vVesCentPointCord,int SiglePatchSize,float spacing111[4],int PatchInfo[4],zxhImageDataT<short>&IntensityImage,float fdivec[3][3],string SavePathname)
{
	
		vector<PointCordProbTypeDef>vPointCord;
		//将正负样本合并，加入到容器中
		vPointCord.insert(vPointCord.end(),vVesCentPointCord.begin(),vVesCentPointCord.end());

		//将所有预备取patch的点打乱
		random_shuffle(vPointCord.begin(), vPointCord.end());
		int nsize=vPointCord.size();
		////将所有的patch中心点分成多个组，每组不超过30000个点
		int ngroup=nsize/30000;

		vector<vector<PointCordProbTypeDef>>vpatchPont;
		vector<PointCordProbTypeDef> tempcurgroupponts;
		for(int ptid=0;ptid<nsize;ptid++)
		{

			PointCordProbTypeDef tempcurpont;
			tempcurpont=vPointCord[ptid];
			tempcurgroupponts.push_back(tempcurpont);
			if (ptid==0) continue;
			if(tempcurgroupponts.size()%30000==0)//取到第30000个点
			{
				vpatchPont.push_back(tempcurgroupponts);
				tempcurgroupponts.clear();
				continue;
			}

		}
		vpatchPont.push_back(tempcurgroupponts);
		int x=0;

		for(int ig=0;ig<vpatchPont.size();ig++)
		{
			//按照组遍历，组编号为ig
			vector<PointCordProbTypeDef> tempcurgroupponts=vpatchPont[ig];
			int npt=tempcurgroupponts.size();
			int P_newsize[] = { SiglePatchSize,npt , 1, 1 };
			zxhImageDataT<short> PatchImage1,PatchImage2,PatchImage3;
			PatchImage1.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
			PatchImage2.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
			PatchImage3.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
			zxhImageData* PatchImageArray[3]={&PatchImage1,&PatchImage2,&PatchImage3};
			int PatchNumIdex = 0;
			string SavePatchinfo=SavePathname+"Patchinfo_"+num2str(ig)+"_"+num2str(numAxis)+".txt";
			ofstream outfile_norm(SavePatchinfo, ios::beg);//output
			outfile_norm << npt<<" "<<0<< "\n";
			for(int ptid=0;ptid<npt;ptid++)
			{

				//-------------------Center--------------------
				//获取当前点的patch
				float InputWorldCoord[4] ={ tempcurgroupponts[ptid].x,  tempcurgroupponts[ptid].y,  tempcurgroupponts[ptid].z, 0 };
				float inputimagecoord[4]={ tempcurgroupponts[ptid].x,  tempcurgroupponts[ptid].y,  tempcurgroupponts[ptid].z, 0 };
				if (LGeneratePatchExtractByWorldWithRandOffset_inxyzTrain(PatchNumIdex,fdivec,InputWorldCoord, PatchInfo,IntensityImage, PatchImageArray) == false)
				{
					continue;//

				}
				PatchNumIdex++;
				if (tempcurgroupponts[ptid].pr>1)
				{
					int x=0;
				}
				outfile_norm << PatchNumIdex << " " << tempcurgroupponts[ptid].pr<< "\n";
					//在mask中标记取过的点
					IntensityImage.GetImageInfo()->WorldToImage(inputimagecoord);
				int nscx = zxh::round(inputimagecoord[0]);
				int nscy = zxh::round(inputimagecoord[1]);
				int nscz = zxh::round(inputimagecoord[2]);
			}
			outfile_norm.close();
			string Str_T1 = SavePathname + "Patch_s1_"+num2str(ig)+"_"+num2str(numAxis)+".nii.gz";
			string Str_T2 = SavePathname + "Patch_s2_"+num2str(ig)+"_"+num2str(numAxis)+".nii.gz";
			string Str_T3 = SavePathname + "Patch_s3_"+num2str(ig)+"_"+num2str(numAxis)+".nii.gz";
			zxh::SaveImage(PatchImageArray[0], Str_T1);
			zxh::SaveImage(PatchImageArray[1], Str_T2);
			zxh::SaveImage(PatchImageArray[2], Str_T3);
		
		}
			return true;
}
int main(int argc, char *argv[])
{
	//这个版本是用来抓取patch，和media18 J.W文章类似
	//但是只有on line 和 offline
	//在jdqPatchExtractForDNNvsls基础上，加入了rotation
	//if( argc < 11)
	//{
	//	cerr << "Usage: " << endl;
	//	cerr << "jdqPatchExtractForTrain	curve	image  labimag pathinfor " << endl;
	//	return -1;
	//}

	//string strintImg =string(argv[1]);
	//string strlabImg =string(argv[2]);
	//string strlabradImg =string(argv[3]);
	//string strROIImg =string(argv[4]); 
	//string strDisDir=string(argv[5]);  
	//string SavePathname = string(argv[6]);  
	//string trainortest = string(argv[7]);  
	//float fonlinerat= atof(argv[8]);  
	//float fofflinerat = atof(argv[9]);  
	//int HalfPatchLength=atoi(argv[10]);  
	//int N =atoi(argv[11]);  
	//int nrottim=atoi(argv[12]);  
	//--------------------------------------
	//if(argc<6) int x=0;
	string strintImg ="J:/work_jdq/for_DNN_vsls_v5/data/whole_GT/dataset00/image.nii.gz";
	string strlabImg ="J:/work_jdq/for_DNN_vsls_v5/data/whole_GT/dataset00/lab_image.nii.gz";
	string strlabradImg ="J:/work_jdq/for_DNN_vsls_v5/data/whole_GT/dataset00/lab_image1_rad.nii.gz";
	string strROIImg ="J:/work_jdq/for_DNN_vsls_v5/data/whole_GT/dataset00/whs_lab_image_ROI.nii.gz";
	string strDisDir="J:/work_jdq/for_DNN_vsls_v5/data/Otherdata/dispersedDirections/dispersedDirectionsSphere500.txt";
	string SavePathname ="J:/work_jdq/for_DNN_vsls_v5/Patch_data_ur/5_5_5/crossValidate/target03/infiles/";
	string trainortest = "-train";  
		float fonlinerat=1;  
	float fofflinerat = 100;  
	int HalfPatchLength=2;  
	int N =2;  
	int nrottim=4;
	//---------------patch大小相关参数设置---------------
	


	srand((unsigned)time(NULL));
	//char * bufferTscar=new char[1048576], * bufferTnormal=new char[1048576], *bufferNlink=new char[1048576] ;
	int SiglePatchSize = (N * 2 + 1)*(N * 2 + 1)*(HalfPatchLength * 2 + 1);
	//---------------读取图像---------------------------
	zxhImageDataT<short> IntensityImage,LabelImage,LabelradImage,MaskImage;
	zxhImageDataT<float>PatchProMask;

	//读取intensity image
	zxh::OpenImageSafe(&IntensityImage,strintImg);
	MaskImage.NewImage(IntensityImage.GetImageInfo());
	PatchProMask.NewImage(IntensityImage.GetImageInfo());
	//读取label image
	zxh::OpenImageSafe(&LabelImage,strlabImg);
	////读取label rad image
	zxh::OpenImageSafe(&LabelradImage,strlabradImg);
	//读取单位球离散化向量
	vector<PointCordTypeDef>vDisDir;
	const char *chDisDir=strDisDir.c_str();
	ReadDireTxt(chDisDir,vDisDir);

	int PatchInfo[4] = { N, HalfPatchLength, 0, 0 };
	float spacing111[] = { 1, 1, 1, 1 }; 
	float fdivec[3][3]={0};
	//x方向
	fdivec[0][0]=1;
	fdivec[0][1]=0;
	fdivec[0][2]=0;

	//y方向
	fdivec[1][0]=0;
	fdivec[1][1]=1;
	fdivec[1][2]=0;

	//z方向
	fdivec[2][0]=0;
	fdivec[2][1]=0;
	fdivec[2][2]=1;

	//将取patch的坐标轴，进行旋转，旋转后，三个坐标轴依然保持互相垂直
int nRottime=nrottim;//旋转 次数
vector<PointCordTypeDef> vxyzori;
PointCordTypeDef xori,yori,zori;
xori.x=1;xori.y=0;xori.z=0;
vxyzori.push_back(xori);
yori.x=0;yori.y=1;yori.z=0;
vxyzori.push_back(yori);
zori.x=0;zori.y=0;zori.z=1;
vxyzori.push_back(zori);
vector<vector<PointCordTypeDef>>vrandDir;
vrandDir.push_back(vxyzori);
Gen_fdivecByRotxyz(vxyzori,vrandDir,nRottime);

	//---------------如果是Train数据,

	if (strcmp(trainortest.c_str(),"-train")==0)
	{
		zxhImageDataT<short>  ROIImage;
		//读取ROI image
		zxh::OpenImageSafe(&ROIImage,strROIImg);
		//-----------------------------------------------读取label image， 并对各个label 的点数统计
		vector<PointCordProbTypeDef> vlab1,vlab2,vlab3;
		//label 1
		int numlab1=SelectPointsFromlabimg(LabelImage,vlab1,1);
		int numlab2=SelectPointsFromlabimg(LabelImage,vlab1,2);
		int numlab3=SelectPointsFromlabimg(LabelImage,vlab1,3);

		vector<PointCordProbTypeDef>vVesCentPointCord;
		vector<PointCordProbTypeDef> vSelecOffSetPointCord;

		//依次挑选正负样本的点
		float fvraPositive=fonlinerat;//选取的正样本占全部正样本的比例
		float flabbal=fofflinerat;//labbal=offline数量/全部正online数量
		cout<<"The ratios of the online samples and offline samples: "<<fvraPositive<<", "<<flabbal<<endl;
		GenOnAndOffPonts_Rand(LabelImage,LabelradImage,ROIImage,vDisDir,vVesCentPointCord,vSelecOffSetPointCord,flabbal, fvraPositive);
		//首先按照原始的坐标轴选取中心线上的点和offcenterline的点
		vector<PointCordTypeDef> vfirstaxis;
		vfirstaxis=vrandDir[0];
		//x方向
		fdivec[0][0]=vfirstaxis[0].x;
		fdivec[0][1]=vfirstaxis[0].y;
		fdivec[0][2]=vfirstaxis[0].z;

		//y方向
		fdivec[1][0]=vfirstaxis[1].x;
		fdivec[1][1]=vfirstaxis[1].y;
		fdivec[1][2]=vfirstaxis[1].z;

		//z方向
		fdivec[2][0]=vfirstaxis[2].x;
		fdivec[2][1]=vfirstaxis[2].y;
		fdivec[2][2]=vfirstaxis[2].z;
		Gen_OnCentandOffCentPoints(0,vVesCentPointCord,vSelecOffSetPointCord,SiglePatchSize,spacing111,PatchInfo,IntensityImage,MaskImage,PatchProMask,fdivec,SavePathname);
	//其次，对oncenterline上的点，按照旋转生成的坐标轴选取patch
		for(int k=1;k<vrandDir.size();k++)
		{
			vector<PointCordTypeDef> vlefttaxis;
		vlefttaxis=vrandDir[k];
		//x方向
		fdivec[0][0]=vlefttaxis[0].x;
		fdivec[0][1]=vlefttaxis[0].y;
		fdivec[0][2]=vlefttaxis[0].z;

		//y方向
		fdivec[1][0]=vlefttaxis[1].x;
		fdivec[1][1]=vlefttaxis[1].y;
		fdivec[1][2]=vlefttaxis[1].z;

		//z方向
		fdivec[2][0]=vlefttaxis[2].x;
		fdivec[2][1]=vlefttaxis[2].y;
		fdivec[2][2]=vlefttaxis[2].z;
				
		cout<<"Number of online Points: "<<vVesCentPointCord.size()<<", Using "<<num2str(k)<<"/"<<num2str(vrandDir.size()-1)<<" coordinate system. "<<endl;
		Gen_OnCentPointsUsingRot(k,vVesCentPointCord,SiglePatchSize,spacing111,PatchInfo,IntensityImage,fdivec,SavePathname);
		}
	}
	if (strcmp(trainortest.c_str(),"-test")==0)
	{


		zxhImageDataT<short>  ROIImage;
		//读取ROI image
		zxh::OpenImageSafe(&ROIImage,strROIImg);
		//读取中心点的坐标
		vector<PointCordProbTypeDef> vVesCentPointCord,vVesCentPointCordrad;//金标准点，概率为1,online
		int nvcp=GenerateVesCentPontsRadset(LabelImage,LabelradImage,vVesCentPointCordrad);
		//这是用来在testing data中取所有的patch
		//为了加速，挑选一定的层，并将层内划分为4*4个分块
		int ImgSizeraw[]={0,0,0,0};
		IntensityImage.GetImageSize(ImgSizeraw[0],ImgSizeraw[1],ImgSizeraw[2],ImgSizeraw[3]);
		int numofPatch=ImgSizeraw[2];
		//
		int numfenkuai=4;
		int slecnumofPatch=1;
		int PatchNumIdex=0;
		for (int i=0;i<numofPatch;i=i+slecnumofPatch)
			//for (int i=53;i<54;i=i+slecnumofPatch)
		{
			vector<PointCordTypeDef>vPointCord;

			for (int XN=0;XN<numfenkuai;XN++)
				for (int YN=0;YN<numfenkuai;YN++)
				{
					vPointCord.clear();
					int XNStart=XN*ImgSizeraw[0]/numfenkuai;
					int YNStart=YN*ImgSizeraw[1]/numfenkuai;
					PatchNumIdex=0;

					for (int j=YNStart;j<YNStart+ImgSizeraw[1]/numfenkuai;j++)
						for (int k=XNStart;k<XNStart+ImgSizeraw[0]/numfenkuai;k++)
						{

							short shintenROI=ROIImage.GetPixelGreyscale( k,j,i,0);
							if (shintenROI==0)continue;
							float curvpointsWolrd[] = { k,j,i,0 };
							IntensityImage.GetImageInfo()->ImageToWorld(curvpointsWolrd);//世界坐标转成图像坐标
							PointCordTypeDef TempPoint;
							TempPoint.x =curvpointsWolrd[0];
							TempPoint.y =curvpointsWolrd[1];
							TempPoint.z =curvpointsWolrd[2];
							vPointCord.push_back(TempPoint);

						}
						int npt=vPointCord.size();
						if (npt==0)continue;

						string SavePatchinfo=SavePathname+"Patchinfo_"+num2str(i)+"_" +num2str(YN)+"_"+num2str(XN)+".txt";
						ofstream outfile_norm(SavePatchinfo, ios::beg);//output
						outfile_norm << npt<<" "<<0<< "\n";
						zxhImageDataT<short> PatchImage1,PatchImage2,PatchImage3;
						int P_newsize[] = { SiglePatchSize,npt , 1, 1 };
						PatchImage1.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
						PatchImage2.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
						PatchImage3.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
						zxhImageData* PatchImageArray[3]={&PatchImage1,&PatchImage2,&PatchImage3};


						//生成一个patch
						for(int ptid=0;ptid<vPointCord.size();ptid++)
						{
							//获取当前点的patch
							float InputWorldCoord[4] ={ vPointCord[ptid].x,  vPointCord[ptid].y,  vPointCord[ptid].z, 0 };
							float inputimagecoord[4]={ vPointCord[ptid].x,  vPointCord[ptid].y,  vPointCord[ptid].z, 0 };
							LabelImage.GetImageInfo()->WorldToImage(inputimagecoord);
							int scx = zxh::round(inputimagecoord[0]);
							int scy = zxh::round(inputimagecoord[1]);
							int scz = zxh::round(inputimagecoord[2]);
							float finten=LabelImage.GetPixelGreyscale(scx,scy,scz,0);
							MaskImage.SetPixelByGreyscale(scx,scy,scz,0,1);
						
							if (LGeneratePatchExtractByWorldWithRandOffset_inxyzTest(PatchNumIdex,fdivec,InputWorldCoord, PatchInfo,IntensityImage, PatchImageArray) == false)
							{
								continue;//

							}
							short shlabnum =LabelImage.GetPixelGreyscale(scx,scy,scz);
							//根据train的函数设定Test金标准

							float fbprob=Calc_probOfPont_MedIA18(InputWorldCoord,vVesCentPointCordrad);
							PatchNumIdex++;
							outfile_norm << PatchNumIdex << " " <<fbprob<< "\n";
								PatchProMask.SetPixelByGreyscale(scx,scy,scz,0,fbprob);
						}
						outfile_norm.close();

						string Str_T1 = SavePathname + "Patch_s1_"+num2str(i)+"_" +num2str(YN)+"_"+num2str(XN)+".nii.gz";
						string Str_T2 = SavePathname + "Patch_s2_"+num2str(i)+"_" +num2str(YN)+"_"+num2str(XN)+".nii.gz";
						string Str_T3 = SavePathname + "Patch_s3_"+num2str(i)+"_" +num2str(YN)+"_"+num2str(XN)+".nii.gz";
						zxh::SaveImage(PatchImageArray[0], Str_T1);
						zxh::SaveImage(PatchImageArray[1], Str_T2);
						zxh::SaveImage(PatchImageArray[2], Str_T3);
				}
		}
			string Str_Mask = SavePathname + "P_Mask.nii.gz";
		    zxh::SaveImage(&MaskImage, Str_Mask);
			 string Str_ProbMask = SavePathname + "P_ProbMask.nii.gz";
		    zxh::SaveImage(&PatchProMask, Str_ProbMask);
	}
	if (strcmp(trainortest.c_str(),"-test_v1")==0)
	{
		//	//这是用来在testing data中取所有的patch
		//	//为了加速，挑选一定的层，并将层内划分为4*4个分块
		//	int ImgSizeraw[]={0,0,0,0};
		//	IntensityImage.GetImageSize(ImgSizeraw[0],ImgSizeraw[1],ImgSizeraw[2],ImgSizeraw[3]);
		//	int numofPatch=ImgSizeraw[2];
		//	//
		//	int numfenkuai=4;
		//	int slecnumofPatch=1;
		//	int PatchNumIdex=0;
		//	for (int i=0;i<numofPatch;i=i+slecnumofPatch)
		//	{
		//		vector<PointCordTypeDef>vPointCord;

		//		for (int XN=0;XN<numfenkuai;XN++)
		//			for (int YN=0;YN<numfenkuai;YN++)
		//			{
		//				vPointCord.clear();
		//				int XNStart=XN*ImgSizeraw[0]/numfenkuai;
		//				int YNStart=YN*ImgSizeraw[1]/numfenkuai;
		//				PatchNumIdex=0;
		//				string SavePatchinfo=SavePathname+"Patchinfo_"+num2str(i)+"_" +num2str(YN)+"_"+num2str(XN)+".txt";
		//				ofstream outfile_norm(SavePatchinfo, ios::beg);//output
		//				for (int j=YNStart;j<YNStart+ImgSizeraw[1]/numfenkuai;j++)
		//					for (int k=XNStart;k<XNStart+ImgSizeraw[0]/numfenkuai;k++)
		//					{

		//						float curvpointsWolrd[] = { k,j,i, 0 };
		//						IntensityImage.GetImageInfo()->ImageToWorld(curvpointsWolrd);//物理坐标转成图像坐标
		//						PointCordTypeDef TempPoint;
		//						TempPoint.x =curvpointsWolrd[0];
		//						TempPoint.y =curvpointsWolrd[1];
		//						TempPoint.z =curvpointsWolrd[2];
		//						vPointCord.push_back(TempPoint);

		//					}
		//					int npt=vPointCord.size();
		//					outfile_norm << npt<<" "<<0<< "\n";
		//					zxhImageDataT<short> PatchImage1,PatchImage2,PatchImage3;
		//					int P_newsize[] = { SiglePatchSize,npt , 1, 1 };
		//					PatchImage1.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
		//					PatchImage2.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
		//					PatchImage3.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
		//					zxhImageData* PatchImageArray[3]={&PatchImage1,&PatchImage2,&PatchImage3};


		//					//生成一个patch
		//					for(int ptid=0;ptid<vPointCord.size();ptid++)
		//					{
		//						//获取当前点的patch
		//						float InputWorldCoord[4] ={ vPointCord[ptid].x,  vPointCord[ptid].y,  vPointCord[ptid].z, 0 };
		//						float inputimagecoord[4]={ vPointCord[ptid].x,  vPointCord[ptid].y,  vPointCord[ptid].z, 0 };
		//						LabelImage.GetImageInfo()->WorldToImage(inputimagecoord);
		//						int scx = zxh::round(inputimagecoord[0]);
		//						int scy = zxh::round(inputimagecoord[1]);
		//						int scz = zxh::round(inputimagecoord[2]);
		//						float finten=LabelImage.GetPixelGreyscale(scx,scy,scz,0);


		//						if (LGeneratePatchExtractByWorldWithRandOffset_inxyzTest(PatchNumIdex,fdivec,InputWorldCoord, PatchInfo,IntensityImage, PatchImageArray) == false)
		//						{
		//							continue;//

		//						}
		//						short shlabnum =LabelImage.GetPixelGreyscale(scx,scy,scz);
		//						float fprob=0;
		//						if(shlabnum==0)
		//						{
		//							fprob=0;

		//						}
		//						if(shlabnum==2)
		//						{
		//							fprob=0.5;

		//						}
		//						if(shlabnum==1)
		//						{
		//							fprob=1;

		//						}
		//						PatchNumIdex++;
		//						outfile_norm << PatchNumIdex << " " <<fprob<< "\n";

		//					}
		//					outfile_norm.close();

		//					string Str_T1 = SavePathname + "Patch_s1_"+num2str(i)+"_" +num2str(YN)+"_"+num2str(XN)+".nii.gz";
		//					string Str_T2 = SavePathname + "Patch_s2_"+num2str(i)+"_" +num2str(YN)+"_"+num2str(XN)+".nii.gz";
		//					string Str_T3 = SavePathname + "Patch_s3_"+num2str(i)+"_" +num2str(YN)+"_"+num2str(XN)+".nii.gz";
		//					zxh::SaveImage(PatchImageArray[0], Str_T1);
		//					zxh::SaveImage(PatchImageArray[1], Str_T2);
		//					zxh::SaveImage(PatchImageArray[2], Str_T3);
		//			}
		//	}
	}
	return 0;
}

