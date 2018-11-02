

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



#include "zxhImageGipl.h" 
#include "zxhImageData.h"
#include "zxhImageNifti.h"

#include "PatchExtractByWorldCoordinate.h"
#include <io.h>  
#include <direct.h> //创建文件夹
using namespace std;

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
bool BoundaryCorrect(int *PointPos,int ImgNewVslsSize[4])
{
	for (int i=0;i<3;i++)
	{
		PointPos[i]=zxh::maxf(0,PointPos[i]);
		PointPos[i]=zxh::minf(ImgNewVslsSize[i]-1,PointPos[i]);
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
					/*	if(intlabinte==3)
					{
					Ponttemp.pr=0.2;

					}*/
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
int GenerateRandomNonVesPontsPrbset(zxhImageData &LabelImage,zxhImageData &ROIImage,vector<PointCordProbTypeDef>&vNonVesPointCord)
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
int main(int argc, char *argv[])
{
	if( argc < 4 )
	{
		cerr << "Usage: " << endl;
		cerr << "jdqPatchExtractForTrain	curve	image  labimag pathinfor " << endl;
		return -1;
	}

	string strintImg =string(argv[1]);
	string strlabImg =string(argv[2]);
	string strROIImg =string(argv[3]); 
	string SavePathname = string(argv[4]);  
	string trainortest = string(argv[5]);  


	//--------------------------------------
	//string strintImg ="J:/work_jdq/for_DNN_vsls_v2/data_DNN/train/dataset01/image.nii.gz";
	//string strlabImg ="J:/work_jdq/for_DNN_vsls_v2/data_DNN/train/dataset01/lab_image.nii.gz";
	//string strROIImg ="J:/work_jdq/for_DNN_vsls_v2/data_DNN/train/dataset01/whs_lab_image_ROI.nii.gz";
	//string SavePathname ="J:/work_jdq/infiles/for_DNN_vsls_v2/Patch_data/dataset01/";
	//string trainortest = "-train";  

	//---------------相关参数设置---------------
	int HalfPatchLength = 2;
	int N =2;


	srand((unsigned)time(NULL));
	//char * bufferTscar=new char[1048576], * bufferTnormal=new char[1048576], *bufferNlink=new char[1048576] ;
	int SiglePatchSize = (N * 2 + 1)*(N * 2 + 1)*(HalfPatchLength * 2 + 1);
	//---------------读取图像---------------------------
	zxhImageDataT<short> IntensityImage,LabelImage, ROIImage;

	//读取intensity image
	zxh::OpenImageSafe(&IntensityImage,strintImg);
	//读取label image
	zxh::OpenImageSafe(&LabelImage,strlabImg);

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
	//---------------如果是Train数据,

	if (strcmp(trainortest.c_str(),"-train")==0)
	{
		zxhImageDataT<short>  ROIImage;
	//读取ROUI image
	zxh::OpenImageSafe(&ROIImage,strROIImg);
		//-----------------读取label image 并记录label1,label2,label3的所有点,并随机打乱---------------------------
		vector<PointCordProbTypeDef> vVesPointCord,vNonVesPointCord,vPointCord;//patch的中心点
		int labbal=100;//labbal=负样本数量/正样本数量
		int nvptsize=GenerateVesPontsPrbset(LabelImage,vVesPointCord);
		int nnonvptsize=GenerateRandomNonVesPontsPrbset(LabelImage,ROIImage,vNonVesPointCord);
		random_shuffle(vNonVesPointCord.begin(), vNonVesPointCord.end());
		vPointCord.insert(vPointCord.end(),vVesPointCord.begin(),vVesPointCord.end());
		vPointCord.insert(vPointCord.end(),vNonVesPointCord.begin(),min(vNonVesPointCord.begin()+nvptsize*labbal,vNonVesPointCord.end()));
		random_shuffle(vPointCord.begin(), vPointCord.end());
		int nsize=vPointCord.size();
		//将所有的patch中心点分成多个组，每组不超过30000个点
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
			string SavePatchinfo=SavePathname+"Patchinfo_"+num2str(ig)+".txt";
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
				outfile_norm << PatchNumIdex << " " << tempcurgroupponts[ptid].pr<< "\n";
			}
			outfile_norm.close();
			string Str_T1 = SavePathname + "Patch_s1_"+num2str(ig)+".nii.gz";
			string Str_T2 = SavePathname + "Patch_s2_"+num2str(ig)+".nii.gz";
			string Str_T3 = SavePathname + "Patch_s3_"+num2str(ig)+".nii.gz";
			zxh::SaveImage(PatchImageArray[0], Str_T1);
			zxh::SaveImage(PatchImageArray[1], Str_T2);
			zxh::SaveImage(PatchImageArray[2], Str_T3);
		}
	}
		if (strcmp(trainortest.c_str(),"-test")==0)
	{
		//这是用来在testing data中取所有的patch
		//为了加速，挑选一定的层，并将层内划分为4*4个分块
		int ImgSizeraw[]={0,0,0,0};
		IntensityImage.GetImageSize(ImgSizeraw[0],ImgSizeraw[1],ImgSizeraw[2],ImgSizeraw[3]);
		int numofPatch=ImgSizeraw[2];
		//
		int numfenkuai=1;
		int slecnumofPatch=1;
		int PatchNumIdex=0;
		for (int i=0;i<numofPatch;i=i+slecnumofPatch)
		{
			vector<PointCordTypeDef>vPointCord;

			for (int XN=0;XN<numfenkuai;XN++)
				for (int YN=0;YN<numfenkuai;YN++)
				{
					vPointCord.clear();
					int XNStart=XN*ImgSizeraw[0]/numfenkuai;
					int YNStart=YN*ImgSizeraw[1]/numfenkuai;
					PatchNumIdex=0;
					string SavePatchinfo=SavePathname+"Patchinfo_"+num2str(i)+"_" +num2str(YN)+"_"+num2str(XN)+".txt";
					ofstream outfile_norm(SavePatchinfo, ios::beg);//output
					for (int j=YNStart;j<YNStart+ImgSizeraw[1]/numfenkuai;j++)
						for (int k=XNStart;k<XNStart+ImgSizeraw[0]/numfenkuai;k++)
						{

							float curvpointsWolrd[] = { k,j,i, 0 };
							IntensityImage.GetImageInfo()->ImageToWorld(curvpointsWolrd);//物理坐标转成图像坐标
							PointCordTypeDef TempPoint;
							TempPoint.x =curvpointsWolrd[0];
							TempPoint.y =curvpointsWolrd[1];
							TempPoint.z =curvpointsWolrd[2];
							vPointCord.push_back(TempPoint);

						}
						int npt=vPointCord.size();
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


							if (LGeneratePatchExtractByWorldWithRandOffset_inxyzTest(PatchNumIdex,fdivec,InputWorldCoord, PatchInfo,IntensityImage, PatchImageArray) == false)
							{
								continue;//

							}
							short shlabnum =LabelImage.GetPixelGreyscale(scx,scy,scz);
							float fprob=0;
							if(shlabnum==0)
							{
								fprob=0;

							}
							if(shlabnum==2)
							{
								fprob=0.5;

							}
							if(shlabnum==1)
							{
								fprob=1;

							}
							PatchNumIdex++;
							outfile_norm << PatchNumIdex << " " <<fprob<< "\n";

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
	}
	return 0;
}

