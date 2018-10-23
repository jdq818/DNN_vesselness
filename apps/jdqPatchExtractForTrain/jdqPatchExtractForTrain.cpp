

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
int CalculatePatchsize(zxhImageData &IntensityImage,std::vector<jdq2017::point3D>&vcurveponts,std::vector<jdq2017::point3D>&vpoints)
{
	int n=0;
	int ImgSizeraw[]={0,0,0,0};
	IntensityImage.GetImageSize(ImgSizeraw[0],ImgSizeraw[1],ImgSizeraw[2],ImgSizeraw[3]);
	short *shMask=new short[ImgSizeraw[0]*ImgSizeraw[1]*ImgSizeraw[2]-1];
	for (int i = 0; i < ImgSizeraw[0]*ImgSizeraw[1]*ImgSizeraw[2]-1; i ++)shMask[i]=0;
	for (int i=0;i<vcurveponts.size();i++)
	{
		float curvpointsWolrd[] = { vcurveponts[i]._x,  vcurveponts[i]._y,  vcurveponts[i]._z, 0 };
		IntensityImage.GetImageInfo()->WorldToImage(curvpointsWolrd);//物理坐标转成图像坐标
		int scx = zxh::round(curvpointsWolrd[0]);
		int scy = zxh::round(curvpointsWolrd[1]);
		int scz = zxh::round(curvpointsWolrd[2]);
		bool bIsInsideImage = IntensityImage.InsideImage(scx, scy, scz, 0); // 超过图像边界的，不给予考虑，也就是说，默认为normal myo
		if (!bIsInsideImage)
		{
			std::cout << "warning: node " << i << " not inside image " << "\n"; //-----------------
			continue;
		}
		int scxyz[]={scx, scy, scz};
		int scxyznp=0;
		TransDatanum2Greynum(scxyz,scxyznp,ImgSizeraw);
		if(shMask[scxyznp]==1)
		{
			continue;
		}
		else
		{

			vpoints.push_back(vcurveponts[i]);
			shMask[scxyznp]=1;
			n++;
		}
	}
	return n;
}
bool Calc3vectorsFromTangetVecAndPoints(float  fdivec1[3],float  fPointWord[3],float fdivec[3][3])
{

	//通过一个向量和其通过的点求解三个垂直向量


	//get the second direction of the plane
	float fdivec2[3]={0,0,0};
	Get_dirvec_byPrj(fdivec1,fPointWord,fdivec2);
	//get the third direction of the plane
	float fdivec3[3]={0,0,0};
	zxh::VectorOP_CrossProduct3D(fdivec1,fdivec2,fdivec3);
	zxh::VectorOP_Normalise(fdivec3,3);
	zxh::VectorOP_Normalise(fdivec2,3);
	zxh::VectorOP_Normalise(fdivec1,3) ;
	float x1=zxh::absf(zxh::VectorOP_Cosine(fdivec1,fdivec2,2));
	float x2=zxh::absf(zxh::VectorOP_Cosine(fdivec1,fdivec2,3));
	float x3=zxh::absf(zxh::VectorOP_Cosine(fdivec2,fdivec2,3));
	if( zxh::absf(zxh::VectorOP_Cosine(fdivec1,fdivec2,3))>1.0e-3 ||zxh::absf(zxh::VectorOP_Cosine(fdivec1,fdivec3,3))>1.0e-3 ||zxh::absf(zxh::VectorOP_Cosine(fdivec1,fdivec2,3))>1.0e-3 )
	{
		std::cerr<<"error: axises should be perpendicular\n";
		return false ;
	}
	if( zxh::absf(zxh::VectorOP_Magnitude(fdivec1,3)-1.0 )>1.0e-3 || zxh::absf(zxh::VectorOP_Magnitude(fdivec3,3)-1.0 )>1.0e-3 || zxh::absf(zxh::VectorOP_Magnitude(fdivec2,3)-1.0 )>1.0e-3 )
	{
		std::cerr<<"error: axis vectors should be normalized\n";
		return false ;
	}
	fdivec[0][0]=fdivec1[0];
	fdivec[0][1]=fdivec1[1];
	fdivec[0][2]=fdivec1[2];

	fdivec[1][0]=fdivec2[0];
	fdivec[1][1]=fdivec2[1];
	fdivec[1][2]=fdivec2[2];

	fdivec[2][0]=fdivec3[0];
	fdivec[2][1]=fdivec3[1];
	fdivec[2][2]=fdivec3[2];
}
bool LGeneratePatchExtractByWorldWithRandOffset(int PatchNumIdex,float fdivec[3][3],float InputOutputWorldCoord[], int PatchSize[],zxhImageData &IntensityImage,zxhImageData &PatchImage)
{

	int N = PatchSize[0];
	int HalfPatchLength = PatchSize[1]; 



	int PatchPointWorldCoord[3] = { 0 };
	//Computing the intensity of points in the tangent plane using interpolation	
	//Interpolation
	zxhImageModelingLinear InterpolationMod;
	InterpolationMod.SetImage(&IntensityImage);
	for (int i = -N; i < N + 1; i++)
	{
		for (int j = -N; j < N + 1; j++)
		{
			for (int k = -HalfPatchLength; k < HalfPatchLength + 1; k++)//extend the patch along the gradient orientation
			{ 
				PatchPointWorldCoord[0] =InputOutputWorldCoord[0] + i*fdivec[0][0] + j*fdivec[1][0] + k*fdivec[2][0];
				PatchPointWorldCoord[1] = InputOutputWorldCoord[1] + i*fdivec[0][1] + j*fdivec[1][1] + k*fdivec[2][1];
				PatchPointWorldCoord[2] =InputOutputWorldCoord[2]+  i*fdivec[0][2] + j*fdivec[1][2] + k*fdivec[2][2];

				//InterpolationMod.SetImage(pLALabel);//check the orientation of gradient				
				float IntensityValue = InterpolationMod.GetPixelFloatValueWithCheckByWorld(PatchPointWorldCoord[0], PatchPointWorldCoord[1], PatchPointWorldCoord[2], 0);

				int W2I_Coor_SetPatchPixel = (i + N) + (j + N) * (2 * N + 1) + (k + HalfPatchLength)* (2 * N + 1)* (2 * N + 1);
				PatchImage.SetPixelByGreyscale(W2I_Coor_SetPatchPixel, PatchNumIdex, 0, 0, IntensityValue);

			}
		}
	}

	return true ;

}


int main(int argc, char *argv[])
{
	//if( argc < 5 )
	//{
	//	cerr << "Usage: " << endl;
	//	cerr << "jdqPatchExtractForTrain.cpp	imageRaw(.nii)	labelimage.nii Line or points(.vtk) results -NorR" << endl;
	//	return -1;
	//}


	//输入一条曲线，一个intensity image和 一个label image
	string mainfold="E:\\work_jdq\\DL_train\\ForDNNvsls\\dataset00\\";
	string strCurvename =mainfold+"centerlines\\vessel0.txt";
	string strintImg =mainfold + "image\\image.nii.gz";
	string strlabImg =mainfold + "image\\lab_image.nii.gz";
	string strResultfold =mainfold;

	string SavePatchinfo = mainfold + "Patch_info.txt";  
	string Str_T = mainfold + "Patch.nii.gz";

	//---------------相关参数设置---------------
	int HalfPatchLength = 2;
	int N =2;
	int Offset =8;
	int PatchNumIdex = 0;

	srand((unsigned)time(NULL));
	//char * bufferTscar=new char[1048576], * bufferTnormal=new char[1048576], *bufferNlink=new char[1048576] ;
	int SiglePatchSize = (N * 2 + 1)*(N * 2 + 1)*(HalfPatchLength * 2 + 1);
	//---------------读取图像---------------------------
	zxhImageDataT<short> IntensityImage, LabelImage, PatchImage;

	//读取intensity image
	zxh::OpenImageSafe(&IntensityImage,strintImg);
	//读取label image
	zxh::OpenImageSafe(&LabelImage,strlabImg);
	//---------------读取中心线,并resample---------------------------

	std::vector<jdq2017::point3D>vcurveponts;
	vcurveponts.clear();
	ReadCurve(strCurvename,vcurveponts);

	//
	cout<<"Resample the curve by image spacing!"<<endl;
	float fImgSpacing[]={1,1,1,1};//Add by JDQ
	IntensityImage.GetImageSpacing(fImgSpacing[0],fImgSpacing[1],fImgSpacing[2],fImgSpacing[3] );//Add by JDQ
	float fresamle=zxh::minf(fImgSpacing[1],fImgSpacing[2])*0.5;
	cout<<"Old curve points: "<<vcurveponts.size()-1<<endl;
	ResampleCurve(vcurveponts,fresamle);
	cout<<"New curve points: "<<vcurveponts.size()-1<<endl;
	//-----------生成 Patch Image ----------------------
	//计算patch的数量；size(patch)=size(imagepoints(vcurveponts))

	std::vector<jdq2017::point3D>vponts;
	vponts.clear();
	int nptsize=CalculatePatchsize(IntensityImage,vcurveponts,vponts);
	vector<PointCordTypeDef> vPointCord;
	for (int i = 0; i < vponts.size(); i++)
	{
		PointCordTypeDef TempPoint;
		TempPoint.x =vponts[i]._x;
		TempPoint.y =vponts[i]._y;
		TempPoint.z =vponts[i]._z;
		vPointCord.push_back(TempPoint);
	}
	float spacing111[] = { 1, 1, 1, 1 }; 
	ofstream outfile_norm(SavePatchinfo, ios::beg);//output
	outfile_norm << nptsize << "\n";
	//定义patch library的size

	int P_newsize[] = { SiglePatchSize, 4*nptsize, 1, 1 };
	PatchImage.NewImage(2, P_newsize, spacing111, IntensityImage.GetImageInfo());
	zxhImageData* ImageArray[4] = { &IntensityImage, &LabelImage, &PatchImage };
	//定义一个mask，避免重复取点
	int PatchInfo[4] = { N, HalfPatchLength, 0, 0 };

	for(int ptid=0;ptid<vPointCord.size();ptid++)
	{
		if(ptid==157)
			int xxx=0;
		//-------------------Center--------------------
		//获取当前点的patch
		float InputWorldCoord[4] ={ vPointCord[ptid].x,  vPointCord[ptid].y,  vPointCord[ptid].z, 0 };
		float inputimagecoord[4]={ vPointCord[ptid].x,  vPointCord[ptid].y,  vPointCord[ptid].z, 0 };
		LabelImage.GetImageInfo()->WorldToImage(inputimagecoord);
		int scx = zxh::round(inputimagecoord[0]);
		int scy = zxh::round(inputimagecoord[1]);
		int scz = zxh::round(inputimagecoord[2]);
		float finten=LabelImage.GetPixelGreyscale(scx,scy,scz,0);

		float fdivec[3][3]={0};//{fdivec[0][0],fdivec[0][1],fdivec[0][2]}是切线方向
		float ftdivec[3]={0};
		//计算切向量（第一个向量）

		Calc_divec(ptid,ftdivec,vPointCord,fresamle);
		if (zxh::VectorOP_Magnitude(ftdivec,3) <ZXH_FloatInfinitesimal)
		{
			std::cout << "error: magnitude too small for node " << ptid << "\n";
			return false ;
		}
		//通过切向量和通过点计算其他两个相互垂直的向量
		Calc3vectorsFromTangetVecAndPoints(ftdivec,InputWorldCoord, fdivec);
		if (LGeneratePatchExtractByWorldWithRandOffset(PatchNumIdex,fdivec,InputWorldCoord, PatchInfo,IntensityImage, PatchImage) == false)
		{
			continue;//

		}
		//-------------------附近三个点--------------------
		short labnum[4]={-1,-1,-1,-1};
		labnum[1]=1;
		int numNP=1;
		PatchNumIdex++;
		outfile_norm << PatchNumIdex << " " <<2<< "\n";//金标
		//法平面随机取点，直到取出三个点，都是不同的label值，如2，3，0；
		srand((unsigned)time(NULL));  


		while(PatchNumIdex<3000)
		{

		float randoffset1 = float((rand() % (200 * Offset))) / 100.0;	
		float randoffset2 = float((rand() % (200 * Offset))) / 100.0;	

		float fvec1[3]={randoffset1*fdivec[1][0],randoffset1*fdivec[1][1],randoffset1*fdivec[1][2]};
		float fvec2[3]={randoffset2*fdivec[2][0],randoffset2*fdivec[2][1],randoffset2*fdivec[2][2]};
		float fsumvec[3]={fvec1[0]+fvec2[0],fvec1[1]+fvec2[1],fvec1[2]+fvec2[2]};

		//Rand shift
		zxh::VectorOP_Normalise(fsumvec,3);
		float randradius= float((rand() % (200 * Offset))) / 100.0;

		float InputWorldCoord_x = InputWorldCoord[0] + randradius*fsumvec[0];
		float InputWorldCoord_y = InputWorldCoord[1] + randradius*fsumvec[1];
		float InputWorldCoord_z = InputWorldCoord[2] + randradius*fsumvec[2];
		//
		//

		float InputNeiPWorldCoord[4]={InputWorldCoord_x,InputWorldCoord_y,InputWorldCoord_y,0};
		//cout<< "warning: radius"<<randradius <<endl;
		//
		float InputNeiNodeWorldCoord[4]={InputWorldCoord_x,InputWorldCoord_y,InputWorldCoord_z,0 };
		IntensityImage.GetImageInfo()->WorldToImage(InputNeiNodeWorldCoord);
		int nscx = zxh::round(InputNeiNodeWorldCoord[0]);
		int nscy = zxh::round(InputNeiNodeWorldCoord[1]);
		int nscz = zxh::round(InputNeiNodeWorldCoord[2]);

		bool bIsInsideImage = IntensityImage.InsideImage(nscx, nscy, nscz, 0); // 超过图像边界的，不给予考虑，也就是说，默认为normal myo
		if (!bIsInsideImage)
		{
			std::cout << "warning: niebour node of point"<< ptid<< "is not inside image " << "\n"; //-----------------
			continue;
		}
		short shlabnum=LabelImage.GetPixelGreyscale(nscx, nscy, nscz, 0);
		bool bexitlab=false;
		for (int k=0;k<4;k++)
		{
			int ndinten=labnum[k];
			if (ndinten==shlabnum)
			{
				bexitlab=true;
				break;
			}
		}

		if(bexitlab)
		{
			continue;
		}
		else
		{
			//计算取patch主方向
			float fNormsumvec[3]={fsumvec[0],fsumvec[1],fsumvec[2],};
			zxh::VectorOP_Normalise(fNormsumvec,3);

			if (zxh::VectorOP_Magnitude(fNormsumvec,3) <ZXH_FloatInfinitesimal)
			{
				std::cout << "error: magnitude too small for node " << ptid << "\n";
				return false ;
			}
			//通过切向量和通过点计算其他两个相互垂直的向量
			float fdivecN[3][3]={0};
			Calc3vectorsFromTangetVecAndPoints(fNormsumvec,InputNeiNodeWorldCoord, fdivecN);

			if (LGeneratePatchExtractByWorldWithRandOffset(PatchNumIdex,fdivecN,InputNeiPWorldCoord, PatchInfo, IntensityImage, PatchImage) == false)
			{    
				continue;//
			}
			PatchNumIdex++;
			numNP++;


		}
		float fprob=0;
		if(shlabnum==0)
		{
			fprob=0;
		}
		if(shlabnum==2)
		{
			fprob=1;
		}
		if(shlabnum==3)
		{
			fprob=0.5;
		}
		labnum[shlabnum]=shlabnum;	
		outfile_norm << PatchNumIdex << " " <<fprob<< "\n";
		if(numNP==4) 
			break;	

		}//while	

	}
	outfile_norm.close();
	short x=PatchImage.GetPixelGreyscale(0,0,0,0);
	zxh::SaveImage(&PatchImage, Str_T);
	return 0;
}

