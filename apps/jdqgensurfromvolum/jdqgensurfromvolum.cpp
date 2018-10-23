

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
using namespace std;


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
bool Generate_sur_from_vol(short *smaskData,int nlabnum,zxhImageDataT<short> &imglabsur)
{
	int gNbr[6][3] = { {-1, 0, 0}, \
	{ 1, 0, 0}, \
	{ 0,-1, 0}, \
	{ 0, 1, 0}, \
	{ 0, 0,-1}, \
	{ 0, 0, 1} };
	int ImgNewSize[4]={1};
	imglabsur.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	int NPointPos[3]={0,0,0};
	short NPointinten[6]={0,0,0,0,0,0};
	for(int it=0;it<ImgNewSize[3];++it)
		for(int iz=0;iz<ImgNewSize[2];++iz)
			for(int iy=0;iy<ImgNewSize[1];++iy)
				for(int ix=0;ix<ImgNewSize[0];++ix)
				{
					int snpd=0;
					int np[3]={ix,iy,iz};
					TransDatanum2Greynum(np,snpd,ImgNewSize);
					bool bsamlab=true;
					float fG[3]={0,0,0};
					short curpinten=smaskData[snpd];
					float fnormfG=0;
					if(ix==119&&iy==77&&iz==54)
						int x119=0;
					if(curpinten!=0)
					{
						for (int i = 0; i < 6; i++)
						{
							NPointPos[0]= ix + gNbr[i][0];
							NPointPos[1]= iy + gNbr[i][1];
							NPointPos[2]= iz + gNbr[i][2];
							BoundaryCorrect(NPointPos,ImgNewSize);
							int nNpd=0;
							TransDatanum2Greynum(NPointPos,nNpd,ImgNewSize);
							NPointinten[i]=smaskData[nNpd];
						}
						fG[0]=0.5*(NPointinten[1]-NPointinten[0]);
						fG[1]=0.5*(NPointinten[3]-NPointinten[2]);
						fG[2]=0.5*(NPointinten[5]-NPointinten[4]);
						zxh::VectorOP_Normalise(fG,3);
						fnormfG=zxh::VectorOP_Magnitude(fG,3);	
						if (fnormfG!=0)
						{
							int nG=1;
						imglabsur.SetPixelByGreyscale(ix,iy,iz,it,nG);
						}
					}
				}
				return true;
}

int main(int argc, char *argv[])
{
	//if( argc < 4 )
	//{
	//	cerr << "Usage: " << endl;
	//	cerr << "jdqMapPont2Img.cpp	imageRaw(.nii)	Line or points(.vtk) results -NorR" << endl;
	//	return -1;
	//}
	string strFileNameRaw ="E:/work_jdq/RCAAEF/train/rcaaef_32_img00/whs_lab_image.nii.gz";
	int nlabnum=500;
	char *chResultName="E:/work_jdq/RCAAEF/train/rcaaef_32_img00/sur_whs_lab_image500.nii.gz";
	zxhImageDataT<short> imglabvol;
	zxhImageDataT<short> imglabsur;
	if( zxh::OpenImage( &imglabvol, strFileNameRaw ) == false )
	{
		std::cerr << "Raw image(nifti-file) is not found!"<<endl; 
		return -1;
	}
	imglabsur.NewImage( imglabvol.GetImageInfo() );
	int ImgNewSize[4]={1};
	imglabsur.GetImageSize(ImgNewSize[0],ImgNewSize[1],ImgNewSize[2],ImgNewSize[3]);
	const short *sintData = imglabvol.GetImageData();	
	short *smaskData = new short[ImgNewSize[0]*ImgNewSize[1]*ImgNewSize[2]];
	//sintData[nz * ImgNewSize[1] * ImgNewSize[0] + ny * ImgNewSize[0] +nx];
	int np1[3]={54,111,55};
	int np2[3]={0,0,0};
	int npd=0;
	//TransDatanum2Greynum(np1,npd,ImgNewSize);
	//TransGreynum2Datanum(npd,np2,ImgNewSize);
	//short curpinten1=sintData[npd];
	//short curpinten2=imglabvol.GetPixelGreyscale(np1[0],np1[1],np1[2]);
	int vnum=imglabsur.GetNumberOfPixels();
	for (int i = 0; i < vnum; i++)
	{
		short curpinten1=sintData[i];
		if(curpinten1==nlabnum)
			smaskData[i]=1;
		else
			smaskData[i]=0;
	}

	cout<<"initialization: success!"<<endl;
	Generate_sur_from_vol(smaskData,nlabnum,imglabsur);
	zxh::SaveImage(&imglabsur,chResultName);
	cout << "Generate surface: successfully" << endl;


}

