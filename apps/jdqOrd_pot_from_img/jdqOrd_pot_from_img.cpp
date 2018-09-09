

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
	int nx,ny,nz;	
	for (int i = 0; i < nPointNum; i++)
	{
		nx= vPBr[i].x+1;
		ny= vPBr[i].y+1;
		nz= vPBr[i].z+1;
		WriteFileTxt <<right<<fixed<<setfill('0')<<setprecision(4)<< nx<<' '<<ny<<' '<<nz<<'\n';

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
bool InitBra(BraTypeDef &BParPont,PointCordTypeDef initPont)
{
	BParPont.p=initPont;
	BParPont.lchild=initPont;
	BParPont.rchild=initPont;
	BParPont.parent=initPont;
	return true;
}
bool FindEndPinvParP(PointCordTypeDef PEndPonts,vector<BraTypeDef> vBraParPonts,int &nPosi)
{
	for(int i=0;i<vBraParPonts.size();i++)
	{
		if(PEndPonts.x==vBraParPonts[i].p.x&&PEndPonts.y==vBraParPonts[i].p.y&&PEndPonts.z==vBraParPonts[i].p.z)
		{
			nPosi=i;
			return true;
		}
	}
	return false;
}
bool FindEndPinvParP1(PointCordTypeDef PEndPonts,vector<PointCordTypeDef> vBraParPonts,int *flag)
{
	bool bfind=false;
	for(int i=0;i<vBraParPonts.size();i++)
	{
		if(PEndPonts.x==vBraParPonts[i].x&&PEndPonts.y==vBraParPonts[i].y&&PEndPonts.z==vBraParPonts[i].z)
		{

			flag[i]=1;
			bfind= true;
		}
	}
	return false;
}
float CalcDist(PointCordTypeDef PBPont,PointCordTypeDef PBPontNext,zxhImageDataT<short>&imgReadRaws)
{
	float nP[3]={PBPont.x,PBPont.y,PBPont.z};
	float nPN[3]={PBPontNext.x,PBPontNext.y,PBPontNext.z};
	imgReadRaws.GetImageInfo()->ImageToWorld(nP);
	imgReadRaws.GetImageInfo()->ImageToWorld(nPN);
	float fDist=zxh::VectorOP_Distance(nP,nPN,3);
	return fDist;
}
//if (DeteChlPonts(vPParPLR,vPChlPLR,imgReadRaws,consfMindist,PparPont,vPchlPont))//find the two children points
bool DeteChlPonts(vector<PointCordTypeDef>vPParPLR,vector<PointCordTypeDef> vPChlPLR,zxhImageDataT<short>& imgReadRaws,const float consfMindist,PointCordTypeDef &PparPont,vector<PointCordTypeDef>&vPChlPont)
{

	vector<PointCordTypeDef> vPChldPont,vPChldOrderPont;
	vPChldPont.clear();
	vPChldOrderPont.clear();
	PointCordTypeDef PpbP,PpbE,PpbNS;
	for(int i=0;i<vPParPLR.size();i++)
	{
		PpbP=vPParPLR[i];
		for(int j=0;j<vPChlPLR.size();j++)
		{
			PointCordTypeDef PCPont=vPChlPLR[j];
			float fDist=CalcDist(PpbP,PCPont,imgReadRaws);
			if(fDist<consfMindist)
			{
				vPChlPont.push_back(PCPont);
			}
		}
		if (vPChlPont.size()==2)
		{
			if(vPChlPont[0].b<vPChlPont[1].b)
			{
				vPChldOrderPont.push_back(vPChlPont[0]);
				vPChldOrderPont.push_back(vPChlPont[1]);
			}
			else
			{
				vPChldOrderPont.push_back(vPChlPont[1]);
				vPChldOrderPont.push_back(vPChlPont[0]);
			}
			vPChlPont.clear();
			vPChlPont.assign(vPChldOrderPont.begin(), vPChldOrderPont.end());  
			PparPont=PpbP;
			return true;
		}
		else
		{
			vPChlPont.clear();
		}
	}


	return false;
}
bool IsNormTri(PointCordTypeDef Pp,PointCordTypeDef PCp,PointCordTypeDef PCNextp,zxhImageDataT<short>& imgReadRaws)
{
	float fPp[3]={Pp.x,Pp.y,Pp.z};
	float fPCp[3]={PCp.x,PCp.y,PCp.z};
	float fPCNextp[3]={PCNextp.x,PCNextp.y,PCNextp.z};
	imgReadRaws.GetImageInfo()->ImageToWorld(fPp);
	imgReadRaws.GetImageInfo()->ImageToWorld(fPCp);
	imgReadRaws.GetImageInfo()->ImageToWorld(fPCNextp);
	float fVec1[3]={fPp[0]-fPCp[0],fPp[1]-fPCp[1],fPp[2]-fPCp[2]};
	float fVec2[3]={fPp[0]-fPCNextp[0],fPp[1]-fPCNextp[1],fPp[2]-fPCNextp[2]};
	float fVec3[3]={fPCp[0]-fPCNextp[0],fPCp[1]-fPCNextp[1],fPCp[2]-fPCNextp[2]};
	float fcosDing=zxh::VectorOP_Cosine(fVec1,fVec2,3);
	float fcosDi1=zxh::VectorOP_Cosine(fVec1,fVec3,3);
	float fcosDi2=zxh::VectorOP_Cosine(fVec2,fVec3,3);
	float fNormDing=zxh::VectorOP_Magnitude(fVec1,3);
	float fNormDi1=zxh::VectorOP_Magnitude(fVec2,3);
	float fNormDi2=zxh::VectorOP_Magnitude(fVec3,3);
	if (abs(fcosDing)<=0.960&&abs(fcosDi1)<=0.960&&abs(fcosDi2)<=0.960&&fNormDing<6&&fNormDi1<6&&fNormDi2<6)
	{
		return true;
	}
	else
		return false;
}
bool DeteChlPonts_Tri(vector<PointCordTypeDef>vPParPLR,vector<PointCordTypeDef> vPChlPLR,zxhImageDataT<short>& imgReadRaws,const float consfMindist,PointCordTypeDef &PparPont,vector<PointCordTypeDef>&vPChlPont)
{

	vector<PointCordTypeDef> vPChldPont,vPChldOrderPont;
	vPChldPont.clear();
	vPChldOrderPont.clear();
	PointCordTypeDef PPbP,PpbE,PpbNS;
	for(int i=0;i<vPParPLR.size();i++)
	{
		PointCordTypeDef PpbP=vPParPLR[i];
		for(int j=0;j<vPChlPLR.size()-1;j++)
		{
			PointCordTypeDef PCPont=vPChlPLR[j];
			PointCordTypeDef PCNextPont=vPChlPLR[j+1];

			if(IsNormTri(PpbP,PCPont,PCNextPont,imgReadRaws))
			{
				PPbP=PpbP;
				vPChlPont.push_back(PCPont);
				vPChlPont.push_back(PCNextPont);
			}
		}
		if (vPChlPont.size()==2)//order left and right
		{
			if(vPChlPont[0].b<vPChlPont[1].b)
			{
				vPChldOrderPont.push_back(vPChlPont[0]);
				vPChldOrderPont.push_back(vPChlPont[1]);
			}
			else
			{
				vPChldOrderPont.push_back(vPChlPont[1]);
				vPChldOrderPont.push_back(vPChlPont[0]);
			}
			vPChlPont.clear();
			vPChlPont.assign(vPChldOrderPont.begin(), vPChldOrderPont.end());  
			PparPont=PpbP;
			return true;
		}
		else
		{
			vPChlPont.clear();
		}
	}


	return false;
}
bool EraseOvPonts(vector<PointCordTypeDef>&vPParPLR )//remove the overlap points
{
	vector<PointCordTypeDef> vPParPLRCopy,vPParNew;
	vPParPLRCopy.assign(vPParPLR.begin(), vPParPLR.end());  
	int nPosi;
	int nnsize=vPParPLR.size();
	int *nflag=new int[nnsize];
	for(int i=0;i<nnsize;i++)
	{
		nflag[i] =0;
	}
	for(int i=0;i<vPParPLR.size();i++)
	{
		if (nflag[i]==1)
		{
			continue;
		}
		else
		{
			PointCordTypeDef PEndPonts=vPParPLR[i];
			vPParNew.push_back(PEndPonts);
			FindEndPinvParP1(PEndPonts,vPParPLRCopy,nflag);

		}
	}

	delete[] nflag;
	vPParPLR.clear();
	vPParPLR.assign(vPParNew.begin(), vPParNew.end());  
	return true;
}
bool InitPoint(PointCordTypeDef &Ppar)
{
	Ppar.b=0;
	Ppar.x=0;
	Ppar.y=0;
	Ppar.z=0;
	return true;
}
bool BackTrack(PointCordTypeDef PEndPont,vector<vector<PointCordTypeDef>>vvBraOrderPonts,vector<BraTypeDef>vBraParPonts,vector<BraTypeDef>vBraChlPonts,PointCordTypeDef PInitPont,vector<PointCordTypeDef> &vPBranch)
{
	int nBraNUM=PEndPont.b;
	int nBraNUMInvvBra=0;//find the right branch
	for(int nvvNUM=0;nvvNUM<vvBraOrderPonts.size();nvvNUM++)
	{
		vector<PointCordTypeDef> Pv=vvBraOrderPonts[nvvNUM];
		int nPvBra=Pv[0].b;
		if(nPvBra==nBraNUM)
		{
			nBraNUMInvvBra=nvvNUM;
		}
	}
	vector<PointCordTypeDef> PLastv=vvBraOrderPonts[nBraNUMInvvBra];

	PointCordTypeDef Ppar;
	InitPoint(Ppar);
	BraTypeDef BraCur,BraPar;
	InitBra(BraCur,PInitPont);
	InitBra(BraPar,PInitPont);
	for(int i=PLastv.size()-1;i>=0;i--)//find the children points in children points vector
	{
		int nPosi;
		PointCordTypeDef PCur=PLastv[i];
		if(FindEndPinvParP(PCur,vBraChlPonts,nPosi))
		{
			BraCur=vBraChlPonts[nPosi];
			vPBranch.push_back(PCur);
			break;
		}
		vPBranch.push_back(PCur);
	}
	if(BraCur.parent.b!=0)
	{
		Ppar=BraCur.parent;
		BackTrack(Ppar,vvBraOrderPonts,vBraParPonts,vBraChlPonts,PInitPont,vPBranch);
	}
	else
	{
		return false;
	}
}
bool TransToWorld(vector<PointCordTypeDef> vPBranch,vector<PointCordTypeDef> &vPBranchWorld,zxhImageDataT<short>&imgReadRaws)
{
	vPBranchWorld.clear();
	for(int i=vPBranch.size()-1;i>=0;i--)
	{
		float fp[3]={vPBranch[i].x,vPBranch[i].y,vPBranch[i].z};
		imgReadRaws.GetImageInfo()->ImageToWorld(fp);
		PointCordTypeDef pontWorld;
		pontWorld.x=fp[0];
		pontWorld.y=fp[1];
		pontWorld.z=fp[2];
		pontWorld.b=vPBranch[i].b;
		vPBranchWorld.push_back(pontWorld);
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

int main(int argc, char *argv[])
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

	string strFileNameRaw ="E:/work_jdq/CTA_data/RCAAEF/training_jdq/01/image01_R0_1.nii.gz";
	char *chResFilefold ="K:/JDQ/CCTA_CAR/RCAA_32/training/dataset00";


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

				for(int i=0;i<vPKPont.size();i++)//找出距离当前点最近的点
				{
					float nP[3]={vPKPont[i].x,vPKPont[i].y,vPKPont[i].z};
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
				
				PointCordTypeDef pmindstP;
				pmindstP=vPKPont[nMinIndex];
				pmindstP.b=fMinDist;
				//cout<<"current point"<<PpbE.x<<","<<PpbE.y<<","<<PpbE.z<<" To-"<<endl;
				cout<<pmindstP.x<<","<<pmindstP.y<<","<<pmindstP.z<<":"<<pmindstP.b<<endl;
				vOrderPonts.push_back(pmindstP);
				vPKPont.erase(vPKPont.begin()+nMinIndex);      
				PpbE=pmindstP;
			}

}




