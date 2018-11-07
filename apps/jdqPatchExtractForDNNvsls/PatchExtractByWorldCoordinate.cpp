#include "PatchExtractByWorldCoordinate.h"
//

bool PatchExtractByWorldCoordinate::GeneratePatchExtractByWorldWithRandOffset(int PatchNumIdex,float fdivec[3][3],float InputOutputWorldCoord[], int PatchSize[], zxhImageData* ImageArray[])
{

	zxhImageData * pIntensityImg = ImageArray[0];
	zxhImageData * pLabelImg =  ImageArray[1];
	zxhImageData * pPatchImage = ImageArray[2];

	int N = PatchSize[0];
	int HalfPatchLength = PatchSize[1]; 



	int PatchPointWorldCoord[3] = { 0 };
	//Computing the intensity of points in the tangent plane using interpolation	
	//Interpolation
	zxhImageModelingLinear InterpolationMod;
	InterpolationMod.SetImage(pIntensityImg);
	for (int i = -N; i < N + 1; i++)
	{
		for (int j = -N; j < N + 1; j++)
		{
			for (int k = -HalfPatchLength; k < HalfPatchLength + 1; k++)//extend the patch along the gradient orientation
			{ 
				PatchPointWorldCoord[0] = InputOutputWorldCoord[0] + i*fdivec[0][0] + j*fdivec[1][0] + k*fdivec[2][0];
				PatchPointWorldCoord[1] = InputOutputWorldCoord[1] + i*fdivec[0][1] + j*fdivec[1][1] + k*fdivec[2][1];
				PatchPointWorldCoord[2] =InputOutputWorldCoord[2]+  i*fdivec[0][2] + j*fdivec[1][2] + k*fdivec[2][2];

				//InterpolationMod.SetImage(pLALabel);//check the orientation of gradient				
				float IntensityValue = InterpolationMod.GetPixelFloatValueWithCheckByWorld(PatchPointWorldCoord[0], PatchPointWorldCoord[1], PatchPointWorldCoord[2], 0);

				int W2I_Coor_SetPatchPixel = (i + N) + (j + N) * (2 * N + 1) + (k + HalfPatchLength)* (2 * N + 1)* (2 * N + 1);
				pPatchImage->SetPixelByGreyscale(W2I_Coor_SetPatchPixel, PatchNumIdex+1, 1, 0, IntensityValue);

			}
		}
	}

	return true ;

}


PatchExtractByWorldCoordinate::~PatchExtractByWorldCoordinate()
{

}
