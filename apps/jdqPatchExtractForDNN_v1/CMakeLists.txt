#
# Add the executable
#

ADD_EXECUTABLE(jdqPatchExtractForDNN_v1 jdqPatchExtractForDNN_v1.cpp PatchExtractByWorldCoordinate.cpp PatchExtractByWorldCoordinate.h)
TARGET_LINK_LIBRARIES(jdqPatchExtractForDNN_v1 zxhImageModelings zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqPatchExtractForDNN_v1)
