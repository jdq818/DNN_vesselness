#
# Add the executable
#

ADD_EXECUTABLE(jdqPatchExtractForDNN_v3 jdqPatchExtractForDNN_v3.cpp PatchExtractByWorldCoordinate.cpp PatchExtractByWorldCoordinate.h)
TARGET_LINK_LIBRARIES(jdqPatchExtractForDNN_v3 zxhImageModelings zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqPatchExtractForDNN_v3)
