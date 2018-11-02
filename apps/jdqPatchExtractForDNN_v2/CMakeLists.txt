#
# Add the executable
#

ADD_EXECUTABLE(jdqPatchExtractForDNN_v2 jdqPatchExtractForDNN_v2.cpp PatchExtractByWorldCoordinate.cpp PatchExtractByWorldCoordinate.h)
TARGET_LINK_LIBRARIES(jdqPatchExtractForDNN_v2 zxhImageModelings zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqPatchExtractForDNN_v2)
