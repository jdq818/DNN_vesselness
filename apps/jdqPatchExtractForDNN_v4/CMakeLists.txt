#
# Add the executable
#

ADD_EXECUTABLE(jdqPatchExtractForDNN_v4 jdqPatchExtractForDNN_v4.cpp PatchExtractByWorldCoordinate.cpp PatchExtractByWorldCoordinate.h)
TARGET_LINK_LIBRARIES(jdqPatchExtractForDNN_v4 zxhImageModelings zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqPatchExtractForDNN_v4)
