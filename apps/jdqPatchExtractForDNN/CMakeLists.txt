#
# Add the executable
#

ADD_EXECUTABLE(jdqPatchExtractForDNN jdqPatchExtractForDNN.cpp PatchExtractByWorldCoordinate.cpp PatchExtractByWorldCoordinate.h)
TARGET_LINK_LIBRARIES(jdqPatchExtractForDNN zxhImageModelings zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqPatchExtractForDNN)
