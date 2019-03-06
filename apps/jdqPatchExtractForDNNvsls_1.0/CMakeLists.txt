#
# Add the executable
#

ADD_EXECUTABLE(jdqPatchExtractForDNNvsls_1.0 jdqPatchExtractForDNNvsls_1.0.cpp PatchExtractByWorldCoordinate.cpp PatchExtractByWorldCoordinate.h)
TARGET_LINK_LIBRARIES(jdqPatchExtractForDNNvsls_1.0 zxhImageModelings zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqPatchExtractForDNNvsls_1.0)
