#
# Add the executable
#

ADD_EXECUTABLE(jdqPatchExtractForDNNvsls_2.0 jdqPatchExtractForDNNvsls_2.0.cpp PatchExtractByWorldCoordinate.cpp PatchExtractByWorldCoordinate.h)
TARGET_LINK_LIBRARIES(jdqPatchExtractForDNNvsls_2.0 zxhImageModelings zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqPatchExtractForDNNvsls_2.0)
