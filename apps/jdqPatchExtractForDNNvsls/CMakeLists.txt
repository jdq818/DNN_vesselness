#
# Add the executable
#

ADD_EXECUTABLE(jdqPatchExtractForDNNvsls jdqPatchExtractForDNNvsls.cpp PatchExtractByWorldCoordinate.cpp PatchExtractByWorldCoordinate.h)
TARGET_LINK_LIBRARIES(jdqPatchExtractForDNNvsls zxhImageModelings zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqPatchExtractForDNNvsls)
