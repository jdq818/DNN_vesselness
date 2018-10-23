#
# Add the executable
#

ADD_EXECUTABLE(jdqresampleCurve jdqresampleCurve.cpp)
TARGET_LINK_LIBRARIES(jdqresampleCurve zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqresampleCurve)
