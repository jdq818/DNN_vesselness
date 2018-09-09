#
# Add the executable
#

ADD_EXECUTABLE(jdqBesselCurve jdqBesselCurve.cpp)
TARGET_LINK_LIBRARIES(jdqBesselCurve zxhImages zxhCommon znz zlib niftiio vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqBesselCurve)
