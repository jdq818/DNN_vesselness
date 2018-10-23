#
# Add the executable
#

ADD_EXECUTABLE(jdqDilPoiWR jdqDilPoiWR.cpp)
TARGET_LINK_LIBRARIES(jdqDilPoiWR zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqDilPoiWR)
