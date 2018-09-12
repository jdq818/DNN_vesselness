#
# Add the executable
#

ADD_EXECUTABLE(jdqMapPont2Img jdqMapPont2Img.cpp)
TARGET_LINK_LIBRARIES(jdqMapPont2Img zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqMapPont2Img)
