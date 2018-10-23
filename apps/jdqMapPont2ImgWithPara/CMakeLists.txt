#
# Add the executable
#

ADD_EXECUTABLE(jdqMapPont2ImgWithPara jdqMapPont2ImgWithPara.cpp)
TARGET_LINK_LIBRARIES(jdqMapPont2ImgWithPara zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqMapPont2ImgWithPara)
