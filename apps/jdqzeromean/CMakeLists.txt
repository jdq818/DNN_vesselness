#
# Add the executable
#

ADD_EXECUTABLE(jdqzeromean jdqzeromean.cpp)
TARGET_LINK_LIBRARIES(jdqzeromean zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqzeromean)
