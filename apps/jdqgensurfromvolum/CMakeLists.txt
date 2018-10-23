#
# Add the executable
#

ADD_EXECUTABLE(jdqgensurfromvolum jdqgensurfromvolum.cpp)
TARGET_LINK_LIBRARIES(jdqgensurfromvolum zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqgensurfromvolum)
