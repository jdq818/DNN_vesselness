#
# Add the executable
#

ADD_EXECUTABLE(jdqGenDirectionOfCur jdqGenDirectionOfCur.cpp)
TARGET_LINK_LIBRARIES(jdqGenDirectionOfCur zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqGenDirectionOfCur)
