#
# Add the executable
#

ADD_EXECUTABLE(jdqOrd_pot_from_img jdqOrd_pot_from_img.cpp)
TARGET_LINK_LIBRARIES(jdqOrd_pot_from_img zxhcaeMinPathLib zxhImages zxhCommon znz zlib niftiio  vtkRendering vtkIO)

INSTALL_TARGETS(/bin jdqOrd_pot_from_img)
