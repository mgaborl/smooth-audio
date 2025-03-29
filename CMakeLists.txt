cmake_minimum_required(VERSION 3.19)
project(smooth-audio)

find_library(sndfile libsndfile)
set(CMAKE_CXX_STANDARD 20)

add_executable(smooth-audio main.cpp)

target_link_libraries(
        smooth-audio 
        PRIVATE
                sndfile
)
install(TARGETS smooth-audio DESTINATION bin)

set(CPACK_GENERATOR "DEB")
set(CPACK_PACKAGE_NAME "smooth-audio")
set(CPACK_PACKAGE_VERSION "1.0.0")
set(CPACK_PACKAGE_RELEASE "1")
set(CPACK_PACKAGE_CONTACT "mgaborl@mgaborl.com")
set(CPACK_DEBIAN_PACKAGE_MAINTAINER "mgaborl@mgaborl.com")
set(CPACK_DEBIAN_PACKAGE_DEPENDS "libsndfile1")
set(CPACK_DEBIAN_FILE_NAME DEB-DEFAULT)
include(CPack)
