cmake_minimum_required(VERSION 3.19)
project(smooth-audio)

find_library(sndfile libsndfile)

add_executable(smooth-audio main.cpp)

target_link_libraries(
	smooth-audio 
	PRIVATE
		sndfile
)
include(CPack)
