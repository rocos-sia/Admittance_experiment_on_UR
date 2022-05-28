
cmake_minimum_required(VERSION 3.10)

project(ur_examples)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_BUILD_TYPE Release)

find_package(catkin REQUIRED COMPONENTS
  roscpp
)

find_package(ur_rtde REQUIRED)

catkin_package(
   CATKIN_DEPENDS   roscpp
)

INCLUDE_DIRECTORIES(3_party_library/kdl/include /usr/include/eigen3   3_party_library/track_Ik/include   include  ${catkin_INCLUDE_DIRS})
LINK_DIRECTORIES(3_party_library/track_Ik/lib  3_party_library/kdl/lib    )



add_executable(admittance_Control src/UR_control/admittance_Control/admittance_Control.cpp   )
target_link_libraries(admittance_Control  PUBLIC ur_rtde::rtde    kdl-1-5  orocos-kdl-models   trac_ik   nlopt  ${catkin_LIBRARIES}  )



add_executable(moveL_admittance_Control src/UR_control/admittance_Control/moveL_admittance_Control.cpp   3_party_library/velocity_profile/interpolate.cc )
target_link_libraries(moveL_admittance_Control  PUBLIC ur_rtde::rtde    kdl-1-5  orocos-kdl-models   trac_ik   nlopt  ${catkin_LIBRARIES}  )

