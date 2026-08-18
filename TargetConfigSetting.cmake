##########################################################################
# CMake requirement
##########################################################################
cmake_minimum_required( VERSION 3.20 )

set( ${TARGET_NAME}_RequiredLibsPublic
    Eigen3::Eigen
    SMRobotSpray::WorkpieceCore
    SMRobotSpray::SprayCore
    SMRobotSpray::SprayTrajectoryCore
)

set( ${TARGET_NAME}_RequiredLibsPrivate
)
