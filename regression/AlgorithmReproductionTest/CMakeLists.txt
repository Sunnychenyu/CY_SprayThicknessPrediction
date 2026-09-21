cmake_minimum_required(VERSION 3.20)

add_executable(AlgorithmReproductionTest main.cpp)
target_link_libraries(
    AlgorithmReproductionTest PRIVATE SMRobotSpray::SprayThicknessPrediction)
target_compile_features(AlgorithmReproductionTest PRIVATE cxx_std_17)
set_target_properties(
    AlgorithmReproductionTest PROPERTIES FOLDER tests/SMRobotSpray)
add_test(
    NAME AlgorithmReproductionTest
    COMMAND $<TARGET_FILE:AlgorithmReproductionTest>)
