add_executable (simulator simulator.cpp)
target_link_libraries (simulator ${PROJECT_NAME}-shared)
set_target_properties (simulator PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
