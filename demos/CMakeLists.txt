option(OMPL_BUILD_DEMOS "Build OMPL demos" ON)

if (OMPL_BUILD_DEMOS)

  macro(add_ompl_demo demo_name)
    add_executable(${ARGV})
    target_link_libraries(${demo_name}
      ompl
      ${Boost_FILESYSTEM_LIBRARY}
      ${Boost_SYSTEM_LIBRARY}
      ${Boost_THREAD_LIBRARY}
      ${Boost_DATE_TIME_LIBRARY}
      ${Boost_PROGRAM_OPTIONS_LIBRARY})
  endmacro(add_ompl_demo)

  add_ompl_demo(demo_FMT3 FMTtest.cpp)

endif(OMPL_BUILD_DEMOS)

file(GLOB OMPL_DEMO_CXX_FILES "*.cpp")
file(GLOB OMPL_DEMO_PY_FILES "*.py")
install(FILES ${OMPL_DEMO_CXX_FILES} ${OMPL_DEMO_PY_FILES}
    DESTINATION "${OMPL_DEMO_INSTALL_DIR}"
    COMPONENT ompl)
install(DIRECTORY Koules
    DESTINATION "${OMPL_DEMO_INSTALL_DIR}"
    COMPONENT ompl)
