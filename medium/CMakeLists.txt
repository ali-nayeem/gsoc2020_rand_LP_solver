# VolEsti (volume computation and sampling library)
# Copyright (c) 20012-2018 Vissarion Fisikopoulos
# Copyright (c) 2018 Apostolos Chalkis
# Licensed under GNU LGPL.3, see LICENCE file

project( VolEsti )

enable_testing()

CMAKE_MINIMUM_REQUIRED(VERSION 2.4.5)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
       cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

#if(NOT DEFINED BOOST)
#  message(FATAL_ERROR "This program requires the boost library, and will not be compiled.")
#else()
#  message(STATUS "BOOST is DEFINED")
#endif()

find_library(LP_SOLVE NAMES liblpsolve55.so PATHS /usr/lib/lp_solve)

if (NOT LP_SOLVE)

  message(FATAL_ERROR "This program requires the lp_solve library, and will not be compiled.")

else ()

  message(STATUS "Library lp_solve found: ${LP_SOLVE}")

  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  #include_directories (BEFORE ../external/Eigen)
  include_directories (BEFORE ../external)
  include_directories (BEFORE ../external/minimum_ellipsoid)
  #include_directories (BEFORE ../include/cartesian_geom)
  #include_directories (BEFORE ../include/convex_bodies)
  include_directories (BEFORE ../external/LPsolve_src/run_headers)
  include_directories (BEFORE ../external/boost)
  #include_directories (BEFORE BOOST)
  include_directories (BEFORE ../include/generators)
  include_directories (BEFORE ../include/volume)
  include_directories (BEFORE ../include)
  include_directories (BEFORE ../include/convex_bodies)
  include_directories (BEFORE ../include/annealing)
  include_directories (BEFORE ../include/samplers)
  include_directories (BEFORE ../include/lp_oracles)
  include_directories (BEFORE ../include/misc)


  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard
  #add_definitions(${CMAKE_CXX_FLAGS} "-O3")  # optimization of the compiler
  add_definitions(${CMAKE_CXX_FLAGS} "-g")  # debug
  add_definitions(${CMAKE_CXX_FLAGS} "-gdwarf-3")  # debug
  #add_definitions(${CMAKE_CXX_FLAGS} "-fpermissive")  # debug
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgsl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lm")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-ldl")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DBOOST_NO_AUTO_PTR")
  add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-DVOLESTI_DEBUG")
  #add_definitions(${CXX_COVERAGE_COMPILE_FLAGS} "-lgslcblas")
  #add_definitions( "-O3 -lgsl -lm -ldl -lgslcblas" )

  add_executable (vol vol.cpp)
  #add_executable (volume volume_example.cpp)
  add_executable (generate generator.cpp)
  add_executable (lp_prob lp_problem.cpp)


  add_library(test_main OBJECT test_main.cpp)

  add_executable (volume_test volume_test.cpp $<TARGET_OBJECTS:test_main>)
  add_executable (cheb_test chebychev_test.cpp $<TARGET_OBJECTS:test_main>)
  #add_executable (rounding_test rounding_test.cpp $<TARGET_OBJECTS:test_main>)
  add_executable (volumeCG_test volumeCG_test.cpp $<TARGET_OBJECTS:test_main>)
  add_executable (VpolyCG_test VpolyCG_test.cpp $<TARGET_OBJECTS:test_main>)
  add_executable (VpolyVol_test VpolyVol_test.cpp $<TARGET_OBJECTS:test_main>)
  add_executable (ZonotopeVol_test ZonotopeVol_test.cpp $<TARGET_OBJECTS:test_main>)
  add_executable (cool_bodies_bill_test cooling_bodies_bill_test.cpp $<TARGET_OBJECTS:test_main>)
  #add_executable (ZonotopeVolCG_test ZonotopeVolCG_test.cpp $<TARGET_OBJECTS:test_main>)
  
  add_test(NAME volume_cube COMMAND volume_test -tc=cube)
  add_test(NAME volume_cross COMMAND volume_test -tc=cross)
  add_test(NAME volume_birkhoff COMMAND volume_test -tc=birk)
  #add_test(NAME volume_prod_simplex COMMAND volume_test -tc=prod_simplex)
  add_test(NAME volume_simplex COMMAND volume_test -tc=simplex)
  add_test(NAME volume_skinny_cube COMMAND volume_test -tc=skinny_cube)

  add_test(NAME volumeCG_cube COMMAND volumeCG_test -tc=cube)
  add_test(NAME volumeCG_cross COMMAND volumeCG_test -tc=cross)
  add_test(NAME volumeCG_birkhoff COMMAND volumeCG_test -tc=birk)
  add_test(NAME volumeCG_prod_simplex COMMAND volumeCG_test -tc=prod_simplex)
  add_test(NAME volumeCG_simplex COMMAND volumeCG_test -tc=simplex)

  #add_test(NAME VpolyCG_cube COMMAND VpolyCG_test -tc=cube)
  add_test(NAME VpolyCG_cross COMMAND VpolyCG_test -tc=cross)
  add_test(NAME VpolyCG_simplex COMMAND VpolyCG_test -tc=simplex)
  
  add_test(NAME VpolyVol_cube COMMAND VpolyVol_test -tc=cube)
  add_test(NAME VpolyVol_cross COMMAND VpolyVol_test -tc=cross)
  add_test(NAME VpolyVol_simplex COMMAND VpolyVol_test -tc=simplex)

  add_test(NAME ZonotopeVol4 COMMAND ZonotopeVol_test -tc=4_dimensional)

  #add_test(NAME ZonotopeVolCG4 COMMAND ZonotopeVolCG_test -tc=4_dimensional)
  #add_test(NAME ZonotopeVolCG5 COMMAND ZonotopeVolCG_test -tc=5_dimensional)

  add_test(NAME cheb_cube COMMAND cheb_test -tc=cheb_cube)
  add_test(NAME cheb_cross COMMAND cheb_test -tc=cheb_cross)
  add_test(NAME cheb_birkhoff COMMAND cheb_test -tc=cheb_birkhoff)
  add_test(NAME cheb_prod_simplex COMMAND cheb_test -tc=cheb_prod_simplex)
  add_test(NAME cheb_simplex COMMAND cheb_test -tc=cheb_simplex)
  add_test(NAME cheb_skinny_cube COMMAND cheb_test -tc=cheb_skinny_cube)

  add_test(NAME cool_bodies_cube COMMAND cool_bodies_bill_test -tc=cube)
  add_test(NAME cool_bodies_cross COMMAND cool_bodies_bill_test -tc=cross)
  add_test(NAME cool_bodies_birkhoff COMMAND cool_bodies_bill_test -tc=birk)
  add_test(NAME cool_bodies_prod_simplex COMMAND cool_bodies_bill_test -tc=prod_simplex)
  add_test(NAME cool_bodies_simplex COMMAND cool_bodies_bill_test -tc=simplex)
  add_test(NAME cool_bodies_skinny_cube COMMAND cool_bodies_bill_test -tc=skinny_cube)

  #add_test(NAME round_skinny_cube COMMAND rounding_test -tc=round_skinny_cube)
  #add_test(NAME round_rot_skinny_cube COMMAND rounding_test -tc=round_rot_skinny_cube)

  TARGET_LINK_LIBRARIES(vol ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(volume ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(generate ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(lp_prob ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volume_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(cheb_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(rounding_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(volumeCG_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(VpolyCG_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(VpolyVol_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(ZonotopeVol_test ${LP_SOLVE})
  TARGET_LINK_LIBRARIES(cool_bodies_bill_test ${LP_SOLVE})
  #TARGET_LINK_LIBRARIES(ZonotopeVolCG_test ${LP_SOLVE})

endif()


