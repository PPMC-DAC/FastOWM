# add CXX the executable
# add_executable( cesgaOWM cesgaOWM.c )

# target_link_libraries( cesgaOWM PUBLIC src )
# target_include_directories( cesgaOWM PUBLIC "${PRJ_INC}")

# add CXX the executable
# add_executable( owm_tfm owm_tfm.cpp )

# target_link_libraries( owm_tfm PUBLIC src )
# target_include_directories( owm_tfm PUBLIC "${PRJ_INC}")

add_executable( owm_gpu owm_gpu.cpp )

target_link_libraries( owm_gpu PUBLIC src )
target_include_directories( owm_gpu PUBLIC "${PRJ_INC}")

# add CXX the executable
# add_executable( prueba prueba.cpp )
#
# target_link_libraries( prueba PUBLIC src )
# target_include_directories( prueba PUBLIC "${PRJ_INC}")


# file( GLOB APP_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/cesgaOWM.c
#                         ${CMAKE_CURRENT_SOURCE_DIR}/owm_tfm.cpp )

# foreach( sourcefile ${APP_SOURCES} )

#     file(RELATIVE_PATH filename ${CMAKE_CURRENT_SOURCE_DIR} ${sourcefile})

#     string( REPLACE ".cpp" "" file ${filename} )

#     # add CXX the executable
#     add_executable( ${file} ${sourcefile} )

#     target_link_libraries( ${file} PUBLIC src )
#     # target_include_directories( ${file} PUBLIC "${PRJ_INC}")

# endforeach( sourcefile ${APP_SOURCES} )
