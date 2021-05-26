cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project3B)

SET(OpenGL_GL_PREFERENCE LEGACY)

message("GL include dir: ${OPENGL_INCLUDE_DIR}")
message("GL link libraries: ${OPENGL_gl_LIBRARY}")

add_executable(project3B project3B.cxx)
if(APPLE)
  target_link_libraries(project3B ${OPENGL_gl_LIBRARY} GLEW::glew_s glfw)
else()
  target_link_libraries(project3B ${OPENGL_gl_LIBRARY} GLEW glfw)
endif()
