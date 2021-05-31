cmake_minimum_required(VERSION 2.8.12.1)

PROJECT(project3C)

SET(OpenGL_GL_PREFERENCE LEGACY)


add_executable(project3C project3C.cxx)
if(APPLE)
  target_link_libraries(project3C ${OPENGL_gl_LIBRARY} GLEW::glew_s glfw)
else()
  target_link_libraries(project3C ${OPENGL_gl_LIBRARY} GLEW glfw)
endif()
