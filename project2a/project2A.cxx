/*
* Project 2A
* 
* It looks like my output matches what is in the PDF. It might have been nice to have a larger image or even a separate screenshot for easier comparison.
* 
* Thomas Mitchell
*/

#include <GL/glew.h>    // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library
#include <stdio.h>
#include <stdlib.h>
#include <iostream>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "proj2_data.h"

#define PHASE3
#define PHASE4
#define PHASE5

void _print_shader_info_log(GLuint shader_index) 
{
  int max_length = 2048;
  int actual_length = 0;
  char shader_log[2048];
  glGetShaderInfoLog(shader_index, max_length, &actual_length, shader_log);
  printf("shader info log for GL index %u:\n%s\n", shader_index, shader_log);
}

GLuint SetupPhase2DataForRendering()
{
  printf("Getting data for Phase 2\n");

  float points[] = {0.5f, 0.0f, 0.0f,
                    0.0f, 0.0f, 0.0f,
                    0.0f, 0.5f, 0.0f,
                   -0.5f, 0.0f, 0.0f};

  float colors[] = {1.0f, 0.0f, 0.0f,
                    0.0f, 1.0f, 0.0f,
                    0.0f, 0.0f, 1.0f,
                    1.0f, 0.0f, 0.0f};

  GLuint indices[] = {0 ,1, 2,
                      1, 2, 3};

  GLuint points_vbo = 0;
  glGenBuffers(1, &points_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(points), points, GL_STATIC_DRAW);

  GLuint colors_vbo = 0;
  glGenBuffers(1, &colors_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(colors), colors, GL_STATIC_DRAW);

  GLuint index_vbo;    // Index buffer object
  glGenBuffers( 1, &index_vbo);
  glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, index_vbo );
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(indices), indices, GL_STATIC_DRAW);

  GLuint vao = 0;
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ARRAY_BUFFER, colors_vbo);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer( GL_ELEMENT_ARRAY_BUFFER, index_vbo );

  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);

  return vao;
}

const char *phase2VertexShader =
  "#version 400\n"
  "layout (location = 0) in vec3 vertex_position;\n"
  "layout (location = 1) in vec3 vertex_color;\n"
  "out vec3 color;\n"
  "void main() {\n"
  "  color = vertex_color;\n"
  "  gl_Position = vec4(vertex_position, 1.0);\n"
  "}\n";

const char *phase2FragmentShader =
  "#version 400\n"
  "in vec3 color;\n"
  "out vec4 frag_color;\n"
  "void main() {\n"
  "  frag_color = vec4(color, 1.0);\n"
  "}\n";


GLuint SetupPhase345DataForRendering()
{
  printf("Getting data for Phase 3\n");

  GLuint points_vbo = 0;
  glGenBuffers(1, &points_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(tri_points), tri_points, GL_STATIC_DRAW);

  GLuint data_vbo = 0;
  glGenBuffers(1, &data_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, data_vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(tri_data), tri_data, GL_STATIC_DRAW);

  GLuint normals_vbo = 0;
  glGenBuffers(1, &normals_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
  glBufferData(GL_ARRAY_BUFFER, sizeof(tri_normals), tri_normals, GL_STATIC_DRAW);

  GLuint index_vbo = 0;
  glGenBuffers(1, &index_vbo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_vbo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(tri_indices), tri_indices, GL_STATIC_DRAW);

  GLuint vao = 0;
  glGenVertexArrays(1, &vao);
  glBindVertexArray(vao);
  glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ARRAY_BUFFER, data_vbo);
  glVertexAttribPointer(1, 1, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
  glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_vbo);
  
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
  glEnableVertexAttribArray(2);

  return vao;
}

const char *phase345VertexShader =
  "#version 400\n"
  "layout (location = 0) in vec3 vertex_position;\n"
  "layout (location = 1) in float vertex_data;\n"
  "layout (location = 2) in vec3 vertex_normal;\n"
  "uniform mat4 MVP;\n"
  "uniform vec3 cameraloc;  // Camera position \n"
  "uniform vec3 lightdir;   // Lighting direction \n"
  "uniform vec4 lightcoeff; // Lighting coeff, Ka, Kd, Ks, alpha\n"
  "out float data;\n"
  "out float shading_amount;\n"
  "void main() {\n"
  "  vec4 position = vec4(vertex_position, 1.0);\n"
  "  gl_Position = MVP*position;\n"
  "  data = vertex_data;\n"

#ifdef PHASE5
  // Assign shading_amount a value by calculating phong shading
  // camaraloc  : is the location of the camera
  // lightdir   : is the direction of the light
  // lightcoeff : represents a vec4(Ka, Kd, Ks, alpha) from LightingParams of 1F

    "float diffuse = max(0.f, (lightdir[0] * vertex_normal[0] + lightdir[1] * vertex_normal[1] + lightdir[2] * vertex_normal[2])) * lightcoeff[1];\n"

    "float dotProdLightNormal =  lightdir[0] * vertex_normal[0] + lightdir[1] * vertex_normal[1] + lightdir[2] * vertex_normal[2];\n"
    "float viewDirection[3];\n"
    "viewDirection[0] = cameraloc[0] - vertex_position[0];\n"
    "viewDirection[1] = cameraloc[1] - vertex_position[1];\n"
    "viewDirection[2] = cameraloc[2] - vertex_position[2];\n"

    "float magnitude = sqrt(pow(viewDirection[0], 2) + pow(viewDirection[1], 2) + pow(viewDirection[2], 2));\n"

    "for (int i = 0; i < 3; i++)\n"
    "  viewDirection[i] /= magnitude; \n"

    "vec3 R = vec3(2.f * dotProdLightNormal * vertex_normal[0] - lightdir[0],\n"
    "  2.f * dotProdLightNormal * vertex_normal[1] - lightdir[1], \n"
    "  2.f * dotProdLightNormal * vertex_normal[2] - lightdir[2]);\n"
    "magnitude = sqrt(pow(R[0], 2) + pow(R[1], 2) + pow(R[2], 2));\n"

    "for (int i = 0; i < 3; i++)\n"
    "  R[i] /= magnitude; \n"

    "float specular =  pow(max(0.f, viewDirection[0] * R[0] + viewDirection[1] * R[1] + viewDirection[2] * R[2]), lightcoeff[3]) * lightcoeff[2];\n"

    "shading_amount = lightcoeff[0] + diffuse + specular;\n"
#endif

  "}\n";

const char *phase345FragmentShader =
  "#version 400\n"
  "in float data;\n"
  "in float shading_amount;\n"
  "out vec4 frag_color;\n"
  "void main() {\n"
  "  frag_color = vec4(1.0, 0.0, 0.0, 1.0);\n"

#ifdef PHASE4
  // Update frag_color by color based on data
    "if (1.0 <= data && data < 4.5)\n"
    "    frag_color = vec4(.25 + .75 * (data - 1) / 3.5, .25 + .75 * (data - 1) / 3.5, 1.0, 1.0);\n"
    "else\n"
    "    frag_color = vec4(1.0, 1.0 - .75 * (data - 4.5) / 1.5, 1.0 - .75 * (data - 4.5) / 1.5, 1.0);\n"
    

#endif

#ifdef PHASE5
  // Update frag_color by mixing the shading factor
#endif

    "frag_color[0] = min(1.f, frag_color[0] * shading_amount);\n"
    "frag_color[1] = min(1.f, frag_color[1] * shading_amount);\n"
    "frag_color[2] = min(1.f, frag_color[2] * shading_amount);\n"

  "}\n";

int main() {
    
  // start GL context and O/S window using the GLFW helper library
  if (!glfwInit()) {
    fprintf(stderr, "ERROR: could not start GLFW3\n");
    return 1;
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  GLFWwindow *window = glfwCreateWindow(700, 700, "CIS 441", NULL, NULL);
  if (!window) {
    fprintf(stderr, "ERROR: could not open window with GLFW3\n");
    glfwTerminate();
    return 1;
  }
  glfwMakeContextCurrent(window);
  // start GLEW extension handler
  glewExperimental = GL_TRUE;
  glewInit();

  // get version info
  const GLubyte *renderer = glGetString(GL_RENDERER); // get renderer string
  const GLubyte *version = glGetString(GL_VERSION);   // version as a string
  printf("Renderer: %s\n", renderer);
  printf("OpenGL version supported %s\n", version);

  // tell GL to only draw onto a pixel if the shape is closer to the viewer
  glEnable(GL_DEPTH_TEST); // enable depth-testing
  glDepthFunc(GL_LESS); // depth-testing interprets a smaller value as "closer"

  GLuint vao = 0;

#ifndef PHASE3
  vao = SetupPhase2DataForRendering();
  const char* vertex_shader = phase2VertexShader;
  const char* fragment_shader = phase2FragmentShader;
#else
  vao = SetupPhase345DataForRendering();
  const char* vertex_shader = phase345VertexShader;
  const char* fragment_shader = phase345FragmentShader;
#endif

  GLuint vs = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vs, 1, &vertex_shader, NULL);
  glCompileShader(vs);
  int params = -1;
  glGetShaderiv(vs, GL_COMPILE_STATUS, &params);
  if (GL_TRUE != params) {
    fprintf(stderr, "ERROR: GL shader index %i did not compile\n", vs);
    _print_shader_info_log(vs);
    exit(EXIT_FAILURE);
  }

  GLuint fs = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fs, 1, &fragment_shader, NULL);
  glCompileShader(fs);
  glGetShaderiv(fs, GL_COMPILE_STATUS, &params);
  if (GL_TRUE != params) {
    fprintf(stderr, "ERROR: GL shader index %i did not compile\n", fs);
    _print_shader_info_log(fs);
    exit(EXIT_FAILURE);
  }

  GLuint shader_programme = glCreateProgram();
  glAttachShader(shader_programme, fs);
  glAttachShader(shader_programme, vs);
  glLinkProgram(shader_programme);

  glUseProgram(shader_programme);


#ifdef PHASE3  // Code block for camera transforms
  // Projection matrix : 30?? Field of View
  // display size  : 1000x1000
  // display range : 5 unit <-> 200 units
  glm::mat4 Projection = glm::perspective(
      glm::radians(30.0f), (float)1000 / (float)1000,  5.0f, 200.0f);
  glm::vec3 camera(0, 40, 40);
  glm::vec3 origin(0, 0, 0);
  glm::vec3 up(0, 1, 0);
  // Camera matrix
  glm::mat4 View = glm::lookAt(
    camera, // Camera in world space
    origin, // looks at the origin
    up      // and the head is up
  );
  // Model matrix : an identity matrix (model will be at the origin)
  glm::mat4 Model = glm::mat4(1.0f);
  // Our ModelViewProjection : multiplication of our 3 matrices
  glm::mat4 mvp = Projection * View * Model;

  // Get a handle for our "MVP" uniform
  // Only during the initialisation
  GLuint mvploc = glGetUniformLocation(shader_programme, "MVP");
  // Send our transformation to the currently bound shader, in the "MVP" uniform
  // This is done in the main loop since each model will have a different MVP matrix
  // (At least for the M part)
  glUniformMatrix4fv(mvploc, 1, GL_FALSE, &mvp[0][0]);
#endif

#ifdef PHASE5 // Code block for shading parameters
  GLuint camloc = glGetUniformLocation(shader_programme, "cameraloc");
  glUniform3fv(camloc, 1, &camera[0]);
  glm::vec3 lightdir = glm::normalize(camera - origin);   // Direction of light
  GLuint ldirloc = glGetUniformLocation(shader_programme, "lightdir");
  glUniform3fv(ldirloc, 1, &lightdir[0]);
  glm::vec4 lightcoeff(0.3, 0.7, 2.8, 50.5); // Lighting coeff, Ka, Kd, Ks, alpha
  GLuint lcoeloc = glGetUniformLocation(shader_programme, "lightcoeff");
  glUniform4fv(lcoeloc, 1, &lightcoeff[0]);
#endif

  while (!glfwWindowShouldClose(window)) {
    // wipe the drawing surface clear
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glBindVertexArray(vao);
    // Draw triangles

#ifndef PHASE3
    glDrawElements( GL_TRIANGLES, 6, GL_UNSIGNED_INT, NULL );
#else
    // Add correct number of indices
    glDrawElements( GL_TRIANGLES, sizeof(tri_indices) / sizeof(tri_indices[0]), GL_UNSIGNED_INT, NULL );
#endif

    // update other events like input handling
    glfwPollEvents();
    // put the stuff we've been drawing onto the display
    glfwSwapBuffers(window);
  }

  // close GL context and any other GLFW resources
  glfwTerminate();
  
  return 0;
}
