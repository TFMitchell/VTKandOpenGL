/*
* Project 3B
* 
* Although the output is not consistent with that of the specs on Windows, I have been testing concurrently on the provided virtual Linux box, which has the correct output for each of the steps.
* 
* Thomas Mitchell
* 
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>

using std::endl;
using std::cerr;

#include <GL/glew.h>    // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper library

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/vec3.hpp>   // glm::vec3
#include <glm/vec4.hpp>   // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>  // glm::translate, glm::rotate, glm::scale

class RenderManager;

const char *GetVertexShader();
const char *GetFragmentShader();

class Triangle
{
  public:
    glm::vec3 v0;
    glm::vec3 v1;
    glm::vec3 v2;
};

std::vector<Triangle> SplitTriangle(std::vector<Triangle> &list)
{
    std::vector<Triangle> output(4*list.size());
    output.resize(4*list.size());
    for (unsigned int i = 0 ; i < list.size() ; i++)
    {
        Triangle t = list[i];
        glm::vec3 vmid1, vmid2, vmid3;
        vmid1 = (t.v0 + t.v1) / 2.0f;
        vmid2 = (t.v1 + t.v2) / 2.0f;
        vmid3 = (t.v0 + t.v2) / 2.0f;
        output[4*i+0].v0 = t.v0;
        output[4*i+0].v1 = vmid1;
        output[4*i+0].v2 = vmid3;
        output[4*i+1].v0 = t.v1;
        output[4*i+1].v1 = vmid2;
        output[4*i+1].v2 = vmid1;
        output[4*i+2].v0 = t.v2;
        output[4*i+2].v1 = vmid3;
        output[4*i+2].v2 = vmid2;
        output[4*i+3].v0 = vmid1;
        output[4*i+3].v1 = vmid2;
        output[4*i+3].v2 = vmid3;
    }
    return output;
}

void PushVertex(std::vector<float>& coords,
                const glm::vec3& v)
{
  coords.push_back(v.x);
  coords.push_back(v.y);
  coords.push_back(v.z);
}

// Sets up a sphere with equation x^2+y^2+z^2=1
void GetSphereData(std::vector<float>& coords, std::vector<float>& normals)
{
  int recursionLevel = 3;
  std::vector<Triangle> list;
  {
    Triangle t;
    t.v0 = glm::vec3(1.0f,0.0f,0.0f);
    t.v1 = glm::vec3(0.0f,1.0f,0.0f);
    t.v2 = glm::vec3(0.0f,0.0f,1.0f);
    list.push_back(t);
  }
  for (int r = 0 ; r < recursionLevel ; r++)
  {
      list = SplitTriangle(list);
  }

  for (int octant = 0 ; octant < 8 ; octant++)
  {
    glm::mat4 view(1.0f);
    float angle = 90.0f*(octant%4);
    if(angle != 0.0f)
      view = glm::rotate(view, glm::radians(angle), glm::vec3(1, 0, 0));
    if (octant >= 4)
      view = glm::rotate(view, glm::radians(180.0f), glm::vec3(0, 0, 1));
    for(int i = 0; i < list.size(); i++)
    {
      Triangle t = list[i];
      float mag_reci;
      glm::vec3 v0 = glm::vec3(view*glm::vec4(t.v0, 1.0f));
      glm::vec3 v1 = glm::vec3(view*glm::vec4(t.v1, 1.0f));
      glm::vec3 v2 = glm::vec3(view*glm::vec4(t.v2, 1.0f));
      mag_reci = 1.0f / glm::length(v0);
      v0 = glm::vec3(v0.x * mag_reci, v0.y * mag_reci, v0.z * mag_reci);
      mag_reci = 1.0f / glm::length(v1);
      v1 = glm::vec3(v1.x * mag_reci, v1.y * mag_reci, v1.z * mag_reci);
      mag_reci = 1.0f / glm::length(v2);
      v2 = glm::vec3(v2.x * mag_reci, v2.y * mag_reci, v2.z * mag_reci);
      PushVertex(coords, v0);
      PushVertex(coords, v1);
      PushVertex(coords, v2);
      PushVertex(normals, v0);
      PushVertex(normals, v1);
      PushVertex(normals, v2);
    }
  }
}

void _print_shader_info_log(GLuint shader_index) {
  int max_length = 2048;
  int actual_length = 0;
  char shader_log[2048];
  glGetShaderInfoLog(shader_index, max_length, &actual_length, shader_log);
  printf("shader info log for GL index %u:\n%s\n", shader_index, shader_log);
}

class RenderManager
{
  public:
   enum ShapeType
   {
      SPHERE,
      CYLINDER
   };

                 RenderManager();
   void          SetView(glm::vec3 &c, glm::vec3 &, glm::vec3 &);
   void          SetUpGeometry();
   void          SetColor(double r, double g, double b);
   void          Render(ShapeType, glm::mat4 model);
   GLFWwindow   *GetWindow() { return window; };

  private:
   glm::vec3 color;
   GLuint sphereVAO;
   GLuint sphereNumPrimitives;
   GLuint cylinderVAO;
   GLuint cylinderNumPrimitives;
   GLuint mvploc;
   GLuint colorloc;
   GLuint camloc;
   GLuint ldirloc;
   glm::mat4 projection;
   glm::mat4 view;
   GLuint shaderProgram;
   GLFWwindow *window;

   void SetUpWindowAndShaders();
   void MakeModelView(glm::mat4 &);
};

RenderManager::RenderManager()
{
  SetUpWindowAndShaders();
  SetUpGeometry();
  projection = glm::perspective(
        glm::radians(45.0f), (float)1000 / (float)1000,  5.0f, 100.0f);

  // Get a handle for our MVP and color uniforms
  mvploc = glGetUniformLocation(shaderProgram, "MVP");
  colorloc = glGetUniformLocation(shaderProgram, "color");
  camloc = glGetUniformLocation(shaderProgram, "cameraloc");
  ldirloc = glGetUniformLocation(shaderProgram, "lightdir");

  glm::vec4 lightcoeff(0.3, 0.7, 2.8, 50.5); // Lighting coeff, Ka, Kd, Ks, alpha
  GLuint lcoeloc = glGetUniformLocation(shaderProgram, "lightcoeff");
  glUniform4fv(lcoeloc, 1, &lightcoeff[0]);
}

void
RenderManager::SetView(glm::vec3 &camera, glm::vec3 &origin, glm::vec3 &up)
{ 
   glm::mat4 v = glm::lookAt(
                       camera, // Camera in world space
                       origin, // looks at the origin
                       up      // and the head is up
                 );
   view = v; 
   glUniform3fv(camloc, 1, &camera[0]);
   // Direction of light
   glm::vec3 lightdir = glm::normalize(camera - origin);   
   glUniform3fv(ldirloc, 1, &lightdir[0]);
};

void
RenderManager::SetUpWindowAndShaders()
{
  // start GL context and O/S window using the GLFW helper library
  if (!glfwInit()) {
    fprintf(stderr, "ERROR: could not start GLFW3\n");
    exit(EXIT_FAILURE);
  }

  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

  window = glfwCreateWindow(700, 700, "CIS 441", NULL, NULL);
  if (!window) {
    fprintf(stderr, "ERROR: could not open window with GLFW3\n");
    glfwTerminate();
    exit(EXIT_FAILURE);
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

  const char* vertex_shader = GetVertexShader();
  const char* fragment_shader = GetFragmentShader();

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

  shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, fs);
  glAttachShader(shaderProgram, vs);
  glLinkProgram(shaderProgram);
  glUseProgram(shaderProgram);
}

void RenderManager::SetColor(double r, double g, double b)
{
   color[0] = r;
   color[1] = g;
   color[2] = b;
}

void RenderManager::MakeModelView(glm::mat4 &model)
{
   glm::mat4 modelview = projection * view * model;
   glUniformMatrix4fv(mvploc, 1, GL_FALSE, &modelview[0][0]);
}

void RenderManager::Render(ShapeType st, glm::mat4 model)
{
   int numPrimitives = 0;
   if (st == SPHERE)
   {
      glBindVertexArray(sphereVAO);
      numPrimitives = sphereNumPrimitives;
   }
   else if (st == CYLINDER)
   {
      glBindVertexArray(cylinderVAO);
      numPrimitives = cylinderNumPrimitives;
   }
   MakeModelView(model);
   glUniform3fv(colorloc, 1, &color[0]);
   glDrawElements(GL_TRIANGLES, numPrimitives, GL_UNSIGNED_INT, NULL);
}

void SetUpVBOs(std::vector<float> &coords, std::vector<float> &normals,
               GLuint &points_vbo, GLuint &normals_vbo, GLuint &index_vbo)
{
  int numIndices = coords.size()/3;
  std::vector<GLuint> indices(numIndices);
  for(int i = 0; i < numIndices; i++)
    indices[i] = i;

  points_vbo = 0;
  glGenBuffers(1, &points_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, points_vbo);
  glBufferData(GL_ARRAY_BUFFER, coords.size() * sizeof(float), coords.data(), GL_STATIC_DRAW);

  normals_vbo = 0;
  glGenBuffers(1, &normals_vbo);
  glBindBuffer(GL_ARRAY_BUFFER, normals_vbo);
  glBufferData(GL_ARRAY_BUFFER, normals.size() * sizeof(float), normals.data(), GL_STATIC_DRAW);

  index_vbo = 0;    // Index buffer object
  glGenBuffers(1, &index_vbo);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_vbo);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);
}

void RenderManager::SetUpGeometry()
{
  std::vector<float> sphereCoords;
  std::vector<float> sphereNormals;
  GetSphereData(sphereCoords, sphereNormals);
  sphereNumPrimitives = sphereCoords.size() / 3;
  GLuint sphere_points_vbo, sphere_normals_vbo, sphere_indices_vbo;
  SetUpVBOs(sphereCoords, sphereNormals, 
            sphere_points_vbo, sphere_normals_vbo, sphere_indices_vbo);

  GLuint vao[1];
  glGenVertexArrays(1, vao);

  glBindVertexArray(vao[SPHERE]);
  sphereVAO = vao[SPHERE];
  glBindBuffer(GL_ARRAY_BUFFER, sphere_points_vbo);
  glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ARRAY_BUFFER, sphere_normals_vbo);
  glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphere_indices_vbo);
  glEnableVertexAttribArray(0);
  glEnableVertexAttribArray(1);
}

// Rotation matrix
glm::mat4 RotateMatrix(float degrees, float x, float y, float z)
{
   glm::mat4 identity(1.0f);
   glm::mat4 rotation = glm::rotate(identity, 
                                    glm::radians(degrees), 
                                    glm::vec3(x, y, z));
   return rotation;
}
// Scale matrix
glm::mat4 ScaleMatrix(double x, double y, double z)
{
   glm::mat4 identity(1.0f);
   glm::vec3 scale(x, y, z);
   return glm::scale(identity, scale);
}
// Translation matrix
glm::mat4 TranslateMatrix(double x, double y, double z)
{
   glm::mat4 identity(1.0f);
   glm::vec3 translate(x, y, z);
   return glm::translate(identity, translate);
}

// Vertex shader
const char *GetVertexShader()
{
   static char vertexShader[4096];
   strcpy(vertexShader, 
           "#version 400\n"
           "layout (location = 0) in vec3 vertex_position;\n"
           "layout (location = 1) in vec3 vertex_normal;\n"
           "uniform mat4 MVP;\n"
           "uniform vec3 cameraloc;  // Camera position \n"
           "uniform vec3 lightdir;   // Lighting direction \n"
           "uniform vec4 lightcoeff; // Lighting coeff, Ka, Kd, Ks, alpha\n"
           "out float shading_amount;\n"
           "void main() {\n"
           "  gl_Position = MVP*vec4(vertex_position, 1.0);\n"
           "  float diffuse = max(0.f, (lightdir[0] * vertex_normal[0] + lightdir[1] * vertex_normal[1] + lightdir[2] * vertex_normal[2])) * lightcoeff[1];\n"

           "  float dotProdLightNormal =  lightdir[0] * vertex_normal[0] + lightdir[1] * vertex_normal[1] + lightdir[2] * vertex_normal[2];\n"
           "  float viewDirection[3];\n"
           "  viewDirection[0] = cameraloc[0] - vertex_position[0];\n"
           "  viewDirection[1] = cameraloc[1] - vertex_position[1];\n"
           "  viewDirection[2] = cameraloc[2] - vertex_position[2];\n"

           "  float magnitude = sqrt(pow(viewDirection[0], 2) + pow(viewDirection[1], 2) + pow(viewDirection[2], 2));\n"

           "  for (int i = 0; i < 3; i++)\n"
           "    viewDirection[i] /= magnitude; \n"

           "  vec3 R = vec3(2.f * dotProdLightNormal * vertex_normal[0] - lightdir[0],\n"
           "  2.f * dotProdLightNormal * vertex_normal[1] - lightdir[1], \n"
           "  2.f * dotProdLightNormal * vertex_normal[2] - lightdir[2]);\n"

           "  magnitude = sqrt(pow(R[0], 2) + pow(R[1], 2) + pow(R[2], 2));\n"

           "  for (int i = 0; i < 3; i++)\n"
           "    R[i] /= magnitude; \n"

           "  float specular =  pow(max(0.f, viewDirection[0] * R[0] + viewDirection[1] * R[1] + viewDirection[2] * R[2]), lightcoeff[3]) * lightcoeff[2];\n"

           "  shading_amount = lightcoeff[0] + diffuse + specular;\n"
           "}\n"
         );
   return vertexShader;
}
// Fragment shader
const char *GetFragmentShader()
{
   static char fragmentShader[1024];
   strcpy(fragmentShader, 
           "#version 400\n"
           "in float shading_amount;\n"
           "uniform vec3 color;\n"
           "out vec4 frag_color;\n"
           "void main() {\n"
           "  frag_color = vec4(color, 1.0);\n"
           "  frag_color[0] = min(1.0,frag_color[0]*shading_amount);\n"
           "  frag_color[1] = min(1.0,frag_color[1]*shading_amount);\n"
           "  frag_color[2] = min(1.0,frag_color[2]*shading_amount);\n"
           "}\n"
         );
   return fragmentShader;
}


const double boundingBox[] = {-5,5, -5,1, -5,5};
class Ball
{
  public:
    // ==================
    // DON'T CHANGE THESE
    // ==================
    Ball() 
    {
  	// No collisions yet
  	tSinceAccident = 0;
  	color[0] = 0;
  	color[1] = 1;
  	color[2] = 1;
  	radius = 0.3;
    }; 
    // Initializer
    void initialize(int id)
    {
    	srand(id*2);
    	// Set position and dir to random values
    	pos[0] = (double)rand()/RAND_MAX*6.0-3.0;
  	    pos[1] = (double)rand()/RAND_MAX*6.0-5.0;
  	    pos[2] = (double)rand()/RAND_MAX*6.0-3.0;
  	    dir[0] = (double)rand()/RAND_MAX*0.4-0.2;
  	    dir[1] = (double)rand()/RAND_MAX*0.4-0.2;
  	    dir[2] = (double)rand()/RAND_MAX*0.4-0.2;
    };
    
    void draw(RenderManager &rm)
    {
      glm::mat4 identity(1.0f);
      glm::mat4 translate = TranslateMatrix(pos[0], pos[1], pos[2]);
      glm::mat4 scale = ScaleMatrix(radius, radius, radius);
      rm.SetColor(color[0],color[1],color[2]);
      rm.Render(RenderManager::SPHERE,identity*translate*scale);
    } 
    // =================
    
    // ============
    // CHANGE THESE
    // ============
    void  UpdatePosition()
    {
       // STEP 3: Uncomment and run
       // Keep moving in the direction we should go
       
       pos[0] += dir[0]*0.1;
       pos[1] += dir[1]*0.1;
       pos[2] += dir[2]*0.1;
       
       // Don't escape the box
       if (pos[0] < boundingBox[0])
           dir[0] *= -1;
       if (pos[0] > boundingBox[1])
           dir[0] *= -1;
       if (pos[1] < boundingBox[2])
           dir[1] *= -1;
       if (pos[1] > boundingBox[3])
           dir[1] *= -1;
       if (pos[2] < boundingBox[4])
           dir[2] *= -1;
       if (pos[2] > boundingBox[5])
           dir[2] *= -1;        
    }
    
    double getDistance(Ball &s)
    {
    	// Return the distance between these two balls

        double distanceBetweenCenters = sqrt(
                                            pow(pos[0] - s.pos[0], 2)
                                            + pow(pos[1] - s.pos[1], 2)
                                            + pow(pos[2] - s.pos[2], 2)
                                            );

	    return distanceBetweenCenters - s.radius - radius;
    }

    void LookForCollision(Ball &s)
    {
       // Check if we overlap with another sphere
       // If we do, call CollisionDetected

        if (getDistance(s) < 0.f)
            CollisionDetected(s);
    }
    
    void CollisionDetected(Ball &s)
    {        
        collideDir(pos, dir, s.pos, s.dir);

        color[0] += 0.7 * (1 - color[0]);
        color[1] -= 0.1 * color[1];
        color[2] -= 0.1 * color[2];
        radius -= 0.25 * radius;
        tSinceAccident = 0;

        s.color[0] += 0.7 * (1 - s.color[0]);
        s.color[1] -= 0.1 * s.color[1];
        s.color[2] -= 0.1 * s.color[2];
        s.radius -= 0.25 * s.radius;
        s.tSinceAccident = 0;
    }
    
    void collideDir(double pos1[3], double dir1[3], double pos2[3], double dir2[3])
    {
    	// Calculate the new collision speed and direction
    	// by editing the direction vector for each sphere
    	
    	// 1. Find the normal vector between the two balls
    	//     n = (pos1-pos2)/(|pos1-pos2|)
    	// 2. Calculate the relative velocity
    	//     vRelative = dir1-dir2
    	// 3. Calculate the relative along the normal
    	//     vNormal = (vRelative (.) n ) * n  [where (.) is the dot operator]
    	// 4. Modify dir using the relative along the normal
    	//     add or subtract vNormal from dir1 and dir2 as appropriate
    	
    	// Any collision between two spheres can be modeled as a collision between two circles on a 2d plane bisecting them both. 
    	// This is the intuition for deriving the above vector math. 
    	// However, all of the components (finding normals, dot products, normalizing, etc.) are things we have done before. 
    	// Feel free to write/copy your own helper methods to implement this clearly.	

        double n[] = 
        {
            (pos1[0] - pos2[0]),
            (pos1[1] - pos2[1]),
            (pos1[2] - pos2[2])
        };

        double magnitudeN = sqrt(pow(n[0], 2) + pow(n[1], 2) + pow(n[2], 2));

        for (int i = 0; i < 3; i++)
            n[i] /= magnitudeN;

        double vRelative[] =
        {
            dir1[0] - dir2[0],
            dir1[1] - dir2[1],
            dir1[2] - dir2[2]
        };

        double dotProdVelN = vRelative[0] * n[0] + vRelative[1] * n[1] + vRelative[2] * n[2];

        double vNormal[] =
        {
            dotProdVelN * n[0],
            dotProdVelN * n[1],
            dotProdVelN * n[2]
        };

        for (int i = 0; i < 3; i++)
        {
            dir1[i] -= vNormal[i];
            dir2[i] += vNormal[i];
        }
    }
    
    // VARIABLES
    double pos[3];
    double dir[3]; 
    int tSinceAccident;
    double color[3];
    double radius;
};

int main() 
{
  // Create the balls
  const int numBalls = 50;
  Ball *ballList = new Ball[numBalls];

  // Set this to -1 to run forever
  int TICK_LIMIT = 1000;
  
  // Create the render manager, window, etc.
  RenderManager rm;
  GLFWwindow *window = rm.GetWindow();
  glm::vec3 origin(0, 0, 0);
  glm::vec3 up(0, 1, 0);
  
  // Initialize the balls (necessary to seed randomness)
  for (int i = 0; i < numBalls; i++) {ballList[i].initialize(i);}

  // MAIN LOOP
  int tick=0;
  while (!glfwWindowShouldClose(window)) 
  {
      if (TICK_LIMIT != 0)
      {
          // Set up the camera
          glm::vec3 camera(10, -10, -15);
          rm.SetView(camera, origin, up);
          // wipe the drawing surface clear
          glClearColor(0, 0, 0, 1.0);
          glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

          // Each tick, check all the balls for collisions 
          for (int i = 0; i < numBalls; i++)
          {
              // Look for a collision with every other ball
              for (int j = i + 1; j < numBalls; j++)
                  ballList[i].LookForCollision(ballList[j]);

              // Move the ball
              ballList[i].UpdatePosition();
              // Draw the ball
              ballList[i].draw(rm);
          }
          // put the stuff we've been drawing onto the display
          glfwSwapBuffers(window);

          //loop again to see which ones have a tSinceAccident of 0
          for (int i = 0; i < numBalls; i++)
          {
              //check if we need to make it bigger (if it hasn't run into anything in a while)
              if (++ballList[i].tSinceAccident % 51 == 0) //incrementing before comparison saves us a comparison to zero, but requires modulo by 51 instead of 50
              {
                  ballList[i].color[0] -= 0.05 * ballList[i].color[0];
                  ballList[i].color[1] += 0.35 * (1 - ballList[i].color[1]);
                  ballList[i].color[2] += 0.35 * (1 - ballList[i].color[2]);
                  ballList[i].radius += 0.01 * (1 - ballList[i].radius);
              }
          }

          // Make a "Tick" pass
          tick++;
          if (TICK_LIMIT > 0) { TICK_LIMIT--; }
      }

    // update other events like input handling
    glfwPollEvents();
  }

  // close GL context and any other GLFW resources
  glfwTerminate();
  return 0;
}