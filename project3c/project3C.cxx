/*
* CIS 441 - Project 3C
* 
* Please see Findings.txt for more information.
* 
* Thomas Mitchell
* 
*/

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <cmath>
#include <random>
#include <algorithm>

using std::endl;
using std::cerr;
#include <GL/glew.h>    // include GLEW and new version of GL on Windows
#include <GLFW/glfw3.h> // GLFW helper libray

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/vec3.hpp>   // glm::vec3
#include <glm/vec4.hpp>   // glm::vec4
#include <glm/mat4x4.hpp> // glm::mat4
#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtc/matrix_transform.hpp>  // glm::translate, glm::rotate, glm::scale

#define M_PI 3.14159265358979323846

// Macros to define the level of recursion for each sphere
#define L1 5
#define L2 6
#define L3 7

// How many number of balls to render
#define numBalls 200

class RenderManager;

struct Ball_physics {
    double h0;         // initial height
    double v = 0;          // m/s, current velocity
    double g = 10;         // m/s/s
    double t = 0;          // starting time
    double dt = 0.02;     // time step
    double rho = 0.75;     // coefficient of restitution
    double tau = 0.10;     // contact time for bounce
    double hmax;      // keep track of the maximum height
    double h;
    double hstop = 0.01;   // stop when bounce is less than 1 cm
    bool freefall = true;
    double t_last;
    double vmax;
};

// Ball has a x,z coordinate, distance to camera and a bunch of paramters to track its movement
struct Ball{
    double x; // x coord of ball movement
    double z; // z coord of ball movement
    double distToCamera; 
    Ball_physics BP;
};


void BounceBall(std::vector<Ball> &balls, RenderManager &rm, glm::vec3 camPos);
const char *GetVertexShader();
const char *GetFragmentShader();

glm::mat4 RotateMatrix(float degrees, float x, float y, float z)
{
    glm::mat4 identity(1.0f);
    glm::mat4 rotation = glm::rotate(identity, 
            glm::radians(degrees), 
            glm::vec3(x, y, z));
    return rotation;
}

glm::mat4 ScaleMatrix(double x, double y, double z)
{
    glm::mat4 identity(1.0f);
    glm::vec3 scale(x, y, z);
    return glm::scale(identity, scale);
}

glm::mat4 TranslateMatrix(double x, double y, double z)
{
    glm::mat4 identity(1.0f);
    glm::vec3 translate(x, y, z);
    return glm::translate(identity, translate);
}


//
//
// PART 1: code to set up spheres and plane
//
//

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

void GetPlaneData(std::vector<float>& coords, std::vector<float>& normals){
    glm::vec3 v1 = glm::vec3(-15, -1, -15);
    glm::vec3 v2 = glm::vec3(0, -1,  -15);
    glm::vec3 v3 = glm::vec3(0.0, -1, 0.0);
    glm::vec3 v4 = glm::vec3(-15, -1, 0);
    
    glm::vec3 norm = glm::vec3(0.0f, 1.0f, 0.0f);
    
    PushVertex(coords, v1);
    PushVertex(coords, v2);
    PushVertex(coords, v4);

    PushVertex(coords, v2);
    PushVertex(coords, v3);
    PushVertex(coords, v4);

    PushVertex(normals, norm);
    PushVertex(normals, norm);
    PushVertex(normals, norm);
    
    PushVertex(normals, norm);
    PushVertex(normals, norm);
    PushVertex(normals, norm);
}
//
// Sets up a sphere with equation x^2+y^2+z^2=1
//
void
GetSphereData(std::vector<float>& coords, std::vector<float>& normals, int recursionLevel)
{
    std::vector<Triangle> list;
    {
        Triangle t;
        t.v0 = glm::vec3(1.0f,0.0f,0.0f);
        t.v1 = glm::vec3(0.0f,1.0f,0.0f);
        t.v2 = glm::vec3(0.0f,0.0f,1.0f);
        list.push_back(t);
    }
    for (int r = 0 ; r < recursionLevel ; r++)
        list = SplitTriangle(list);

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


//
//
// PART 2: RenderManager module
//
//

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
            SPHERE1,
            SPHERE2,
            SPHERE3,
            PLANE
        };

        RenderManager();
        void          SetView(glm::vec3 &c, glm::vec3 &, glm::vec3 &);
        void          SetUpGeometry();
        void          SetColor(double r, double g, double b);
        void          Render(ShapeType, glm::mat4 model);
        GLFWwindow   *GetWindow() { return window; };

    private:
        glm::vec3 color;
        GLuint sphere1VAO;
        GLuint sphere1NumPrimitives;
        GLuint sphere2VAO;
        GLuint sphere2NumPrimitives;
        GLuint sphere3VAO;
        GLuint sphere3NumPrimitives;
        GLuint planeVAO;
        GLuint planeNumPrimitives;
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
};

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
    if (st == SPHERE1)
    {
        glBindVertexArray(sphere1VAO);
        numPrimitives = sphere1NumPrimitives;
    }
    else if (st == SPHERE2)
    {
        glBindVertexArray(sphere2VAO);
        numPrimitives = sphere2NumPrimitives;
    }
    else if (st == SPHERE3)
    {
        glBindVertexArray(sphere3VAO);
        numPrimitives = sphere3NumPrimitives;
    }
    else if (st == PLANE)
    {
        glBindVertexArray(planeVAO);
        numPrimitives = planeNumPrimitives;
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
    std::vector<float> sphere1Coords;
    std::vector<float> sphere1Normals;
    GetSphereData(sphere1Coords, sphere1Normals, L1);
    sphere1NumPrimitives = sphere1Coords.size() / 3;
    GLuint sphere1_points_vbo, sphere1_normals_vbo, sphere1_indices_vbo;
    SetUpVBOs(sphere1Coords, sphere1Normals, 
            sphere1_points_vbo, sphere1_normals_vbo, sphere1_indices_vbo);
    
    std::vector<float> sphere2Coords;
    std::vector<float> sphere2Normals;
    GetSphereData(sphere2Coords, sphere2Normals, L2);
    sphere2NumPrimitives = sphere2Coords.size() / 3;
    GLuint sphere2_points_vbo, sphere2_normals_vbo, sphere2_indices_vbo;
    SetUpVBOs(sphere2Coords, sphere2Normals, 
            sphere2_points_vbo, sphere2_normals_vbo, sphere2_indices_vbo);
    
    std::vector<float> sphere3Coords;
    std::vector<float> sphere3Normals;
    GetSphereData(sphere3Coords, sphere3Normals, L3);
    sphere3NumPrimitives = sphere3Coords.size() / 3;
    GLuint sphere3_points_vbo, sphere3_normals_vbo, sphere3_indices_vbo;
    SetUpVBOs(sphere3Coords, sphere3Normals, 
            sphere3_points_vbo, sphere3_normals_vbo, sphere3_indices_vbo);
    
    std::vector<float> planeCoords;
    std::vector<float> planeNormals;
    GetPlaneData(planeCoords, planeNormals);
    planeNumPrimitives = planeCoords.size() / 3;
    GLuint plane_points_vbo, plane_normals_vbo, plane_indices_vbo;
    SetUpVBOs(planeCoords, planeNormals, 
            plane_points_vbo, plane_normals_vbo, plane_indices_vbo);

    GLuint vao[4];
    glGenVertexArrays(4, vao);

    glBindVertexArray(vao[SPHERE1]);
    sphere1VAO = vao[SPHERE1];
    glBindBuffer(GL_ARRAY_BUFFER, sphere1_points_vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindBuffer(GL_ARRAY_BUFFER, sphere1_normals_vbo);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphere1_indices_vbo);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);

    glBindVertexArray(vao[SPHERE2]);
    sphere2VAO = vao[SPHERE2];
    glBindBuffer(GL_ARRAY_BUFFER, sphere2_points_vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindBuffer(GL_ARRAY_BUFFER, sphere2_normals_vbo);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphere2_indices_vbo);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    
    glBindVertexArray(vao[SPHERE3]);
    sphere3VAO = vao[SPHERE3];
    glBindBuffer(GL_ARRAY_BUFFER, sphere3_points_vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindBuffer(GL_ARRAY_BUFFER, sphere3_normals_vbo);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphere3_indices_vbo);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
   
    glBindVertexArray(vao[PLANE]);
    planeVAO = vao[PLANE];
    glBindBuffer(GL_ARRAY_BUFFER, plane_points_vbo);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindBuffer(GL_ARRAY_BUFFER, plane_normals_vbo);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 0, NULL);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, plane_indices_vbo);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
}


void DrawPlane(RenderManager &rm){
    glm::mat4 identity(1.0f);
    rm.SetColor(1,1,1);
    rm.Render(RenderManager::PLANE, identity);
    glm::mat4 transform = TranslateMatrix(15, 0, 0);
    rm.SetColor(192/255.0, 192/255.0, 192/255.0);
    rm.Render(RenderManager::PLANE, transform);
    transform = TranslateMatrix(0,0, 15);
    rm.Render(RenderManager::PLANE, transform);
    transform = TranslateMatrix(15, 0, 15);
    rm.SetColor(1,1,1);
    rm.Render(RenderManager::PLANE, transform);
}


//
// PART3: main function
//
int main() 
{
    std::random_device r;
    RenderManager rm;
    GLFWwindow *window = rm.GetWindow();

    glm::vec3 origin(0, 0, 0);
    glm::vec3 up(0, 1, 0);

    // Vector that holds all the balls
    std::vector<Ball> balls(numBalls);

    for (int i = 0; i < numBalls; i++)
    {
        balls.at(i).x = (double)rand() / RAND_MAX * 20 - 10;
        balls.at(i).z = (double)rand() / RAND_MAX * 20 - 10;
        double y = (double)rand() / RAND_MAX * 30 + 10;
        balls.at(i).BP.h0 = y;
        balls.at(i).BP.h = y;
        balls.at(i).BP.hmax = y;
        balls.at(i).BP.t_last = -1.f * sqrt(2.f * y / balls.at(i).BP.g);
        balls.at(i).BP.vmax = sqrt(2.f * y * balls.at(i).BP.g);
    }
    
    glm::mat4 identity(1.0f);

    int counter=0;
    double t, t0;
    t0 = glfwGetTime();
    
    while (!glfwWindowShouldClose(window)) 
    {
        t = glfwGetTime();

        double angle=counter/1000.0*2*M_PI;
        counter++;
        // print FPS to console
        if (counter % 100 == 99)
        {
            double time = t-t0;
            cerr << "Frame rate = " << 100.0 / time << " FPS" << endl;
            t0 = t;
        } 
        // adjust the camera
        glm::vec3 cameraPos = glm::vec3(15.0f*sin(angle), 10.0f, 15.0f*cos(angle));
        glm::vec3 cameraDirection = glm::normalize(cameraPos - origin);
        glm::vec3 cameraRight = glm::normalize(glm::cross(up, cameraDirection));
        glm::vec3 cameraUp = glm::cross(cameraDirection, cameraRight);

        rm.SetView(cameraPos, origin, cameraUp);

        // wipe the drawing surface clear
        glClearColor(0.3, 0.3, 0.8, 1.0);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // draw the plane at the bottom
        DrawPlane(rm);

        // render the bouncing balls!!
        BounceBall(balls, rm, cameraPos); 

        // update other events like input handling
        glfwPollEvents();
        // put the stuff we've been drawing onto the display
        glfwSwapBuffers(window);
    }

    // close GL context and any other GLFW resources
    glfwTerminate();
    return 0;
}


void
UpdateBallPhysics(Ball& ball){
    auto& ball_physics = ball.BP;
    // If the ball is still bouncing, calculate the new height
    if (ball_physics.hmax > ball_physics.hstop)
    {
        if (ball_physics.freefall){
            double hnew = ball_physics.h + ball_physics.v*ball_physics.dt - 0.5*ball_physics.g*ball_physics.dt*ball_physics.dt;
            if (hnew<0){
                ball_physics.t = ball_physics.t_last + 2*sqrt(2*ball_physics.hmax/ball_physics.g);
                ball_physics.freefall = false;
                ball_physics.t_last = ball_physics.t + ball_physics.tau;
                ball_physics.h = 0;
            }
            else{
                ball_physics.t += ball_physics.dt;
                ball_physics.v -=  ball_physics.g*ball_physics.dt;
                ball_physics.h = hnew;
            }
        }
        else{
            ball_physics.t +=  ball_physics.tau;
            ball_physics.vmax *= ball_physics.rho;
            ball_physics.v = ball_physics.vmax;
            ball_physics.freefall = true;
            ball_physics.h = 0;
        }
        ball_physics.hmax = 0.5*ball_physics.vmax*ball_physics.vmax/ball_physics.g;
    } 
    else // When its done, make a new ball
    { 
        ball.BP.h = ball.BP.h0;
        ball.BP.hmax = ball.BP.h0;
        ball.BP.t_last = -1.f * sqrt(2.f * ball.BP.h0 / ball.BP.g);
        ball.BP.vmax = sqrt(2.f * ball.BP.hmax * ball.BP.g);
    }
}


void
BounceBall(std::vector<Ball> &balls, RenderManager &rm, glm::vec3 camPos)
{
    for (int i = 0; i < numBalls; i++)
    {
        UpdateBallPhysics(balls.at(i));
        balls.at(i).distToCamera = sqrt(pow(balls.at(i).x - camPos[0], 2)
                                        + pow(balls.at(i).BP.h - camPos[1], 2)
                                        + pow(balls.at(i).z - camPos[2], 2));
    }
    
    sort(balls.begin(), balls.end(), [](Ball ball1, Ball ball2) {return ball1.distToCamera < ball2.distToCamera;});


    for (int i = 0; i < numBalls / 3; i++)
    {
        glm::mat4 translate = TranslateMatrix(balls.at(i).x, balls.at(i).BP.h, balls.at(i).z);
        glm::mat4 scale = ScaleMatrix(.3, .3, .3);
        rm.SetColor(1.f, 0.f, 0.f);
        rm.Render(RenderManager::SPHERE3, translate * scale);
    }
    for (int i = numBalls / 3; i < 2 * numBalls / 3; i++)
    {
        glm::mat4 translate = TranslateMatrix(balls.at(i).x, balls.at(i).BP.h, balls.at(i).z);
        glm::mat4 scale = ScaleMatrix(.3, .3, .3);
        rm.SetColor(0.f, 1.f, 0.f);
        rm.Render(RenderManager::SPHERE2, translate * scale);
    }
    for (int i = 2 * numBalls / 3; i < numBalls; i++)
    {
        glm::mat4 translate = TranslateMatrix(balls.at(i).x, balls.at(i).BP.h, balls.at(i).z);
        glm::mat4 scale = ScaleMatrix(.3, .3, .3);
        rm.SetColor(0.f, 0.f, 1.f);
        rm.Render(RenderManager::SPHERE1, translate * scale);
    }
}


const char *GetVertexShader()
{
    static char vertexShader[2048];
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
                "gl_Position = MVP*vec4(vertex_position, 1.0);\n"
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
            "}\n"
            );
    return vertexShader;
}

const char *GetFragmentShader()
{
    static char fragmentShader[1024];
    strcpy(fragmentShader, 
            "#version 400\n"
            "in float shading_amount;\n"
            "uniform vec3 color;\n"
            "out vec4 frag_color;\n"
            "void main() {\n"
            "  frag_color = vec4(min(1.0, shading_amount*color[0]),min(1.0, shading_amount*color[1]), min(1.0, shading_amount*color[2]), 1.0);\n"
            "}\n"
          );
    return fragmentShader;
}