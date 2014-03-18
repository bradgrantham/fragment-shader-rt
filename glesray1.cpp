//
// Copyright 2013-2014, Bradley A. Grantham
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//      http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// 

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <string>
#include <ctime>
#include <cerrno>
#include <sys/time.h>
#include <vector>
#include <limits>

const bool debug_gl = false;

#define GL3_PROTOTYPES
#include <GL3/gl3.h>
#define __gl_h_

void query_log(void)
{
    GLint totalMessages, len;

    glGetIntegerv(GL_DEBUG_LOGGED_MESSAGES_ARB, &totalMessages);
    printf("Number of messages in the log:%d\n", totalMessages);

    for(int i = 0; i < totalMessages; i++) {
        GLenum source, type, id, severity;

        glGetIntegerv(GL_DEBUG_NEXT_LOGGED_MESSAGE_LENGTH_ARB, &len);
        char *message = new char[len];

        glGetDebugMessageLogARB(1, len, &source, &type, &id, &severity, NULL, message);

        printf("GL: \"%s\"\n", message);

        delete[] message;
    }
}

#include <GL/freeglut.h>

#include "world.h"

int gWindowWidth, gWindowHeight;
unsigned char *gImage;
world *gWorld;
timeval prev_frame_time;
float zoom = 0.0f;
float object_rotation[4];
float light_rotation[4];
vec3 light_dir(0, 0, -1);
vec3 object_position(0, 0, 0);
int which;

void drag_to_rotation(float dx, float dy, float rotation[4])
{
    float dist;

    /* XXX grantham 990825 - this "dist" doesn't make me confident. */
    /* but I put in the *10000 to decrease chance of underflow  (???) */
    dist = sqrt(dx * 10000 * dx * 10000 + dy * 10000 * dy * 10000) / 10000;
    /* dist = sqrt(dx * dx + dy * dy); */

    rotation[0] = M_PI * dist;
    rotation[1] = dy / dist;
    rotation[2] = dx / dist;
    rotation[3] = 0.0f;
}

void trackball_motion(float prevrotation[4], float dx, float dy, float newrotation[4])
{
    if(dx != 0 || dy != 0) {
        float rotation[4];
        drag_to_rotation(dx, dy, rotation);
        rotation_mult_rotation(prevrotation, rotation, newrotation);
    }
}

void create_camera_matrix(const vec3& viewpoint, float yaw, float pitch, float roll, float matrix[16], float normal_matrix[16])
{
    mat4_make_identity(matrix);

    /* rotate around Z to roll */
    float roll_matrix[16];
    mat4_make_rotation(roll, 0, 0, 1, roll_matrix);
    mat4_mult(roll_matrix, matrix, matrix);

    /* rotate around X to pitch */
    float pitch_matrix[16];
    mat4_make_rotation(pitch, 1, 0, 0, pitch_matrix);
    mat4_mult(pitch_matrix, matrix, matrix);

    /* rotate around Y to yaw */
    float yaw_matrix[16];
    mat4_make_rotation(yaw, 0, 1, 0, yaw_matrix);
    mat4_mult(yaw_matrix, matrix, matrix);

    float viewpoint_matrix[16];
    mat4_make_translation(-viewpoint.x, -viewpoint.y, -viewpoint.z, viewpoint_matrix);
    mat4_mult(viewpoint_matrix, matrix, matrix);

    mat4_invert(matrix, normal_matrix);
    mat4_transpose(normal_matrix, normal_matrix);
    normal_matrix[3] = 0.0;
    normal_matrix[7] = 0.0;
    normal_matrix[11] = 0.0;
}

void create_object_matrix(const vec3& center, const float rotation[4], const vec3& position, float matrix[16], float inverse[16], float normal[16], float normal_inverse[16])
{
    mat4_make_rotation(rotation[0], rotation[1], rotation[2], rotation[3], matrix);
    float m2[16];
    mat4_make_translation(center.x + position.x, center.y + position.y, center.z + position.z, m2);
    mat4_mult(matrix, m2, matrix);

    mat4_invert(matrix, inverse);
    mat4_transpose(matrix, normal);
    mat4_invert(normal, normal);
    normal[3] = 0.0;
    normal[7] = 0.0;
    normal[11] = 0.0;
    mat4_transpose(matrix, normal_inverse);
    normal_inverse[3] = 0.0;
    normal_inverse[7] = 0.0;
    normal_inverse[11] = 0.0;
}

void update_light()
{
    float light_matrix[16];
    float light_normal[16];
    float light_transpose[16];

    mat4_make_rotation(light_rotation[0], light_rotation[1], light_rotation[2], light_rotation[3], light_matrix);
    mat4_transpose(light_matrix, light_transpose);
    mat4_invert(light_transpose, light_normal);
    light_normal[3] = 0.0;
    light_normal[7] = 0.0;
    light_normal[11] = 0.0;

    vec4 l1(0, 0, -1, 0), l2;
    l2 = mat4_mult_vec4(light_normal, l1);
    light_dir.x = l2.x;
    light_dir.y = l2.y;
    light_dir.z = l2.z;
}

void update_view_params(world *world, float zoom)
{
    vec3 offset;

    offset.x = 0;
    offset.y = 0;
    offset.z = zoom;

    world->cam.eye = vec3_subtract(world->scene_center, offset);

    create_camera_matrix(offset, 0, 0, 0, world->camera_matrix, world->camera_normal_matrix);

    create_object_matrix(world->scene_center, object_rotation, object_position, world->object_matrix, world->object_inverse, world->object_normal_matrix, world->object_normal_inverse);
}

enum render_mode_enum {CPU, OPENGL_GEOMETRY, FRAGMENT_SHADER} render_mode = FRAGMENT_SHADER;

static void check_opengl(const char *filename, int line)
{
    int glerr;

    if((glerr = glGetError()) != GL_NO_ERROR) {
        if(debug_gl) {
            printf("GL Error discovered at %s:%d\n", filename, line);
            query_log();
        }
        else
            printf("GL Error: %s(%04X) at %s:%d\n", gluErrorString(glerr), glerr,
                filename, line);
    }
}

bool gPrintShaderErrorInfo = true;

static bool CheckShaderCompile(GLhandleARB shader, std::string shader_name)
{
    if(gPrintShaderErrorInfo) {
        int length;
        glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &length);
        query_log();

        if (length > 0) {
            char log[length];
            glGetShaderInfoLog(shader, length, NULL, log);
            fprintf(stderr, "%s shader error log:\n%s\n", shader_name.c_str(), log);
        }

        // fprintf(stderr, "shader text:\n");
        // glGetShaderiv(shader, GL_SHADER_SOURCE_LENGTH, &length);
        // char source[length];
        // glGetShaderSource(shader, length, NULL, source);
        // fprintf(stderr, "%s\n", source);
    }

    int status;
    glGetShaderiv(shader, GL_COMPILE_STATUS, &status);
    if(status == GL_TRUE)
	return true;

    fprintf(stderr, "%s compile failure.\n", shader_name.c_str());
    return false;
}

static bool CheckProgramLink(GLhandleARB program)
{
    if(gPrintShaderErrorInfo) {
        int log_length;
        glGetProgramiv(program, GL_INFO_LOG_LENGTH, &log_length);

        if (log_length > 0) {
            char log[log_length];
            glGetProgramInfoLog(program, log_length, NULL, log);
            fprintf(stderr, "program link log: %s\n",log);
        }
    }

    int status;
    glGetProgramiv(program, GL_LINK_STATUS, &status);
    return status == GL_TRUE;
}

template <class T>
inline T round_up(T v, unsigned int r)
{
    return ((v + r - 1) / r) * r;
}

static char *load_text(FILE *fp)
{
    char *text = NULL;
    char line[1024];
    int newsize;

    text = strdup("");
    while(fgets(line, 1023, fp) != NULL) {
        newsize = round_up(strlen(text) + strlen(line) + 1, 65536);
        text = (char *)realloc(text,  newsize);
        if(text == NULL) {
            fprintf(stderr, "loadShader: Couldn't realloc program string to"
                " add more characters\n");
            exit(EXIT_FAILURE);
        }
        strcat(text, line);
    }
    
    return text;
}

GLint modelview_uniform = -1;
GLint texture_uniform = -1;
GLint pos_attrib = 0, texcoord_attrib = 1;
GLuint decalprogram;

static const char *gDecalVertexShaderText =
"uniform mat4 modelview;\n"
"attribute vec4 pos;\n"
"attribute vec2 vtex;\n"
"varying vec2 ftex;\n"
"\n"
"void main()\n"
"{\n"
"\n"
"    gl_Position = modelview * pos;\n"
"    ftex = vtex;\n"
"}\n";

static const char *gDecalFragmentShaderText =
"uniform sampler2D texture_image;\n"
"\n"
"varying vec2 ftex;\n"
"\n"
"void main()\n"
"{\n"
"    gl_FragColor = texture2D(texture_image, ftex);\n"
"}\n";

static char *gRayTracingFragmentShaderText;
static char *gRayTracingVertexShaderText;

GLuint vert_buffer;
GLuint texcoord_buffer;
GLuint screenquad_vao;
GLuint decal_vertex_shader;

struct raytracer_gl_binding
{
    GLint sphere_geometry_uniform;
    GLint sphere_colors_uniform;
    GLint group_children_uniform;
    GLint group_objects_uniform;
    GLint group_hitmiss_uniform;
    GLint group_directions_uniform;
    GLint group_boxmax_uniform;
    GLint group_boxmin_uniform;

    GLuint sphere_geometry_texture;
    GLuint sphere_colors_texture;
    GLuint group_children_texture;
    GLuint group_objects_texture;
    GLuint group_hitmiss_texture;
    GLuint group_directions_texture;
    GLuint group_boxmin_texture;
    GLuint group_boxmax_texture;

    GLint which_uniform;
    GLint tree_root_uniform;

    GLint sphere_data_row_size_uniform;
    GLint group_data_row_size_uniform;
    GLint sphere_data_rows_uniform;
    GLint group_data_rows_uniform;

    GLint modelview_uniform;
    GLint eye_uniform;
    GLint camera_matrix_uniform;
    GLint camera_normal_matrix_uniform;
    GLint object_matrix_uniform;
    GLint object_inverse_uniform;
    GLint object_normal_matrix_uniform;
    GLint object_normal_inverse_uniform;
    GLint fov_uniform;
    GLint aspect_uniform;
    GLint background_color_uniform;
    GLint light_dir_uniform;

    GLuint vertex_shader;
    GLuint fragment_shader;
    GLuint program;
};

unsigned int scene_data_row_pad = 512;
scene_shader_data scene_data;
raytracer_gl_binding raytracer_gl;


void load_scene_data(world *w, raytracer_gl_binding &binding)
{
    const char *filename;
    if(getenv("SHADER") != NULL)
        filename = getenv("SHADER");
    else
        filename = "raytracer.es.fs";
    FILE *fp = fopen(filename, "r");
    if(fp == NULL) {
        fprintf(stderr, "couldn't open raytracer.fs\n");
        exit(EXIT_FAILURE);
    }
    gRayTracingFragmentShaderText = load_text(fp);
    fclose(fp);

    fp = fopen("raytracer.vs", "r");
    if(fp == NULL) {
        fprintf(stderr, "couldn't open raytracer.vs\n");
        exit(EXIT_FAILURE);
    }
    gRayTracingVertexShaderText = load_text(fp);
    fclose(fp);

    unsigned int spheresqrt = (unsigned int)ceil(sqrt(w->sphere_count));
    scene_data_row_pad = 1;
    while(scene_data_row_pad < spheresqrt)
        scene_data_row_pad *= 2;
    printf("scene data row pad = %d\n", scene_data_row_pad);

    get_shader_data(w, scene_data, scene_data_row_pad);

    if(false) {
        printf("memory use by scene data\n");
        printf("%f megabytes in sphere geometry\n", scene_data.sphere_count * sizeof(float) * 4 / 1000000.0);
        printf("%f megabytes in sphere colors\n", scene_data.sphere_count * sizeof(float) * 3 / 1000000.0);
        printf("%d groups\n", scene_data.group_count);
        printf("%f megabytes in group bounds\n", scene_data.group_count * sizeof(float) * 4 / 1000000.0);
        printf("%f megabytes in group children\n", scene_data.group_count * sizeof(int) * 2 / 1000000.0);
        printf("total %f megabytes\n", (scene_data.sphere_count * (sizeof(float) * 4 + sizeof(float) * 3) + scene_data.group_count * (sizeof(float) * 4 + sizeof(int) * 2)) / 1000000.0);
    }

    raytracer_gl.fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    char *string = gRayTracingFragmentShaderText;
    glShaderSource(raytracer_gl.fragment_shader, 1, &string, NULL);
    glCompileShader(raytracer_gl.fragment_shader);
    if(!CheckShaderCompile(raytracer_gl.fragment_shader, "ray tracer fragment shader"))
	exit(1);

    raytracer_gl.vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    string = gRayTracingVertexShaderText;
    glShaderSource(raytracer_gl.vertex_shader, 1, &string, NULL);
    glCompileShader(raytracer_gl.vertex_shader);
    if(!CheckShaderCompile(raytracer_gl.vertex_shader, "ray tracer vertex shader"))
	exit(1);

    raytracer_gl.program = glCreateProgram();
    glAttachShader(raytracer_gl.program, raytracer_gl.vertex_shader);
    glAttachShader(raytracer_gl.program, raytracer_gl.fragment_shader);
    glBindAttribLocation(raytracer_gl.program, pos_attrib, "pos");
    glBindAttribLocation(raytracer_gl.program, texcoord_attrib, "vtex");
    glLinkProgram(raytracer_gl.program);
    if(!CheckProgramLink(raytracer_gl.program))
	exit(1);

    glUseProgram(raytracer_gl.program);
    check_opengl(__FILE__, __LINE__);

    raytracer_gl.light_dir_uniform = glGetUniformLocation(raytracer_gl.program, "light_dir");
    raytracer_gl.modelview_uniform = glGetUniformLocation(raytracer_gl.program, "modelview");
    raytracer_gl.sphere_geometry_uniform = glGetUniformLocation(raytracer_gl.program, "sphere_geometry");
    raytracer_gl.sphere_colors_uniform = glGetUniformLocation(raytracer_gl.program, "sphere_colors");

    raytracer_gl.group_children_uniform = glGetUniformLocation(raytracer_gl.program, "group_children");
    raytracer_gl.group_objects_uniform = glGetUniformLocation(raytracer_gl.program, "group_objects");
    raytracer_gl.group_hitmiss_uniform = glGetUniformLocation(raytracer_gl.program, "group_hitmiss");
    raytracer_gl.group_directions_uniform = glGetUniformLocation(raytracer_gl.program, "group_directions");
    raytracer_gl.group_boxmin_uniform = glGetUniformLocation(raytracer_gl.program, "group_boxmin");
    raytracer_gl.group_boxmax_uniform = glGetUniformLocation(raytracer_gl.program, "group_boxmax");

    raytracer_gl.sphere_data_row_size_uniform = glGetUniformLocation(raytracer_gl.program, "sphere_data_row_size");
    raytracer_gl.sphere_data_rows_uniform = glGetUniformLocation(raytracer_gl.program, "sphere_data_rows");
    raytracer_gl.group_data_row_size_uniform = glGetUniformLocation(raytracer_gl.program, "group_data_row_size");
    raytracer_gl.group_data_rows_uniform = glGetUniformLocation(raytracer_gl.program, "group_data_rows");
    raytracer_gl.which_uniform = glGetUniformLocation(raytracer_gl.program, "which");
    raytracer_gl.tree_root_uniform = glGetUniformLocation(raytracer_gl.program, "tree_root");
    raytracer_gl.eye_uniform = glGetUniformLocation(raytracer_gl.program, "eye");
    raytracer_gl.camera_matrix_uniform = glGetUniformLocation(raytracer_gl.program, "camera_matrix");
    raytracer_gl.camera_normal_matrix_uniform = glGetUniformLocation(raytracer_gl.program, "camera_normal_matrix");
    raytracer_gl.object_matrix_uniform = glGetUniformLocation(raytracer_gl.program, "object_matrix");
    raytracer_gl.object_inverse_uniform = glGetUniformLocation(raytracer_gl.program, "object_inverse");
    raytracer_gl.object_normal_matrix_uniform = glGetUniformLocation(raytracer_gl.program, "object_normal_matrix");
    raytracer_gl.object_normal_inverse_uniform = glGetUniformLocation(raytracer_gl.program, "object_normal_inverse");
    raytracer_gl.fov_uniform = glGetUniformLocation(raytracer_gl.program, "fov");
    raytracer_gl.aspect_uniform = glGetUniformLocation(raytracer_gl.program, "aspect");
    raytracer_gl.background_color_uniform = glGetUniformLocation(raytracer_gl.program, "background_color");
    check_opengl(__FILE__, __LINE__);

    printf("sphere data textures are %dx%d.\n", scene_data_row_pad, scene_data.sphere_data_rows);
    printf("group data textures are %dx%d.\n", scene_data_row_pad, scene_data.group_data_rows);

    glGenTextures(1, &raytracer_gl.sphere_geometry_texture);
    glBindTexture(GL_TEXTURE_2D, raytracer_gl.sphere_geometry_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, scene_data_row_pad, scene_data.sphere_data_rows, 0, GL_RGBA, GL_FLOAT, scene_data.sphere_geometries);

    glGenTextures(1, &raytracer_gl.sphere_colors_texture);
    glBindTexture(GL_TEXTURE_2D, raytracer_gl.sphere_colors_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, scene_data_row_pad, scene_data.sphere_data_rows, 0, GL_RGB, GL_FLOAT, scene_data.sphere_colors);

    glGenTextures(1, &raytracer_gl.group_children_texture);
    glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_children_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    check_opengl(__FILE__, __LINE__);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, scene_data_row_pad, scene_data.group_data_rows, 0, GL_RG, GL_FLOAT, scene_data.group_children);
    check_opengl(__FILE__, __LINE__);

    glGenTextures(1, &raytracer_gl.group_objects_texture);
    glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_objects_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    check_opengl(__FILE__, __LINE__);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, scene_data_row_pad, scene_data.group_data_rows, 0, GL_RG, GL_FLOAT, scene_data.group_objects);
    check_opengl(__FILE__, __LINE__);

    glGenTextures(1, &raytracer_gl.group_hitmiss_texture);
    glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_hitmiss_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    check_opengl(__FILE__, __LINE__);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RG32F, scene_data_row_pad, scene_data.group_data_rows, 0, GL_RG, GL_FLOAT, scene_data.group_hitmiss);
    check_opengl(__FILE__, __LINE__);

    glGenTextures(1, &raytracer_gl.group_directions_texture);
    glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_directions_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    check_opengl(__FILE__, __LINE__);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, scene_data_row_pad, scene_data.group_data_rows, 0, GL_RGB, GL_FLOAT, scene_data.group_directions);
    check_opengl(__FILE__, __LINE__);

    glGenTextures(1, &raytracer_gl.group_boxmin_texture);
    glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_boxmin_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    check_opengl(__FILE__, __LINE__);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, scene_data_row_pad, scene_data.group_data_rows, 0, GL_RGB, GL_FLOAT, scene_data.group_boxmin);
    check_opengl(__FILE__, __LINE__);

    glGenTextures(1, &raytracer_gl.group_boxmax_texture);
    glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_boxmax_texture);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    check_opengl(__FILE__, __LINE__);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, scene_data_row_pad, scene_data.group_data_rows, 0, GL_RGB, GL_FLOAT, scene_data.group_boxmax);
    check_opengl(__FILE__, __LINE__);

    glBindTexture(GL_TEXTURE_2D, 0);
}

void init_decal_program()
{
    check_opengl(__FILE__, __LINE__);

    decal_vertex_shader = glCreateShader(GL_VERTEX_SHADER);
    const char *string = gDecalVertexShaderText;
    glShaderSource(decal_vertex_shader, 1, &string, NULL);
    glCompileShader(decal_vertex_shader);
    if(!CheckShaderCompile(decal_vertex_shader, "vertex shader"))
	exit(1);

    GLuint fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
    string = gDecalFragmentShaderText;
    glShaderSource(fragment_shader, 1, &string, NULL);
    glCompileShader(fragment_shader);
    if(!CheckShaderCompile(fragment_shader, "fragment shader"))
	exit(1);

    decalprogram = glCreateProgram();
    glAttachShader(decalprogram, decal_vertex_shader);
    glAttachShader(decalprogram, fragment_shader);
    glBindAttribLocation(decalprogram, pos_attrib, "pos");
    glBindAttribLocation(decalprogram, texcoord_attrib, "vtex");
    glLinkProgram(decalprogram);
    if(!CheckProgramLink(decalprogram))
	exit(1);

    modelview_uniform = glGetUniformLocation(decalprogram, "modelview");
    texture_uniform = glGetUniformLocation(decalprogram, "texture_image");
    glUseProgram(decalprogram);
}

void init_decal_geometry(void)
{
    float verts[4][4];
    float texcoords[4][2];

    verts[0][0] = -1;
    verts[0][1] = -1;
    verts[0][2] = 0;
    verts[0][3] = 1;
    verts[1][0] = 1; // gWindowWidth;
    verts[1][1] = -1;
    verts[1][2] = 0;
    verts[1][3] = 1;
    verts[2][0] = -1;
    verts[2][1] = 1; // gWindowHeight;
    verts[2][2] = 0;
    verts[2][3] = 1;
    verts[3][0] = 1; // gWindowWidth;
    verts[3][1] = 1; // gWindowHeight;
    verts[3][2] = 0;
    verts[3][3] = 1;

    texcoords[0][0] = 0;
    texcoords[0][1] = 1;
    texcoords[1][0] = 1;
    texcoords[1][1] = 1;
    texcoords[2][0] = 0;
    texcoords[2][1] = 0;
    texcoords[3][0] = 1;
    texcoords[3][1] = 0;

    glGenVertexArrays(1, &screenquad_vao);
    glBindVertexArray(screenquad_vao);
    glGenBuffers(1, &vert_buffer);
    glGenBuffers(1, &texcoord_buffer);

    check_opengl(__FILE__, __LINE__);

    glBindBuffer(GL_ARRAY_BUFFER, vert_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(verts), verts, GL_STATIC_DRAW);
    glVertexAttribPointer(pos_attrib, 4, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(pos_attrib);

    glBindBuffer(GL_ARRAY_BUFFER, texcoord_buffer);
    glBufferData(GL_ARRAY_BUFFER, sizeof(texcoords), texcoords, GL_STATIC_DRAW);
    glVertexAttribPointer(texcoord_attrib, 2, GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(texcoord_attrib);
}

void init()
{
    if(false) {
        printf("GL_RENDERER: \"%s\"\n", glGetString(GL_RENDERER));
        printf("GL_VERSION: \"%s\"\n", glGetString(GL_VERSION));
        printf("GL_VENDOR: \"%s\"\n", glGetString(GL_VENDOR));
        GLint num_extensions;
        glGetIntegerv(GL_NUM_EXTENSIONS, &num_extensions);
        printf("GL_EXTENSIONS: \n");
        for(int i = 0; i < num_extensions; i++)
            printf("    \"%s\"\n", glGetStringi(GL_EXTENSIONS, i));

        GLint max_tex;
        glGetIntegerv(GL_MAX_TEXTURE_SIZE, &max_tex);
        printf("max texture size: %d\n", max_tex);
        glGetIntegerv(GL_MAX_RECTANGLE_TEXTURE_SIZE, &max_tex);
        printf("max texture rectangle size: %d\n", max_tex);
    }

    init_decal_program();
    init_decal_geometry();
    load_scene_data(gWorld, raytracer_gl);
}

void noredraw(void) {}

void redraw(void)
{
    glClearColor(1, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    if(render_mode == CPU) {

        timeval t1, t2;
        gettimeofday(&t1, NULL);

        trace_image(gWindowWidth, gWindowHeight, gWindowHeight / (1.0f * gWindowWidth), gImage, gWorld, light_dir);
        gettimeofday(&t2, NULL);
        long long micros = (t2.tv_sec - t1.tv_sec) * 1000000 + t2.tv_usec - t1.tv_usec;
        printf("trace: %lld millis\n", micros / 1000);

        glUseProgram(decalprogram);

        glActiveTexture(GL_TEXTURE0 + 0);
        glBindTexture(GL_TEXTURE_2D, 0);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, gWindowWidth, gWindowHeight, 0, GL_RGB, GL_UNSIGNED_BYTE, gImage);

        glBindBuffer(GL_ARRAY_BUFFER, vert_buffer);
        glVertexAttribPointer(pos_attrib, 4, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(pos_attrib);

        glBindBuffer(GL_ARRAY_BUFFER, texcoord_buffer);
        glVertexAttribPointer(texcoord_attrib, 2, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(texcoord_attrib);

        glUniformMatrix4fv(modelview_uniform, 1, GL_FALSE, mat4_identity);

        glBindVertexArray(screenquad_vao);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    } else if(render_mode == FRAGMENT_SHADER) {

        glUseProgram(raytracer_gl.program);
        check_opengl(__FILE__, __LINE__);

        int which_texture = 0;

        glActiveTexture(GL_TEXTURE0 + which_texture);
        glBindTexture(GL_TEXTURE_2D, raytracer_gl.sphere_geometry_texture);
        glUniform1i(raytracer_gl.sphere_geometry_uniform, which_texture);
        which_texture++;

        glActiveTexture(GL_TEXTURE0 + which_texture);
        glBindTexture(GL_TEXTURE_2D, raytracer_gl.sphere_colors_texture);
        glUniform1i(raytracer_gl.sphere_colors_uniform, which_texture);
        which_texture++;

        glActiveTexture(GL_TEXTURE0 + which_texture);
        glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_objects_texture);
        glUniform1i(raytracer_gl.group_objects_uniform, which_texture);
        which_texture++;

        glActiveTexture(GL_TEXTURE0 + which_texture);
        glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_hitmiss_texture);
        glUniform1i(raytracer_gl.group_hitmiss_uniform, which_texture);
        which_texture++;

        glActiveTexture(GL_TEXTURE0 + which_texture);
        glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_children_texture);
        glUniform1i(raytracer_gl.group_children_uniform, which_texture);
        which_texture++;

        glActiveTexture(GL_TEXTURE0 + which_texture);
        glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_directions_texture);
        glUniform1i(raytracer_gl.group_directions_uniform, which_texture);
        which_texture++;

        glActiveTexture(GL_TEXTURE0 + which_texture);
        glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_boxmin_texture);
        glUniform1i(raytracer_gl.group_boxmin_uniform, which_texture);
        which_texture++;

        glActiveTexture(GL_TEXTURE0 + which_texture);
        glBindTexture(GL_TEXTURE_2D, raytracer_gl.group_boxmax_texture);
        glUniform1i(raytracer_gl.group_boxmax_uniform, which_texture);
        which_texture++;

        glUniform1i(raytracer_gl.tree_root_uniform, scene_data.tree_root);
        glUniform1i(raytracer_gl.sphere_data_row_size_uniform, scene_data_row_pad);
        glUniform1i(raytracer_gl.sphere_data_rows_uniform, scene_data.sphere_data_rows);
        glUniform1i(raytracer_gl.group_data_row_size_uniform, scene_data_row_pad);
        glUniform1i(raytracer_gl.group_data_rows_uniform, scene_data.group_data_rows);
        glUniform3f(raytracer_gl.eye_uniform, gWorld->cam.eye.x, gWorld->cam.eye.y, gWorld->cam.eye.z);

        glUniformMatrix4fv(raytracer_gl.camera_matrix_uniform, 1, GL_FALSE, gWorld->camera_matrix);
        glUniformMatrix4fv(raytracer_gl.camera_normal_matrix_uniform, 1, GL_FALSE, gWorld->camera_normal_matrix);
        glUniformMatrix4fv(raytracer_gl.object_matrix_uniform, 1, GL_FALSE, gWorld->object_matrix);
        glUniformMatrix4fv(raytracer_gl.object_inverse_uniform, 1, GL_FALSE, gWorld->object_inverse);
        glUniformMatrix4fv(raytracer_gl.object_normal_matrix_uniform, 1, GL_FALSE, gWorld->object_normal_matrix);
        glUniformMatrix4fv(raytracer_gl.object_normal_inverse_uniform, 1, GL_FALSE, gWorld->object_normal_inverse);

        glUniform1f(raytracer_gl.fov_uniform, gWorld->cam.fov);
        glUniform1f(raytracer_gl.aspect_uniform, gWindowHeight / (1.0f * gWindowWidth));
        glUniform3f(raytracer_gl.background_color_uniform, gWorld->background.x, gWorld->background.y, gWorld->background.z);

        glBindBuffer(GL_ARRAY_BUFFER, vert_buffer);
        glVertexAttribPointer(pos_attrib, 4, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(pos_attrib);

        glBindBuffer(GL_ARRAY_BUFFER, texcoord_buffer);
        glVertexAttribPointer(texcoord_attrib, 2, GL_FLOAT, GL_FALSE, 0, 0);
        glEnableVertexAttribArray(texcoord_attrib);

        glUniformMatrix4fv(raytracer_gl.modelview_uniform, 1, GL_FALSE, mat4_identity);
        vec4 light_dir4(light_dir.x, light_dir.y, light_dir.z, 0.0);

        glUniform4fv(raytracer_gl.light_dir_uniform, 1, (GLfloat*)&light_dir4);

        glBindVertexArray(screenquad_vao);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
        check_opengl(__FILE__, __LINE__);

    } else /* OPENGL_GEOMETRY */ {

#if 0
        set projection from fov
        construct modelview
            viewpoint and rotation
        need vertex shader - pass transformed vertex, normal
        need fragment shader - normalize normal, do lighting 
        need sphere geometry
        for all spheres
            set transformation matrix
            set color
            draw sphere
#endif

    }

    glutSwapBuffers();
    check_opengl(__FILE__, __LINE__);

    timeval now;
    gettimeofday(&now, NULL);
    long long micros = (now.tv_sec - prev_frame_time.tv_sec) * 1000000 + now.tv_usec - prev_frame_time.tv_usec;
    printf("fps: %f (estimated %lld millis)\n", 1000000.0 / micros, micros / 1000);
    prev_frame_time = now;
}

void keyup(unsigned char key, int x, int y)
{
}

enum {
    MOVE_OBJECT,
    MOVE_LIGHT,
} motion_target = MOVE_OBJECT;

/* snapshot the whole frontbuffer.  Return -1 on error, 0 on success. */
int screenshot(const char *colorName, const char *alphaName)
{
    static unsigned char *pixels;
    static size_t pixelsSize;
    GLint viewport[4];
    int i;
    FILE *fp;
    GLint prevReadBuf;

    glPixelStorei(GL_PACK_ALIGNMENT, 1);
    glGetIntegerv(GL_READ_BUFFER, &prevReadBuf);
    glReadBuffer(GL_FRONT);

    glGetIntegerv(GL_VIEWPORT, viewport);
    if(pixelsSize < (unsigned int)(viewport[2] * viewport[3] * 3)) {
	pixelsSize = (unsigned int)(viewport[2] * viewport[3] * 3);
	pixels = (unsigned char *)realloc(pixels, pixelsSize);
	if(pixels == NULL) {
	    fprintf(stderr, "snapshot: couldn't allocate %zd bytes for"
		" screenshot.\n", pixelsSize);
	    return -1;
	}
    }

	/* color ---------------------------------------- */
    if(colorName != NULL) {
	if((fp = fopen(colorName, "wb")) == NULL) {
	    fprintf(stderr, "snapshot: couldn't open \"%s\".\n", colorName);
	    return -1;
	}
	glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3],
	    GL_RGB, GL_UNSIGNED_BYTE, pixels);
	fprintf(fp, "P6 %d %d 255\n", viewport[2], viewport[3]);
	for(i = viewport[3] - 1; i >= 0; i--)
	    fwrite(pixels + viewport[2] * 3 * i, 3, viewport[2], fp);
	fclose(fp);
    }

	/* alpha ---------------------------------------- */
    if(alphaName != NULL) {
	if((fp = fopen(alphaName, "wb")) == NULL) {
	    fprintf(stderr, "snapshot: couldn't open \"%s\".\n", alphaName);
	    return -1;
	}
	glReadPixels(viewport[0], viewport[1], viewport[2], viewport[3],
	    GL_ALPHA, GL_UNSIGNED_BYTE, pixels);
	fprintf(fp, "P5 %d %d 255\n", viewport[2], viewport[3]);
	for(i = viewport[3] - 1; i >= 0; i--)
	    fwrite(pixels + viewport[2] * i, 1, viewport[2], fp);
	fclose(fp);
    }

    glReadBuffer(prevReadBuf);

    return 0;
}


void keydown(unsigned char key, int x, int y)
{
    switch(key) {
	case 'q': case 'Q': case '\033':
	    exit(0);
	    break;
        case 'o':
            motion_target = MOVE_OBJECT;
            break;
        case 'l':
            motion_target = MOVE_LIGHT;
            break;
        case 'c':
            render_mode = CPU;
            glutPostRedisplay();
            break;
        case 'g':
            render_mode = FRAGMENT_SHADER;
            glutPostRedisplay();
            break;

        case 's':
            screenshot("color.ppm", NULL);
            break;

        case 'p':
            printf("%f;%f;%f\n", gWorld->cam.eye.x, gWorld->cam.eye.y, gWorld->cam.eye.z);
            printf("%f;%f;%f\n",
                to_degrees(gWorld->cam.yaw),
                to_degrees(gWorld->cam.pitch),
                to_degrees(gWorld->cam.roll));
    }
}

static int ox, oy;

void button(int b, int state, int x, int y)
{
    if(state == 1) {
    } else {
	ox = x;
	oy = y;
    }
    glutPostRedisplay();
}


void motion(int x, int y)
{
    int dx, dy;

    dx = x - ox;
    dy = y - oy;

    if(glutGetModifiers() & GLUT_ACTIVE_SHIFT) {
        zoom *= exp(log(5.0) / gWindowHeight / 2 * -dy);
    } else {
        if(motion_target == MOVE_OBJECT) {
            trackball_motion(object_rotation, (dx / (float)gWindowWidth), (dy / (float)gWindowHeight), object_rotation);
        } else {
            trackball_motion(light_rotation, (dx / (float)gWindowWidth), (dy / (float)gWindowHeight), light_rotation);
        }
    }
    update_view_params(gWorld, zoom);
    update_light();
    glutPostRedisplay();

    ox = x;
    oy = y;
}

void reshape(int width, int height)
{
check_opengl(__FILE__, __LINE__);
    static bool initialized = false;

    if(!initialized) {
        init();
        initialized = true;
        glutDisplayFunc(redraw);
        // glutIdleFunc(redraw);
        gettimeofday(&prev_frame_time, NULL);
    }

    glViewport(0, 0, width, height);

    delete[] gImage;
    gImage = new unsigned char[(((width * 3) + 3) & ~3) * height];

    gWindowWidth = width;
    gWindowHeight = height;

check_opengl(__FILE__, __LINE__);
    glutPostRedisplay();
}

void usage(char *progname)
{
    fprintf(stderr, "usage: %s <inputfilename>\n", progname);
}

int main(int argc, char *argv[])
{
    glutInitContextVersion (3, 1);
#if 0
    if(debug_gl)
        glutInitContextFlags (GLUT_FORWARD_COMPATIBLE | GLUT_DEBUG);
    else
        glutInitContextFlags (GLUT_FORWARD_COMPATIBLE);
#endif

    glutInit(&argc, argv);
    glutInitWindowSize(512,512);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
    glutCreateWindow("ray1 interactive program");
    glutDisplayFunc(noredraw);
    // glutIdleFunc(noredraw);
    glutKeyboardFunc(keydown);
    glutKeyboardUpFunc(keyup);
    glutMotionFunc(motion);
    glutMouseFunc(button);
    glutReshapeFunc(reshape);

    if(argc < 2) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if((!strcmp(argv[1], "-h")) || (!strcmp(argv[1], "--help"))) {
        usage(argv[0]);
        exit(EXIT_FAILURE);
    }

    if((gWorld = load_world(argv[1])) == NULL) {
        fprintf(stderr, "Cannot set up world.\n");
        exit(EXIT_FAILURE);
    }
    fprintf(stderr, "loaded\n");

    gWorld->xsub = gWorld->ysub = 1;
    gWorld->cam.fov = to_radians(40.0);
    zoom = gWorld->scene_extent / 2 / sinf(gWorld->cam.fov / 2);

    light_rotation[0] = to_radians(20.0);
    light_rotation[1] = .707;
    light_rotation[2] = -.707;
    light_rotation[3] = 0;

    update_view_params(gWorld, zoom);
    update_light();

    glutMainLoop();

    exit(EXIT_SUCCESS);
}
