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

#include "vectormath.h"

void set_keep_stats(bool keep);
void get_stats(
    unsigned long long *sphere_tests_,
    unsigned long long *sphere_shadings_,
    unsigned long long *sphere_intersections_,
    unsigned long long *group_intersections_);

struct camera { /* Viewpoint specification. */
    vec3 eye; /* Origin of rays */
    float yaw, pitch, roll; /* Also known as Heading, Pitch, and Bank. */
    float fov; /* Entire View angle, left to right. */
};

struct surface_hit
{
    float t;
    vec3 normal;
    vec3 point;
    vec3 color;
};

struct range
{
    float t0, t1;
    range() :
        t0(-std::numeric_limits<float>::max()),
        t1(std::numeric_limits<float>::max())
    {}
    range(const range &r) :
        t0(r.t0),
        t1(r.t1)
    {}
    range(float t0_, float t1_) :
        t0(t0_),
        t1(t1_)
    {}
    range &intersect(const range &r2)
    {
        t0 = std::max(t0, r2.t0);
        t1 = std::min(t1, r2.t1);
        return *this;
    }
    operator bool() { return t0 < t1; }
};

struct sphere
{
    vec3 center; /* Center of sphere */
    float radius; /* (radius) */
    vec3 color;
    bool intersect(const ray& ray, const range& r, surface_hit *hit);
    sphere() {};
    sphere(const vec3& center_, float radius_, const vec3& color_)
    {
        center = center_;
        radius = radius_;
        color = color_;
    }
    ~sphere() {}
};

struct group
{
    vec3 D; /* split direction */

    vec3 boxmin, boxmax;

    group *negative, *positive;

    sphere *spheres;
    int start;
    unsigned int count;

    bool intersect(const ray& ray, const range& r, surface_hit *hit);

    group(sphere *spheres_, group *neg, group *pos, const vec3& direction, const vec3& boxmin_, const vec3& boxmax_);
    group(sphere *spheres_, int start_, unsigned int count_);
    ~group();

    int my_index;
};

struct world { /* world specification */
    char *Title;

    int sphere_count;
    sphere *spheres; // base spheres, only traced through "root"

    group *root;
    vec3 background; /* Background color */

    vec3 scene_center;
    float scene_extent;

    camera cam; /* Viewpoint */
    int xsub, ysub; /* Yaw and Pitch divisions per pixel */

    world() :
        Title(0),
        spheres(0),
        root(0)
    {};
    ~world()
    {
        delete root;
        delete[] spheres;
        delete[] Title;
    }

    float camera_matrix[16];
    float camera_normal_matrix[16];

    float object_matrix[16];
    float object_inverse[16];
    float object_normal_matrix[16];
    float object_normal_inverse[16];

};

world *load_world(char *fname);
void trace_image(int width, int height, float aspect, unsigned char *image, const world* Wd, const vec3& light_dir);

struct scene_shader_data {
    unsigned int row_pad;
    unsigned int sphere_count;
    float *sphere_geometries; // array of float4: {x, y, z, radius}
    float *sphere_colors; // array of float3: {r, g, b}
    unsigned int sphere_data_rows;

    int group_count;
    int tree_root;
    float *group_boxmin; // array of float3: {x, y, z}
    float *group_boxmax; // array of float3: {x, y, z}
    float *group_directions; // array of float3: {x, y, z}
    float *group_children; // array of float, [0]>65535.0  if leaf

    float *group_hitmiss; // array of int2
    // [0] = node to go to on hit
    // [1] = node to go to on miss
    // > 65535.0 on terminate

    float *group_objects; // array of {start, count}, count==0 if not leaf
    unsigned int group_data_rows;

    scene_shader_data();
    ~scene_shader_data();
};

void get_shader_data(world *w, scene_shader_data &data, unsigned int row_pad);
