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

#include <cerrno>
#include <ctime>
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cassert>
#include <memory>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>
#ifdef USE_OPENMP
#include <omp.h>
#endif
#include "world.h"

const float NO_t = -1.0;

bool keep_stats = false;
static unsigned long long sphere_tests = 0;
static unsigned long long sphere_shadings = 0;
static unsigned long long sphere_intersections = 0;
static unsigned long long group_intersections = 0;

void set_keep_stats(bool keep)
{
    keep_stats = keep;
}

void get_stats(
    unsigned long long *sphere_tests_,
    unsigned long long *sphere_shadings_,
    unsigned long long *sphere_intersections_,
    unsigned long long *group_intersections_)
{
    *sphere_tests_ = sphere_tests;
    *sphere_shadings_ = sphere_shadings;
    *sphere_intersections_ = sphere_intersections;
    *group_intersections_ = group_intersections;
}

bool sphere_possible(vec3 center, float radius, const ray& ray)
{
    if(keep_stats) sphere_tests++;
    vec3 diff = vec3_subtract(ray.o, center);

    float b = 2 * vec3_dot(ray.d, diff);

    float radicand = b * b - 4 * (vec3_dot(diff, diff) - radius * radius);
    if(radicand < 0)
        return false;

    return true;
}

range box_intersect(const vec3& boxmin, const vec3& boxmax, const ray& theray)
{
    range r = range(-100000000.0, 100000000.0);

    float t0, t1;

    t0 = (boxmin.x - theray.o.x) / theray.d.x;
    t1 = (boxmax.x - theray.o.x) / theray.d.x;
    if(theray.d.x >= 0.0)
        r.intersect(range(t0, t1));
    else
        r.intersect(range(t1, t0));

    t0 = (boxmin.y - theray.o.y) / theray.d.y;
    t1 = (boxmax.y - theray.o.y) / theray.d.y;
    if(theray.d.y >= 0.0)
        r.intersect(range(t0, t1));
    else
        r.intersect(range(t1, t0));

    t0 = (boxmin.z - theray.o.z) / theray.d.z;
    t1 = (boxmax.z - theray.o.z) / theray.d.z;
    if(theray.d.z >= 0.0)
        r.intersect(range(t0, t1));
    else
        r.intersect(range(t1, t0));

    return r;
}

range sphere_intersect(vec3 center, float radius, const ray& ray)
{
    if(keep_stats) sphere_intersections++;
    vec3 diff = vec3_subtract(ray.o, center);

    float b = 2 * vec3_dot(ray.d, diff);

    float radicand = b * b - 4 * (vec3_dot(diff, diff) - radius * radius);
    if(radicand < 0)
        return range(std::numeric_limits<float>::max(), -std::numeric_limits<float>::max());

    float t0 = (-b - sqrtf(radicand)) / 2;
    float t1 = (-b + sqrtf(radicand)) / 2;
    return range(t0, t1);
}

float sphere_intersect(vec3 center, float radius, const ray& ray, const range &r)
{
    if(keep_stats) sphere_intersections++;
    vec3 diff = vec3_subtract(ray.o, center);

    float b = 2 * vec3_dot(ray.d, diff);

    float radicand = b * b - 4 * (vec3_dot(diff, diff) - radius * radius);
    if(radicand < 0)
        return NO_t;

    float t1 = (-b + sqrtf(radicand)) / 2;

    if(t1 < r.t0)
         return NO_t;

    float t0 = (-b - sqrtf(radicand)) / 2;

    if(t0 > r.t1)
        return NO_t;

    if(t0 < r.t0)
        return t1;

    return t0;
}

ray transform_ray(const ray& r, const float matrix[16], const float normal_matrix[16])
{
    ray result;

    vec4 o1(r.o.x, r.o.y, r.o.z, 1.0), o2;
    vec4 d1(r.d.x, r.d.y, r.d.z, 0.0), d2;

    o2 = mat4_mult_vec4(matrix, o1);
    d2 = mat4_mult_vec4(normal_matrix, d1);

    // XXX ignoring w
    result.o = vec3(o2.x, o2.y, o2.z);
    result.d = vec3(d2.x, d2.y, d2.z);

    return result;
}

bool sphere::intersect(const ray& ray, const range& r, surface_hit* hit)
{
    float t = sphere_intersect(center, radius, ray, r);

    if(t == NO_t)
        return false;

    if(t > hit->t)
        return false;

    if(keep_stats) sphere_shadings++;

    vec3 point = vec3_add(ray.o, vec3_scale(ray.d, t));
    // snap to sphere surface
    vec3 to_surface = vec3_subtract(point, center);
    float distance = sqrtf(vec3_dot(to_surface, to_surface));
    hit->point = vec3_add(center, vec3_scale(to_surface, radius / distance));
    hit->normal = vec3_divide(to_surface, radius);
    hit->color = color;
    hit->t = t;

    return true;
}

group::group(sphere *spheres_, group *neg, group *pos, const vec3& direction, const vec3& boxmin_, const vec3& boxmax_) :
    D(direction),
    boxmin(boxmin_),
    boxmax(boxmax_),
    negative(neg),
    positive(pos),
    spheres(spheres_),
    start(0),
    count(0)
{
}

group::group(sphere *spheres_, int start_, unsigned int count_) :
    negative(NULL),
    positive(NULL),
    spheres(spheres_),
    start(start_),
    count(count_)
{
    box_init(boxmin, boxmax);
    for(unsigned int i = 0; i < count; i++) {
        add_sphere(&boxmin, &boxmax, spheres[start + i].center, spheres[start + i].radius);
    }
}

group::~group()
{
    delete negative;
    delete positive;
}

bool tree_intersect_stack(group *root, const ray& ray, const range& r, surface_hit *hit)
{
    bool have_hit = false;
    group* stack_groups[25];
    int stack_top;

    range r2 = box_intersect(root->boxmin, root->boxmax, ray);
    r2.intersect(r);
    if(!r2 || (r2.t0 > hit->t))
        return false;

    stack_groups[0] = root;
    stack_top = 0;

    while(stack_top >= 0) {
        group *g = stack_groups[stack_top--];

        if(g->negative != NULL) {
            group *g1, *g2;

            if(vec3_dot(ray.d, g->D) > 0) {
                g1 = g->negative;
                g2 = g->positive;
            } else {
                g1 = g->positive;
                g2 = g->negative;
            }

            range r3 = box_intersect(g2->boxmin, g2->boxmax, ray);
            r3.intersect(r);
            if(r3 && (r3.t0 < hit->t))
                stack_groups[++stack_top] = g2;

            r3 = box_intersect(g1->boxmin, g1->boxmax, ray);
            r3.intersect(r);
            if(r3 && (r3.t0 < hit->t))
                stack_groups[++stack_top] = g1;

        } else {

            for(unsigned int i = 0; i < g->count; i++) {
                if(g->spheres[g->start + i].intersect(ray, r, hit))
                    have_hit = true;
            }
        }
    }

    return have_hit;
}

bool set_bad_hit(surface_hit *hit, float r, float g, float b)
{
    hit->color = vec3(r, g, b);
    hit->normal = vec3(0, 0, -1);
    hit->t = 1.0f;
    return true;
}

bool group::intersect(const ray& ray, const range &r, surface_hit* hit)
{
    return tree_intersect_stack(this, ray, r, hit);
}

static const bool cast_shadows = true;

vec3 trace(const ray& theray, const world* world, const vec3& light_dir)
{
    range r(0, std::numeric_limits<float>::max());
    surface_hit hit;
    hit.t = std::numeric_limits<float>::max();

    ray object_ray = transform_ray(theray, world->object_matrix, world->object_normal_matrix);
    bool have_hit = world->root->intersect(object_ray, r, &hit);
    if(!have_hit)
        return world->background;

    vec4 n1(hit.normal.x, hit.normal.y, hit.normal.z, 0.0), n2;
    vec4 p1(hit.point.x, hit.point.y, hit.point.z, 1.0), p2;
    n2 = mat4_mult_vec4(world->object_normal_inverse, n1);
    p2 = mat4_mult_vec4(world->object_inverse, p1);
    vec3 normal = vec3(n2.x, n2.y, n2.z);
    vec3 point = vec3(p2.x, p2.y, p2.z);

    float brightness;

    float shadowed = false;

    if(cast_shadows)
    {
        struct ray shadowray;
        shadowray.o = vec3_add(point, vec3_scale(normal, .0001));
        shadowray.d = light_dir;
        surface_hit shadow_hit;
        shadow_hit.t = std::numeric_limits<float>::max();
        ray object_ray = transform_ray(shadowray, world->object_matrix, world->object_normal_matrix);
        shadowed = world->root->intersect(object_ray, range(0, std::numeric_limits<float>::max()), &shadow_hit);
    }

    if(shadowed) {

        brightness = .1;

    } else {

        brightness = vec3_dot(normal, light_dir);
        if(brightness < 0)
            brightness = 0;
        if(brightness > 1)
            brightness = 1;
    }

    return vec3_scale(hit.color, brightness);
}

vec3 make_eye_ray(float u, float v, float aspect, float fov)
{
    float eye_alpha = fov * (u - .5);
    float eye_beta = fov * (v - .5) * aspect;

    // eye space ray
    vec3 eye;
    eye.x = sinf(eye_alpha) * cosf(eye_beta);
    eye.y = sinf(eye_beta);
    eye.z = cosf(eye_alpha) * cosf(eye_beta);

    return eye;
}

void trace_image(int width, int height, float aspect, unsigned char *image, const world* world, const vec3& light_dir)
{
#pragma omp parallel for schedule(dynamic)
    for(int yloop = 0; yloop < height; yloop++) {
        unsigned char *row = image + (height - yloop - 1) * (((width * 3) + 3) & ~3);
        for(int xloop = 0; xloop < width; xloop++) {
            // printf("%d, %d\n", xloop, yloop);
            int cols = 0;
            vec3 color(0, 0, 0);
            for(int qky = 0; qky < world->ysub; qky++) {
                for(int qkx = 0; qkx < world->xsub; qkx++) {
                    float u = ((xloop + qkx / (float)world->xsub) + .5) / width;
                    float v = ((yloop + qky / (float)world->ysub) + .5) / height;
                    ray eye_ray;
                    eye_ray.d = make_eye_ray(u, v, aspect, world->cam.fov);
                    eye_ray.o = vec3(0, 0, 0);
                    ray world_ray = transform_ray(eye_ray, world->camera_matrix, world->camera_normal_matrix);
                    vec3 sample = trace(world_ray, world, light_dir);
                    color = vec3_add(color, sample);
                    ++cols;
                }
            }
            vec3 final_color = vec3_divide(color, cols);
            unsigned char *pixel = row + xloop * 3;
            pixel[0] = final_color.x * 255;
            pixel[1] = final_color.y * 255;
            pixel[2] = final_color.z * 255;
        }
    }
}

// max-depth has to be less than 24 to work with fragment shader short
// stack implementation
static int bvh_max_depth = 23;
static unsigned int bvh_leaf_max = 8;
static bool bvh_split_median = false;

static void initialize_BVH_parameters() __attribute__((constructor));
static void initialize_BVH_parameters()
{
    if(getenv("BVH_MAX_DEPTH") != 0) {
        bvh_max_depth = atoi(getenv("BVH_MAX_DEPTH"));
        fprintf(stderr, "BVH max depth set to %d\n", bvh_max_depth);
    }
    if(getenv("BVH_LEAF_MAX") != 0) {
        bvh_leaf_max = atoi(getenv("BVH_LEAF_MAX"));
        fprintf(stderr, "BVH max spheres per leaf set to %d\n", bvh_leaf_max);
    }
    if(getenv("BVH_SPLIT_MEDIAN") != 0) {
        bvh_split_median = true;
        fprintf(stderr, "BVH will split at median\n");
    }
}

int total_treed = 0;
time_t previous;
int bvh_level_counts[64];
int bvh_leaf_size_counts[64];
int bvh_node_count = 0;
int bvh_leaf_count = 0;

FILE *bvh_spheres_debug_output = NULL;

void print_tree_stats()
{
    fprintf(stderr, "%d bvh nodes\n", bvh_node_count);
    fprintf(stderr, "%d of those are leaves\n", bvh_leaf_count);
    for(int i = 0; i < bvh_max_depth + 1; i++) {
        fprintf(stderr, "bvh level %2d: %6d nodes\n", i, bvh_level_counts[i]);
    }
    int largest_leaf_count = 63;
    while((largest_leaf_count > 0) && (bvh_leaf_size_counts[largest_leaf_count]) == 0)
        largest_leaf_count--;

    for(int i = 0; i <= largest_leaf_count; i++) {
        fprintf(stderr, "%2d spheres in %6d leaves\n", i, bvh_leaf_size_counts[i]);
    }
    if(bvh_leaf_size_counts[63] > 0)
        fprintf(stderr, "63 or more spheres in %6d leaves\n", bvh_leaf_size_counts[63]);
}

struct sphere_sorter
{
    vec3 split_plane_normal;
    sphere_sorter(vec3 split_plane_normal_) :
        split_plane_normal(split_plane_normal_)
    {}
    bool operator() (const sphere& s1, const sphere& s2)
    {
        vec3 diff = vec3_subtract(s2.center, s1.center);
        return vec3_dot(diff, split_plane_normal) > 0;
    }
};

group* make_tree(sphere* spheres, int start, unsigned int count, int level = 0)
{
    if(level == 0) {
        previous = time(NULL);
        if(getenv("DEBUG_BVH") != NULL) {
            bvh_spheres_debug_output = fopen("bvh.r9", "w");
        }
    }
    if(time(NULL) > previous) {
        fprintf(stderr, "total treed = %d\n", total_treed);
        previous = time(NULL);
    }

    if((level >= bvh_max_depth) || count <= bvh_leaf_max) {
        total_treed += count;
        group* g = new group(spheres, start, count);
        if(bvh_spheres_debug_output != NULL) {
            // fprintf(bvh_spheres_debug_output, "sphere 7 %f %f %f %f 1 1 1 # %d\n", g->radius, g->center.x, g->center.y, g->center.z, level );
            if(level == 0) {
                fclose(bvh_spheres_debug_output);
                bvh_spheres_debug_output = NULL;
            }
        }
        bvh_leaf_size_counts[std::min(63U, count)]++;
        bvh_leaf_count++;
        bvh_level_counts[level]++;
        bvh_node_count++;
        return g;
    }

    vec3 split_pivot;
    vec3 split_plane_normal;

    // find bounding box
    vec3 boxmin, boxmax;
    box_init(boxmin, boxmax);
    for(unsigned int i = 0; i < count; i++) {
        sphere &s = spheres[start + i];
        vec3 spheremin = vec3_subtract(s.center, vec3(s.radius));
        vec3 spheremax = vec3_add(s.center, vec3(s.radius));
        box_extend(boxmin, boxmax, spheremin, spheremax);
    }

    split_pivot = vec3_scale(vec3_add(boxmin, boxmax), .5);

    vec3 boxdim = vec3_subtract(boxmax, boxmin);
    if(boxdim.x > boxdim.y && boxdim.x > boxdim.z) {
        split_plane_normal = vec3(1, 0, 0);
    } else if(boxdim.y > boxdim.z) {
        split_plane_normal = vec3(0, 1, 0);
    } else {
        split_plane_normal = vec3(0, 0, 1);
    }

    // XXX output split plane to BVH debug file

    int startA;
    int countA;
    int startB;
    int countB;

    if(!bvh_split_median) {

        int s1 = start - 1;
        int s2 = start + count;

        do {
            // from start to s1, not including s1, is negative
            // from s2 to start + count - 1 is positive
            do {
                s1 += 1;
            } while((s1 < s2) && vec3_dot(vec3_subtract(spheres[s1].center, split_pivot), split_plane_normal) < 0);

            // If there wasn't a positive sphere before s2, done.
            if(s1 >= s2)
                break;

            // s1 is now location of lowest positive sphere 

            do {
                s2 -= 1;
            } while((s1 < s2) && vec3_dot(vec3_subtract(spheres[s2].center, split_pivot), split_plane_normal) >= 0);


            // If there wasn't a negative sphere between s1 and s2, done
            if(s1 >= s2)
                break;

            // s2 is now location of highest negative sphere 
            std::swap(spheres[s1], spheres[s2]);
        } while(true);

        // s1 is the first of the positive spheres
        startA = start;
        countA = s1 - startA;
        startB = s1;
        countB = start + count - s1;

    } else {

        sphere_sorter sorter(split_plane_normal);
        std::sort(spheres + start, spheres + start + count - 1, sorter);

        startA = start;
        countA = count / 2;
        startB = startA + countA;
        countB = count - countA;

    }

    group *g;

    if(countA > 0 && countB > 0) {

        // get a tighter bound around children than hierarchical bounds
        vec3 boxmin, boxmax;
        box_init(boxmin, boxmax);
        for(unsigned int i = 0; i < count; i++) {
            add_sphere(&boxmin, &boxmax, spheres[start + i].center, spheres[start + i].radius);
        }

        // construct children
        group *g1 = make_tree(spheres, startA, countA, level + 1);
        g = new group(spheres, g1, NULL, split_plane_normal, boxmin, boxmax);
        group *g2 = make_tree(spheres, startB, countB, level + 1);
        g->positive = g2;
        bvh_level_counts[level]++;
        bvh_node_count++;

    } else {

        total_treed += count;
        fprintf(stderr, "Leaf node at %d, %u spheres, total %d\n", level, count, total_treed);
        g = new group(spheres, start, count);
        bvh_leaf_size_counts[std::min(63U, count)]++;
        bvh_leaf_count++;
        bvh_level_counts[level]++;
        bvh_node_count++;
    }

    if(bvh_spheres_debug_output != NULL) {
        // fprintf(bvh_spheres_debug_output, "sphere 7 %f %f %f %f 1 1 1 # %d\n", g->radius, g->center.x, g->center.y, g->center.z, level );
        if(level == 0) {
            fclose(bvh_spheres_debug_output);
            bvh_spheres_debug_output = NULL;
        }
    }

    return g;
}

char *getstr(FILE *fp)
{
     static char str[256];
     char *s = str;
     char c;

     while(1) {
	c = fgetc(fp);
	if(c == ';' || c == '\n')
	     break;
	*s++ = c;
     }
     *s = '\0';

     return str;
}

int getint(FILE *fp, int *num) /* Call getstr() and */
                               /*  convert result to integer */
{
    char *str; 

    if((str = getstr(fp)) == NULL)
        return 0;
    *num = atoi(str);
    return 1;
}

float fltget(FILE *fp, float *num) /* Call getstr() and convert result to*/
                                   /*  floating point */
{
    char *str;

    if((str = getstr(fp)) == NULL)
        return 0;
    *num = atof(str);
    return 1;
}

struct scoped_FILE
{
    FILE *fp;
    scoped_FILE(FILE *fp_) : fp(fp_) {}
    ~scoped_FILE() {if(fp) fclose(fp);}
    operator FILE*() { return fp; }
};

world *load_world(char *fname) // Get world and return pointer.
{
    char *inpstr;
    int wkd;

    std::auto_ptr<world> w(new world);

    time_t prev = time(NULL);

    scoped_FILE fp(fopen(fname, "r"));
    if(fp == NULL) {
        fprintf(stderr, "Cannot open file %s for input.\nE#%d\n", fname, errno);
        return NULL;
    }
    if((inpstr = getstr(fp)) == NULL) {
        fprintf(stderr, "Cannot read title.\n");
        return NULL;
    }
    if(strcmp(inpstr, ".") != 0) {
        w->Title = NULL;
    } else {
        w->Title = new char[strlen(inpstr) + 1];
        strcpy(w->Title, inpstr);
    }

    // lace now ignored
    int lace;
    if(!getint(fp, &lace)) {
        fprintf(stderr, "*!LACE\n");
        return NULL;
    }

    int wdth, lnth; // now ignored
    if(!getint(fp, &wdth)) {
        fprintf(stderr, "*!WDTH\n");
        return NULL;
    }
    if(!getint(fp, &lnth)) {
        fprintf(stderr, "*!LNTH\n");
        return NULL;
    }

    int xsub, ysub; // now ignored
    if(!getint(fp, &xsub)) {
        fprintf(stderr, "*!xsub\n");
        return NULL;
    }
    if(!getint(fp, &ysub)) {
        fprintf(stderr, "*!ysub\n");
        return NULL;
    }
    w->xsub = 1;
    w->ysub = 1;
    
    int r, g, b, brt;
    if(!getint(fp, &r)) {
        fprintf(stderr, "*!backgroundr\n");
        return NULL;
    }
    if(!getint(fp, &g)) {
        fprintf(stderr, "*!backgroundg\n");
        return NULL;
    }
    if(!getint(fp, &b)) {
        fprintf(stderr, "*!backgroundb\n");
        return NULL;
    }
    if(!getint(fp, &brt)) {
        fprintf(stderr, "*!backgroundb\n");
        return NULL;
    }
    w->background.x = r / 16.0 * brt / 100.0;
    w->background.y = g / 16.0 * brt / 100.0;
    w->background.z = b / 16.0 * brt / 100.0;

    int df;
    if(!getint(fp, &df)) {
        fprintf(stderr, "*!DIFFUSION\n");
        return NULL;
    }
    if(!getint(fp, &w->sphere_count)) {
        fprintf(stderr, "*!#Sphere\n");
        return NULL;
    }

    prev = time(NULL);

    w->spheres = new sphere[w->sphere_count];
    for(int i = 0; i < w->sphere_count; i++) {
        if(time(NULL)  > prev) {
            prev = time(NULL);
            fprintf(stderr, "loaded %u spheres\n", i);
        }
        float x, y, z, radius;
        wkd = fltget(fp, &x) && fltget(fp, &y);
        wkd = (wkd && fltget(fp, &z) && fltget(fp, &radius));

        int r, g, b, brt;
        wkd = (wkd && getint(fp, &r) && getint(fp, &g));
        wkd = (wkd && getint(fp, &b) && getint(fp, &brt));

        if(!wkd) {
            fprintf(stderr, "*!Sphere #%d\n", i);
            return NULL;
        }

        vec3 color(r / 16.0 * brt / 100.0,
            g / 16.0 * brt / 100.0, b / 16.0 * brt / 100.0);

        w->spheres[i] = sphere(vec3(x, y, z), radius, color);
    }

    w->scene_center = w->spheres[0].center;
    for(int i = 1; i < w->sphere_count; i++) {
        w->scene_center = vec3_add(w->spheres[i].center, w->scene_center);
    }
    w->scene_center = vec3_divide(w->scene_center, w->sphere_count);

    w->scene_extent = 0;
    for(int i = 0; i < w->sphere_count; i++) {
        vec3 to_center = vec3_subtract(w->scene_center, w->spheres[i].center);
        float distance = sqrtf(vec3_dot(to_center, to_center)) + w->spheres[i].radius;
        w->scene_extent = std::max(w->scene_extent, distance);
    }
    w->scene_extent *= 2;

    w->root = make_tree(w->spheres, 0, w->sphere_count);
    print_tree_stats();

    wkd = fltget(fp, &w->cam.eye.x) && fltget(fp, &w->cam.eye.y);
    wkd = (wkd && fltget(fp, &w->cam.eye.z) && fltget(fp, &w->cam.yaw));
    wkd = (wkd && fltget(fp, &w->cam.pitch) && fltget(fp, &w->cam.roll));
    wkd = (wkd && fltget(fp, &w->cam.fov));
    if(!wkd) {
        fprintf(stderr, "*!Viewpoint\n");
        return NULL;
    }

    w->cam.pitch = to_radians(w->cam.pitch);
    w->cam.yaw = to_radians(w->cam.yaw);
    w->cam.roll = to_radians(w->cam.roll);
    w->cam.fov = to_radians(w->cam.fov);

    return w.release();
}

int get_node_count(group *g)
{
    int count = 1;
    if(g->negative != NULL) {
        count += get_node_count(g->negative);
        count += get_node_count(g->positive);
    }
    return count;
}

void generate_group_indices(group *g, int starting, int *used_, int max, int rowsize)
{
    assert(starting < max);

    int mine;
    int used;

    if(g->negative != NULL) {

        int neg_used;
        int pos_used;

        generate_group_indices(g->negative, starting, &neg_used, max, rowsize);
        assert(starting + neg_used <= max);

        mine = starting + neg_used;

        generate_group_indices(g->positive, mine + 1, &pos_used, max, rowsize);
        assert(mine + 1 + pos_used <= max);

        used = neg_used + 1 + pos_used;

    } else {

        mine = starting;
        used = 1;
    }

    assert(mine < max);
    // g->my_index = renumber_to_improve_bvh_coherence(mine, rowsize);
    g->my_index = mine;
    *used_ = used;
}

void store_group_data(group *g, scene_shader_data &data)
{
    int mine = g->my_index;

    if(g->negative != NULL) {

        store_group_data(g->negative, data);
        store_group_data(g->positive, data);

        data.group_boxmin[mine * 3 + 0] = g->boxmin.x;
        data.group_boxmin[mine * 3 + 1] = g->boxmin.y;
        data.group_boxmin[mine * 3 + 2] = g->boxmin.z;
        data.group_boxmax[mine * 3 + 0] = g->boxmax.x;
        data.group_boxmax[mine * 3 + 1] = g->boxmax.y;
        data.group_boxmax[mine * 3 + 2] = g->boxmax.z;
        data.group_directions[mine * 3 + 0] = g->D.x;
        data.group_directions[mine * 3 + 1] = g->D.y;
        data.group_directions[mine * 3 + 2] = g->D.z;
        data.group_children[mine * 2 + 0] = g->negative->my_index;
        data.group_children[mine * 2 + 1] = g->positive->my_index;
        data.group_objects[mine * 2 + 0] = 0;
        data.group_objects[mine * 2 + 1] = 0;

    } else {

        data.group_boxmin[mine * 3 + 0] = g->boxmin.x;
        data.group_boxmin[mine * 3 + 1] = g->boxmin.y;
        data.group_boxmin[mine * 3 + 2] = g->boxmin.z;
        data.group_boxmax[mine * 3 + 0] = g->boxmax.x;
        data.group_boxmax[mine * 3 + 1] = g->boxmax.y;
        data.group_boxmax[mine * 3 + 2] = g->boxmax.z;
        data.group_children[mine * 2 + 0] = 0x7fffffff;
        data.group_children[mine * 2 + 1] = 0x7fffffff;
        data.group_objects[mine * 2 + 0] = g->start;
        data.group_objects[mine * 2 + 1] = g->count;
    }
}

void create_hitmiss(group *root, scene_shader_data& data)
{
    group *stack[25];
    group *g = root;
    int stack_top = -1;

    while(g != NULL) {

        unsigned int miss_index;
        if(stack_top == -1)
            miss_index = 0x7fffffffU;
        else
            miss_index = stack[stack_top]->my_index;

        if(g->negative == NULL) {

            data.group_hitmiss[g->my_index * 2 + 0] = miss_index;
            data.group_hitmiss[g->my_index * 2 + 1] = miss_index;
            if(stack_top > -1)
                g = stack[stack_top--];
            else
                g = NULL;

        } else {

            group* g1;
            group* g2;
            if(drand48() > .5) {
                g1 = g->positive;
                g2 = g->negative;
            } else {
                g2 = g->positive;
                g1 = g->negative;
            }

            data.group_hitmiss[g->my_index * 2 + 0] = g1->my_index;
            data.group_hitmiss[g->my_index * 2 + 1] = miss_index;
            stack[++stack_top] = g2;
            g = g1;
        }
    }
}

template <class T>
inline T round_up(T v, unsigned int r)
{
    return ((v + r - 1) / r) * r;
}

void xxx_dump_group_code(group *g, std::string name)
{
    if(g->negative == NULL)
        printf("    group *%s = new group(\"%s\");\n", name.c_str(), name.c_str());
    else {
        std::string s1 = name + "0";
        std::string s2 = name + "1";
        xxx_dump_group_code(g->negative, s1);
        xxx_dump_group_code(g->positive, s2);
        printf("    group *%s = new group(\"%s\", %s, %s);\n", name.c_str(), name.c_str(), s1.c_str(), s2.c_str());
        fflush(stdout);
    }
}

void get_shader_data(world *w, scene_shader_data &data, unsigned int row_pad)
{
    data.row_pad = row_pad;
    data.sphere_count = w->sphere_count;
    data.sphere_data_rows = (data.sphere_count + row_pad - 1) / row_pad;
    data.sphere_geometries = new float[4 * row_pad * data.sphere_data_rows];
    data.sphere_colors = new float[3 * row_pad * data.sphere_data_rows];

    for(int i = 0; i < w->sphere_count; i++) {
        const sphere &s = w->spheres[i];
        data.sphere_geometries[i * 4 + 0] = s.center.x;
        data.sphere_geometries[i * 4 + 1] = s.center.y;
        data.sphere_geometries[i * 4 + 2] = s.center.z;
        data.sphere_geometries[i * 4 + 3] = s.radius;
        data.sphere_colors[i * 3 + 0] = s.color.x;
        data.sphere_colors[i * 3 + 1] = s.color.y;
        data.sphere_colors[i * 3 + 2] = s.color.z;
    }

    data.group_count = get_node_count(w->root);
    data.group_data_rows = (data.group_count + row_pad - 1) / row_pad;
    data.group_directions = new float[3 * row_pad * data.group_data_rows];
    data.group_boxmin = new float[3 * row_pad * data.group_data_rows];
    data.group_boxmax = new float[3 * row_pad * data.group_data_rows];
    data.group_children = new float[2 * row_pad * data.group_data_rows];
    data.group_objects = new float[2 * row_pad * data.group_data_rows];
    data.group_hitmiss = new float[2 * row_pad * data.group_data_rows];

    int used;
    generate_group_indices(w->root, 0, &used, data.group_count, data.row_pad);
    assert(used == data.group_count);
    data.tree_root = w->root->my_index;

    store_group_data(w->root, data);

    create_hitmiss(w->root, data);

    // xxx_dump_group_code(w->root, "g");
}

scene_shader_data::scene_shader_data() :
    sphere_geometries(NULL),
    sphere_colors(NULL),
    group_boxmin(NULL),
    group_boxmax(NULL),
    group_directions(NULL),
    group_children(NULL),
    group_objects(NULL)
{
}

scene_shader_data::~scene_shader_data()
{
    delete[] sphere_colors;
    delete[] sphere_geometries;
    delete[] group_directions;
    delete[] group_children;
    delete[] group_objects;
    delete[] group_boxmin;
    delete[] group_boxmax;
}
